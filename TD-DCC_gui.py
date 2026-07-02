import os
import re
import sys
import logging
import threading
import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from datetime import datetime

import matplotlib
matplotlib.use("TkAgg", force=True)
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import matplotlib.colors as mcolors
from matplotlib import cm
# NEW: Forces the Matplotlib toolbar to default to the current running directory
matplotlib.rcParams['savefig.directory'] = os.getcwd()

# MDAnalysis for trajectory I/O (replaces MDTraj / lib.* dependencies)
try:
    import MDAnalysis as mda
    _MDA_AVAILABLE = True
except ImportError:
    _MDA_AVAILABLE = False

import concurrent.futures

print("Matplotlib backend:", matplotlib.get_backend())

# ---------------------------------------------------------------------------
# Internal logger (replaces lib.utils.Logger)
# ---------------------------------------------------------------------------
logging.basicConfig(level=logging.INFO,
                    format="%(levelname)s::%(message)s",
                    stream=sys.stdout)
_log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Helper: natural sort for corr_dcc_1.txt < corr_dcc_2.txt < corr_dcc_10.txt
# ---------------------------------------------------------------------------
def _natural_sort_key(filename):
    return [int(c) if c.isdigit() else c.lower()
            for c in re.split(r'(\d+)', filename)]


# ---------------------------------------------------------------------------
# Modern button helper 
# ---------------------------------------------------------------------------
def create_modern_button(parent, text, bg_color, command):
    btn = tk.Button(
        parent, text=text, bg=bg_color, fg="white",
        font=("Helvetica", 11, "bold"), relief="flat",
        cursor="hand2", pady=2, command=command
    )
    return btn


# ===========================================================================
# DCC computation engine  (dcc_no_mdtraj logic, embedded directly)
# ===========================================================================

_DCC_CMAP = mcolors.LinearSegmentedColormap.from_list(
    'dcc_cmap',
    [('white')] + [(cm.jet(i)) for i in range(40,250)]
   # [(0.00, 'darkblue'), (0.25, 'white'), (0.75, 'white'), (1.00, 'darkred')]
)


def _correlate_vectorized(coords_array: np.ndarray) -> np.ndarray:
    """
    DCC for one time window. Input: (T, N, 3) float32. Returns (N, N) float32.
    NaN-free: zero-variance atom pairs are set to 0.
    """
    T, N, _ = coords_array.shape
    chunk = 256
    delta = coords_array - coords_array.mean(axis=0, keepdims=True)

    cov = np.zeros((N, N), dtype=np.float64)
    for start in range(0, N, chunk):
        end = min(start + chunk, N)
        cov[:, start:end] = (
            np.einsum('tix,tjx->ij', delta, delta[:, start:end, :],
                      optimize=True) / T
        )

    var      = np.diag(cov)
    mag      = np.sqrt(np.maximum(var, 0.0))
    mag_prod = np.outer(mag, mag)
    safe     = mag_prod > 0.0
    corr     = np.where(safe, cov / np.where(safe, mag_prod, 1.0), 0.0)
    return corr.astype(np.float32)


def _process_single_window(coords_array: np.ndarray, window_end: int,
                            out_dir: str) -> int:
    """Worker: compute DCC, save PNG + TXT, return window_end."""
    corr   = _correlate_vectorized(coords_array)
    prefix = os.path.join(out_dir, f"corr_dcc_{window_end}")

    # Save text matrix (consumed by run_analysis via np.loadtxt)
    np.savetxt(f"{prefix}.txt", corr, fmt="%.6f")

    # Save heatmap PNG
    fig, ax = plt.subplots()
    hm = ax.pcolor(corr, cmap=_DCC_CMAP, vmin=-1, vmax=1)
    ax.set_frame_on(False)
    ax.grid(False)
    plt.xticks(rotation=90, fontsize=20)
    for axis in (ax.xaxis, ax.yaxis):
        for t in axis.get_major_ticks():
            t.tick1line.set_visible(False)
            t.tick2line.set_visible(False)
    plt.yticks(fontsize=20)
    plt.xlabel('Residue Index', fontsize=30)
    plt.ylabel('Residue Index', fontsize=30)
    cbar = plt.colorbar(hm, orientation='vertical')
    cbar.ax.tick_params(labelsize=20)
    cbar.set_label('DCC', fontsize=20)
    plt.savefig(f'{prefix}.png', dpi=150, format='png', bbox_inches='tight')
    plt.close('all')

    return window_end


def run_dcc_pipeline(trajectory: str, topology: str, step: int,
                     window_size: int, out_dir: str, n_cores: int, 
                     log_callback, done_callback):
    """
    Full DCC pipeline run in a background thread.
    Streams frames via MDAnalysis, dispatches windows to a ProcessPoolExecutor.
    Calls log_callback(msg: str) for progress and done_callback(out_dir) on finish.
    """
    if not _MDA_AVAILABLE:
        log_callback("ERROR: MDAnalysis not installed.\n"
                     "Run:  pip install MDAnalysis")
        done_callback(None)
        return

    os.makedirs(out_dir, exist_ok=True)

    try:
        u  = mda.Universe(topology, trajectory)
        ca = u.select_atoms("name CA")
        N  = len(ca)
        total_frames = len(u.trajectory[::step])
    except Exception as e:
        log_callback(f"ERROR loading files:\n{e}")
        done_callback(None)
        return

    if N == 0:
        log_callback("ERROR: No CA atoms found. Check topology file.")
        done_callback(None)
        return

    log_callback(f"Found {N} CA atoms  |  ~{total_frames} frames (step={step})")

    window_buffer     = np.zeros((window_size, N, 3), dtype=np.float32)
    frames_in_window  = 0
    current_window_end = window_size
    futures: set      = set()

    try:
        with concurrent.futures.ProcessPoolExecutor(max_workers=n_cores) as executor:
            for ts in u.trajectory[::step]:
                window_buffer[frames_in_window] = ca.positions.astype(np.float32)
                frames_in_window += 1

                if frames_in_window == window_size:
                    start = current_window_end - window_size
                    log_callback(f"Dispatching window {start} --- {current_window_end} ...")
                    coords_copy = window_buffer.copy()
                    future = executor.submit(
                        _process_single_window, coords_copy,
                        current_window_end, out_dir
                    )
                    futures.add(future)
                    frames_in_window   = 0
                    current_window_end += window_size

                    while len(futures) >= n_cores:
                        done, futures = concurrent.futures.wait(
                            futures,
                            return_when=concurrent.futures.FIRST_COMPLETED
                        )
                        for f in done:
                            try:
                                we = f.result()
                                log_callback(f" Window ending at frame {we} done.")
                            except Exception as e:
                                log_callback(f" Worker error: {e}")

            # Final partial window
            if frames_in_window >= 2:
                final_end = (current_window_end - window_size) + frames_in_window
                log_callback(f"Dispatching final partial window "
                             f"({frames_in_window} frames → {final_end}) ...")
                coords_copy = window_buffer[:frames_in_window].copy()
                future = executor.submit(
                    _process_single_window, coords_copy,
                    final_end, out_dir
                )
                futures.add(future)
            elif frames_in_window == 1:
                log_callback("Skipping 1-frame remainder (DCC undefined).")

            for future in concurrent.futures.as_completed(futures):
                try:
                    we = future.result()
                    log_callback(f" Window ending at frame {we} done.")
                except Exception as e:
                    log_callback(f" Worker error: {e}")

    except Exception as e:
        log_callback(f"PIPELINE ERROR: {e}")
        done_callback(None)
        return

    log_callback(f"\nAll windows complete.  Output → {out_dir}")
    done_callback(out_dir)


# ===========================================================================
# Intro Page  (unchanged)
# ===========================================================================

class IntroPage:
    def __init__(self, root):
        self.root = root
        self.root.title("TD-DCC")
        self.root.geometry("800x600")
        self.root.configure(bg="#93ba91")

        frame = tk.Frame(root, bg="#93ba91")
        frame.place(relx=0.5, rely=0.5, anchor="center")

        tk.Label(
            frame, text="Dynamic Cross-Correlation Analysis", 
            font=("Helvetica", 28, "bold"), fg="#87511f", bg="#93ba91" 
            ).pack(pady=20)

        tk.Label(
            frame, text="Explore correlated motions in biomolecular systems",
            font=("Helvetica", 18, "bold"), fg="#1f6387", bg="#93ba91"
        ).pack(pady=10)

        tk.Button(
            frame, text="Launch Analysis Dashboard",
            font=("Helvetica", 14, "bold"), bg="#1f6fd6", fg="white",
            relief="flat", cursor="hand2", padx=40, pady=12,
            command=self.launch
        ).pack(pady=40)

    def launch(self):
        for w in self.root.winfo_children():
            w.destroy()
        CorrelationGUI(self.root)


# ===========================================================================
# Main GUI
# ===========================================================================

class CorrelationGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("TDDCC")
        self.root.geometry("1400x900")
        self.root.configure(bg="#f4f6f9")

        # ── Original variables (unchanged) ──────────────────────────────────
        self.input_dir = tk.StringVar()
        self.i_min = tk.IntVar(value=1)
        self.i_max = tk.IntVar(value=60)
        self.j_min = tk.IntVar(value=1)
        self.j_max = tk.IntVar(value=60)
        self.res_dist = tk.IntVar(value=4)

        self.c_pos_min = tk.DoubleVar(value=0.2)
        self.c_pos_max = tk.DoubleVar(value=1.0)
        self.c_neg_min = tk.DoubleVar(value=-1.0)
        self.c_neg_max = tk.DoubleVar(value=-0.2)
        self.p_pos_min = tk.DoubleVar(value=0.1)
        self.p_pos_max = tk.DoubleVar(value=1.0)
        self.p_neg_min = tk.DoubleVar(value=0.1)
        self.p_neg_max = tk.DoubleVar(value=1.0)


        self.df_all   = None
        self.raw_data = None
        self.toolbar  = None

        # ── New DCC-compute variables ────────────────────────────────────────
        self.traj_path    = tk.StringVar()
        self.topo_path    = tk.StringVar()
        self.dcc_out_dir  = tk.StringVar()
        self.dcc_window   = tk.IntVar(value=100)
        self.dcc_step     = tk.IntVar(value=1)
        self.dcc_cores    = tk.IntVar(value=max(1, (os.cpu_count() or 2) - 1))

        self.build_gui()

    # -----------------------------------------------------------------------
    # GUI construction
    # -----------------------------------------------------------------------
    def build_gui(self):
        # ── Left panel (scrollable) ─────────────────────────────────────────
        control_container = tk.Frame(self.root, width=400, bg="#ffffff",
                                     highlightthickness=1,
                                     highlightbackground="#d1d5db")
        control_container.pack(side=tk.LEFT, fill=tk.Y)
        control_container.pack_propagate(False)

        canvas   = tk.Canvas(control_container, bg="#ffffff", highlightthickness=0)
        scrollbar = ttk.Scrollbar(control_container, orient="vertical",
                                  command=canvas.yview)
        
        self.left_panel = tk.Frame(canvas, bg="#ffffff", padx=15, pady=2)

        self.left_panel.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        canvas.create_window((0, 0), window=self.left_panel, anchor="nw", width=380)
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        section_font = ("Helvetica", 12, "bold")
        label_font   = ("Helvetica", 10)
        entry_font   = ("Helvetica", 10)

        # ── SECTION 0: Compute DCC from Trajectory (NEW) ────────────────────
        f_dcc = tk.LabelFrame(
            self.left_panel,
            text="Compute DCC from Trajectory ",
            font=section_font, bg="#ffffff", pady=4, padx=10
        )
        f_dcc.pack(fill="x", pady=2)

        # Trajectory file
        tk.Label(f_dcc, text="Trajectory (.xtc / .dcd / .trr):",
                 font=label_font, bg="#ffffff", anchor="w").pack(fill="x")
        tf = tk.Frame(f_dcc, bg="#ffffff")
        tf.pack(fill="x", pady=(0, 2))
        tk.Entry(tf, textvariable=self.traj_path, font=entry_font).pack(
            side="left", fill="x", expand=True, padx=(0, 5))
        tk.Button(tf, text="Browse", font=label_font,
                  command=self._browse_trajectory).pack(side="left")

        # Topology file
        tk.Label(f_dcc, text="Topology (.pdb / .tpr / .gro):",
                 font=label_font, bg="#ffffff", anchor="w").pack(fill="x")
        topof = tk.Frame(f_dcc, bg="#ffffff")
        topof.pack(fill="x", pady=(0, 2))
        tk.Entry(topof, textvariable=self.topo_path, font=entry_font).pack(
            side="left", fill="x", expand=True, padx=(0, 5))
        tk.Button(topof, text="Browse", font=label_font,
                  command=self._browse_topology).pack(side="left")

        # Output directory
        tk.Label(f_dcc, text="Output directory:",
                 font=label_font, bg="#ffffff", anchor="w").pack(fill="x")
        outf = tk.Frame(f_dcc, bg="#ffffff")
        outf.pack(fill="x", pady=(0, 2))
        tk.Entry(outf, textvariable=self.dcc_out_dir, font=entry_font).pack(
            side="left", fill="x", expand=True, padx=(0, 5))
        tk.Button(outf, text="Browse", font=label_font,
                  command=self._browse_out_dir).pack(side="left")

        # Parameters row
        param_f = tk.Frame(f_dcc, bg="#ffffff")
        param_f.pack(fill="x", pady=1)

        for col, (lbl, var, w) in enumerate([
            ("Window",  self.dcc_window, 7),
            ("Step",    self.dcc_step,   4),
            ("Cores",   self.dcc_cores,  4),
        ]):
            tk.Label(param_f, text=lbl, font=label_font,
                     bg="#ffffff").grid(row=0, column=col*2, sticky="w", padx=(4,1))
            tk.Entry(param_f, textvariable=var, width=w,
                     font=entry_font).grid(row=0, column=col*2+1, padx=(0,6))

        # Log box
        self.dcc_log = tk.Text(f_dcc, height=2, font=("Courier", 9),
                               bg="#f0f4f8", relief="flat", state="disabled",
                               wrap="word")
        self.dcc_log.pack(fill="x", pady=(2, 1))

        # Compute button + progress bar
        self.dcc_progress = ttk.Progressbar(f_dcc, mode="indeterminate")
        self.dcc_progress.pack(fill="x", pady=(0, 2))

        self.dcc_btn = create_modern_button(
            f_dcc, "COMPUTE DCC", "#8e44ad", self._start_dcc_compute
        )
        self.dcc_btn.pack(fill="x", pady=2)

        # Auto-load hint label
        self.dcc_hint = tk.Label(
            f_dcc, text="", font=("Helvetica", 9, "italic"),
            fg="#27ae60", bg="#ffffff", wraplength=340, justify="left"
        )
        self.dcc_hint.pack(fill="x")

        # ── SECTION 1: Input Data ──────────────────────
        f_input = tk.LabelFrame(self.left_panel, text="Input Data ",
                                font=section_font, bg="#ffffff", pady=4, padx=10)
        f_input.pack(fill="x", pady=2)
        entry_frame = tk.Frame(f_input, bg="#ffffff")
        entry_frame.pack(fill="x")
        tk.Entry(entry_frame, textvariable=self.input_dir,
                 font=entry_font).pack(side="left", fill="x", expand=True, padx=(0, 5))
        tk.Button(entry_frame, text="Browse", font=label_font,
                  command=self.select_input).pack(side="left")

        # ── SECTION 2: Residue Selection ───────────────
        f_res = tk.LabelFrame(self.left_panel, text="Residue Limits ",
                              font=section_font, bg="#ffffff", pady=4, padx=10)
        f_res.pack(fill="x", pady=2)
        self.range_entry(f_res, "Residue i range:", self.i_min, self.i_max,
                         label_font, entry_font)
        self.range_entry(f_res, "Residue j range:", self.j_min, self.j_max,
                         label_font, entry_font)
        self.single_entry(f_res, "Min distance |j - i| >=", self.res_dist,
                          label_font, entry_font)

        # ── SECTION 3: Correlation Filters ─────────────
        f_adv = tk.LabelFrame(self.left_panel, text="Correlation Filters ",
                              font=section_font, bg="#ffffff", pady=4, padx=10)
        f_adv.pack(fill="x", pady=2)
        self.range_entry(f_adv, "Positive correlation:", self.c_pos_min,
                         self.c_pos_max, label_font, entry_font)
        self.range_entry(f_adv, "Negative correlation:", self.c_neg_min,
                         self.c_neg_max, label_font, entry_font)
        self.range_entry(f_adv, "Positive corr prob:",   self.p_pos_min,
                         self.p_pos_max, label_font, entry_font)
        self.range_entry(f_adv, "Negative corr prob:",   self.p_neg_min,
                         self.p_neg_max, label_font, entry_font)

        # Status
        self.status_label = tk.Label(self.left_panel, text="",
                                     font=("Helvetica", 10, "italic"),
                                     fg="#e67e22", bg="#ffffff")
        self.status_label.pack(pady=2)

        # Actions
        create_modern_button(self.left_panel, "RUN ANALYSIS", "#1f6a87",
                             self.run_analysis).pack(fill="x", pady=1)

        f_dash = tk.Frame(self.left_panel, bg="#ffffff")
        f_dash.pack(fill="x", pady=2)
        create_modern_button(f_dash, "View Correlation Dashboard", "#61871f",
                             self.open_correlation_dashboard).pack(fill="x", pady=1)

        # Downloads
        f_csv = tk.LabelFrame(self.left_panel, text=" Export Data ",
                              font=section_font, bg="#ffffff", pady=4, padx=10)
        f_csv.pack(fill="x", pady=2)
        create_modern_button(f_csv, "Download Combined CSV", "#7f8047",
                             self.download_combined_csv).pack(fill="x", pady=1)
        create_modern_button(f_csv, "Download Positive CSV", "#301f87",
                             self.download_positive_csv).pack(fill="x", pady=1)
        create_modern_button(f_csv, "Download Negative CSV", "#c0392b",
                             self.download_negative_csv).pack(fill="x", pady=1)

        # ── Right plot area  (original, unchanged) ───────────────────────────
        self.plot_frame = tk.Frame(self.root, bg="#f4f6f9")
        self.plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True,
                             padx=10, pady=10)

        self.toolbar_frame = tk.Frame(self.plot_frame, bg="#f4f6f9")
        self.toolbar_frame.pack(side=tk.BOTTOM, fill=tk.X)

        self.fig_dash = plt.Figure(figsize=(10, 8), dpi=100)
        self.fig_dash.patch.set_facecolor('#f4f6f9')

        self.canvas_dash = FigureCanvasTkAgg(self.fig_dash, master=self.plot_frame)
        self.canvas_dash.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.canvas_dash.mpl_connect("pick_event", self._on_heatmap_pick)
        self.canvas_dash.draw()

    # -----------------------------------------------------------------------
    # New: DCC file browsers
    # -----------------------------------------------------------------------
    def _browse_trajectory(self):
        path = filedialog.askopenfilename(
            title="Select trajectory file",
            filetypes=[
                ("Trajectory files", "*.xtc *.dcd *.trr *.nc *.ncdf"),
                ("All files", "*.*")
            ]
        )
        if path:
            self.traj_path.set(path)

    def _browse_topology(self):
        path = filedialog.askopenfilename(
            title="Select topology file",
            filetypes=[
                ("Topology files", "*.pdb *.tpr *.gro *.psf *.prmtop"),
                ("All files", "*.*")
            ]
        )
        if path:
            self.topo_path.set(path)

    def _browse_out_dir(self):
        folder = filedialog.askdirectory(title="Select output directory for DCC files")
        if folder:
            self.dcc_out_dir.set(folder)

    # -----------------------------------------------------------------------
    # New: DCC log helper (thread-safe via root.after)
    # -----------------------------------------------------------------------
    def _dcc_log(self, msg: str):
        def _append():
            self.dcc_log.configure(state="normal")
            self.dcc_log.insert("end", msg + "\n")
            self.dcc_log.see("end")
            self.dcc_log.configure(state="disabled")
        self.root.after(0, _append)

    # -----------------------------------------------------------------------
    # New: Start DCC computation in background thread
    # -----------------------------------------------------------------------
    def _start_dcc_compute(self):
        if not _MDA_AVAILABLE:
            messagebox.showerror(
                "Missing dependency",
                "MDAnalysis is not installed.\n\nRun:\n  pip install MDAnalysis"
            )
            return

        traj    = self.traj_path.get().strip()
        topo    = self.topo_path.get().strip()
        out_dir = self.dcc_out_dir.get().strip()

        if not traj or not os.path.isfile(traj):
            messagebox.showerror("Error", "Please select a valid trajectory file.")
            return
        if not topo or not os.path.isfile(topo):
            messagebox.showerror("Error", "Please select a valid topology file.")
            return
        if not out_dir:
            messagebox.showerror("Error", "Please select an output directory.")
            return

        window  = self.dcc_window.get()
        step    = self.dcc_step.get()
        cores   = self.dcc_cores.get()

        if window < 2:
            messagebox.showerror("Error", "Window size must be ≥ 2.")
            return
        if step < 1:
            messagebox.showerror("Error", "Step must be ≥ 1.")
            return
        if cores < 1:
            messagebox.showerror("Error", "Cores must be ≥ 1.")
            return

        # Clear log and start spinner
        self.dcc_log.configure(state="normal")
        self.dcc_log.delete("1.0", "end")
        self.dcc_log.configure(state="disabled")
        self.dcc_hint.config(text="")
        self.dcc_btn.config(state="disabled")
        self.dcc_progress.start(12)
        self._dcc_log(f"Started at {datetime.now().strftime('%H:%M:%S')}")

        def _done(result_dir):
            def _ui():
                self.dcc_progress.stop()
                self.dcc_btn.config(state="normal")
                if result_dir:
                    # Auto-populate Section 1 input_dir so user can click Run Analysis
                    self.input_dir.set(result_dir)
                    self._auto_set_residue_limits(result_dir)
                    self.dcc_hint.config(
                        text=f" Done! Input Data (Section 1) auto-filled with:\n{result_dir}"
                    )
                    messagebox.showinfo(
                        "DCC Complete",
                        f"DCC matrices saved to:\n{result_dir}\n\n"
                        #"Section 1 'Input Data' has been filled automatically.\n"
                        "Click  RUN ANALYSIS  to visualise."
                    )
                else:
                    self.dcc_hint.config(text="")
                    messagebox.showerror("DCC Failed",
                                         "Computation failed. Check the log for details.")
            self.root.after(0, _ui)

        threading.Thread(
            target=run_dcc_pipeline,
            args=(traj, topo, step, window, out_dir, cores,
                  self._dcc_log, _done),
            daemon=True
        ).start()

    # -----------------------------------------------------------------------
    # UI Helpers 
    # -----------------------------------------------------------------------
    def range_entry(self, parent, text, v1, v2, lf, ef):
        f = tk.Frame(parent, bg="#ffffff")
        f.pack(anchor="w", fill="x", pady=1)
        tk.Label(f, text=text, font=lf, bg="#ffffff", anchor="w").grid(
            row=0, column=0, sticky="w")
        tk.Entry(f, textvariable=v1, width=6, font=ef).grid(row=0, column=1, padx=(5, 2))
        tk.Label(f, text="to", font=lf, bg="#ffffff").grid(row=0, column=2, padx=2)
        tk.Entry(f, textvariable=v2, width=6, font=ef).grid(row=0, column=3, padx=2)
        f.grid_columnconfigure(0, minsize=165)

    def single_entry(self, parent, text, var, lf, ef):
        f = tk.Frame(parent, bg="#ffffff")
        f.pack(anchor="w", fill="x", pady=1)
        tk.Label(f, text=text, font=lf, bg="#ffffff", anchor="w").grid(
            row=0, column=0, sticky="w")
        tk.Entry(f, textvariable=var, width=6, font=ef).grid(row=0, column=1, padx=(5, 2))
        f.grid_columnconfigure(0, minsize=165)

    def select_input(self):
        folder = filedialog.askdirectory(
            title="Select folder containing corr_dcc_*.txt files")
        if folder:
            self.input_dir.set(folder)
            self._auto_set_residue_limits(folder)
            # NEW: Tell the Matplotlib toolbars to use this new folder for saving PNGs
            matplotlib.rcParams['savefig.directory'] = folder

    def _auto_set_residue_limits(self, folder):
        """Peeks at the first txt file to detect the number of CA atoms."""
        try:
            files = [f for f in os.listdir(folder) if f.endswith(".txt")]
            if not files:
                return

            # Read just the first line of the first file (very fast, no full memory load)
            first_file = os.path.join(folder, files[0])
            with open(first_file, 'r') as f:
                first_line = f.readline()
                if first_line:
                    n_residues = len(first_line.split())

                    # Update the GUI boxes automatically
                    self.i_max.set(n_residues)
                    self.j_max.set(n_residues)

                    # Optional: Give the user visual feedback
                    self.status_label.config(
                        text=f"Auto-detected {n_residues} residues from files.",
                        fg="#2980b9"
                    )
        except Exception:
            pass # If anything fails, just keep the default numbers

    # -----------------------------------------------------------------------
    # Core Logic  (original, unchanged)
    # -----------------------------------------------------------------------
    def run_analysis(self):
        folder = self.input_dir.get()
        if not os.path.isdir(folder):
            messagebox.showerror("Error", "Invalid input directory")
            return

        raw_files = [f for f in os.listdir(folder) if f.endswith(".txt")]
        files     = sorted(raw_files, key=_natural_sort_key)

        if not files:
            messagebox.showerror("Error", "No .txt files found in selected folder")
            return

        # -------------------------------------------------------------------
        # NEW SAFETY CHECK: Verify limits against actual data size and minimums
        # -------------------------------------------------------------------
        i_min_val, i_max_val = self.i_min.get(), self.i_max.get()
        j_min_val, j_max_val = self.j_min.get(), self.j_max.get()

        # 1. Enforce the floor (Must start at 1)
        if i_min_val < 1 or j_min_val < 1:
            messagebox.showerror("Invalid Input", "Residue limits must start at 1. Zero or negative numbers are not allowed.")
            return

        # 2. Prevent backward ranges
        if i_min_val > i_max_val or j_min_val > j_max_val:
            messagebox.showerror("Invalid Input", "Minimum residue cannot be greater than the maximum residue.")
            return

        # 3. Enforce the ceiling (Must not exceed file limits)
        try:
            first_file_path = os.path.join(folder, files[0])
            with open(first_file_path, 'r') as f:
                actual_n_residues = len(f.readline().split())

            if i_max_val > actual_n_residues or j_max_val > actual_n_residues:
                messagebox.showerror(
                    "Limit Exceeded",
                    f"Your requested residue range exceeds the actual data.\n\n"
                    f"The loaded files only contain {actual_n_residues} residues. "
                    f"Please lower your i/j maximums."
                )
                return
        except Exception:
            pass # If reading fails for some reason, the main block will catch the error later
        # -------------------------------------------------------------------

        self.status_label.config(text="Processing data, please wait...", fg="#e67e22")
        self.root.update_idletasks()

        try:
            data = np.array([np.loadtxt(os.path.join(folder, f)) for f in files])

            rows = []
            for i in range(self.i_min.get(), self.i_max.get() + 1):
                j_start = max(self.j_min.get(), i + self.res_dist.get())
                for j in range(j_start, self.j_max.get() + 1):
                    if i - 1 < data.shape[1] and j - 1 < data.shape[2]:
                        vals = data[:, i - 1, j - 1]
                        rows.append([
                            i, j, np.mean(vals), np.std(vals),
                            np.mean(vals > 0), np.mean(vals < 0)
                        ])

            if not rows:
                messagebox.showerror("Error",
                                     "No valid matrix data within selected range.")
                self.status_label.config(text="")
                return

            rows.sort(key=lambda x: (x[0], x[1]))
            self.df_all = pd.DataFrame(rows, columns=[
                "Residue_i", "Residue_j", "Mean_Correlation",
                "Std_Correlation", "Positive_Probability", "Negative_Probability"
            ])
            self.raw_data   = data
            self.frame_files = files

            self.status_label.config(text="Analysis complete!", fg="#27ae60")
            messagebox.showinfo("Done", "Data processing complete.")

        except Exception as e:
            messagebox.showerror("Processing Error", f"An error occurred: {str(e)}")
            self.status_label.config(text="")

    def compute_masks(self, corr, pos_prob, neg_prob):
        mask_pos = ((corr >= self.c_pos_min.get()) & (corr <= self.c_pos_max.get()) &
                    (pos_prob >= self.p_pos_min.get()) & (pos_prob <= self.p_pos_max.get()))
        mask_neg = ((corr >= self.c_neg_min.get()) & (corr <= self.c_neg_max.get()) &
                    (neg_prob >= self.p_neg_min.get()) & (neg_prob <= self.p_neg_max.get()))
        return mask_pos, mask_neg, "Range-based Cutoff"

    def _get_pivot_data(self):
        i_min, i_max = self.i_min.get(), self.i_max.get()
        j_min, j_max = self.j_min.get(), self.j_max.get()
        i_range = list(range(i_min, i_max + 1))
        j_range = list(range(j_min, j_max + 1))

        corr = (self.df_all
                .pivot(index="Residue_i", columns="Residue_j", values="Mean_Correlation")
                .reindex(index=i_range, columns=j_range)
                .clip(-1, 1))
        pos_prob = (self.df_all
                    .pivot(index="Residue_i", columns="Residue_j",
                           values="Positive_Probability")
                    .reindex(index=i_range, columns=j_range)
                    .clip(0, 1))
        neg_prob = (self.df_all
                    .pivot(index="Residue_i", columns="Residue_j",
                           values="Negative_Probability")
                    .reindex(index=i_range, columns=j_range)
                    .clip(0, 1))

        grid_extent = [j_min - 0.5, j_max + 0.5, i_min - 0.5, i_max + 0.5]
        return corr, pos_prob, neg_prob, grid_extent

    def _decorate_axes(self, ax):
        ax.set_xlabel("Residue j", fontsize=11, fontweight="bold")
        ax.set_ylabel("Residue i", fontsize=11, fontweight="bold")
        ax.tick_params(labelsize=10)
        # NEW LINE: Hides x, y, and z values on hover for the heatmaps
        ax.format_coord = lambda x, y: ""

    def _reset_toolbar(self):
        if self.toolbar:
            self.toolbar.destroy()
        self.toolbar = NavigationToolbar2Tk(self.canvas_dash, self.toolbar_frame)
        # NEW LINE: Absolutely forces the toolbar to never display hover messages
        self.toolbar.set_message = lambda msg: None
        self.toolbar.update()

    # -----------------------------------------------------------------------
    # Dashboards  (original, unchanged)
    # -----------------------------------------------------------------------
    def open_correlation_dashboard(self):
        if self.df_all is None or self.df_all.empty:
            messagebox.showerror("Error", "Run analysis first.")
            return

        self.fig_dash.clear()
        self.fig_dash.patch.set_facecolor('#f4f6f9')

        ax1 = self.fig_dash.add_subplot(3, 1, 1)
        ax2 = self.fig_dash.add_subplot(3, 1, 2)
        ax3 = self.fig_dash.add_subplot(3, 1, 3)

        corr, pos_prob, neg_prob, grid_extent = self._get_pivot_data()
        mask_pos, mask_neg, _ = self.compute_masks(corr, pos_prob, neg_prob)

        cmap_bwr1  = mcolors.LinearSegmentedColormap.from_list(
            "bwr", ["white", "blue"])
        im1 = ax1.imshow(corr.where(mask_pos), cmap=cmap_bwr1, vmin=0, vmax=1,
                         origin="lower", extent=grid_extent, aspect="auto", picker=True)
        self.fig_dash.colorbar(im1, ax=ax1, fraction=0.035,
                               pad=0.02).set_label("Positive Corr", fontsize=10)
        ax1.set_title("Positive Correlation", fontweight="bold")

        cmap_bwr2  = mcolors.LinearSegmentedColormap.from_list(
            "bwr", ["red", "white"])
        im2 = ax2.imshow(corr.where(mask_neg), cmap=cmap_bwr2, vmin=-1, vmax=0,
                         origin="lower", extent=grid_extent, aspect="auto", picker=True)
        self.fig_dash.colorbar(im2, ax=ax2, fraction=0.035,
                               pad=0.02).set_label("Negative Corr", fontsize=10)
        ax2.set_title("Negative Correlation", fontweight="bold")

        combined  = corr.where(mask_pos | mask_neg)
        cmap_bwr3  = mcolors.LinearSegmentedColormap.from_list(
            "bwr", ["red", "white", "blue"])
        im3 = ax3.imshow(combined.values, cmap=cmap_bwr3, vmin=-1, vmax=1,
                         origin="lower", extent=grid_extent, aspect="auto", picker=True)
        self.fig_dash.colorbar(im3, ax=ax3, fraction=0.035,
                               pad=0.02).set_label("Combined Corr", fontsize=10)
        ax3.set_title("Combined Correlation", fontweight="bold")

        for ax in (ax1, ax2, ax3):
            self._decorate_axes(ax)

        self.fig_dash.tight_layout()
        self.canvas_dash.draw()
        self._reset_toolbar()

    # -----------------------------------------------------------------------
    # Mouse Events & Time Series  (original, unchanged)
    # -----------------------------------------------------------------------
    def _on_heatmap_pick(self, event):
        mouse = event.mouseevent
        if not mouse.dblclick or mouse.xdata is None or mouse.ydata is None:
            return
        i, j = int(round(mouse.ydata)), int(round(mouse.xdata))
        # Add these two lines to force the index to never drop below 1
        i = max(1, i)
        j = max(1, j)
        if self.raw_data is None:
            return
        if not (self.i_min.get() <= i <= self.i_max.get()):
            return
        if not (self.j_min.get() <= j <= self.j_max.get()):
            return
        try:
            ts = self.raw_data[:, i - 1, j - 1]
            self.open_time_series_window(i, j, ts)
        except IndexError:
            pass

    def open_time_series_window(self, i, j, ts):
        win = tk.Toplevel(self.root)
        win.title(f"Time Series — Residue ({i}, {j})")
        win.geometry("750x520")
        win.configure(bg="#f4f6f9")
        
        # NEW: Calculate actual frame numbers
        window = self.dcc_window.get()
        step = self.dcc_step.get()
        x_values = [k * window * step for k in range(1, len(ts) + 1)]

        fig = plt.Figure(figsize=(8, 4), dpi=110)
        ax  = fig.add_subplot(111)
        # MODIFIED: Pass x_values to the plot
        ax.plot(x_values, ts, lw=2, color="#2980b9")
        #ax.plot(ts, lw=2, color="#2980b9")
        ax.set_title(f"Frame-wise Correlation  i={i}, j={j}",
                     fontsize=13, fontweight="bold")
        ax.set_xlabel("Frame", fontsize=11, fontweight="bold")
        ax.set_ylabel("Correlation", fontsize=11, fontweight="bold")
        ax.grid(True, linestyle="--", linewidth=0.8, alpha=0.6)
        
        # NEW LINE: Hides x, y values on hover for the time series line graph
        ax.format_coord = lambda x, y: ""

        canvas = FigureCanvasTkAgg(fig, master=win)
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        toolbar = NavigationToolbar2Tk(canvas, win)
        toolbar.update()
        canvas.draw()

        btn_frame = tk.Frame(win, bg="#f4f6f9")
        btn_frame.pack(fill=tk.X, padx=10, pady=6)

        def save_csv():
            path = filedialog.asksaveasfilename(
                parent=win, defaultextension=".csv",
                initialfile=f"timeseries_i{i}_j{j}.csv",
                filetypes=[("CSV file", "*.csv"), ("All files", "*.*")],
                title="Save data as CSV"
            )
            if not path:
                return
            try:
                frame_labels = (self.frame_files
                                if hasattr(self, "frame_files") and self.frame_files
                                else [str(k) for k in range(len(ts))])
                df_ts = pd.DataFrame({
                    "Frame":       x_values, #range(1, len(ts) + 1),
                 #   "File":        frame_labels,
                    "Correlation": ts
                })
                df_ts.to_csv(path, index=False)
                messagebox.showinfo("Saved", f"CSV saved:\n{path}", parent=win)
            except Exception as e:
                messagebox.showerror("Error", f"Could not save CSV:\n{e}", parent=win)

        tk.Button(btn_frame, text="Save CSV", bg="#27ae60", fg="white",
                  font=("Helvetica", 10, "bold"), relief="flat", cursor="hand2",
                  padx=16, pady=5, command=save_csv).pack(side=tk.LEFT, padx=5)
        tk.Label(btn_frame,
                 text=f"i={i}     j={j}    "# |  frames={len(ts)}  |  "
                      f"mean={np.mean(ts):.4f}    std={np.std(ts):.4f}",
                 font=("Helvetica", 9), fg="#555555", bg="#f4f6f9").pack(
                 side=tk.LEFT, padx=15)

    # -----------------------------------------------------------------------
    # CSV Exports  (original, unchanged)
    # -----------------------------------------------------------------------
    def download_positive_csv(self):
        if self.df_all is None or self.df_all.empty:
            return messagebox.showerror("Error", "Run analysis first.")
        c_pos_min, c_pos_max = self.c_pos_min.get(), self.c_pos_max.get()
        p_pos_min, p_pos_max = self.p_pos_min.get(), self.p_pos_max.get()
        df = self.df_all[
            (self.df_all["Mean_Correlation"] >= c_pos_min) &
            (self.df_all["Mean_Correlation"] <= c_pos_max) &
            (self.df_all["Positive_Probability"] >= p_pos_min) &
            (self.df_all["Positive_Probability"] <= p_pos_max)
        ]
        self._save_csv(df[["Residue_i", "Residue_j",
                            "Mean_Correlation", "Positive_Probability"]],
                       "positive_correlation.csv")

    def download_negative_csv(self):
        if self.df_all is None or self.df_all.empty:
            return messagebox.showerror("Error", "Run analysis first.")
        c_neg_min, c_neg_max = self.c_neg_min.get(), self.c_neg_max.get()
        p_neg_min, p_neg_max = self.p_neg_min.get(), self.p_neg_max.get()
        df = self.df_all[
            (self.df_all["Mean_Correlation"] >= c_neg_min) &
            (self.df_all["Mean_Correlation"] <= c_neg_max) &
            (self.df_all["Negative_Probability"] >= p_neg_min) &
            (self.df_all["Negative_Probability"] <= p_neg_max)
        ]
        self._save_csv(df[["Residue_i", "Residue_j",
                            "Mean_Correlation", "Negative_Probability"]],
                       "negative_correlation.csv")

    def download_combined_csv(self):
        if self.df_all is None or self.df_all.empty:
            return messagebox.showerror("Error", "Run analysis first.")
        c_pos_min, c_pos_max = self.c_pos_min.get(), self.c_pos_max.get()
        p_pos_min, p_pos_max = self.p_pos_min.get(), self.p_pos_max.get()
        c_neg_min, c_neg_max = self.c_neg_min.get(), self.c_neg_max.get()
        p_neg_min, p_neg_max = self.p_neg_min.get(), self.p_neg_max.get()
        cond_pos = (
            (self.df_all["Mean_Correlation"] >= c_pos_min) &
            (self.df_all["Mean_Correlation"] <= c_pos_max) &
            (self.df_all["Positive_Probability"] >= p_pos_min) &
            (self.df_all["Positive_Probability"] <= p_pos_max)
        )
        cond_neg = (
            (self.df_all["Mean_Correlation"] >= c_neg_min) &
            (self.df_all["Mean_Correlation"] <= c_neg_max) &
            (self.df_all["Negative_Probability"] >= p_neg_min) &
            (self.df_all["Negative_Probability"] <= p_neg_max)
        )
        df = self.df_all[cond_pos | cond_neg]
        self._save_csv(df[["Residue_i", "Residue_j", "Mean_Correlation",
                            "Positive_Probability", "Negative_Probability"]],
                       "combined_correlation.csv")

    def _save_csv(self, df, default_name):
        if df.empty:
            messagebox.showwarning("No Data",
                                   "No entries satisfy the selected thresholds.")
            return
        file_path = filedialog.asksaveasfilename(
            defaultextension=".csv", initialfile=default_name,
            filetypes=[("CSV files", "*.csv")], title="Save CSV"
        )
        if not file_path:
            return
        try:
            df.to_csv(file_path, index=False)
            messagebox.showinfo("Success", f"CSV saved:\n{file_path}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save CSV:\n{e}")


# ===========================================================================
if __name__ == "__main__":
    root = tk.Tk()
    IntroPage(root)
    root.mainloop()
