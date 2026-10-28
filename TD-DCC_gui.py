from pathlib import Path
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
class BlueProgressBar(tk.Frame):
    """OS-independent indeterminate progress bar.

    Avoids ttk theme engines changing the trough/indicator to black on
    macOS/Windows/Linux. Only standard Tk Frame widgets are used.
    """
    def __init__(self, parent, height=12, **kwargs):
        super().__init__(
            parent, height=height, bg="#e2e6eb",
            highlightbackground="#aeb7c2",
            highlightcolor="#1f6fd6",
            highlightthickness=1, bd=0, relief="flat", **kwargs
        )
        self._height = height
        self._running = False
        self._after_id = None
        self._position = 0
        self._segment_width = 0.28
        self._canvas = tk.Canvas(
            self, height=height, bg="#e2e6eb",
            highlightthickness=0, bd=0, relief="flat"
        )
        self._canvas.pack(fill="both", expand=True)
        self.bind("<Configure>", self._redraw)

    def _redraw(self, event=None):
        self._canvas.delete("all")
        w = max(1, self.winfo_width())
        h = max(1, self.winfo_height())
        self._canvas.create_rectangle(
            0, 0, w, h, fill="#e2e6eb", outline=""
        )
        if self._running:
            seg = max(20, int(w * self._segment_width))
            x = int((w + seg) * self._position / 100.0) - seg
            self._canvas.create_rectangle(
                x, 1, x + seg, h - 1,
                fill="#1f6fd6", outline=""
            )

    def _animate(self):
        if not self._running:
            return
        self._position = (self._position + 3) % 110
        self._redraw()
        self._after_id = self.after(25, self._animate)

    def start(self, interval=12):
        self.stop()
        self._running = True
        self._position = 0
        self._animate()

    def stop(self):
        self._running = False
        if self._after_id is not None:
            try:
                self.after_cancel(self._after_id)
            except tk.TclError:
                pass
            self._after_id = None
        self._redraw()

class ColoredButton(tk.Frame):
    """Theme-independent coloured button for macOS Light/Dark mode.

    A Frame + Label is used instead of native tk.Button so macOS Aqua
    cannot replace the requested background colour in Dark Mode.
    """

    def __init__(self, parent, text, bg_color, command,
                 font=("Helvetica", 11, "bold"), padx=10, pady=7):
        super().__init__(
            parent, bg=bg_color, highlightbackground=bg_color,
            highlightcolor=bg_color, highlightthickness=0,
            bd=0, relief="flat", cursor="hand2"
        )
        self._command = command
        self._enabled = True

        self.label = tk.Label(
            self, text=text, bg=bg_color, fg="#ffffff",
            activebackground=bg_color, activeforeground="#ffffff",
            font=font, padx=padx, pady=pady, bd=0, relief="flat",
            cursor="hand2", anchor="center"
        )
        self.label.pack(fill="both", expand=True)

        self.bind("<Button-1>", self._on_click)
        self.label.bind("<Button-1>", self._on_click)

    def _on_click(self, event=None):
        if self._enabled and callable(self._command):
            self._command()

    def config(self, cnf=None, **kwargs):
        """Keep compatibility with existing .config(state=...) calls."""
        if cnf:
            kwargs.update(cnf)
        state = kwargs.pop("state", None)
        if state is not None:
            self._enabled = (state != "disabled")
            cursor = "hand2" if self._enabled else "arrow"
            self.configure(cursor=cursor)
            self.label.configure(cursor=cursor, fg="#ffffff")
        if kwargs:
            return super().config(**kwargs)
        return None

    configure = config


def create_modern_button(parent, text, bg_color, command,
                         font=("Helvetica", 11, "bold"), padx=10, pady=7):
    return ColoredButton(
        parent, text, bg_color, command,
        font=font, padx=padx, pady=pady
    )


def create_readable_entry(parent, **kwargs):
    """Explicit entry colours so text remains readable in macOS Dark Mode."""
    options = {
        "bg": "#ffffff",
        "fg": "#1a1a1a",
        "insertbackground": "#1a1a1a",
        "selectbackground": "#1f6fd6",
        "selectforeground": "#ffffff",
        "disabledbackground": "#f2f2f2",
        "disabledforeground": "#555555",
        "relief": "solid",
        "bd": 1,
        "highlightthickness": 1,
        "highlightbackground": "#b8b8b8",
        "highlightcolor": "#1f6fd6",
    }
    options.update(kwargs)
    return tk.Entry(parent, **options)


def create_browse_button(parent, command):
    """Blue Browse button that stays blue in macOS Light/Dark mode."""
    return ColoredButton(
        parent, "Browse", "#1f6fd6", command,
        font=("Helvetica", 10, "bold"), padx=10, pady=3
    )


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
                     selection_str: str, log_callback, done_callback):
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
        # Handle multiple files or wildcard support
        if ";" in trajectory:
            traj_list = [t.strip() for t in trajectory.split(";") if t.strip()]
            u = mda.Universe(topology, *traj_list)
        else:
            u = mda.Universe(topology, trajectory)
            
        ca = u.select_atoms(selection_str)
        N  = len(ca)
        total_frames = len(u.trajectory[::step])
    except Exception as e:
        log_callback(f"ERROR loading files or selection:\n{e}")
        done_callback(None)
        return

    if N == 0:
        log_callback(f"ERROR: No atoms found for selection '{selection_str}'. Check topology file.")
        done_callback(None)
        return

    log_callback(f"Found {N} atoms  |  ~{total_frames} frames (step={step})")

    window_buffer     = np.zeros((window_size, N, 3), dtype=np.float32)
    frames_in_window  = 0
    futures: set      = set()

    try:
        with concurrent.futures.ProcessPoolExecutor(max_workers=n_cores) as executor:
            for ts in u.trajectory[::step]:
                window_buffer[frames_in_window] = ca.positions.astype(np.float32)
                frames_in_window += 1

                if frames_in_window == window_size:
                    # Capture exact frame index from MDAnalysis (1-indexed for clarity)
                    actual_frame = ts.frame + 1 
                    log_callback(f"Dispatching window ending at frame {actual_frame} ...")
                    coords_copy = window_buffer.copy()
                    future = executor.submit(
                        _process_single_window, coords_copy,
                        actual_frame, out_dir
                    )
                    futures.add(future)
                    frames_in_window   = 0

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
                actual_frame = ts.frame + 1
                log_callback(f"Dispatching final partial window "
                             f"({frames_in_window} frames → {actual_frame}) ...")
                coords_copy = window_buffer[:frames_in_window].copy()
                future = executor.submit(
                    _process_single_window, coords_copy,
                    actual_frame, out_dir
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

        create_modern_button(
            frame,
            "Launch Analysis Dashboard",
            "#1f6fd6",
            self.launch,
            font=("Helvetica", 14, "bold"),
            padx=40,
            pady=12,
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
        self.dcc_selection = tk.StringVar(value="name CA")
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
            font=section_font, bg="#ffffff", fg="black", pady=4, padx=10
        )
        f_dcc.pack(fill="x", pady=2)

        # Trajectory file
        tk.Label(f_dcc, text="Trajectory (.xtc/.dcd/.trr):",# [Use ';' for multiple]:",
                 font=label_font, bg="#ffffff", fg="black", anchor="w").pack(fill="x")
        tf = tk.Frame(f_dcc, bg="#ffffff")
        tf.pack(fill="x", pady=(0, 2))
        create_readable_entry(tf, textvariable=self.traj_path, font=entry_font).pack(
            side="left", fill="x", expand=True, padx=(0, 5))
        create_browse_button(tf, self._browse_trajectory).pack(side="left")

        # Topology file
        tk.Label(f_dcc, text="Topology (.pdb / .tpr / .gro):",
                 font=label_font, bg="#ffffff", fg="black", anchor="w").pack(fill="x")
        topof = tk.Frame(f_dcc, bg="#ffffff")
        topof.pack(fill="x", pady=(0, 2))
        create_readable_entry(topof, textvariable=self.topo_path, font=entry_font).pack(
            side="left", fill="x", expand=True, padx=(0, 5))
        create_browse_button(topof, self._browse_topology).pack(side="left")

        # Output directory
        tk.Label(f_dcc, text="Output directory:",
                 font=label_font, bg="#ffffff", fg="black", anchor="w").pack(fill="x")
        outf = tk.Frame(f_dcc, bg="#ffffff")
        outf.pack(fill="x", pady=(0, 2))
        create_readable_entry(outf, textvariable=self.dcc_out_dir, font=entry_font).pack(
            side="left", fill="x", expand=True, padx=(0, 5))
        create_browse_button(outf, self._browse_out_dir).pack(side="left")

        # Atom Selection
        tk.Label(f_dcc, text="Atom Selection (MDAnalysis format):",
                 font=label_font, bg="#ffffff", fg="black", anchor="w").pack(fill="x")
        create_readable_entry(
            f_dcc,
            textvariable=self.dcc_selection,
            font=entry_font
        ).pack(fill="x", pady=(0, 2))

        # Parameters row
        param_f = tk.Frame(f_dcc, bg="#ffffff")
        param_f.pack(fill="x", pady=1)

        for col, (lbl, var, w) in enumerate([
            ("Window",  self.dcc_window, 7),
            ("Step",    self.dcc_step,   4),
            ("Cores",   self.dcc_cores,  4),
        ]):
            tk.Label(param_f, text=lbl, font=label_font,
                     bg="#ffffff", fg="black").grid(row=0, column=col*2, sticky="w", padx=(4,1))
            create_readable_entry(
                param_f,
                textvariable=var,
                width=w,
                font=entry_font
            ).grid(row=0, column=col*2+1, padx=(0,6))

        # Log box
        self.dcc_log = tk.Text(
            f_dcc, height=2, font=("Courier", 9),
            bg="#f0f4f8", fg="#1a1a1a",
            insertbackground="#1a1a1a",
            selectbackground="#1f6fd6",
            selectforeground="#ffffff",
            relief="flat", state="disabled", wrap="word"
        )
        self.dcc_log.pack(fill="x", pady=(2, 1))

        # Compute button + OS-independent progress bar
        # Do not use ttk.Progressbar here: macOS Aqua / Windows themes can
        # recolour ttk widgets and produce the black/white glitch.
        self.dcc_progress = BlueProgressBar(f_dcc, height=12)
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
                                font=section_font, bg="#ffffff", fg="black", pady=4, padx=10)
        f_input.pack(fill="x", pady=2)
        entry_frame = tk.Frame(f_input, bg="#ffffff")
        entry_frame.pack(fill="x")
        create_readable_entry(entry_frame, textvariable=self.input_dir,
                 font=entry_font).pack(side="left", fill="x", expand=True, padx=(0, 5))
        create_browse_button(entry_frame, self.select_input).pack(side="left")

        # ── SECTION 2: Residue Selection ───────────────
        f_res = tk.LabelFrame(self.left_panel, text="Residue Limits ",
                              font=section_font, bg="#ffffff", fg="black", pady=4, padx=10)
        f_res.pack(fill="x", pady=2)
        self.range_entry(f_res, "Residue i range:", self.i_min, self.i_max,
                         label_font, entry_font)
        self.range_entry(f_res, "Residue j range:", self.j_min, self.j_max,
                         label_font, entry_font)
        self.single_entry(f_res, "Min distance |j - i| >=", self.res_dist,
                          label_font, entry_font)

        # ── SECTION 3: Correlation Filters ─────────────
        f_adv = tk.LabelFrame(self.left_panel, text="Correlation Filters ",
                              font=section_font, bg="#ffffff", fg="black", pady=4, padx=10)
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
                              font=section_font, bg="#ffffff", fg="black", pady=4, padx=10)
        f_csv.pack(fill="x", pady=2)
        create_modern_button(f_csv, "Download Combined CSV", "#7f8047",
                             self.download_combined_csv).pack(fill="x", pady=1)
        create_modern_button(f_csv, "Download Positive CSV", "#301f87",
                             self.download_positive_csv).pack(fill="x", pady=1)
        create_modern_button(f_csv, "Download Negative CSV", "#c0392b",
                             self.download_negative_csv).pack(fill="x", pady=1)

        # ── Right plot area  (original, unchanged) ───────────────────────────
        self.plot_frame = tk.Frame(self.root, bg="#f4f6f9", highlightthickness=0, bd=0)
        self.plot_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True,
                             padx=10, pady=10)

        self.toolbar_frame = tk.Frame(
            self.plot_frame,
            bg="#d9dde3",
            highlightbackground="#b8bec7",
            highlightthickness=1,
            bd=0,
            relief="flat",
        )
        self.toolbar_frame.pack(side=tk.BOTTOM, fill=tk.X)

        self.fig_dash = plt.Figure(figsize=(10, 8), dpi=100)
        self.fig_dash.patch.set_facecolor("#ffffff")

        self.canvas_dash = FigureCanvasTkAgg(self.fig_dash, master=self.plot_frame)
        self.canvas_dash.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)
        self.canvas_dash.mpl_connect("pick_event", self._on_heatmap_pick)
        self.canvas_dash.draw()

    # -----------------------------------------------------------------------
    # New: DCC file browsers
    # -----------------------------------------------------------------------
    def _browse_trajectory(self):
        paths = filedialog.askopenfilenames(
            title="Select trajectory file(s)",
            filetypes=[
                ("Trajectory files", "*.xtc *.dcd *.trr *.nc *.ncdf"),
                ("All files", "*.*")
            ]
        )
        if paths:
            self.traj_path.set("; ".join(paths))

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
        sel_str = self.dcc_selection.get().strip()

        if not traj:
            messagebox.showerror("Error", "Please select a valid trajectory file.")
            return
        if not topo or not os.path.isfile(topo):
            messagebox.showerror("Error", "Please select a valid topology file.")
            return
        if not out_dir:
            messagebox.showerror("Error", "Please select an output directory.")
            return
        if not sel_str:
            messagebox.showerror("Error", "Please enter a valid atom selection.")
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
            args=(traj, topo, step, window, out_dir, cores, sel_str,
                  self._dcc_log, _done),
            daemon=True
        ).start()

    # -----------------------------------------------------------------------
    # UI Helpers 
    # -----------------------------------------------------------------------
    def range_entry(self, parent, text, v1, v2, lf, ef):
        f = tk.Frame(parent, bg="#ffffff")
        f.pack(anchor="w", fill="x", pady=1)
        tk.Label(f, text=text, font=lf, bg="#ffffff", fg="black", anchor="w").grid(
            row=0, column=0, sticky="w")
        create_readable_entry(f, textvariable=v1, width=6, font=ef).grid(row=0, column=1, padx=(5, 2))
        tk.Label(f, text="to", font=lf, bg="#ffffff", fg="black").grid(row=0, column=2, padx=2)
        create_readable_entry(f, textvariable=v2, width=6, font=ef).grid(row=0, column=3, padx=2)
        f.grid_columnconfigure(0, minsize=165)

    def single_entry(self, parent, text, var, lf, ef):
        f = tk.Frame(parent, bg="#ffffff")
        f.pack(anchor="w", fill="x", pady=1)
        tk.Label(f, text=text, font=lf, bg="#ffffff", fg="black", anchor="w").grid(
            row=0, column=0, sticky="w")
        create_readable_entry(f, textvariable=var, width=6, font=ef).grid(row=0, column=1, padx=(5, 2))
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
        # ax.format_coord = lambda x, y: ""

    # -----------------------------------------------------------------------
    # Cross-platform Matplotlib toolbar
    # -----------------------------------------------------------------------
    def _build_cross_platform_mpl_toolbar(self, parent, canvas):
        """Use real Matplotlib navigation with theme-independent blue controls.

        The visible controls use Matplotlib's official PNG icons, while the
        hidden NavigationToolbar2Tk object remains the functional controller.
        This avoids macOS Light/Dark Aqua recolouring of native Tk buttons.
        """
        try:
            from PIL import Image, ImageTk
            _PIL_AVAILABLE = True
        except ImportError:
            Image = ImageTk = None
            _PIL_AVAILABLE = False

        TOOLBAR_BG = "#d9dde3"
        BUTTON_BLUE = "#1f6fd6"
        BUTTON_HOVER = "#185db5"
        BUTTON_ACTIVE = "#124f9d"
        BORDER = "#164f96"

        controller = NavigationToolbar2Tk(
            canvas, parent, pack_toolbar=False
        )
        controller.update()
        controller.set_message = lambda msg: None

        try:
            controller.pack_forget()
        except tk.TclError:
            pass

        bar = tk.Frame(parent, bg=TOOLBAR_BG, bd=0, relief="flat")
        bar.pack(side=tk.LEFT, padx=3, pady=2)

        image_dir = Path(matplotlib.get_data_path()) / "images"
        icon_files = {
            "home": "home.png",
            "back": "back.png",
            "forward": "forward.png",
            "pan": "move.png",
            "zoom": "zoom_to_rect.png",
            "subplots": "subplots.png",
            "save_figure": "filesave.png",
        }

        icon_refs = {}

        def load_white_icon(filename):
            if not _PIL_AVAILABLE:
                return None
            path = image_dir / filename
            with Image.open(path).convert("RGBA") as source:
                alpha = source.getchannel("A")
                white = Image.new("RGBA", source.size, (255, 255, 255, 255))
                white.putalpha(alpha)
                return ImageTk.PhotoImage(white, master=bar)

        for key, filename in icon_files.items():
            try:
                icon_refs[key] = load_white_icon(filename)
            except Exception:
                icon_refs[key] = None

        buttons = {}
        active_key = [None]

        def repaint():
            for key, holder in buttons.items():
                bg = BUTTON_ACTIVE if key == active_key[0] else BUTTON_BLUE
                try:
                    holder.configure(bg=bg, highlightbackground=bg)
                    holder._icon_label.configure(bg=bg)
                except tk.TclError:
                    pass

        def add_icon_button(key, command):
            holder = tk.Frame(
                bar, bg=BUTTON_BLUE, width=36, height=32,
                bd=0, relief="flat", highlightthickness=1,
                highlightbackground=BUTTON_BLUE,
                highlightcolor=BUTTON_BLUE, cursor="hand2"
            )
            holder.pack(side=tk.LEFT, padx=2, pady=1)
            holder.pack_propagate(False)

            photo = icon_refs.get(key)
            if photo is not None:
                label = tk.Label(
                    holder, image=photo, bg=BUTTON_BLUE, bd=0,
                    relief="flat", highlightthickness=0, cursor="hand2"
                )
                label._toolbar_photo = photo
            else:
                fallback = {
                    "home": "⌂", "back": "←", "forward": "→",
                    "pan": "✚", "zoom": "⌕", "subplots": "≡",
                    "save_figure": "▣",
                }.get(key, "")
                label = tk.Label(
                    holder, text=fallback, font=("Helvetica", 16, "bold"),
                    fg="#ffffff", bg=BUTTON_BLUE, bd=0,
                    relief="flat", highlightthickness=0, cursor="hand2"
                )

            label.pack(fill=tk.BOTH, expand=True)
            holder._icon_label = label

            def enter(_event=None, holder=holder, label=label, key=key):
                if active_key[0] != key:
                    try:
                        holder.configure(
                            bg=BUTTON_HOVER,
                            highlightbackground=BORDER
                        )
                        label.configure(bg=BUTTON_HOVER)
                    except tk.TclError:
                        pass

            def leave(_event=None, holder=holder, label=label, key=key):
                try:
                    bg = BUTTON_ACTIVE if active_key[0] == key else BUTTON_BLUE
                    holder.configure(bg=bg, highlightbackground=bg)
                    label.configure(bg=bg)
                except tk.TclError:
                    pass

            def click(_event=None, command=command):
                command()

            for widget in (holder, label):
                widget.bind("<Enter>", enter)
                widget.bind("<Leave>", leave)
                widget.bind("<Button-1>", click)

            buttons[key] = holder

        def do_home():
            active_key[0] = None
            repaint()
            controller.home()

        def do_back():
            active_key[0] = None
            repaint()
            controller.back()

        def do_forward():
            active_key[0] = None
            repaint()
            controller.forward()

        def do_pan():
            controller.pan()
            mode = str(getattr(controller, "mode", "")).lower()
            active_key[0] = "pan" if "pan" in mode else None
            repaint()

        def do_zoom():
            controller.zoom()
            mode = str(getattr(controller, "mode", "")).lower()
            active_key[0] = "zoom" if "zoom" in mode else None
            repaint()

        def do_subplots():
            active_key[0] = None
            repaint()
            controller.configure_subplots()

        def do_save():
            active_key[0] = None
            repaint()
            controller.save_figure()

        for key, fn in [
            ("home", do_home),
            ("back", do_back),
            ("forward", do_forward),
            ("pan", do_pan),
            ("zoom", do_zoom),
            ("subplots", do_subplots),
            ("save_figure", do_save),
        ]:
            add_icon_button(key, fn)

        parent._mpl_controller = controller
        parent._mpl_icon_refs = icon_refs
        parent._mpl_icon_buttons = buttons
        parent._mpl_visible_bar = bar

        tk.Frame(parent, bg=TOOLBAR_BG, bd=0,
                 highlightthickness=0).pack(
            side=tk.LEFT, fill=tk.X, expand=True
        )

        return controller

    def _reset_toolbar(self):
        if self.toolbar:
            try:
                self.toolbar.destroy()
            except tk.TclError:
                pass

        for child in self.toolbar_frame.winfo_children():
            try:
                child.destroy()
            except tk.TclError:
                pass

        self.toolbar_frame.configure(
            bg="#d9dde3",
            highlightbackground="#b8bec7",
            highlightthickness=1,
            bd=0,
            relief="flat",
        )

        self.toolbar = self._build_cross_platform_mpl_toolbar(
            self.toolbar_frame, self.canvas_dash
        )

    # -----------------------------------------------------------------------
    # Dashboards  (original, unchanged)
    # -----------------------------------------------------------------------
    def open_correlation_dashboard(self):
        if self.df_all is None or self.df_all.empty:
            messagebox.showerror("Error", "Run analysis first.")
            return

        self.fig_dash.clear()
        self.fig_dash.patch.set_facecolor("#ffffff")

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
        
        # Extract correct X-axis frames from filenames (formerly: Calculate actual frame numbers)
        x_values = []
        for f in self.frame_files:
            match = re.search(r'corr_dcc_(\d+)\.txt', f)
            if match:
                x_values.append(int(match.group(1)))
            else:
                x_values.append(len(x_values) + 1)
        
        # Calculate 95% Confidence Intervals using Fisher z-transformation
        window_size = self.dcc_window.get()
        r_safe = np.clip(ts, -0.9999, 0.9999)
        z = np.arctanh(r_safe)
        
        if window_size > 3:
            se = 1.0 / np.sqrt(window_size - 3)
            ci_upper = np.tanh(z + 1.96 * se)
            ci_lower = np.tanh(z - 1.96 * se)
        else:
            ci_upper = ts
            ci_lower = ts

        fig = plt.Figure(figsize=(8, 4), dpi=110)
        ax  = fig.add_subplot(111)
        
        # MODIFIED: Pass x_values to the plot
        ax.plot(x_values, ts, lw=2, color="#2980b9", label="Correlation")
        #ax.plot(ts, lw=2, color="#2980b9")
        ax.fill_between(x_values, ci_lower, ci_upper, color="#2980b9", alpha=0.3, label="95% CI")
        
        ax.set_title(f"Frame-wise Correlation  i={i}, j={j}",
                     fontsize=13, fontweight="bold")
        ax.set_xlabel("Frame", fontsize=11, fontweight="bold")
        ax.set_ylabel("Correlation", fontsize=11, fontweight="bold")
        ax.grid(True, linestyle="--", linewidth=0.8, alpha=0.6)
        ax.legend(loc="best")
        
        # NEW LINE: Hides x, y values on hover for the time series line graph
        # ax.format_coord = lambda x, y: ""

        canvas = FigureCanvasTkAgg(fig, master=win)
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        toolbar_frame = tk.Frame(
            win,
            bg="#d9dde3",
            highlightbackground="#b8bec7",
            highlightthickness=1,
            bd=0,
            relief="flat",
        )
        toolbar_frame.pack(side=tk.BOTTOM, fill=tk.X, padx=8, pady=(0, 4))
        self._build_cross_platform_mpl_toolbar(toolbar_frame, canvas)

        canvas.draw()

        btn_frame = tk.Frame(win, bg="#d9dde3")
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
                    "Correlation": ts,
                    "CI_Lower":    ci_lower,
                    "CI_Upper":    ci_upper
                })
                df_ts.to_csv(path, index=False)
                messagebox.showinfo("Saved", f"CSV saved:\n{path}", parent=win)
            except Exception as e:
                messagebox.showerror("Error", f"Could not save CSV:\n{e}", parent=win)

        create_modern_button(
            btn_frame,
            "Save CSV",
            "#27ae60",
            save_csv,
            font=("Helvetica", 10, "bold"),
            padx=16,
            pady=5,
        ).pack(side=tk.LEFT, padx=5)
        tk.Label(btn_frame,
                 text=f"i={i}     j={j}    "# |  frames={len(ts)}  |  "
                      f"mean={np.mean(ts):.4f}    std={np.std(ts):.4f}",
                 font=("Helvetica", 9), fg="#333333", bg="#d9dde3").pack(
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
