# TD-DCC — Time-Dependent Dynamic Cross-Correlation Analysis

**TD-DCC** is a standalone desktop application for computing and visualizing time-resolved dynamic cross-correlation matrices (DCCMs) from molecular dynamics (MD) trajectories.

Unlike a conventional DCCM — which collapses an entire trajectory into a single static matrix — TD-DCC computes correlations over a **sliding window of frames** and tracks how each residue pair's correlation evolves across the simulation: whether it stays consistently positive, consistently negative, or fluctuates between the two.

> For full details, see the [User Manual](TD-DCC_User_Manual.md).

---

## Features

- **Built-in compute engine** — load a trajectory + topology directly; correlation matrices are computed in parallel across sliding windows
- **Statistics dashboard** — mean correlation, standard deviation, and positive/negative occurrence probability for every residue pair
- **Interactive heatmaps** — filtered Positive / Negative / Combined correlation maps with pan, zoom, and save
- **Per-pair time series** — double-click any residue pair to inspect its correlation trajectory window-by-window
- **CSV export** — download filtered results for downstream analysis in Python, R, or a spreadsheet
- **Fully local** — no internet connection or API key needed at any point

---

## System Requirements

- Python 3.9 or later (3.10–3.12 recommended)
- Windows, MacOs or Linux
- 8 GB RAM or more recommended for large proteins
- A graphical desktop session (Tkinter GUI; X11 forwarding supported for remote HPC use)

---

## Installation

### 1. Install system dependencies

**Ubuntu / Debian**
```bash
sudo apt update
sudo apt install python3 python3-pip python3-venv python3-tk
```

**Fedora / RHEL**
```bash
sudo dnf install python3 python3-pip python3-tkinter
```

**Windows** — download from [python.org/downloads](https://python.org/downloads) and check *"Add python.exe to PATH"* during install. Tkinter is included automatically.

**MacOs** — download the latest macOS universal installer (`.pkg`) from [python.org/downloads/macos](https://www.python.org/downloads/macos/).

---

### 2. Create a virtual environment (recommended)

**Conda**
```bash
conda create -n TD-DCC_env python=3.12
conda activate TD-DCC_env
```

**venv**
```bash
python3 -m venv TD-DCC_env
source TD-DCC_env/bin/activate      # Linux
TD-DCC_env\Scripts\activate         # Windows
```

---

### 3. Install Python dependencies

```bash
pip install numpy pandas matplotlib
```

To use the built-in trajectory compute engine, also install MDAnalysis:

```bash
pip install MDAnalysis
```

---

### 4. Clone the repository

```bash
git clone https://github.com/sandipandrils/TD-DCC.git
cd TD-DCC
```

---

## Quick Start

```bash
python3 TD-DCC_gui.py
```

Click **"Launch Analysis Dashboard"** on the intro screen to open the main workspace.

---

## Workflows


1. Browse to your **Trajectory** file (`.xtc`, `.dcd`, `.trr`, or any MDAnalysis-supported format)
2. Browse to your **Topology** file (`.pdb`, `.tpr`, `.gro`, etc.)
3. Browse to an **Output directory** for the computed matrices
4. Set **Window**, **Step**, and **Cores**, then click **COMPUTE DCC**
5. TD-DCC fills the Input Data path automatically when done — proceed to analysis
6. Point **Input Data** at a folder of pre-computed `.txt` matrix files
7. Set **Residue Limits** and **Correlation Filters**
8. Click **RUN ANALYSIS** → **View Correlation Dashboard**

> **TIP:** Double-click any cell in a heatmap to open the per-pair time series inspector and confirm whether a correlation is persistent or transient.

> **NOTE:** To calculate DCC of the entire trajectory, instead of the time windowed ones, put the total number of frames in the window. This will calculate a single time-averaged DCC plot from a .txt matrix file.

---

## Interface Overview

The main workspace is a 1400×900 window with a scrollable control panel on the left and a plotting area on the right. Control panel sections from top to bottom:

| Section | Purpose |
|---|---|
| **Compute DCC from Trajectory** | Built-in engine to generate matrices from a trajectory |
| **Input Data** | Folder of `.txt` correlation matrices to analyze |
| **Residue Limits** | Restrict the residue-pair range; exclude sequence-adjacent trivial correlations |
| **Correlation Filters** | Set mean correlation thresholds and occurrence probability cutoffs |
| **Run Analysis** | Load matrices and compute per-pair statistics |
| **View Correlation Dashboard** | Display filtered Positive / Negative / Combined heatmaps |
| **Export Data** | Save filtered results to CSV |

> **NOTE:** Scroll the left panel down to reach Run Analysis, View Correlation Dashboard, and Export Data.

---

## Output Files

| File | Description |
|---|---|
| `corr_dcc_<frame>.txt` | Plain-text N×N correlation matrix for one sliding window |
| `corr_dcc_<frame>.png` | Heatmap image of the same window (blue–white–red colormap) |
| `*_positive.csv` / `*_negative.csv` / `*_combined.csv` | Filtered residue-pair statistics: mean correlation, SD, occurrence probability |
| Time-series PNG / CSV | Per-pair correlation trajectory exported from the inspector |

---

## Sample Data

The `example/trajectory_and_topology/` directory contains a 500 ns AAV2 MD trajectory split into 5 parts to comply with GitHub file size limits. 

You can load all segments as a single, continuous trajectory in one of two ways:
* **Multi-file Selection:** Open the trajectory browser and hold `Ctrl` (or `Cmd` on macOS) to select multiple `.xtc` files at once.
* **Wildcard Pattern:** Type a glob pattern directly into the trajectory field (e.g., `traj_AAV2_500ns_part*of5.xtc`). 

Matched files are automatically expanded, naturally sorted, and passed as a list to `MDAnalysis.Universe()`, which chains the segments together seamlessly for analysis in a single run.


---

## Troubleshooting

| Symptom | Fix |
|---|---|
| `No module named 'tkinter'` | `sudo apt install python3-tk` (Linux) |
| COMPUTE DCC fails / MDAnalysis unavailable | `pip install MDAnalysis` |
| "Please select a valid trajectory / topology file" | Browse to the correct file before clicking COMPUTE DCC |
| "No .txt files found in selected folder" | Point Input Data at the matrix output folder from the Compute DCC step |
| "No valid matrix data with selected range" | Ensure residue j ≥ residue i in Residue Limits |
| Buttons cut off at bottom of screen | Scroll the left control panel downward |

For the full troubleshooting list, see [Section 8 of the User Manual](TD-DCC_User_Manual.md#8-troubleshooting).

---

## Citation

If you use TD-DCC in your research, please cite:

> Prasun Pal, Chaitanya Chintalapati, Sandipan Chakraborty. TD-DCC: An Interactive and Efficient Graphical User Interface-Based Time-Dependent Cross-Correlation Analysis from Molecular Dynamics Trajectories. ChemRxiv. 03 July 2026.
DOI: https://doi.org/10.26434/chemrxiv.15005626/v1

TD-DCC uses **MDAnalysis** for trajectory parsing. If you use the built-in compute engine, please also cite:

> Michaud-Agrawal, N., Denning, E.J., Woolf, T.B. and Beckstein, O. (2011). MDAnalysis: A toolkit for the analysis of molecular dynamics simulations. *Journal of Computational Chemistry*, 32(10): 2319–2327. https://doi.org/10.1002/jcc.21787

---

## License

*License details to be added.*
