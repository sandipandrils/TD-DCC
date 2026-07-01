# TD-DCC User Manual

**Time-Dependent Dynamic Cross-Correlation Analysis**

*User Manual & Installation Guide*

Document version 1.0  
Covers `TD-DCC_gui.py` — standalone desktop edition

---

## Contents

1. [Introduction](#1-introduction)
2. [System Requirements](#2-system-requirements)
3. [Installation Guide](#3-installation-guide)
   - 3.1 [Install Python](#31-install-python)
   - 3.2 [Create a virtual environment (recommended)](#32-create-a-virtual-environment-recommended)
   - 3.3 [Install required packages](#33-install-required-packages)
   - 3.4 [Download TD-DCC](#34-download-td-dcc)
   - 3.5 [Verify the installation](#35-verify-the-installation)
4. [Launching TD-DCC](#4-launching-td-dcc)
5. [Interface Tour](#5-interface-tour)
6. [Step-by-Step Workflows](#6-step-by-step-workflows)
   - 6.1 [Workflow — Starting from an MD trajectory](#61-workflow--starting-from-an-md-trajectory)
   - 6.2 [Setting residue limits and correlation filters](#62-setting-residue-limits-and-correlation-filters)
   - 6.3 [Running the analysis and reading the dashboard](#63-running-the-analysis-and-reading-the-dashboard)
   - 6.4 [Inspecting a residue pair's time series](#64-inspecting-a-residue-pairs-time-series)
   - 6.5 [Exporting results](#65-exporting-results)
7. [Output Files](#7-output-files)
8. [Troubleshooting](#8-troubleshooting)

---

## 1. Introduction

TD-DCC (Time-Dependent Dynamic Cross-Correlation) is a standalone desktop application for computing and exploring dynamical cross-correlation matrices (DCCMs) from molecular dynamics (MD) trajectories. Unlike a conventional DCCM, which averages an entire trajectory into a single static matrix, TD-DCC computes correlations over a sliding window of frames and tracks how each residue pair's correlation behaves across the whole trajectory — whether it stays consistently positive, consistently negative, or fluctuates between the two.

The application provides everything needed for this workflow in a single window: an optional built-in trajectory-to-matrix compute engine, a statistics engine that aggregates per-window matrices into mean correlation, standard deviation, and positive/negative occurrence probability for every residue pair, and an interactive dashboard for filtering, visualizing, and exporting the results.

### TD-DCC System Architecture

The diagram below summarizes the overall pipeline: a trajectory is streamed and split into windows; each window's correlation matrix is computed in parallel; and the resulting matrices are aggregated into statistics that drive the interactive dashboard.

<img width="4500" height="3031" alt="Fig1" src="https://github.com/user-attachments/assets/fbab437f-776b-462b-be74-9fd4d692be0c" />

*Figure 1 — Overview of the TD-DCC compute and visualization pipeline.*

---

## 2. System Requirements

TD-DCC is a Python desktop application that runs on Windows and Linux. It has modest hardware requirements:

- Python 3.9 or later (3.10–3.12 recommended).
- 8 GB RAM or more recommended for proteins larger than a few hundred residues.
- A graphical desktop session. TD-DCC uses Tkinter and will not run over a plain text-only SSH connection without a display (a remote desktop with X11 forwarding works).
- All computation in TD-DCC runs locally; no internet connection or API key is needed at any point.

---

## 3. Installation Guide

This section walks through installing everything TD-DCC needs on Windows and Linux. If you already have a working Python 3.9+ environment, you can skip ahead to [Section 3.3](#33-install-required-packages).

### 3.1 Install Python

#### Windows

1. Download the latest Python 3 installer from [python.org/downloads](https://python.org/downloads).
2. Run the installer and check **"Add python.exe to PATH"** at the bottom of the first screen before clicking **Install Now**.
3. Open Command Prompt and confirm the install:

```bash
python --version
```

> **NOTE:** The standard Windows installer from python.org includes Tkinter automatically — no extra step is needed.

#### Linux (Ubuntu / Debian)

1. Install Python, pip, and the Tkinter system package (Tkinter is not bundled with Python on Debian-based distributions and must be installed separately):

```bash
sudo apt update
sudo apt install python3 python3-pip python3-venv python3-tk
```

#### Linux (Fedora / RHEL)

```bash
sudo dnf install python3 python3-pip python3-tkinter
```

### 3.2 Create a virtual environment (recommended)

A virtual environment keeps TD-DCC's dependencies isolated from other Python projects on your system. This step is optional but strongly recommended. You can use either conda or the standard `venv` module.

**Option A: Using Conda (Anaconda / Miniconda)**

Conda is an excellent choice for isolating environments, especially when deploying on Slurm-based HPC clusters or managing complex scientific pipelines.

```bash
conda create -n TD-DCC_env python=3.12

# Activate it:
conda activate TD-DCC_env
```

**Option B: Using standard venv**

```bash
python3 -m venv TD-DCC_env

# Activate it:
source TD-DCC_env/bin/activate       # Linux
TD-DCC_env\Scripts\activate          # Windows (Command Prompt)
```

Your terminal prompt should now show `(TD-DCC_env)` at the start of the line, confirming the environment is active. Repeat the activation command every time you open a new terminal to work with TD-DCC.

### 3.3 Install required packages

With your environment activated, install the core dependencies:

```bash
pip install numpy pandas matplotlib
```

If you want to use the built-in "Compute DCC from Trajectory" panel to generate correlation matrices directly from a simulation, also install MDAnalysis:

```bash
pip install MDAnalysis
```

> **NOTE:** On recent Ubuntu/Debian systems running outside a virtual environment, pip may refuse to install with an "externally-managed-environment" error. Either use a virtual environment as in Section 3.2 (recommended), or add `--break-system-packages` to the pip command.

Tkinter is part of the Python standard library and does not need to be installed via pip; on Windows (python.org installer) it is included automatically, and on Linux it is installed via your system package manager as shown in Section 3.1.

### 3.4 Download TD-DCC

TD-DCC is distributed as a single Python file, `TD-DCC_gui.py`. Save it anywhere convenient on your computer — for example, a folder named `TD-DCC` in your home directory or Documents folder. No further files are required to run the application.

### 3.5 Verify the installation

From the folder containing `TD-DCC_gui.py`, with your virtual environment activated, run:

```bash
python TD-DCC_gui.py        # Windows
python3 TD-DCC_gui.py       # Linux
```

A small window titled "Dynamic Cross-Correlation Analysis" should appear within a few seconds (Figure 2). If it does, installation is complete and you can proceed to Section 4.

<img width="811" height="638" alt="Screenshot 2026-06-30 104647" src="https://github.com/user-attachments/assets/39d5f957-78f3-4756-9e05-caab5e3f3132" />

*Figure 2 — The TD-DCC intro screen confirms a successful installation.*

> **CAUTION:** If nothing appears and you instead see an error in the terminal, check [Section 8 (Troubleshooting)](#8-troubleshooting) for the most common causes.

---

## 4. Launching TD-DCC

Every time you want to use TD-DCC:

1. Open a terminal (Command Prompt, Terminal, or your shell of choice).
2. If you created a virtual environment, activate it (Section 3.2).
3. Navigate to the folder containing `TD-DCC_gui.py` using `cd`.
4. Run the script, then click **"Launch Analysis Dashboard"** on the screen that appears to open the main workspace.

```bash
cd path/to/your/TD-DCC/folder
python3 TD-DCC_gui.py
```

---

## 5. Interface Tour

Clicking "Launch Analysis Dashboard" opens the main workspace — a single 1400×900 window split into two halves: a scrollable control panel on the left, and a plotting area on the right that stays empty until you run an analysis.

<img width="3798" height="2716" alt="Fig3" src="https://github.com/user-attachments/assets/13d951de-3bc1-4d0d-b914-cd44015d4238" />


*Figure 3 — The top of the control panel: sections related to Run Analysis, View Correlation Dashboard, and Export Data.*

The control panel is organized into different sections from top to bottom:

- **Compute DCC from Trajectory** — runs the built-in compute engine on a trajectory file.
- **Input Data** — the folder of correlation matrices to analyze.
- **Residue Limits** — which residue-pair range to include.
- **Correlation Filters** — the thresholds used to flag a pair as positively or negatively correlated.
- **Run Analysis**, **View Correlation Dashboard**, and **Export Data** — towards the bottom of the panel (scroll down to reach them).

> **NOTE:** The control panel is taller than the window. Use the thin scrollbar at the right edge of the panel, or drag it directly, to reach the lower sections.

---

## 6. Step-by-Step Workflows

### 6.1 Workflow — Starting from an MD trajectory

Use this workflow if you have a raw trajectory and topology file and want TD-DCC to compute the correlation matrices for you. This requires MDAnalysis to be installed (Section 3.3).

1. In the first section, click **Browse** next to **Trajectory** and select your trajectory file (`.xtc`, `.dcd`, `.trr`, or any format MDAnalysis supports).
2. Click **Browse** next to **Topology** and select the matching topology file (`.pdb`, `.tpr`, `.gro`, etc.).
3. Click **Browse** next to **Output directory** and choose where the resulting matrices should be written.
4. Set **Window** (the number of frames per sliding window), **Step** (the frame stride used when streaming the trajectory), and **Cores** (how many windows to compute in parallel). The defaults (Window 100, Step 1) are a reasonable starting point for most trajectories.
5. Click **COMPUTE DCC**. Progress is reported in the log box and the progress bar; this can take anywhere from a few seconds to several minutes, depending on trajectory length and Window/Step settings.
6. When finished, TD-DCC automatically fills the next section (Input Data) with the output directory, so you can go straight to Section 6.2.

<img width="4500" height="1455" alt="Fig4_manual" src="https://github.com/user-attachments/assets/f6b39d7f-7ba7-468f-9862-3ef7cfdb2ba3" />

*Figure 4 — The dashboard, showing DCC calculation for a sample data set.*

> **NOTE:** If MDAnalysis is not installed, it will report an error when you click COMPUTE DCC. Install it with `pip install MDAnalysis`.

### 6.2 Setting residue limits and correlation filters

**Residue Limits** controls which part of the matrix is analyzed:

- **Residue i range / Residue j range** — the rows and columns to include. Narrowing these speeds up analysis for very large proteins.
- **Min distance |j − i| ≥** — excludes residue pairs closer together than this threshold along the sequence, which are almost always trivially correlated and are not usually of biological interest.

**Correlation Filters** controls which residue pairs are highlighted once the statistics are computed:

- **Positive correlation / Negative correlation** — the range of mean correlation values that count as positively or negatively correlated.
- **Positive corr. Prob. / Negative corr. Prob.** — the range of occurrence probability — the fraction of windows in which a pair showed that sign of correlation — required for a pair to be included. A high probability means the correlation is consistent across the trajectory rather than driven by a few outlier windows.

<img width="4388" height="1358" alt="Fig5_manual" src="https://github.com/user-attachments/assets/d3a85ffb-15d0-43df-b1d8-9d61aa4012df" />

*Figure 5 — The populated dashboard, showing residue limit settings and running the analysis for a sample data set.*

> **NOTE:** The default filters (correlation magnitude ≥ 0.2, probability ≥ 0.1) are a reasonable starting point. You can re-enter different values and click Run Analysis again at any time without reloading your data.

### 6.3 Running the analysis and reading the dashboard

1. Scroll down and click **RUN ANALYSIS**. TD-DCC loads every matrix in the Input Data folder and computes, for every qualifying residue pair, the mean correlation, standard deviation, and positive/negative occurrence probability across all windows.
2. When the status line reads "Analysis complete!", click **View Correlation Dashboard**.
3. Three heatmaps appear in the right-hand panel: **Positive Correlation**, **Negative Correlation**, and **Combined Correlation**, each showing only the residue pairs that satisfy the filters from Section 6.2.

<img width="4398" height="1386" alt="Fig6_manual" src="https://github.com/user-attachments/assets/48eaf38d-30b9-4c6f-9c02-f8423ee02055" />

*Figure 6 — The populated dashboard, showing positive, negative, and combined correlation heatmaps for a sample data set.*

Each heatmap supports the standard Matplotlib navigation toolbar beneath it, so you can pan, zoom, and save any panel as an image independently of the rest of the dashboard.

### 6.4 Inspecting a residue pair's time series

To examine how a specific residue pair's correlation behaves across the trajectory rather than just its summary statistic, double-click that cell in any of the three heatmaps. A new window opens showing the pair's correlation value across the frames, along with its mean and standard deviation. From this window you can save the plot as a PNG or export the underlying values as a CSV.

<img width="4398" height="1394" alt="Fig7_manual" src="https://github.com/user-attachments/assets/04920bdf-8fb9-4052-b6a3-9c2b1ca4c08d" />

*Figure 7 — The dashboard, showing the residue pair time series for a sample data set.*

> **TIP:** This view is the best way to confirm that a filtered hit reflects a persistent, biologically meaningful correlation rather than one driven by a handful of unusual windows.

### 6.5 Exporting results

The **Export Data** section provides three buttons:

- **Download Combined CSV** — every residue pair satisfying either the positive or the negative filter.
- **Download Positive CSV** — pairs satisfying only the positive filter.
- **Download Negative CSV** — pairs satisfying only the negative filter.

Each CSV contains one row per residue pair with columns for residue i, residue j, mean correlation, standard deviation, and positive/negative occurrence probability — ready for further analysis in a spreadsheet, R, or Python.

---

## 7. Output Files

| File | Description |
|---|---|
| `corr_dcc_<frame>.txt` | Plain-text N×N correlation matrix for one window, named by the index of the window's last frame. Produced by "Compute DCC from Trajectory" and consumed by "Input Data". |
| `corr_dcc_<frame>.png` | A heatmap image of the same window's matrix, using TD-DCC's diverging blue–white–red colormap, saved alongside the text file for quick visual inspection without reopening the dashboard. |
| Combined / Positive / Negative CSV | Exported from the dashboard via the Export Data buttons; one row per qualifying residue pair with mean correlation, standard deviation, and occurrence probability. |
| Time-series PNG / CSV | Exported from the double-click time-series inspector window for a single residue pair. |

---

## 8. Troubleshooting

| Symptom | Likely cause and fix |
|---|---|
| `ModuleNotFoundError: No module named 'tkinter'` | Tkinter is missing from your Python installation. On Linux, install it with your system package manager (e.g. `sudo apt install python3-tk`). See Section 3.1. |
| "Please select a valid trajectory file" when clicking COMPUTE DCC | Browse to a folder containing the MD simulation trajectory first. |
| "Please select a valid topology file" when clicking COMPUTE DCC | Browse to a folder containing the MD simulation topology. |
| "Please select an output Directory" when clicking COMPUTE DCC | Create and browse to a folder where the time-windowed DCC images and `.txt` files will be generated. |
| COMPUTE DCC fails immediately or reports MDAnalysis is unavailable | MDAnalysis is not installed. Run `pip install MDAnalysis` in the same environment you use to launch TD-DCC. |
| The window opens but appears blank or frozen while computing | Large trajectories or a high Window/Step/Cores combination can take time. Progress is reported in the log box and the main status line; avoid closing the window while a computation is in progress. |
| "Invalid input directory" when clicking Run Analysis | Browse to a folder containing the `.txt` matrices first. |
| "No .txt files found in selected folder" when clicking Run Analysis | The selected Input Data folder does not contain any plain-text matrix files. Confirm you selected the correct output directory from the Compute DCC step. |
| "No valid matrix data with selected range" when clicking Run Analysis | In Residue Limits, j must be greater than or equal to i. Confirm you selected the residues properly. |
| Plots or buttons look cut off at the bottom of the screen | The control panel is scrollable. Drag the thin scrollbar at its right edge downward to reveal Run Analysis, View Correlation Dashboard, and Export Data. |

---

*End of document.*
