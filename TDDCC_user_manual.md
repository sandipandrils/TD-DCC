# TDDCC

**Time-Dependent Dynamical Cross-Correlation Analysis**

A standalone desktop application for computing and exploring dynamical cross-correlation matrices (DCCMs) from molecular dynamics (MD) trajectories.

Unlike a conventional DCCM, which averages an entire trajectory into a single static matrix, TDDCC computes correlations over a sliding window of frames and tracks how each residue pair's correlation behaves across the whole trajectory — whether it stays consistently positive, consistently negative, or fluctuates between the two.

The application provides everything needed for this workflow in a single window: an optional built-in trajectory-to-matrix compute engine, a statistics engine that aggregates per-window matrices into mean correlation, standard deviation, and positive/negative occurrence probability for every residue pair, and an interactive dashboard for filtering, visualizing, and exporting the results.

> **TIP** — You do not need to compute matrices inside TDDCC to use it. If you already have correlation matrices from another pipeline, you can point the dashboard directly at that folder — see [Workflow B](#62-workflow-b--starting-from-existing-correlation-matrices).

Covers `TDDCC_gui.py` — standalone desktop edition · Document version 1.0

---

## Table of Contents

1. [Introduction](#introduction)
2. [System Requirements](#system-requirements)
3. [Installation Guide](#installation-guide)
   - [3.1 Install Python](#31-install-python)
   - [3.2 Create a virtual environment (recommended)](#32-create-a-virtual-environment-recommended)
   - [3.3 Install required packages](#33-install-required-packages)
   - [3.4 Download TDDCC](#34-download-tddcc)
   - [3.5 Verify the installation](#35-verify-the-installation)
4. [Launching TDDCC](#4-launching-tddcc)
5. [Interface Tour](#5-interface-tour)
6. [Step-by-Step Workflows](#6-step-by-step-workflows)
   - [6.1 Workflow A — Starting from an MD trajectory](#61-workflow-a--starting-from-an-md-trajectory)
   - [6.2 Workflow B — Starting from existing correlation matrices](#62-workflow-b--starting-from-existing-correlation-matrices)
   - [6.3 Setting residue limits and correlation filters](#63-setting-residue-limits-and-correlation-filters)
   - [6.4 Running the analysis and reading the dashboard](#64-running-the-analysis-and-reading-the-dashboard)
   - [6.5 Inspecting a residue pair's time series](#65-inspecting-a-residue-pairs-time-series)
   - [6.6 Exporting results](#66-exporting-results)
7. [Output File Reference](#7-output-file-reference)
8. [Troubleshooting](#8-troubleshooting)
9. [Quick Reference](#9-quick-reference)

---

## Introduction

The diagram below summarizes the overall pipeline: a trajectory is streamed and split into windows, each window's correlation matrix is computed in parallel, and the resulting matrices are aggregated into the statistics that drive the interactive dashboard.

> *Figure 1 — Overview of the TDDCC compute and visualization pipeline.*

## System Requirements

TDDCC is a pure-Python desktop application and runs on Windows, macOS, and Linux. It has modest hardware requirements:

- Python 3.9 or later (3.10–3.12 recommended).
- 4 GB RAM minimum; 8 GB or more recommended for proteins larger than a few hundred residues.
- A graphical desktop session. TDDCC uses Tkinter and will not run over a plain text-only SSH connection without a display (a remote desktop with X11 forwarding works).
- All computation in TDDCC runs locally; no internet connection or API key is needed at any point.

> **NOTE** — MDAnalysis (used only by the optional "Compute DCC from Trajectory" panel) requires a small amount of additional disk space for its own dependencies. If you only plan to visualize matrices computed elsewhere, you can skip installing it entirely.

## Installation Guide

This section walks through installing everything TDDCC needs, on Windows, macOS, and Linux. If you already have a working Python 3.9+ environment, you can skip ahead to [Section 3.3](#33-install-required-packages).

### 3.1 Install Python

**Windows**

1. Download the latest Python 3 installer from [python.org/downloads](https://python.org/downloads).
2. Run the installer and check "Add python.exe to PATH" at the bottom of the first screen before clicking Install Now.
3. Open Command Prompt and confirm the install:

   ```bash
   python --version
   ```

   > **NOTE** — The standard Windows installer from python.org includes Tkinter automatically — no extra step is needed.

**macOS**

1. Install Python 3 either from [python.org/downloads](https://python.org/downloads) (recommended for the simplest Tkinter setup) or via Homebrew:

   ```bash
   brew install python-tk@3.12
   ```

2. Confirm the install in Terminal:

   ```bash
   python3 --version
   ```

   > **CAUTION** — If you installed Python via Homebrew without the matching `python-tk` formula, Tkinter will be missing and TDDCC will fail to start with a "No module named `_tkinter`" error. Installing `python-tk@<version>` for the same version you use resolves this.

**Linux (Ubuntu / Debian)**

Install Python, pip, and the Tkinter system package (Tkinter is not bundled with Python on Debian-based distributions and must be installed separately):

```bash
sudo apt update
sudo apt install python3 python3-pip python3-venv python3-tk
```

**Linux (Fedora / RHEL)**

```bash
sudo dnf install python3 python3-pip python3-tkinter
```

### 3.2 Create a virtual environment (recommended)

A virtual environment keeps TDDCC's dependencies isolated from other Python projects on your system. This step is optional but strongly recommended. You can use either conda or the standard `venv` module.

**Option A: Using Conda (Anaconda / Miniconda)**

Conda is an excellent choice for isolating environments, especially when deploying on Slurm-based HPC clusters or managing complex scientific pipelines.

Create a new conda environment specifically for TDDCC, specifying a compatible Python version, and activate it:

```bash
conda create -n tddcc_env python=3.12

# Activate it:
conda activate tddcc_env
```

**Option B: Using standard venv**

```bash
python3 -m venv tddcc_env

# Activate it:
source tddcc_env/bin/activate        # macOS / Linux
tddcc_env\Scripts\activate           # Windows (Command Prompt)
```

Your terminal prompt should now show `(tddcc_env)` at the start of the line, confirming the environment is active. Repeat the activation command every time you open a new terminal to work with TDDCC.

> **NOTE** — Anaconda or Miniconda can also be used similarly, and one can create separate environments accordingly.

### 3.3 Install required packages

With your environment activated, install the core dependencies:

```bash
pip install numpy pandas matplotlib
```

If you want to use the built-in "Compute DCC from Trajectory" panel to generate correlation matrices directly from a simulation, also install MDAnalysis:

```bash
pip install MDAnalysis
```

> **NOTE** — On recent Ubuntu/Debian systems running outside a virtual environment, pip may refuse to install with an "externally-managed-environment" error. Either use a virtual environment as in [Section 3.2](#32-create-a-virtual-environment-recommended) (recommended), or add `--break-system-packages` to the pip command.

Tkinter is part of the Python standard library and does not need to be installed via pip; on Windows and macOS (python.org installer) it is included automatically, and on Linux it is installed via your system package manager as shown in [Section 3.1](#31-install-python).

### 3.4 Download TDDCC

TDDCC is distributed as a single Python file, `TDDCC_gui.py`. Save it anywhere convenient on your computer, for example a folder named `TDDCC` in your home directory or Documents folder. No further files are required to run the application.

### 3.5 Verify the installation

From the folder containing `TDDCC_gui.py`, with your virtual environment activated, run:

```bash
python TDDCC_gui.py        # Windows
python3 TDDCC_gui.py       # macOS / Linux
```

A small window titled "Dynamic Cross-Correlation Analysis" should appear within a few seconds (Figure 2). If it does, installation is complete and you can proceed to [Section 4](#4-launching-tddcc).

> *Figure 2 — The TDDCC intro screen confirms a successful installation.*

> **CAUTION** — If nothing appears and you instead see an error in the terminal, check [Section 8 (Troubleshooting)](#8-troubleshooting) for the most common causes.

## 4. Launching TDDCC

Every time you want to use TDDCC:

1. Open a terminal (Command Prompt, Terminal, or your shell of choice).
2. If you created a virtual environment, activate it ([Section 3.2](#32-create-a-virtual-environment-recommended)).
3. Navigate to the folder containing `TDDCC_gui.py` using `cd`.
4. Run the script, then click "Launch Analysis Dashboard" on the screen that appears to open the main workspace.

```bash
cd path/to/your/TDDCC/folder
python3 TDDCC_gui.py
```

## 5. Interface Tour

Clicking "Launch Analysis Dashboard" opens the main workspace, a single 1400×900 window split into two halves: a scrollable control panel on the left, and a plotting area on the right that stays empty until you run an analysis.

> *Figure 3 — The top of the control panel: Sections 0–3. Scrolling down reveals Run Analysis, View Correlation Dashboard, and Export Data.*

The control panel is organized into numbered sections, from top to bottom:

- **Section 0 — Compute DCC from Trajectory** (optional): runs the built-in compute engine on a trajectory file.
- **Section 1 — Input Data**: the folder of correlation matrices to analyze.
- **Section 2 — Residue Limits**: which residue-pair range to include.
- **Section 3 — Correlation Filters**: the thresholds used to flag a pair as positively or negatively correlated.
- **Run Analysis**, **View Correlation Dashboard**, and **Export Data**, towards the bottom of the panel (scroll down to reach them).

> **NOTE** — The control panel is taller than the window. Use the thin scrollbar at the right edge of the panel, or drag it directly, to reach the lower sections.

## 6. Step-by-Step Workflows

### 6.1 Workflow A — Starting from an MD trajectory

Use this workflow if you have a raw trajectory and topology file and want TDDCC to compute the correlation matrices for you. This requires MDAnalysis to be installed ([Section 3.3](#33-install-required-packages)).

1. In Section 0, click **Browse** next to **Trajectory** and select your trajectory file (`.xtc`, `.dcd`, `.trr`, or any format MDAnalysis supports).
2. Click **Browse** next to **Topology** and select the matching topology file (`.pdb`, `.tpr`, `.gro`, etc.).
3. Click **Browse** next to **Output directory** and choose where the resulting matrices should be written.
4. Set **Window** (the number of frames per sliding window), **Step** (the frame stride used when streaming the trajectory), and **Cores** (how many windows to compute in parallel). The defaults (Window 100, Step 1) are a reasonable starting point for most trajectories.
5. Optionally edit **Plot title**, used as a label on the generated per-window heatmaps.
6. Click **COMPUTE DCC**. Progress is reported in the log box and the progress bar; this can take anywhere from seconds to several minutes depending on trajectory length and Window/Step settings.
7. When finished, TDDCC automatically fills Section 1 (Input Data) with the output directory, so you can go straight to [Section 6.3](#63-setting-residue-limits-and-correlation-filters).

> *Figure 4 — The dashboard, showing DCC calculation for a sample data set.*

> **NOTE** — If MDAnalysis is not installed, Section 0 will report an error when you click COMPUTE DCC. Install it with `pip install MDAnalysis`, or use Workflow B with matrices computed elsewhere.

### 6.2 Workflow B — Starting from existing correlation matrices

Use this workflow if you already have a folder of correlation matrices saved as plain-text files (one file per window, each containing an N×N matrix), whether produced by Workflow A in a previous session or by an external pipeline.

1. In Section 1 (Input Data), click **Browse** and select the folder containing your matrix files.
2. Continue to [Section 6.3](#63-setting-residue-limits-and-correlation-filters) to configure the residue range and filters before running the analysis.

> **TIP** — File names do not need to follow any special convention beyond being plain-text matrices; TDDCC sorts them naturally so that, for example, a file ending in `_2` is read before one ending in `_10`.

### 6.3 Setting residue limits and correlation filters

**Section 2 (Residue Limits)** controls which part of the matrix is analyzed:

- **Residue i range / Residue j range**: the rows and columns to include. Narrowing these speeds up analysis on very large proteins.
- **Min distance |j − i| ≥**: excludes residue pairs closer together than this threshold along the sequence, which are almost always trivially correlated through the backbone and are not usually of biological interest.

**Section 3 (Correlation Filters)** controls which residue pairs are highlighted once the statistics are computed:

- **Positive correlation / Negative correlation**: the range of mean correlation values that count as positively or negatively correlated.
- **Positive corr prob / Negative corr prob**: the range of occurrence probability — the fraction of windows in which a pair showed that sign of correlation — required for a pair to be included. A high probability means the correlation is consistent across the trajectory rather than driven by a few outlier windows.

> *Figure 5 — The populated dashboard, showing setting residue limits and running the analysis for a sample data set.*

> **NOTE** — The default filters (correlation magnitude ≥ 0.2, probability ≥ 0.1) are a reasonable starting point. You can re-enter different values and click Run Analysis again at any time without reloading your data.

### 6.4 Running the analysis and reading the dashboard

1. Scroll down and click **RUN ANALYSIS**. TDDCC loads every matrix in the Input Data folder and computes, for every qualifying residue pair, the mean correlation, standard deviation, and positive/negative occurrence probability across all windows.
2. When the status line reads "Analysis complete!", click **View Correlation Dashboard**.
3. Three heatmaps appear in the right-hand panel: **Positive Correlation**, **Negative Correlation**, and **Combined Correlation**, each showing only the residue pairs that satisfy the filters from Section 3.

> *Figure 6 — The populated dashboard, showing positive, negative, and combined correlation heatmaps for a sample data set.*

Each heatmap supports the standard Matplotlib navigation toolbar beneath it, so you can pan, zoom, and save any panel as an image independently of the rest of the dashboard.

### 6.5 Inspecting a residue pair's time series

To examine how a specific residue pair's correlation behaves across the trajectory rather than just its summary statistic, double-click that cell in any of the three heatmaps. A new window opens showing the pair's correlation value at every individual window plotted against window index, along with its mean and standard deviation. From this window you can save the plot as a PNG or export the underlying values as a CSV.

> *Figure 7 — The populated dashboard, showing the residue pair time series for a sample data set.*

> **TIP** — This view is the best way to confirm that a filtered hit reflects a persistent, biologically meaningful correlation rather than one driven by a handful of unusual windows.

### 6.6 Exporting results

The **Export Data** section provides three buttons:

- **Download Combined CSV** — every residue pair satisfying either the positive or the negative filter.
- **Download Positive CSV** — only pairs satisfying the positive filter.
- **Download Negative CSV** — only pairs satisfying the negative filter.

Each CSV contains one row per residue pair with columns for residue i, residue j, mean correlation, standard deviation, and positive/negative occurrence probability — ready for further analysis in a spreadsheet, R, or Python.

## 7. Output File Reference

| File | Description |
| --- | --- |
| `corr_dcc_<frame>.txt` | Plain-text N×N correlation matrix for one window, named by the index of the window's last frame. Produced by Section 0 (Compute DCC from Trajectory) and consumed by Section 1 (Input Data). |
| `corr_dcc_<frame>.png` | A heatmap image of the same window's matrix, using TDDCC's diverging blue–white–red colormap, saved alongside the text file for quick visual inspection without reopening the dashboard. |
| Combined / Positive / Negative CSV | Exported from the dashboard via the Export Data buttons; one row per qualifying residue pair with mean correlation, standard deviation, and occurrence probability. |
| Time-series PNG / CSV | Exported from the double-click time-series inspector window for a single residue pair. |

## 8. Troubleshooting

| Symptom | Likely cause and fix |
| --- | --- |
| `ModuleNotFoundError: No module named 'tkinter'` | Tkinter is missing from your Python installation. On Linux, install it with your system package manager (e.g. `sudo apt install python3-tk`). On macOS with Homebrew Python, install the matching `python-tk@<version>` formula. See [Section 3.1](#31-install-python). |
| "Invalid input directory" when clicking Run Analysis | Section 1 (Input Data) is empty or points to a folder with no readable correlation matrix files. Browse to a folder containing the `.txt` matrices first. |
| "No .txt files found in selected folder" | The selected Input Data folder does not contain any plain-text matrix files. Confirm you selected the output directory from Workflow A, or the correct external folder in Workflow B. |
| COMPUTE DCC fails immediately or reports MDAnalysis is unavailable | MDAnalysis is not installed. Run `pip install MDAnalysis` in the same environment you use to launch TDDCC, or use Workflow B instead if you already have matrices. |
| The window opens but appears blank or frozen while computing | Large trajectories or a high Window/Step/Cores combination can take time. Progress is reported in the Section 0 log box and the main status line; avoid closing the window while a computation is in progress. |
| Plots or buttons look cut off at the bottom of the screen | The control panel is scrollable. Drag the thin scrollbar at its right edge downward to reveal Run Analysis, View Correlation Dashboard, and Export Data. |
| `pip install` fails with "externally-managed-environment" | Use a virtual environment ([Section 3.2](#32-create-a-virtual-environment-recommended)), or add `--break-system-packages` to the pip command as a quick workaround. |

## 9. Quick Reference

| Feature | What it does |
| --- | --- |
| **Compute DCC from Trajectory** | Streams a trajectory + topology pair through MDAnalysis and computes one DCC matrix per sliding window, in parallel across the chosen number of cores. |
| **Input Data** | Points the dashboard at a folder of correlation matrix files. |
| **Residue Limits** | Restricts the analyzed residue-pair grid and excludes trivially correlated sequence-adjacent neighbors. |
| **Correlation Filters** | Sets the correlation-value and occurrence-probability ranges used to flag positive and negative hits. |
| **Run Analysis** | Loads every matrix and computes mean, standard deviation, and occurrence probability for every qualifying residue pair. |
| **View Correlation Dashboard** | Displays the filtered Positive / Negative / Combined heatmaps; double-click a cell for its time series. |
| **Export Data** | Saves the filtered positive, negative, or combined results to CSV. |

---

*Document version 1.0*
