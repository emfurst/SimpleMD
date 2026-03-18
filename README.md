# SimpleMD

A molecular dynamics code for a Lennard-Jones system.

The purpose is to provide a basic molecular dynamics code base (`CHEG231MD.py`) written in fairly plain Python and standard libraries. The code is used as part of a module on molecular thermodynamics in the University of Delaware Chemical Engineering Thermodynamics I course's honors section. Students use it to explore the basic concepts of molecular simulations and molecular processes in thermodynamics.

*SimpleMD* is not a high performance simulation. Those seeking modern simulation tools should use LAMMPS or other powerful packages.

See the accompanying manual, *Really Simple Molecular Dynamics with Python*, for instructions on getting started and exercises to try with the code.

## Getting started

There are two ways to run SimpleMD: using Google Colab (no installation required) or running locally on your own machine.

### Option A: Google Colab (recommended for getting started quickly)

[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/emfurst/SimpleMD/blob/main/MDsim.ipynb)

Click the badge above to open `MDsim.ipynb` directly in Google Colab. All required packages (`numpy`, `numba`, `matplotlib`) are pre-installed in Colab.

You will need to clone the repository in your first notebook cell so that the simulation code is available:

```python
!git clone https://github.com/emfurst/SimpleMD.git
%cd SimpleMD
```

Then proceed with the notebook as described in the manual.

### Option B: Run locally

You will need the following installed on your computer:

- **Python 3** (3.10 or later recommended) — download from [python.org](https://www.python.org/downloads/)
- **Git** — download from [git-scm.com](https://git-scm.com/downloads)

#### 1. Clone the repository

```bash
git clone https://lem.che.udel.edu/git/furst/SimpleMD.git
cd SimpleMD
```

#### 2. Create a Python virtual environment

Create and activate a virtual environment using Python's built-in `venv` module:

```bash
python3 -m venv .venv
source .venv/bin/activate
```

On Windows, activate with:
```
.venv\Scripts\activate
```

#### 3. Install dependencies

```bash
pip install numpy numba matplotlib jupyter
```

#### 4. Run the simulation

Launch Jupyter and open `MDsim.ipynb`:

```bash
jupyter notebook MDsim.ipynb
```

The notebook walks through creating a simulation, running it, and analyzing the results. See the manual for detailed instructions and exercises.

## Contents

- `CHEG231MD.py` — the molecular dynamics simulation code
- `MDsim.ipynb` — Jupyter notebook for running and analyzing simulations
- `maxwell boltzmann/` — Jupyter notebook for calculating the Maxwell-Boltzmann distribution; students can compare simulation results with the distribution
- `SimpleMD_manual.pdf` — the manual, *Really Simple Molecular Dynamics with Python*
