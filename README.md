# LABS on LR-QAOA Benchmarking

This repository applies **Linear-Ramp QAOA (LR-QAOA)** to the **Low Autocorrelation Binary Sequence (LABS)** problem, with the goal of enabling **reproducible benchmarking** across:

- LR-QAOA schedules
- LABS-trained QAOA schedules (e.g., Shaydulin-style schedules)
- classical heuristics

It is also intended to support turning LABS + these approaches into a **Metriq Task** so algorithm developers and hardware teams can compare methods in a consistent format.

---

## Background

This project tests a strong quantum-advantage candidate—the **LABS** problem—where QAOA has shown evidence consistent with favorable scaling behavior in prior work. The approach here is to evaluate LABS using **scalable LR-QAOA** and to structure the workflow so results are straightforward to reproduce and extend.

Relevant reading:
- LABS as a quantum-advantage candidate: `arXiv:2504.03832`, `arXiv:2308.02342`
- LR-QAOA as a scalable protocol: (see Alejandro Montañez Barrera’s LR-QAOA work)
- LABS-trained QAOA schedules: (see Shaydulin’s LABS schedule work)

---

## Repository structure

Top-level layout:

- `Data/`  
  Input instances and cached outputs used by figure notebooks.
- `Figures/`  
  Generated figures (exported plots).

### LABS-focused notebooks

- `generate_problems_LABS.ipynb`  
  Generates / collects LABS problem instances and related inputs.
- `1D-Chain-Experiments_LABS.ipynb`  
  Runs LABS experiments in the 1D-chain setting.
- `1D-Chain-Experiments_LABS_LR-QAOA.ipynb`  
  Runs LABS experiments using LR-QAOA schedules.
- `1D-Chain-Figures_LABS.ipynb`  
  Produces the LABS figures from saved results (and optionally saves them under `Figures/`).

### Other notebooks (broader benchmarking)

The repo also includes non-LABS-specific experiment/figure notebooks (e.g., `1D-Chain-*`, `FC-*`, `NL-*`) that support the broader LR-QAOA benchmarking workflow.

### Utilities

- `get_objective.py`  
  Objective / scoring utilities used by notebooks.
- `LR-QAOA-Benchmark.md`  
  Notes and background on the LR-QAOA benchmarking protocol.
- `requirements.txt`  
  Python dependencies.

---

## Quickstart

### 1) Clone and install dependencies

```bash
git clone https://github.com/makansij/LABS_on_LR-QAOA_Benchmarking.git
cd LABS_on_LR-QAOA_Benchmarking

python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt

