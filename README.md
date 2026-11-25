# Local Contrast Energy Accounts for Early Neural Signatures of Numerosity

Companion manuscript and analysis code for the *Nature Communications* Matters Arising:

> **Harrison, W. J. & Rideaux, R.** (2026). *Early neural signatures of numerosity are accounted for by local contrast energy.* ([Download Manuscript PDF](Matters_Arising_Karami2026.pdf))  
> Arising from: Karami, A., Castaldi, E., Eger, E., Hebart, M. N. & Piazza, M. (2026). *Numerosity is directly sensed and dynamically transformed in the human brain.* *Nature Communications*, https://doi.org/10.1038/s41467-026-79330-x.

## Study Motivation

Karami et al. (2026) reported that visual numerosity emerges rapidly in human MEG responses (~80–140 ms) independently of four non-numeric geometric properties of dot arrays (average item area, total field area, total surface area, and density), concluding that numerosity is "directly sensed" as a primary visual attribute.

This repository tests a potential low-level sensory limitation of that conclusion: any neural estimate of item count must be computed from local contrast energy in early spatial-frequency-selective visual channels. Applying contrast-sensitivity-function (CSF)-limited ideal-observer theory (Harrison, 2019, *i-Perception*, https://doi.org/10.1177/2041669519842895) solely to condition-matched synthetic dot arrays (including within-array dot-size heterogeneity, 100 exemplars per condition) identifies a single optimal spatial-frequency channel at **2.65 cycles per degree (cpd)**. Across all *n* = 30 MEG datasets, this single channel accounts for **92.6%** of early (80–140 ms) *Number* variance and significantly outperforms *Number* (+18.4% relative variance, *p* = 0.0009), while unique *Number* variance emerges only after ~250 ms.

## Required External Data & Original Analysis Files

This repository contains **only the manuscript PDF and the code required for the new contrast-energy modeling and MEG representational similarity analysis (RSA)**. It does not redistribute Karami et al.'s MEG datasets or original analysis scripts. Because trial-level stimulus images and dot-coordinate metadata were not archived in the public MEG datasets, contrast energy is computed from 100 condition-matched synthetic exemplars per condition.

To reproduce the full MEG analyses, download the datasets and reference files directly from the Karami et al. (2026) team:

- **MEG Datasets (`sub-01.7z` – `sub-30.7z`):** https://osf.io/czswn
- **Original Analysis Repository (`NumerosityMEG-FMRI`):** https://github.com/per2con/NumerosityMEG-FMRI

## Repository Layout

- `Matters_Arising_Karami2026.pdf` — Full manuscript and Supplementary Information (including Extended Analyses, Supplementary Figs. S1–S2, and LLM transparency statement).
- `code/` — All MATLAB and Python analysis scripts:
  - `runAll.m` — Master MATLAB pipeline (MEG preprocessing check, ideal-observer SF derivation, RDM simulation, group MEG RSA, and figure generation).
  - `frf_default_params.m` — Stimulus, CSF, 5-band filter-bank, and RDM parameters (relative paths).
  - `derive_optimal_sf.m` — Five-step CSF-limited ideal-observer derivation (yields 2.646 cpd *a priori*).
  - `generate_karami_stimulus.m` — Synthetic dot-array generator matching Karami et al.'s 32-condition grid with within-array dot-size heterogeneity (±50% area jitter constrained to exact condition-mean *Average Item Area*).
  - `build_frf_filters.m` — Isotropic raised-log-cosine spatial-frequency filter bank builder (5 bands, 1.00-octave spacing).
  - `run_frf_model.m` — Single-stage local contrast-energy model (filter, rectify, spatially pool).
  - `run_frf_rdm_simulation.m` — Builds the primary single-SF (2.65 cpd) RDM, control RDMs, and supplementary 5-band RDMs (100 exemplars per condition).
  - `run_meg_vp_test.m` — Sensor-level MEG RSA, non-negative variance partitioning, matched 5-predictor comparison, multi-band SF tuning, and standalone 300-DPI panel export.
  - `make_manuscript_figures.m` — Generates the 7 standalone 300-DPI main-text panels for Figures 1 and 2.
  - `plot_ideal_observer_walkthrough.m` — Generates diagnostic plots of the five-step ideal-observer derivation.
  - `preprocess_meg_to_rdms.py` — Headless MNE-Python MEG preprocessing matching Karami et al. to build 32 × 32 gradiometer RDMs.
  - `filters/` — Self-contained MATLAB filter helpers (`createCSF.m`, `polarDistFun.m`, `createBandPassFilters_customSizes.m`, `createLowPassFilter.m`, `raisedCos.m`).
  - `requirements.txt` — Python dependencies for MEG preprocessing (`mne`, `numpy`, `scipy`).
- `LICENSE` — MIT License.

## How to Reproduce

1. **Stimulus-only analysis (no MEG download required):** Open MATLAB in the `code/` folder and run `run_frf_rdm_simulation` to reproduce the 2.65 cpd ideal-observer derivation and stimulus RDM collinearity (partial *r* = +0.885).
2. **Full MEG re-analysis (*n* = 30):**
   - Download `sub-01.7z` through `sub-30.7z` from https://osf.io/czswn into `code/MEG_Data/`.
   - Install Python dependencies (`pip install -r code/requirements.txt`) and ensure `7zz`/`7z` (or macOS `tar`) is installed.
   - Run `python3 code/preprocess_meg_to_rdms.py --meg-dir ./code/MEG_Data --out-dir ./code/meg_rdms`.
   - Run `runAll` (and `make_manuscript_figures`) in MATLAB from `code/` to generate `code/results/meg_vp_summary.csv`, `code/results/sf_band_summary.csv`, `code/results/panels_main/`, and `code/results/panels_supplement/`.

## Contact

William J. Harrison (`wharrison@usc.edu.au`, ORCID: 0000-0002-6408-0359), Visual Cognition Lab, School of Health, University of the Sunshine Coast.  
