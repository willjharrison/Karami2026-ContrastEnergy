#!/usr/bin/env python3
"""
Fast headless MEG preprocessing + 32x32 gradiometer RDM builder matching
Karami et al. (2026, Nature Comms; NumerosityMEG-FMRI/Analysis/meg_mne_preprocess.ipynb).

Strips out interactive plotting, PDF saving, unused muscle annotations, and
redundant .copy() calls, and runs filtering with n_jobs=-1.

Usage (in your pyenv Terminal):
    pip install mne scipy numpy
    python3 preprocess_meg_to_rdms.py --subjects sub-01 sub-02

By default reads archives from:
    .../Numerosity versus textures/MEG Data/sub-XX.7z
and saves compact ~2 MB RDM files to:
    .../Numerosity versus textures/FRF_model/meg_rdms/sub-XX_rdm_grad.mat
"""

import argparse
import glob
import os
import shutil
import subprocess
import tempfile
import time
import numpy as np
import scipy.io
import scipy.signal
import mne

mne.set_log_level("WARNING")

DEFAULT_BASE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_MEG_DIR = os.path.join(DEFAULT_BASE, "MEG_Data")
DEFAULT_OUT_DIR = os.path.join(DEFAULT_BASE, "meg_rdms")

# Exact Karami et al. parameters from meg_mne_preprocess.ipynb
TMIN = -0.1
TMAX = 0.8
BASELINE = (-0.1, 0.0)
DESIRED_SFREQ = 200.0
POWERLINE_FREQS = (50, 100, 150, 200)
HP_FREQ = 0.05
LP_FREQ = 330.0
SAVGOL_WINDOW = 25
SAVGOL_POLY = 4

EVENT_SAMPLE_LIST = [
    f"N{n}S{s}TFA{f}"
    for n in range(1, 5)
    for s in range(1, 5)
    for f in range(1, 3)
]
EVENT_SAMPLE_DICT = {lbl: idx + 1 for idx, lbl in enumerate(EVENT_SAMPLE_LIST)}


def extract_archive(archive_path, tmp_dir):
    """Extract .7z archive using 7zz/7z (if available) or macOS bsdtar."""
    seven_z = shutil.which("7zz") or shutil.which("7z")
    if seven_z:
        cmd = [seven_z, "x", "-y", f"-o{tmp_dir}", archive_path]
    else:
        cmd = ["tar", "-xf", archive_path, "-C", tmp_dir]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)


def compute_pearson_rdms(cond_avg):
    """
    cond_avg: (32, n_sensors, n_times)
    Returns:  (32, 32, n_times) matrix of 1 - Pearson correlation.
    """
    n_cond, _, n_times = cond_avg.shape
    rdms = np.zeros((n_cond, n_cond, n_times), dtype=np.float32)
    for t in range(n_times):
        X = cond_avg[:, :, t]
        # Row-wise Pearson correlation (matching cosmo_pdist(..., 'correlation'))
        Xc = X - X.mean(axis=1, keepdims=True)
        norms = np.linalg.norm(Xc, axis=1, keepdims=True)
        norms[norms == 0] = 1.0
        Xn = Xc / norms
        R = Xn @ Xn.T
        D = 1.0 - np.clip(R, -1.0, 1.0)
        np.fill_diagonal(D, 0.0)
        rdms[:, :, t] = D
    return rdms


def process_subject(subj, meg_dir, out_dir, overwrite=False):
    out_path = os.path.join(out_dir, f"{subj}_rdm_grad.mat")
    if os.path.exists(out_path) and not overwrite:
        print(f"[{subj}] Skipping (already exists: {out_path})")
        return

    archive_path = os.path.join(meg_dir, f"{subj}.7z")
    if not os.path.exists(archive_path):
        print(f"[{subj}] Archive not found: {archive_path}")
        return

    t0 = time.time()
    print(f"\n[{subj}] Extracting {os.path.basename(archive_path)}...")

    with tempfile.TemporaryDirectory(prefix=f"meg_{subj}_") as tmp_dir:
        extract_archive(archive_path, tmp_dir)
        fif_files = sorted(glob.glob(os.path.join(tmp_dir, "**", "*.fif"), recursive=True))
        if not fif_files:
            raise RuntimeError(f"No .fif files found inside {archive_path}")

        print(f"[{subj}] Found {len(fif_files)} runs ({time.time() - t0:.1f}s). Checking header...")
        raw0 = mne.io.read_raw_fif(fif_files[0], preload=False, verbose="ERROR")
        proc_hist = [str(p.get("short_name", "")) for p in raw0.info.get("proc_history", [])]
        print(
            f"[{subj}] sfreq={raw0.info['sfreq']} Hz | "
            f"n_grad={len(mne.pick_types(raw0.info, meg='grad'))} | "
            f"bads={raw0.info['bads']} | proc_history={proc_hist}"
        )

        # Process each run in-place and epoch immediately to keep RAM low
        all_epochs = []
        for r_idx, fif in enumerate(fif_files, start=1):
            raw = mne.io.read_raw_fif(fif, preload=True, verbose="ERROR")
            meg_picks = mne.pick_types(raw.info, meg="grad", exclude=[])
            stim_picks = mne.pick_types(raw.info, stim=True)
            keep_picks = np.sort(np.concatenate([meg_picks, stim_picks]))
            raw.pick(keep_picks)

            grad_picks = mne.pick_types(raw.info, meg="grad", exclude=[])

            # 1. Notch filter (50 Hz & harmonics)
            raw.notch_filter(freqs=POWERLINE_FREQS, picks=grad_picks, n_jobs=-1, verbose="ERROR")
            # 2. Band-pass filter (0.05 - 330 Hz)
            raw.filter(HP_FREQ, LP_FREQ, fir_design="firwin", picks=grad_picks, n_jobs=-1, verbose="ERROR")
            # 3. Savitzky-Golay smoothing (window=25, polyorder=4) vectorized over all grad channels
            raw._data[grad_picks] = scipy.signal.savgol_filter(
                raw._data[grad_picks], SAVGOL_WINDOW, SAVGOL_POLY, axis=-1
            )
            # 4. Anti-aliasing low-pass before decimation (matches notebook's obtained_sfreq / 3)
            sfreq = raw.info["sfreq"]
            decim = int(np.round(sfreq / DESIRED_SFREQ))
            aa_lpf = (sfreq / decim) / 3.0
            raw.filter(l_freq=None, h_freq=aa_lpf, picks=grad_picks, n_jobs=-1, verbose="ERROR")

            # 5. Events & epoching
            events = mne.find_events(
                raw, stim_channel="STI101", min_duration=0.005, shortest_event=1,
                output="onset", verbose="ERROR"
            )
            # Keep only sample events (1..32) present in this run
            present_ids = {k: v for k, v in EVENT_SAMPLE_DICT.items() if v in events[:, 2]}
            if not present_ids:
                print(f"  [{subj}] Warning: run {r_idx} has no sample triggers 1..32; skipping.")
                continue

            ep = mne.Epochs(
                raw, events, event_id=present_ids, tmin=TMIN, tmax=TMAX,
                baseline=BASELINE, decim=decim, picks=grad_picks,
                preload=True, verbose="ERROR"
            )
            all_epochs.append(ep)
            if r_idx % 4 == 0 or r_idx == len(fif_files):
                print(f"  [{subj}] Run {r_idx:2d}/{len(fif_files)} done [{time.time() - t0:.1f}s]")

        # Concatenate across runs and compute condition-averaged gradiometer patterns
        epochs_all = mne.concatenate_epochs(all_epochs, verbose="ERROR")
        epochs_all.equalize_event_counts(EVENT_SAMPLE_LIST)

        n_times = len(epochs_all.times)
        n_grad = len(epochs_all.ch_names)
        cond_avg = np.zeros((32, n_grad, n_times), dtype=np.float64)
        trial_counts = np.zeros(32, dtype=int)

        for c_idx, lbl in enumerate(EVENT_SAMPLE_LIST):
            dat = epochs_all[lbl].get_data()  # (n_trials, n_grad, n_times)
            trial_counts[c_idx] = dat.shape[0]
            cond_avg[c_idx] = dat.mean(axis=0)

        rdm_grad = compute_pearson_rdms(cond_avg)

        os.makedirs(out_dir, exist_ok=True)
        scipy.io.savemat(
            out_path,
            {
                "rdm_grad": rdm_grad,                 # 32 x 32 x n_times
                "times_ms": epochs_all.times * 1000.0,
                "trial_counts": trial_counts,
                "n_grad": n_grad,
                "subject": subj,
                "proc_history": proc_hist,
            },
        )
        print(
            f"[{subj}] Saved {out_path} | shape={rdm_grad.shape} | "
            f"trials/cond={trial_counts.min()}..{trial_counts.max()} | "
            f"total={time.time() - t0:.1f}s"
        )


def main():
    parser = argparse.ArgumentParser(description="Fast headless MEG -> 32x32 gradiometer RDM builder.")
    parser.add_argument("--meg-dir", default=DEFAULT_MEG_DIR, help="Directory containing sub-XX.7z archives")
    parser.add_argument("--out-dir", default=DEFAULT_OUT_DIR, help="Output directory for sub-XX_rdm_grad.mat")
    parser.add_argument("--subjects", nargs="*", default=None, help="Specific subjects (e.g., sub-01 sub-02). Default: all .7z in meg-dir.")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing RDM .mat files")
    args = parser.parse_args()

    if args.subjects:
        subjects = args.subjects
    else:
        min_bytes = 1_000_000_000  # strictly complete archives (>= 1 GB)
        archives = sorted(glob.glob(os.path.join(args.meg_dir, "sub-*.7z")))
        complete = [a for a in archives if os.path.getsize(a) >= min_bytes]
        in_progress = sorted(glob.glob(os.path.join(args.meg_dir, "sub-*.7z.part")))
        if in_progress:
            print(f"Ignoring {len(in_progress)} in-progress .part download(s): {[os.path.basename(p) for p in in_progress]}")
        subjects = [os.path.splitext(os.path.basename(a))[0] for a in complete]

    print(f"Found {len(subjects)} complete subject archive(s): {subjects}")
    for subj in subjects:
        process_subject(subj, args.meg_dir, args.out_dir, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
