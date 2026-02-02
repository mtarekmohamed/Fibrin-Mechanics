#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Kelvin–Voigt fit with ONE loading file + MANY consecutive unloading files (continuous time axis).

We compute displacement e(t)=d(t)-d0 (nm) using ONE fixed baseline d0=46.1054 nm,
then fit:
  Loading:    de/dt = -λ_load * e + b,     λ_load = (k_spec + k) / c_spec
  Unloading:  de/dt = -λ_unload * e (+ b), λ_unload =  k_spec      / c_spec

Recovered:
  k_spec = k * λ_unload / (λ_load - λ_unload)
  c_spec = k_spec / λ_unload

Inputs (CSV): columns: step, distance_nm  (comma-separated, one header line)
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Optional, Tuple

# ---------------- Utils ----------------

def load_distance_csv(path: str) -> Tuple[np.ndarray, np.ndarray]:
    steps, dist = np.loadtxt(path, delimiter=",", usecols=(0, 1), unpack=True, skiprows=1)
    return steps.astype(float), dist.astype(float)

def gaussian_smooth(y: np.ndarray, window_pts: int = 21, sigma_frac: float = 0.25) -> np.ndarray:
    window_pts = int(max(5, window_pts) | 1)
    half = window_pts // 2
    x = np.arange(-half, half + 1)
    sigma = max(1.0, sigma_frac * window_pts)
    k = np.exp(-(x**2) / (2.0 * sigma**2)); k /= k.sum()
    y_pad = np.pad(y, (half, half), mode="edge")
    return np.convolve(y_pad, k, mode="valid")

def estimate_derivative(t_ps: np.ndarray, y: np.ndarray, smooth_win: Optional[int]) -> Tuple[np.ndarray, np.ndarray]:
    if smooth_win is None:
        smooth_win = max(9, min(101, (len(y)//100)*2 + 1))
    y_sm = gaussian_smooth(y, smooth_win)
    dy_dt = np.gradient(y_sm, t_ps)
    return dy_dt, y_sm

def fit_line(x: np.ndarray, y: np.ndarray, fit_intercept: bool, min_abs_x: float) -> Tuple[float,float,float,float,int]:
    mask = np.isfinite(x) & np.isfinite(y)
    if min_abs_x > 0:
        mask &= (np.abs(x) >= min_abs_x)
    x = x[mask]; y = y[mask]
    if x.size < 8:
        return np.nan, np.nan, np.nan, np.nan, 0
    if fit_intercept:
        X = np.column_stack([x, np.ones_like(x)])
        theta, *_ = np.linalg.lstsq(X, y, rcond=None)
        m, b = float(theta[0]), float(theta[1])
        yhat = X @ theta
        rss = float(np.sum((y - yhat)**2)); n, p = X.shape
        sigma2 = rss / max(1, n - p)
        cov = sigma2 * np.linalg.inv(X.T @ X)
        se = np.sqrt(np.diag(cov))
        m_se, b_se = float(se[0]), float(se[1])
    else:
        X = x.reshape(-1,1)
        theta, *_ = np.linalg.lstsq(X, y, rcond=None)
        m = float(theta[0]); b = 0.0
        yhat = X @ theta
        rss = float(np.sum((y - yhat)**2)); n, p = X.shape
        sigma2 = rss / max(1, n - p)
        cov = sigma2 * np.linalg.inv(X.T @ X)
        m_se = float(np.sqrt(cov[0,0])); b_se = float("nan")
    return m, b, m_se, b_se, x.size

# ------------- Main --------------

def main():
    ap = argparse.ArgumentParser(description="KV fit with one loading file + consecutive unloading (continuous time axis).")
    ap.add_argument("--load_file", required=True, help="CSV step,distance_nm (loading)")
    ap.add_argument("--unload_files", nargs="+", required=True, help="CSV list for unloading (consecutive order).")
    ap.add_argument("--k", type=float, required=True, help="Actuator spring constant [kJ/mol/nm^2]")
    ap.add_argument("--dt_fs", type=float, default=2.0, help="MD step size [fs] per 'step' in CSV (500k pts in 1 ns => 2 fs)")
    ap.add_argument("--smooth_win", type=int, default=None, help="Smoothing window (odd); default auto")
    ap.add_argument("--min_e", type=float, default=0.0, help="Ignore |e| below this (nm) in fits")
    ap.add_argument("--unload_fit_intercept", action="store_true", help="Fit intercept in unloading (handles small drift)")
    ap.add_argument("--out_prefix", default="kv_fit", help="Output prefix")
    ap.add_argument("--out_dir", default=".", help="Output directory")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    dt_ps = args.dt_fs * 1e-3  # fs -> ps

    # -------- Fixed baseline for ALL strain calculations --------
    d0 = 46.1054  # nm (fixed)
    d0_note = "fixed by user"

    # --- LOAD loading data ---
    print(f"[Loading] {args.load_file}")
    stepsL, dL = load_distance_csv(args.load_file)
    tL = stepsL * dt_ps
    eL = dL - d0

    # Fit loading in (e, de/dt)
    dLdt, eL_sm = estimate_derivative(tL, eL, args.smooth_win)
    mL, bL, mL_se, bL_se, nL = fit_line(eL_sm, dLdt, fit_intercept=True, min_abs_x=args.min_e)
    lam_load = -mL

    # --- LOAD & CONCATENATE unloading files (continuous global time right after loading) ---
    tU_chunks, eU_chunks = [], []
    count = 0
    t_start = tL[-1] + dt_ps  # start unloading immediately after loading ends
    for path in args.unload_files:
        print(f"[Unloading] {path}")
        stepsU, dU = load_distance_csv(path)
        n = len(stepsU)
        tU_global = t_start + dt_ps * np.arange(count, count + n)  # continuous
        count += n
        eU = dU - d0  # SAME fixed baseline as loading
        tU_chunks.append(tU_global)
        eU_chunks.append(eU)
    tU = np.concatenate(tU_chunks)
    eU = np.concatenate(eU_chunks)

    # --- Force continuity: start unloading from the last loading displacement ---
    unload_offset = float(eL[-1] - eU[0])
    eU = eU + unload_offset

    # Fit unloading (after shift)
    dUdt, eU_sm = estimate_derivative(tU, eU, args.smooth_win)
    mU, bU, mU_se, bU_se, nU = fit_line(eU_sm, dUdt, fit_intercept=args.unload_fit_intercept, min_abs_x=args.min_e)
    lam_unload = -mU

    # Check sign; if invalid, auto-correct by flipping e sign in unloading fit
    note = ""
    if not np.isfinite(lam_unload) or lam_unload <= 0:
        mU2, bU2, mU2_se, bU2_se, nU2 = fit_line(-eU_sm, dUdt, fit_intercept=args.unload_fit_intercept, min_abs_x=args.min_e)
        lam_unload2 = -mU2
        if np.isfinite(lam_unload2) and lam_unload2 > 0:
            mU, bU, mU_se, bU_se, nU = mU2, bU2, mU2_se, bU2_se, nU2
            lam_unload = lam_unload2
            note = "NOTE: unloading sign auto-corrected (used x=-e)."

    # Recover k_spec & c_spec
    if not np.isfinite(lam_load) or not np.isfinite(lam_unload) or lam_unload <= 0 or lam_load <= lam_unload:
        k_spec = float("nan"); c_spec = float("nan")
        warn = "WARNING: invalid rates (lam_load <= lam_unload or lam_unload <= 0). Adjust --smooth_win/--min_e."
        note = (note + " " + warn).strip()
    else:
        k_spec = args.k * lam_unload / (lam_load - lam_unload)
        c_spec = k_spec / lam_unload

    # --- Report ---
    files_str = f"LOAD: {args.load_file}\nUNLOAD (consecutive):\n  " + "\n  ".join(args.unload_files)
    report = f"""Kelvin–Voigt Fit (loading + consecutive unloading; continuous time)
==============================================================
Input files:
{files_str}
k_act [kJ/mol/nm^2]  = {args.k:.6g}
Time step Δt         = {args.dt_fs:.6g} fs  ({dt_ps:.6g} ps)
Baseline d0 [nm]     = {d0:.6f}  ({d0_note})
Applied unload shift = {unload_offset:.6g} nm

Fits (units: e in nm, t in ps)
------------------------------
Loading:    de/dt = -λ_load * e + b
  λ_load [1/ps]        = {lam_load:.6g}
  slope (load)         = {mL:.6g} ± {mL_se:.3g}
  intercept b          = {bL:.6g} ± {bL_se:.3g}
  N (load)             = {nL}

Unloading:  de/dt = -λ_unload * e{(' + b' if args.unload_fit_intercept else '')}
  λ_unload [1/ps]      = {lam_unload:.6g}
  slope (unload)       = {mU:.6g} ± {mU_se:.3g}
  intercept b          = {bU:.6g} ± {bU_se:.3g}
  N (unload)           = {nU}

Recovered specimen parameters
-----------------------------
k_spec [kJ/mol/nm^2]    = {k_spec:.6g}
c_spec [kJ/mol/nm^2·ps] = {c_spec:.6g}

{note}
"""
    out_txt = os.path.join(args.out_dir, f"{args.out_prefix}_report.txt")
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write(report)
    print(report)

    # --- Plots ---
    # 1) e(t): loading + continuous unloading (one axis, no gaps)
    plt.figure(figsize=(9,5))
    plt.plot(tL, eL, label="Loading", linewidth=2)
    plt.plot(tU, eU, label="Unloading (continuous)", linewidth=1.5, alpha=0.9)
    plt.xlabel("Time [ps]"); plt.ylabel("Displacement e [nm]")
    plt.title("Displacement vs Time (Loading + Continuous Unloading)")
    plt.grid(True, linestyle="--", alpha=0.5); plt.legend(); plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, f"{args.out_prefix}_e_vs_t.png"), dpi=200)

    # 2) de/dt vs e with linear fits
    plt.figure(figsize=(9,5))
    plt.scatter(eL_sm, dLdt, s=4, alpha=0.4, label="Loading")
    xL = np.linspace(np.min(eL_sm), np.max(eL_sm), 200)
    plt.plot(xL, (-lam_load)*xL + bL, linewidth=2)

    plt.scatter(eU_sm, dUdt, s=3, alpha=0.3, label="Unloading")
    xU = np.linspace(np.min(eU_sm), np.max(eU_sm), 200)
    plt.plot(xU, (-lam_unload)*xU + (bU if np.isfinite(bU) else 0.0), linewidth=2)

    plt.xlabel("Displacement e [nm]"); plt.ylabel("de/dt [nm/ps]")
    plt.title("de/dt vs e with linear fits")
    plt.grid(True, linestyle="--", alpha=0.5); plt.legend(); plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, f"{args.out_prefix}_dedt_vs_e.png"), dpi=200)

    # 3) Overlay: model-predicted displacement vs simulation
    # Loading model: e(t) = e_inf + [e(t0) - e_inf] * exp(-lam_load * (t - t0))
    t0 = float(tL[0])
    e0 = float(eL[0])  # initial loading displacement
    e_inf_load = bL / lam_load if np.isfinite(lam_load) and lam_load != 0 else 0.0
    e_model_load = e_inf_load + (e0 - e_inf_load) * np.exp(-lam_load * (tL - t0))

    # Unloading model: e(t) = e_ss + [e(t_r) - e_ss] * exp(-lam_unload * (t - t_r))
    # Start exactly from the last loading value:
    t_r = float(tU[0])
    e_r = float(eL[-1])
    if args.unload_fit_intercept and np.isfinite(bU) and np.isfinite(lam_unload) and lam_unload > 0:
        e_ss_un = bU / lam_unload
    else:
        e_ss_un = 0.0
    e_model_un = e_ss_un + (e_r - e_ss_un) * np.exp(-lam_unload * (tU - t_r))

    plt.figure(figsize=(9,5))
    plt.plot(tL, eL, label="Loading (sim)", linewidth=2)
    plt.plot(tU, eU, label="Unloading (sim)", linewidth=1.5, alpha=0.9)
    plt.plot(tL, e_model_load, "--", label="Loading (KV fit)", linewidth=2)
    plt.plot(tU, e_model_un,  "--", label="Unloading (KV fit)", linewidth=2)
    plt.xlabel("Time [ps]"); plt.ylabel("Displacement e [nm]")
    plt.title("Displacement vs Time: Simulation vs Kelvin–Voigt Model")
    plt.grid(True, linestyle="--", alpha=0.5); plt.legend(); plt.tight_layout()
    plt.savefig(os.path.join(args.out_dir, f"{args.out_prefix}_e_vs_t_with_fit.png"), dpi=200)

    # Save processed CSVs (note: unloading distance includes the applied shift)
    dfL = pd.DataFrame({
        "time_ps": tL,
        "distance_nm": dL,
        "e_nm": eL,
        "e_smooth_nm": eL_sm,
        "de_dt_nm_per_ps": dLdt
    })
    dfU = pd.DataFrame({
        "time_ps": tU,
        "distance_nm": eU + d0,  # includes shift
        "e_nm": eU,
        "e_smooth_nm": eU_sm,
        "de_dt_nm_per_ps": dUdt
    })
    dfL.to_csv(os.path.join(args.out_dir, f"{args.out_prefix}_loading_processed.csv"), index=False)
    dfU.to_csv(os.path.join(args.out_dir, f"{args.out_prefix}_unloading_processed.csv"), index=False)

if __name__ == "__main__":
    main()

