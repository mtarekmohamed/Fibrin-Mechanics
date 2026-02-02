#!/usr/bin/env python3
"""
Force-by-region analysis over ALL frames using the AVERAGE of TWO replicas.

Inputs (forces):
  - pull_strain_atom_indicies_1000_forces_kjmolnm.npy   (replica 1)
  - pull1_replica2_forces.npy                           (replica 2)

Each forces array must be shaped (n_frames, n_atoms, 3).
Both are converted to nN, then VECTOR-AVERAGED across replicas BEFORE magnitudes.

Outputs:
  - forces_pull1_replica_avg_peratom.csv                  # per-atom forces for all frames (from averaged replica forces)
  - forces_pull1_replica_avg_region_timeseries.csv        # per-frame region mean/std
  - forces_pull1_replica_avg_summary.csv                  # time-averaged region mean/std
  - region_residues_map.csv                               # which residues belong to each region
  - force_allframes_violin_avg.svg                        # violin (all frames × atoms) per region
  - force_region_timeseries_avg.svg                       # per-frame mean force per region
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import MDAnalysis as mda

# ---------------- Config ----------------
TOP_FILE = "fibrin_solv_ions.gro"
DCD_FILE = "replica1/pull_strain_atom_indicies_1000.dcd"

# Two replicas (path, units). If units differ, set individually.
FORCE_FILES = [
    ("replica1/pull_strain_atom_indicies_1000_forces_kjmolnm.npy", "kJ/mol/nm"),  # replica 1
    ("replica2/pull1_replica2_forces_kjmolnm.npy",                           "kJ/mol/nm"),  # replica 2
]

# Interpret region ranges using:
#   - 'resindex' : 0-based contiguous residue ordinal (recommended)
#   - 'resid'    : PDB/topology-style residue ID (often 1-based, non-contiguous)
INDEX_MODE = "resindex"

# Region definitions (ranges are [start, stop), i.e., stop is EXCLUSIVE)
REGION_SLICES = {
    "Gamma1":  (696, 956),
    "Gamma2":  (1652, 1913),
    "Beta1":   (314, 575),
    "Beta2":   (1270, 1531),
    "ABC":     [(20, 174), (193, 314), (583, 696)],
    "DEF":     [(976, 1130), (1149, 1270), (1539, 1652)],
    "Central": [(0, 20), (174, 193), (575, 583), (956, 976), (1130, 1149), (1531, 1539)],
}
ORDERED_REGIONS = ['Gamma1', 'Beta1', 'ABC', 'Central', 'DEF', 'Beta2', 'Gamma2']

# Output
CSV_PER_ATOM_ALL        = "forces_pull1_replica_avg_peratom.csv"
CSV_REGION_TIMESERIES   = "forces_pull1_replica_avg_region_timeseries.csv"
CSV_SUMMARY_TIMEAVG     = "forces_pull1_replica_avg_summary.csv"
CSV_REGION_RESIDUES     = "region_residues_map.csv"
PLOT_VIOLIN_SVG         = "force_allframes_violin_avg.svg"
PLOT_REGION_TS_SVG      = "force_region_timeseries_avg.svg"

# ---------------- Helpers ----------------
def load_forces(path):
    """Load forces array of shape (n_frames, n_atoms, 3) from .npy or .npz."""
    if path.endswith(".npy"):
        F = np.load(path)
    elif path.endswith(".npz"):
        data = np.load(path)
        if "forces" in data.files:
            F = data["forces"]
        else:
            key = list(data.keys())[0]
            F = data[key]
    else:
        raise ValueError(f"Unsupported forces format for {path} (use .npy or .npz)")
    if F.ndim != 3 or F.shape[2] != 3:
        raise ValueError(f"{path}: expected (n_frames, n_atoms, 3), got {F.shape}")
    return F

def to_nN_multiplier(units):
    """Return multiplier to convert given units to nN."""
    u = units.strip().lower().replace(" ", "")
    if u in {"kj/mol/nm", "kjmol/nm", "kjmolnm", "kj/molnm"}:
        # 1 kJ/mol/nm ≈ 1.6605390666 pN = 1.6605390666e-3 nN
        return 1.6605390666e-3
    if u in {"pn"}:
        return 1e-3  # 1 pN = 1e-3 nN
    if u in {"nn"}:
        return 1.0
    raise ValueError(f"Unknown FORCE_UNITS: {units}")

def _normalize_ranges(ranges):
    if isinstance(ranges, tuple):
        return [ranges]
    return list(ranges)

def residues_for_ranges(universe, ranges, index_mode="resindex"):
    """Return list of Residue objects covered by ranges."""
    ranges = _normalize_ranges(ranges)
    if index_mode == "resindex":
        residues = []
        for start, stop in ranges:
            start = max(0, start)
            stop  = min(stop, universe.residues.n_residues)
            residues.extend(list(universe.residues[start:stop]))
        return residues
    elif index_mode == "resid":
        valid_resids = set()
        for start, stop in ranges:
            valid_resids.update(range(start, stop))
        return [r for r in universe.residues if r.resid in valid_resids]
    else:
        raise ValueError("index_mode must be 'resindex' or 'resid'")

def build_region_maps(universe, region_slices, index_mode="resindex"):
    """Return (region_atom_indices, region_residue_info) and print residues."""
    region_atom_indices = {}
    region_residue_info = {}

    print("\n=== Building region → atoms map ===")
    print(f"INDEX_MODE = {index_mode}")
    for region, rngs in region_slices.items():
        residues = residues_for_ranges(universe, rngs, index_mode=index_mode)
        atom_idx = []
        res_info = []
        print(f"\nRegion: {region}")
        if not residues:
            print("  (no residues matched)")
        for r in residues:
            res_index_print = r.resindex if index_mode == "resindex" else r.resid
            atom_ix = r.atoms.ix
            atom_idx.extend(atom_ix.tolist())
            res_info.append({
                "Region": region,
                "ResidueName": r.resname,
                "ResidueIndex": int(res_index_print),
                "ResidueIndexMode": index_mode,
                "N_atoms": int(len(atom_ix)),
                "FirstAtomIndex": int(atom_ix[0]) if len(atom_ix) else -1,
            })
            print(f"  Residue: {r.resname}  {res_index_print}  "
                  f"(atoms: {len(atom_ix)}, first atom ix: "
                  f"{int(atom_ix[0]) if len(atom_ix) else 'NA'})")
        region_atom_indices[region] = np.array(sorted(set(atom_idx)), dtype=int)
        region_residue_info[region] = res_info
        print(f"  -> total unique atoms in region '{region}': {region_atom_indices[region].size}")
    return region_atom_indices, region_residue_info

# ---------------- Load traj ----------------
print("Loading trajectory...")
u = mda.Universe(TOP_FILE, DCD_FILE)
n_atoms = u.atoms.n_atoms
n_frames = len(u.trajectory)
if n_frames < 1:
    raise RuntimeError("No frames found in DCD.")

# ---------------- Load & average replicas ----------------
print("Loading forces from replicas and converting to nN...")
replicas = []
for path, units in FORCE_FILES:
    F = load_forces(path)                        # (T, N, 3)
    if F.shape[0] != n_frames or F.shape[1] != n_atoms:
        raise ValueError(
            f"Shape mismatch for {path}: forces {F.shape} vs traj (frames={n_frames}, atoms={n_atoms})"
        )
    mult = to_nN_multiplier(units)
    replicas.append(F * mult)                    # convert to nN

# Stack -> (R, T, N, 3); average across replicas -> (T, N, 3)
F_nN_avg_vec = np.mean(np.stack(replicas, axis=0), axis=0)

# Magnitudes from averaged vectors -> (T, N)
Fmag = np.linalg.norm(F_nN_avg_vec, axis=2)

# ---------------- Build region maps & residue CSV ----------------
region_atom_indices, region_residue_info = build_region_maps(
    u, REGION_SLICES, index_mode=INDEX_MODE
)
all_res_rows = []
for region, lst in region_residue_info.items():
    all_res_rows.extend(lst)
if all_res_rows:
    pd.DataFrame(all_res_rows).to_csv(CSV_REGION_RESIDUES, index=False)
    print(f"\nWrote residue map: {CSV_REGION_RESIDUES}")

# ---------------- All-frames analysis on averaged forces ----------------
print("\nAnalyzing ALL frames on averaged forces...")

peratom_rows = []
ts_rows = []

for frame_i in range(n_frames):
    F_frame = Fmag[frame_i, :]  # (N,)
    for region in ORDERED_REGIONS:
        idx = region_atom_indices.get(region, np.array([], dtype=int))
        if idx.size == 0:
            ts_rows.append({
                "Frame": frame_i, "Region": region,
                "MeanForce(nN)": np.nan, "StdDev(nN)": 0.0, "N_atoms": 0
            })
            continue
        valid = idx[(idx >= 0) & (idx < n_atoms)]
        vals = F_frame[valid]
        # per-atom rows
        peratom_rows.extend(
            {"Frame": frame_i, "Region": region, "Force_nN": float(v)} for v in vals
        )
        # per-frame region stats
        ts_rows.append({
            "Frame": frame_i, "Region": region,
            "MeanForce(nN)": float(vals.mean()) if vals.size else np.nan,
            "StdDev(nN)": float(vals.std(ddof=1)) if vals.size > 1 else 0.0,
            "N_atoms": int(valid.size)
        })

# Save per-atom all-frames CSV
df_peratom = pd.DataFrame(peratom_rows)
df_peratom.to_csv(CSV_PER_ATOM_ALL, index=False)
print(f"Wrote: {CSV_PER_ATOM_ALL}  (rows: {len(df_peratom)})")

# Save per-frame region time series CSV
df_ts = pd.DataFrame(ts_rows)
df_ts.to_csv(CSV_REGION_TIMESERIES, index=False)
print(f"Wrote: {CSV_REGION_TIMESERIES}  (rows: {len(df_ts)})")

# Time-averaged region summary across all frames × atoms
summary = {"Region": [], "MeanForce(nN)": [], "StdDev(nN)": [], "N_atoms": []}
for region in ORDERED_REGIONS:
    idx = region_atom_indices.get(region, np.array([], dtype=int))
    idx = idx[(idx >= 0) & (idx < n_atoms)]
    if idx.size == 0:
        summary["Region"].append(region)
        summary["MeanForce(nN)"].append(np.nan)
        summary["StdDev(nN)"].append(0.0)
        summary["N_atoms"].append(0)
        continue
    vals = Fmag[:, idx].ravel()
    summary["Region"].append(region)
    summary["MeanForce(nN)"].append(float(vals.mean()))
    summary["StdDev(nN)"].append(float(vals.std(ddof=1)))
    summary["N_atoms"].append(int(idx.size))

df_summary = pd.DataFrame(summary)
df_summary.to_csv(CSV_SUMMARY_TIMEAVG, index=False)
print(f"Wrote: {CSV_SUMMARY_TIMEAVG}")

# ---------------- Plots ----------------
print("Making violin plot (all frames × atoms, averaged replicas)...")
df_violin = df_peratom.rename(columns={"Force_nN": "PerAtomForce_nN"})
plt.figure(figsize=(11, 6))
sns.violinplot(
    data=df_violin,
    x="Region",
    y="PerAtomForce_nN",
    order=ORDERED_REGIONS,
    hue="Region",
    palette="Set2",
    legend=False,
    cut=0,
    inner="box"
)
plt.ylabel("Per-atom |Force| (nN) across frames (replica-averaged)")
plt.title("Forces by Region — all frames × atoms (averaged replicas)")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(PLOT_VIOLIN_SVG, dpi=300)
plt.close()
print(f"Wrote: {PLOT_VIOLIN_SVG}")

print("Making per-region mean-force time series plot (averaged replicas)...")
plt.figure(figsize=(11, 6))
for region in ORDERED_REGIONS:
    sub = df_ts[df_ts["Region"] == region]
    if not sub.empty:
        plt.plot(sub["Frame"].to_numpy(), sub["MeanForce(nN)"].to_numpy(), label=region)
plt.xlabel("Frame")
plt.ylabel("Region mean |Force| (nN)")
plt.title("Per-region mean force over time (averaged replicas)")
plt.legend(ncol=2, frameon=False)
plt.tight_layout()
plt.savefig(PLOT_REGION_TS_SVG, dpi=300)
plt.close()
print(f"Wrote: {PLOT_REGION_TS_SVG}")

print("Done.")
# ===================== Statistical tests across regions =====================
import itertools
from scipy import stats

OMNIBUS_CSV  = "stats_region_omnibus.csv"
PAIRWISE_CSV = "stats_region_pairwise.csv"

# Optional: to speed up with huge data, set a per-group cap (e.g., 50000)
MAX_PER_GROUP = None  # or an integer

# Gather data by region
grouped = {reg: df_peratom.loc[df_peratom["Region"] == reg, "Force_nN"].to_numpy()
           for reg in ORDERED_REGIONS}

# Optional subsampling (without replacement) for speed
if MAX_PER_GROUP is not None:
    rng = np.random.default_rng(123)
    for reg, arr in grouped.items():
        if arr.size > MAX_PER_GROUP:
            grouped[reg] = rng.choice(arr, size=MAX_PER_GROUP, replace=False)

groups_in_order = [grouped[reg] for reg in ORDERED_REGIONS]

# --- Omnibus tests ---
# Kruskal–Wallis (non-parametric)
H, p_kw = stats.kruskal(*groups_in_order)

# One-way ANOVA (parametric; assumes normality/equal variances)
F, p_anova = stats.f_oneway(*groups_in_order)

# Effect sizes
def epsilon_squared_kruskal(H, k, N):
    # ε² = (H - k + 1) / (N - k), where k = number of groups, N = total N
    return float((H - k + 1) / (N - k)) if (N - k) > 0 else np.nan

def eta_squared_anova(F, k, N):
    # η² ≈ (k-1)*F / ((k-1)*F + (N-k))
    num = (k - 1) * F
    den = num + (N - k)
    return float(num / den) if den > 0 else np.nan

k = len(ORDERED_REGIONS)
N = int(sum(arr.size for arr in groups_in_order))
eps2 = epsilon_squared_kruskal(H, k, N)
eta2 = eta_squared_anova(F, k, N)

df_omnibus = pd.DataFrame([{
    "Test": "Kruskal–Wallis",
    "Statistic": H,
    "p_value": p_kw,
    "EffectSize": eps2,
    "EffectSizeName": "epsilon^2",
    "k_groups": k,
    "N_total": N
}, {
    "Test": "One-way ANOVA",
    "Statistic": F,
    "p_value": p_anova,
    "EffectSize": eta2,
    "EffectSizeName": "eta^2",
    "k_groups": k,
    "N_total": N
}])
df_omnibus.to_csv(OMNIBUS_CSV, index=False)
print(f"Wrote: {OMNIBUS_CSV}")

# --- Pairwise tests (Mann–Whitney U with Holm–Bonferroni correction) ---
def holm_correction(pvals):
    """
    Holm–Bonferroni step-down adjustment.
    Returns array of adjusted p-values in the original order.
    """
    m = len(pvals)
    order = np.argsort(pvals)
    adj = np.empty(m, dtype=float)
    prev = 0.0
    for i, idx in enumerate(order):
        rank = i + 1
        adj_p = (m - i) * pvals[idx]
        adj_p = max(adj_p, prev)  # ensure monotonicity
        adj[idx] = min(adj_p, 1.0)
        prev = adj[idx]
    # re-check monotonicity increasing with i in original order
    return adj

pairs = []
p_raw = []

for a, b in itertools.combinations(ORDERED_REGIONS, 2):
    x = grouped[a]
    y = grouped[b]
    # Use two-sided Mann–Whitney; tie correction handled internally
    U, p = stats.mannwhitneyu(x, y, alternative="two-sided", method="asymptotic")
    pairs.append((a, b, U))
    p_raw.append(p)

p_adj = holm_correction(np.array(p_raw, dtype=float))

rows = []
for (a, b, U), p, padj in zip(pairs, p_raw, p_adj):
    # Cliff's delta effect size for ordinal data
    # Efficient approximate via rank-biserial correlation:
    nx, ny = grouped[a].size, grouped[b].size
    # rank-biserial r = 1 - 2*U/(nx*ny); Cliff's delta ≈ r
    rbc = 1.0 - 2.0 * (U / (nx * ny))
    rows.append({
        "RegionA": a, "RegionB": b,
        "U_stat": U,
        "p_value_raw": p,
        "p_value_holm": padj,
        "EffectSize": rbc,
        "EffectSizeName": "CliffsDelta≈RankBiserial",
        "n_A": int(nx), "n_B": int(ny)
    })

df_pairwise = pd.DataFrame(rows).sort_values(["p_value_holm", "RegionA", "RegionB"])
df_pairwise.to_csv(PAIRWISE_CSV, index=False)
print(f"Wrote: {PAIRWISE_CSV}")

# Optional: print a quick verdict
alpha = 0.05
sig_any = (df_pairwise["p_value_holm"] < alpha).any()
if p_kw < alpha:
    print(f"Omnibus Kruskal–Wallis: p={p_kw:.3e} (ε²={eps2:.3f}) → at least one region differs.")
else:
    print(f"Omnibus Kruskal–Wallis: p={p_kw:.3e} (ε²={eps2:.3f}) → no evidence of differences.")

if sig_any:
    n_sig = int((df_pairwise['p_value_holm'] < alpha).sum())
    print(f"{n_sig} pairwise differences remain significant after Holm–Bonferroni (α={alpha}).")
else:
    print(f"No pairwise differences are significant after Holm–Bonferroni (α={alpha}).")
# ==========================================================================#

