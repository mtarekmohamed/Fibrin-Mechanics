Pulling simulations and mechanical analysis of fibrinogen
# Fibrin-Mechanics

Simulation and analysis scripts for fibrin(ogen) mechanics workflows, including:
- OpenMM-based pulling simulations between two atom groups
- Relaxation runs from checkpoints / last pulled structure
- Force-by-region post-processing (two replicas averaged)
- Kelvin–Voigt fitting of loading + multiple unloading segments
- Secondary-structure (DSSP) and Ramachandran-derived metrics
- ANM/ENM stiffness estimates (ProDy)

Repository: Fibrin-Mechanics (GitHub: mtarekmohamed)
---

## Contents

### Simulation
- **`pull_openmm.py`** – pulling simulation using OpenMM `CustomCentroidBondForce` between two centroid groups; logs distance/strain and writes trajectory + checkpoint. :contentReference[oaicite:0]{index=0}  
- **`relax_openmm.py`** – relaxation run that loads from a checkpoint and continues dynamics; also logs distance/strain. :contentReference[oaicite:1]{index=1}  

### Analysis
- **`plot_forces.py`** – force-by-region analysis across **all frames**, using **two replicas averaged as vectors** before taking magnitudes; outputs CSVs + plots + stats. :contentReference[oaicite:2]{index=2}  
- **`fit_kelvin_voigt.py`** – Kelvin–Voigt fit using **one loading file** and **multiple consecutive unloading files**, using a continuous time axis and a fixed baseline distance. :contentReference[oaicite:3]{index=3}  
- **`dssp.py`** – computes DSSP (secondary structure) from a trajectory and makes a simple alpha-helix vs beta-sheet count plot. :contentReference[oaicite:4]{index=4}  
- **`Ramchandranplot_box.py`** – computes per-frame alpha/beta classification from Ramachandran (phi/psi) angles and saves a boxplot + numpy array. :contentReference[oaicite:5]{index=5}  
- **`stiff_ANM.py`** – ProDy ANM stiffness estimate and fluctuations for last-frame PDBs (two examples). :contentReference[oaicite:6]{index=6}  

---

## Expected inputs / file layout

These scripts assume a working directory containing (or pointing to) the following kinds of files:

### OpenMM pulling / relaxation
Typical files used by `pull_openmm.py`:
- `fibrin_solv_ions.gro` (coordinates)
- `topol.top` (topology)
- outputs: `*.dcd`, `*.chk`, and text logs for energies + distance/strain :contentReference[oaicite:7]{index=7}  

Typical files used by `relax_openmm.py`:
- `last_frame_pull1.gro`
- `topol.top`
- a checkpoint such as `relax1_7.chk` (example in the script)
- outputs: `relax1_8.dcd`, `relax1_8.log`, `relax1_8.chk` :contentReference[oaicite:8]{index=8}  

### Force analysis (two replicas)
`plot_forces.py` expects:
- topology: `fibrin_solv_ions.gro`
- trajectory for mapping: e.g. `replica1/pull_strain_atom_indicies_1000.dcd`
- forces arrays shaped `(n_frames, n_atoms, 3)` for each replica, e.g.
  - `replica1/pull_strain_atom_indicies_1000_forces_kjmolnm.npy`
  - `replica2/pull1_replica2_forces_kjmolnm.npy` :contentReference[oaicite:9]{index=9}  
