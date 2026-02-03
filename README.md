
# Fibrin-Mechanics

## Pulling Simulations and Mechanical Analysis of Fibrinogen

This repository contains molecular dynamics simulation and analysis scripts developed
to investigate the **mechanical response of fibrinogen and fibrin under external force**.
The workflows focus on connecting **molecular structure, force transmission, and
viscoelastic behavior** across multiple length scales.

The code supports studies of **force–extension behavior, stiffness, secondary-structure
changes, and viscoelastic modeling**, and is designed to accompany published work on
fibrin molecular mechanics.

**Repository:** Fibrin-Mechanics  
**GitHub:** https://github.com/mtarekmohamed/Fibrin-Mechanics

---

## Scientific Motivation

Fibrin is the load-bearing protein network responsible for blood-clot stability.
Its remarkable mechanical properties—high extensibility, nonlinear elasticity,
strain stiffening, and viscoelastic dissipation—emerge from **hierarchical molecular
rearrangements** within fibrinogen and fibrin fibers.

At the molecular level, these properties are governed by:
- deformation of coiled-coil regions,
- force-induced secondary-structure transitions (α-helix ↔ β-sheet),
- redistribution of load across fibrin domains,
- time-dependent relaxation and hysteresis.

This repository enables **atomistic, force-resolved characterization** of these
mechanisms using molecular dynamics simulations and post-processing analyses.

---

## Analysis Scope Implemented Here

The scripts in this repository enable:

- OpenMM-based **pulling simulations** between two atom or centroid groups
- **Relaxation dynamics** following applied strain
- **Force-by-region decomposition**, averaged across independent replicas
- **Kelvin–Voigt viscoelastic fitting** using loading and multiple unloading segments
- **Secondary-structure analysis** (DSSP and Ramachandran-derived metrics)
- **Elastic-network / ANM stiffness estimation** using ProDy

Together, these analyses connect **local molecular deformation mechanisms** to
**effective mechanical response**.

---

## Repository Structure

Large trajectories and raw simulation outputs are intentionally **not included**.

---

## Simulation Workflow Overview

### 1. Molecular Pulling

Pulling simulations are performed using OpenMM by applying controlled mechanical
strain between two predefined atom or centroid groups (via
`CustomCentroidBondForce`). These simulations generate force–extension trajectories
and deformation pathways.

### 2. Relaxation

Relaxation scripts allow equilibration and recovery analysis starting from checkpoints
or final pulled structures, enabling the study of **hysteresis and time-dependent
response**.

### 3. Force Analysis

Forces are decomposed across residues or structural regions. Two independent replicas
are averaged **as vectors prior to magnitude calculation**, ensuring physically
consistent force estimates.

---

## Structural and Mechanical Analysis

- **Secondary Structure (DSSP)**  
  Tracks helix and β-sheet populations during deformation.

- **Ramachandran Analysis**  
  Quantifies backbone conformational shifts associated with mechanical transitions.

- **Viscoelastic Modeling**  
  Fits Kelvin–Voigt models to loading and multiple unloading segments using a continuous
  time axis and fixed baseline distance.

- **ANM Stiffness Estimation**  
  Provides coarse-grained insight into mechanical heterogeneity across fibrin domains.

---

## Expected Inputs and File Layout

### OpenMM Pulling / Relaxation

Typical inputs for `pull_openmm.py`:
- `fibrin_solv_ions.gro` – coordinates
- `topol.top` – topology  
- outputs: `*.dcd`, `*.chk`, and text logs for energies and distance/strain

Typical inputs for `relax_openmm.py`:
- `last_frame_pull1.gro`
- `topol.top`
- checkpoint file (e.g., `relax1_7.chk`)
- outputs: relaxed trajectories, logs, and checkpoints

---

### Force Analysis (Two Replicas)

`plot_forces.py` expects:
- topology: `fibrin_solv_ions.gro`
- trajectory for atom mapping (e.g.,
  `replica1/pull_strain_atom_indices_1000.dcd`)
- force arrays of shape `(n_frames, n_atoms, 3)` for each replica, e.g.:
  - `replica1/pull_strain_atom_indices_1000_forces_kjmolnm.npy`
  - `replica2/pull1_replica2_forces_kjmolnm.npy`

---

## Relation to Recent Literature

The workflows implemented here are consistent with recent experimental and
computational studies that connect **molecular-scale deformation** to **macroscopic
fibrin mechanics**, including:

- **Biophysical Journal (2025)**  
  *Molecular mechanisms underlying fibrin fiber extensibility and mechanical
  heterogeneity*, demonstrating how force-induced molecular rearrangements govern
  nonlinear elasticity and viscoelastic response.  
  DOI: https://doi.org/10.1016/j.bpj.2025.01.032  
  PDF: https://www.cell.com/action/showPdf?pii=S0006-3495%2825%2900752-0

- **Norouzi et al., Soft Matter (2025)**  
  *Loading causes molecular damage in fibrin fibers*, showing that cyclic mechanical
  loading leads to irreversible molecular damage and viscoelastic dissipation,
  combining experiments with molecular simulations.  
  DOI: https://doi.org/10.1039/D5SM00681C

These studies emphasize the importance of **atomistic, force-resolved simulations**
and **time-dependent viscoelastic modeling**, which form the core of this repository.

---



