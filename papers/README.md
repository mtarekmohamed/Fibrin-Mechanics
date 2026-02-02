# Publications – Fibrin Mechanics

This folder contains **open-access published papers** associated with the
**Fibrin-Mechanics** project and the simulation and analysis code provided
in this GitHub repository.

All papers included here are final published versions and are redistributed
in accordance with their respective open-access licenses.

---

## Paper 1 – Fibrin Mechanics Study I

This paper presents molecular simulations of fibrin(ogen) under mechanical
pulling, characterizing force–extension behavior and region-specific
mechanical responses at the molecular level.

**File included**
- Open-access published PDF (final version)

**Associated code**
- Pulling simulations: `pull_openmm.py`
- Relaxation runs: `relax_openmm.py`
- Force decomposition and analysis: `plot_forces.py`
- Viscoelastic modeling: `fit_kelvin_voigt.py`

**Reproducibility**
- The exact code state used for this study is preserved via a GitHub release
  associated with this repository.

---

## Paper 2 – Fibrin Mechanics Study II

This paper extends the fibrin mechanics analysis by incorporating additional
pulling protocols, structural metrics, and viscoelastic interpretation,
linking molecular deformation pathways to emergent mechanical behavior.

**File included**
- Open-access published PDF (final version)

**Associated code**
- Force-by-region analysis across replicas
- Kelvin–Voigt model fitting
- Secondary structure analysis (DSSP)
- Ramachandran-based alpha/beta characterization
- ANM-based stiffness estimation

**Reproducibility**
- The exact code version corresponding to this paper is preserved via a
  dedicated GitHub release.

---

## Relationship to This Repository

- Simulation drivers and analysis scripts located in the repository root
  reproduce the quantitative results reported in these papers.
- Large trajectories, checkpoints, and raw simulation outputs are not included
  due to size constraints and were generated on HPC resources.
- The PDFs provided here ensure long-term accessibility of the published work
  alongside the code.

---

## Citation

If you use the methods, simulations, or analysis workflows from this
repository, please cite the corresponding publication(s) included in
this folder.

---

## Notes

- This repository is intended to support transparency and reproducibility.
- For exact software versions, parameters, and execution details, refer to
  the corresponding GitHub release and script documentation.

