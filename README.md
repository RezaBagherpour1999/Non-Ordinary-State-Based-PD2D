# NOSB-PD2D
## A 2D MATLAB Implementation of Non-Ordinary State-Based Peridynamics

**NOSB-PD2D** is an extended, two-dimensional peridynamics computational code developed on the foundation of [ORNL/PDMATLAB2D](https://github.com/ORNL/PDMATLAB2D). This code expands the original framework by implementing the Non-Ordinary State-Based (NOSB) formulation, incorporating stress-based fracture criteria, damping for quasi-static problems, and advanced nodal tracking capabilities.

---

## Key Features & Modifications

Compared to the reference `PDMATLAB2D` code, the primary modifications implemented in this repository include:

1. **Non-Ordinary State-Based (NOSB) Formulation:**
   - The force and strain energy density computations (`Forceenergydensity`) have been reformulated to implement the **Non-Ordinary State-Based (NOSB)** peridynamic theory, enabling the use of classical continuum constitutive models within a nonlocal meshfree framework.

2. **Maximum Tensile Stress Failure Criterion:**
   - The bond breaking logic (`bondbreaking`) has been modified to evaluate bond failure based on a **maximum tensile stress criterion** instead of the conventional critical stretch condition, providing physically consistent crack initiation and propagation analysis.

3. **Quasi-Static Analysis Support:**
   - Input files and solver routines have been upgraded to include **damping coefficients (damper variables)**, allowing dynamic relaxation to solve quasi-static problems efficiently.

4. **Nodal History Tracking & Plotting:**
   - The main solver script has been augmented to track, record, and plot response histories (such as kinetic energy density, strain energy density, displacements, and force densities) for selected target nodes across time steps.

---

## How to Run

Simulations are executed by passing the input file name to the top-level driver function:
```matlab
NOSBPD2D('CrackBranching')

- Individual component tests are located in the `Tests/` directory.
- Plotting scripts and post-processing tools are available in the `PlottingExamples/` directory.

---

## Citing PDMATLAB2D

This repository is derived from and builds upon the open-source **PDMATLAB2D** framework. If you use this code in your work, please cite the original reference paper:

bibtex
@article{Seleson2024,
  author = {Seleson, Pablo and Pasetto, Marco and John, Yohan and Trageser, Jeremy and Reeve, Samuel Temple},
  title = {PDMATLAB2D: A Peridynamics MATLAB Two-dimensional Code},
  journal = {Journal of Peridynamics and Nonlocal Modeling},
  volume = {6},
  number = {1},
  pages = {149-205},
  year = {2024},
  url = {https://doi.org/10.1007/s42102-023-00104-w}
}

If you wish to cite the original software release directly, refer to [Zenodo](https://zenodo.org/doi/10.5281/zenodo.7348667).

---

## License

This project is distributed under the terms of the **BSD 3-Clause License**, in compliance with the original upstream software. See the `LICENSE` file for details.
`

---
