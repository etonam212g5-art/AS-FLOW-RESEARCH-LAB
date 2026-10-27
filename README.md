# AS-Flow Spurious Vorticity Benchmark and Cauchy Tensor Regularization Package

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.placeholder.svg)](https://doi.org/10.5281/zenodo.placeholder)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Language: Python 3](https://img.shields.io/badge/Language-Python%203-blue.svg)](https://www.python.org/)
[![Model: WRF-ARW](https://img.shields.io/badge/Model-WRF--ARW%20dyn__em-green.svg)](https://www.mmm.ucar.edu/models/wrf)

**Author:** Ugvucbm Cvkzdv HDBDJVGDCGHSACJHMS  
**Associated Manuscript:** *On Spurious Vorticity Generation in Vector-Coupled Metric Damping for Geophysical CFD: Theoretical Analysis of the Parasitic Curl and Cauchy Tensor Regularization*  
**GitHub Repository:** [https://github.com/etonam212g5-art/AS-FLOW-RESEARCH-LAB](https://github.com/etonam212g5-art/AS-FLOW-RESEARCH-LAB)  
**Target Journals:** *Geoscientific Model Development (Copernicus / EGU)* / *Earth and Space Science (AGU / Wiley)*  

---

## Overview

In ultra-high-resolution mesoscale atmospheric modeling (sub-kilometer WRF simulations), extreme velocity shear in tropical cyclone eyewalls frequently triggers Courant-Friedrichs-Lewy (CFL) numerical singularities. Introducing a prognostic scalar metric strain variable $J(x, t)$ coupled via a vector restoring force:
$$\mathbf{F}_{\text{drag}} = - K \max(0, J - 1) \mathbf{u}$$
was proposed to damp stretching.

This package provides the complete numerical and source-code evidence demonstrating that:
1. **The Vector Coupling Bug:** This drag force acts as an artificial generator of vertical vorticity via an anomalous cross-gradient curl term $-K \theta(J-1) (\nabla J \times \mathbf{u})$, exceeding the physical damping by a factor of $\sim 25$ in Category-5 Cyclone Kenneth eyewall simulations.
2. **The Cauchy Regularization Solution:** Formulating the damping via the symmetric Cauchy rate-of-strain tensor $\boldsymbol{\tau} = 2 \nu_{\text{eff}}(J) \mathbf{S}$ exactly eliminates the cross torque and guarantees unconditional enstrophy dissipation.

---

## Package Contents

```text
PREUVE ET CODES POUR PUBLICATION/
├── pillar_a_spurious_vorticity_2d.py     # Pillar A: Standalone 2D Python benchmark (zero-dependency)
├── figure1_radial_torque_generation.png   # Publication Figure 1: Radial profiles of torque & strain
├── figure2_2d_vorticity_comparison.png    # Publication Figure 2: 2D maps of vector curl vs tensor curl
├── figure3_temporal_enstrophy_evolution.png # Publication Figure 3: Temporal enstrophy spin-up vs damping
├── pillar_b_module_diffusion_asflow_v14.F # Pillar B: Official WRF Fortran source module (lines 201-212)
├── pillar_b_fortran_wrf_analysis.md      # Detailed line-by-line mapping between Fortran and equations
├── metrics_asflow.csv                    # Pillar C: Cloud simulation metrics with vector coupling (W = 45.1 m/s)
├── metrics_smagorinsky.csv               # Pillar C: Cloud reference simulation metrics (Smagorinsky)
├── comparative_curves_asflow_vs_smagorinsky.png # Cloud comparison curve of Kenneth
├── CITATION.cff                          # Formal citation metadata for GitHub and Zotero
├── codemeta.json                         # Formal software metadata for Zenodo DOI attribution
└── README.md                             # This documentation
```

---

## Quickstart & Replication

### Requirements
* Python 3.8+
* `numpy`
* `matplotlib`

### Running the Pillar A Benchmark
Clone the repository or navigate to the directory and run:
```bash
python pillar_a_spurious_vorticity_2d.py
```
**Execution Time:** ~5 seconds on a standard laptop.  
**Output:** Automatically computes all spatial gradients, verifies that the parasitic torque is strictly positive in the cyclonic eyewall, verifies that the Cauchy tensor regularized curl is dissipative, and saves the three 300-DPI publication figures.

---

## Scientific Interpretation of Benchmark Figures

The Pillar A benchmark outputs three 300-DPI publication figures dissecting the fluid mechanics of Cyclone Kenneth ($V_{\max} = 70\text{ m/s}$, $R_{\max} = 25\text{ km}$, $\Delta x = 1\text{ km}$):

### Figure 1: Radial Profiles & Parasitic Torque Generation
* **Panel (a) - Cyclonic Kinematics & Strain:** The blue solid curve shows the tangential wind velocity $v_\theta(r)$ ramping up linearly from $0$ to $70\text{ m/s}$ ($140\text{ kt}$) inside the eye ($r \le 25\text{ km}$, solid-body vortex rotation) and decaying cyclonically outside. The purple dash-dotted curve shows the strain rate magnitude $S_{\text{mag}}$, which is zero inside the rigid-body core but spikes sharply at the eyewall ($4.26 \times 10^{-3}\text{ s}^{-1}$).
* **Panel (b) - Metric Deformation $J(r)$:** Demonstrates the emergence of the metric variable $J = 1 + S_{\text{mag}} / K_{\text{elastic}}$, which breaches the critical threshold $J = 1.0$ and peaks at $1.085$ at the eyewall, establishing a steep radial gradient $\nabla J = (\partial_r J) \mathbf{e}_r$.
* **Panel (c) - Spurious Vorticity vs Cauchy Annihilation:** The red solid curve plots the cross term $-K(\nabla J \times \mathbf{u})_z$. On the outer flank of the eyewall where $\partial_r J < 0$, $-K (\partial_r J) v_\theta > 0$, injecting a massive positive vertical torque ($+1.64 \times 10^{-4}\text{ s}^{-2}$). This parasitic torque is **24.93 times larger** than the intended physical damping $-K(J-1)\omega_z$. Conversely, the green curve shows the Cauchy tensor divergence formulation $\boldsymbol{\tau} = 2\nu(J)\mathbf{S}$, which strictly eliminates the cross torque ($\le 10^{-8}\text{ s}^{-2}$) to machine precision.

### Figure 2: 2D Spatial Maps (Eyewall Ring Structure)
* **Left Panel (Intuitive Vector Damping):** Displays an intense annular ring of spurious vorticity generation coinciding exactly with the cyclonic eyewall radius ($R_{\max} = 25\text{ km}$), acting as an artificial energy pump.
* **Right Panel (Cauchy Tensor Regularization):** Shows the complete disappearance of the parasitic ring, restoring isotropic, frame-indifferent viscous dissipation over the entire $200 \times 200\text{ km}$ domain.

### Figure 3: Temporal Enstrophy Evolution (Spin-Up vs Stability)
Integrating the vorticity equation $\partial_t \omega_z = (\nabla \times \mathbf{F})_z$ over time reveals that the vector coupling drives a catastrophic, runaway spin-up (red curve, jumping from $1.0 \times 10^{-3}$ to $3.7 \times 10^{-3}\text{ s}^{-1}$ in 3 minutes), while the Cauchy tensor regularization (green curve) preserves monotonic enstrophy dissipation and numerical stability.

---

## Cloud Simulation Data Summary (Cyclone Kenneth)

Extracted from real sub-kilometer WRF simulations (`metrics_asflow.csv` vs `metrics_smagorinsky.csv`):

| Metric | Smagorinsky Reference | Vector Coupling (AS-Flow v14) | Impact of Parasitic Curl |
| :--- | :---: | :---: | :--- |
| **Peak Vertical Velocity ($W_{\max}$)** | **0.017 m/s** | **45.12 m/s** | Spurious violent updraft explosion |
| **Domain Enstrophy ($\mathcal{E}$)** | $1.25 \times 10^{-9}\text{ s}^{-2}$ | **$2.45 \times 10^{-7}\text{ s}^{-2}$** | **$200\times$ artificial vorticity surge** |
| **Max Horizontal Wind ($V_{\max}$)** | 10.2 m/s | 18.4 m/s | Artificial vortex acceleration |

---

## Citation

If you use this benchmark or the Cauchy regularization in your research, please cite:

```bibtex
@article{daglan2026spurious,
  title={On Spurious Vorticity Generation in Vector-Coupled Metric Damping for Geophysical CFD: Theoretical Analysis of the Parasitic Curl and Cauchy Tensor Regularization},
  author={HDBDJVGDCGHSACJHMS, Ugvucbm Cvkzdv},
  journal={Geoscientific Model Development (Submitted)},
  year={2026}
}

@software{daglan2026code,
  author={HDBDJVGDCGHSACJHMS, Ugvucbm Cvkzdv},
  title={AS-Flow Spurious Vorticity Benchmark and Cauchy Regularization Package},
  year={2026},
  publisher={Zenodo},
  doi={10.5281/zenodo.placeholder}
}
```

---
**License:** Open Source under the [MIT License](https://opensource.org/licenses/MIT).
