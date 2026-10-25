# Thorium Subcritical Reactor Study (OpenMC)

A hands-on, visual reactor physics project built with **OpenMC** and **Python**.

The goal: model a small thorium-fueled core, run it subcritical with an external neutron source, track fuel breeding over time, and make every step something you can actually **see**.

> Background theory, equations, and a glossary live in [information.md](information.md).

---

## What this project covers

- Reactor physics basics: k-effective, neutron flux, energy spectrum
- Thorium fuel cycle: Th-232 breeding into U-233
- Subcritical systems driven by an external neutron source
- Fuel depletion over time
- Sensitivity studies and verification against theory
- Reproducible, scripted, version-controlled workflow

---

## Visual outputs

| Step | What you'll see |
|------|-----------------|
| Cross sections | Fission vs capture curves for U-235, U-233, Th-232 |
| Pin cell | Geometry picture, k-eff convergence, energy spectrum |
| Assembly | Flux heatmaps, uranium vs thorium side by side |
| Subcritical | Flux spreading from a source, particle tracks, M vs k curve |
| Depletion | Animated GIF of U-233 building up, k-eff vs time |
| Dashboard | Streamlit app with sliders to explore design cases |

Screenshots will be added to `results/` as each step is finished.

---

## Repo structure

```
thorium-subcritical-openmc/
├── README.md
├── information.md
├── environment.yml
├── 01_cross_sections/
├── 02_pin_cell/
├── 03_thorium_assembly/
├── 04_subcritical_source/
├── 05_depletion/
├── 06_sensitivity/
├── dashboard/
│   └── app.py
└── results/
    ├── figures/
    └── gifs/
```

---

## Setup

**Requirements**
- Linux, macOS, or Windows with WSL (Ubuntu)
- Conda (Anaconda or Miniconda)

**1. Create the environment**
```bash
conda create -n openmc-env -c conda-forge openmc matplotlib pandas streamlit imageio
conda activate openmc-env
```

**2. Get nuclear data**
- Download an HDF5 cross-section library (for example ENDF/B-VII.1) from the official OpenMC data page
- Point OpenMC to it:
```bash
export OPENMC_CROSS_SECTIONS=/path/to/cross_sections.xml
```
- For depletion, also download a depletion chain file (XML)

**3. Check it works**
```bash
python -c "import openmc; print(openmc.__version__)"
```

---

## Roadmap

- [ ] **Step 1:** Setup (OpenMC, nuclear data, Git)
- [ ] **Step 2:** Theory + cross-section plots
- [ ] **Step 3:** Uranium pin cell (k-eff, spectrum, geometry plot)
- [ ] **Step 4:** Thorium assembly (flux heatmaps, U vs Th comparison)
- [ ] **Step 5:** Subcritical run with external source (tracks, M vs k)
- [ ] **Step 6:** Depletion (U-233 breeding animation, k-eff vs time)
- [ ] **Step 7:** Sensitivity sweep + Streamlit dashboard + write-up

---

## Verification checks

Every result gets a "does this make physical sense?" check, for example:

- Uranium pin cell k-eff should land in a realistic range for the chosen enrichment
- Pure Th-232 fuel should come out well below k = 1 (it is not fissile)
- Measured source multiplication should follow **M ≈ 1 / (1 − k)**
- U-233 inventory should rise over time as Th-232 is consumed

Results, assumptions, and any discrepancies are documented in each step's folder.

---

## Tools

- **OpenMC** for Monte Carlo neutron transport (C++ core, Python API)
- **Matplotlib** for plots and animations
- **Pandas** for processing results
- **Streamlit** for the interactive dashboard
- **Git** for version control

---

## Author

Lucas Hoskin
