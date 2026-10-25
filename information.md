# Background Information

Plain-English notes on the physics behind this project. Read this before (and during) each step.

---

## 1. The neutron life cycle

A reactor runs on a chain reaction:

1. A neutron hits a fissile nucleus (like U-235 or U-233)
2. The nucleus splits, releasing energy and 2 to 3 new neutrons
3. Those neutrons either cause more fissions, get absorbed without fission, or leak out

Everything in reactor physics is about tracking where neutrons go.

---

## 2. k-effective (k-eff)

The single most important number in a reactor.

**k-eff = neutrons in this generation / neutrons in the previous generation**

| k-eff | State | What happens |
|-------|-------|--------------|
| k < 1 | Subcritical | Chain reaction dies out on its own |
| k = 1 | Critical | Steady, self-sustaining |
| k > 1 | Supercritical | Neutron population grows |

**Example:** k = 0.95 means every 100 neutrons produce 95 in the next generation.

### Six-factor formula

k-eff can be broken into pieces:

**k-eff = η · f · p · ε · P_FNL · P_TNL**

- **η (eta):** neutrons produced per neutron absorbed in fuel
- **f:** fraction of thermal neutrons absorbed in fuel (vs everything else)
- **p:** chance a neutron slows down without getting captured
- **ε (epsilon):** bonus from fast fissions
- **P_FNL / P_TNL:** chance fast / thermal neutrons don't leak out

Useful for explaining *why* k changes when you change the design.

---

## 3. Cross sections

A cross section is basically the "target size" a nucleus shows to a neutron, measured in **barns** (1 barn = 10⁻²⁴ cm²).

- Each reaction (fission, capture, scattering) has its own cross section
- They depend heavily on neutron energy
- Slow (thermal) neutrons are much more likely to cause fission in U-235 and U-233

**Reaction rate = flux × cross section × number density**

R = φ · σ · N

---

## 4. Neutron flux

**Flux (φ)** = how many neutrons pass through an area per second.

- Usually highest in the center of the core, lower near the edges
- Power is roughly proportional to flux times fission cross section
- A flux heatmap shows where the reactor is "working hardest"

---

## 5. Monte Carlo transport (what OpenMC does)

Instead of solving equations directly, OpenMC simulates millions of individual neutrons:

1. Start a neutron somewhere with some energy
2. Randomly sample how far it travels and what it hits
3. Follow it until it's absorbed or escapes
4. Repeat millions of times and average the results

**Key terms**
- **Batches:** groups of neutron histories
- **Inactive batches:** early batches thrown away while the source settles
- **Active batches:** batches used for the final answer
- **Tally:** a quantity you ask OpenMC to count (flux, fission rate, etc.)
- **Statistical uncertainty:** every result comes with a ± error; more particles means smaller error

---

## 6. The thorium fuel cycle

Th-232 is **fertile**, not **fissile**. It can't keep a chain reaction going alone, but it can be turned into fuel:

```
Th-232 + neutron → Th-233 → Pa-233 → U-233
                  (~22 min)  (~27 days)
```

- **U-233** is fissile and very good with thermal neutrons
- Thorium is more abundant than uranium
- Produces less long-lived waste in some respects
- Downside: needs a starter (U-233, U-235, or an external neutron source) to get going

---

## 7. Subcritical systems with an external source

A subcritical core (k < 1) can't sustain itself, but it can **multiply** neutrons from an outside source.

**Source multiplication:**

**M ≈ 1 / (1 − k)**

| k | M (approx.) |
|---|-------------|
| 0.90 | 10 |
| 0.95 | 20 |
| 0.98 | 50 |

**Why this is attractive**
- Turn off the source and the reaction stops on its own
- Can breed fuel (like U-233 from thorium) while running
- Tradeoff: you need a strong, reliable neutron source

---

## 8. Fuel depletion

As a reactor runs, its fuel composition changes:

- Fissile atoms get used up
- Fission products build up (some absorb neutrons and hurt k)
- Fertile atoms (like Th-232) convert into fissile ones (like U-233)

OpenMC handles this by alternating transport runs with "burn" steps that update the material compositions over time.

**Key term: burnup** = energy produced per mass of fuel (usually MWd/kgHM).

---

## 9. Verification and validation (V&V)

- **Verification:** did I solve the equations correctly? (code works as intended)
- **Validation:** do my results match reality? (compare to benchmarks or experiments)

Simple checks used in this project:
- Compare subcritical multiplication to M ≈ 1 / (1 − k)
- Make sure uncertainties shrink as particle count goes up
- Sanity check trends (more fissile fuel → higher k)

---

## 10. Sensitivity and uncertainty

- **Sensitivity study:** change one input, see how much the output moves
- **Uncertainty quantification:** how confident are we in the result given uncertain inputs?

**Example:** sweep fuel rod spacing from tight to loose and plot k-eff. Usually k rises, peaks, then falls (too little vs too much water to slow neutrons down).

---

## Glossary

| Term | Meaning |
|------|---------|
| Fissile | Can fission with slow neutrons (U-235, U-233, Pu-239) |
| Fertile | Can be turned into fissile material (Th-232, U-238) |
| Moderator | Material that slows neutrons down (water, graphite) |
| Thermal neutron | Slow neutron, around 0.025 eV |
| Fast neutron | High-energy neutron, around 1 MeV and up |
| Reactivity (ρ) | How far from critical: ρ = (k − 1) / k |
| Pin cell | One fuel rod with its surrounding coolant |
| Assembly | A grid of fuel rods |
| Shielding | Material that blocks radiation |
| Burnup | Energy produced per mass of fuel |

---

## Resources

- **OpenMC documentation and user guide** (official site), especially the Python API examples
- **OpenMC example notebooks** on GitHub (pin cell, assembly, depletion)
- **MIT OpenCourseWare:** introductory nuclear engineering courses
- **Textbook:** *Nuclear Reactor Analysis* by Duderstadt and Hamilton
- **Textbook:** *Introduction to Nuclear Engineering* by Lamarsh and Baratta
