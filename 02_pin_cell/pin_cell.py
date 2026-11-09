import os
import matplotlib.pyplot as plt
import numpy as np
import openmc

# ---------- Materials ----------
fuel = openmc.Material(name="UO2 fuel")
fuel.add_element("U", 1.0, enrichment=3.0)
fuel.add_element("O", 2.0)
fuel.set_density("g/cm3", 10.4)

clad = openmc.Material(name="Zirconium cladding")
clad.add_element("Zr", 1.0)
clad.set_density("g/cm3", 6.55)

water = openmc.Material(name="Water")
water.add_element("H", 2.0)
water.add_element("O", 1.0)
water.set_density("g/cm3", 1.0)
water.add_s_alpha_beta("c_H_in_H2O")

materials = openmc.Materials([fuel, clad, water])

# ---------- Sizes (cm) ----------
fuel_radius = 0.39
clad_radius = 0.46
pitch = 1.26

# ---------- Geometry ----------
fuel_surface = openmc.ZCylinder(r=fuel_radius)
clad_surface = openmc.ZCylinder(r=clad_radius)
box = openmc.model.RectangularPrism(pitch, pitch, boundary_type="reflective")

fuel_cell = openmc.Cell(fill=fuel, region=-fuel_surface)
clad_cell = openmc.Cell(fill=clad, region=+fuel_surface & -clad_surface)
water_cell = openmc.Cell(fill=water, region=+clad_surface & -box)

geometry = openmc.Geometry([fuel_cell, clad_cell, water_cell])

# ---------- Picture ----------
colors = {fuel: "orange", clad: "gray", water: "skyblue"}
ax = geometry.root_universe.plot(
    width=(pitch, pitch), pixels=(500, 500),
    color_by="material", colors=colors,
)
ax.set_title("Pin cell (top view)")

os.makedirs("results/figures", exist_ok=True)
plt.savefig("results/figures/pin_cell_geometry.png", dpi=150)
print("Saved to results/figures/pin_cell_geometry.png")

# ---------- Simulation settings ----------
settings = openmc.Settings()
settings.run_mode = "eigenvalue"
settings.particles = 10000
settings.batches = 100
settings.inactive = 20
settings.source = openmc.IndependentSource(space=openmc.stats.Point())

# ---------- Energy tally ----------
energy_bins = np.logspace(-5, 7.3, 501)
energy_filter = openmc.EnergyFilter(energy_bins)
tally = openmc.Tally(name="spectrum")
tally.filters = [energy_filter]
tally.scores = ["flux"]
tallies = openmc.Tallies([tally])

model = openmc.Model(geometry=geometry, materials=materials, settings=settings, tallies=tallies)

# ---------- Run ----------
os.makedirs("02_pin_cell/output", exist_ok=True)
sp_path = model.run(threads=8, cwd="02_pin_cell/output")

with openmc.StatePoint(sp_path) as sp:
    k = sp.keff
    print(f"k-eff = {k.nominal_value:.5f} +/- {k.std_dev:.5f}")

    # ---------- k convergence plot ----------
with openmc.StatePoint(sp_path) as sp:
    k_each_round = sp.k_generation
    k_final = sp.keff.nominal_value

rounds = range(1, len(k_each_round) + 1)

plt.figure()
plt.plot(rounds, k_each_round, label="k each round")
plt.axvline(20, color="red", linestyle="--", label="end of warm-up")
plt.axhline(k_final, color="green", label=f"final k = {k_final:.4f}")
plt.xlabel("Round (batch)")
plt.ylabel("k")
plt.title("k convergence")
plt.legend()
plt.savefig("results/figures/k_convergence.png", dpi=150)
print("Saved to results/figures/k_convergence.png")


# ---------- Energy spectrum plot ----------
with openmc.StatePoint(sp_path) as sp:
    t = sp.get_tally(name="spectrum")
    flux = t.mean.flatten()

energy_mid = np.sqrt(energy_bins[:-1] * energy_bins[1:])

plt.figure()
plt.semilogx(energy_mid, flux)
plt.xlabel("Energy (eV)")
plt.ylabel("Flux (per bin)")
plt.title("Neutron energy spectrum")
plt.savefig("results/figures/energy_spectrum.png", dpi=150)
print("Saved to results/figures/energy_spectrum.png")