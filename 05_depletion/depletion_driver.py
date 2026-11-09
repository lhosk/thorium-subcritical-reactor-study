import os
import math
import numpy as np
import matplotlib.pyplot as plt
import openmc
import openmc.deplete

openmc.deplete.pool.USE_MULTIPROCESSING = False
openmc.config["chain_file"] = os.path.expanduser("~/nuclear_data/chain_endfb80_pwr.xml")

# ---------- Paths ----------
ROOT = os.getcwd()
OUT = os.path.join(ROOT, "05_depletion/output/driver")
FIG = os.path.join(ROOT, "results/figures")

# ---------- Materials ----------
def make_th_fuel(u233_frac):
    fuel = openmc.Material(name=f"ThO2 + {u233_frac * 100:.1f}% U-233")
    fuel.add_element("Th", 1.0 - u233_frac)
    fuel.add_nuclide("U233", u233_frac)
    fuel.add_element("O", 2.0)
    fuel.set_density("g/cm3", 10.0)
    return fuel

clad = openmc.Material(name="Zirconium cladding")
clad.add_element("Zr", 1.0)
clad.set_density("g/cm3", 6.55)

water = openmc.Material(name="Water")
water.add_element("H", 2.0)
water.add_element("O", 1.0)
water.set_density("g/cm3", 1.0)
water.add_s_alpha_beta("c_H_in_H2O")

# ---------- Sizes (cm) ----------
fuel_radius = 0.39
clad_radius = 0.46
pitch = 1.26
half = 17 * pitch / 2
u233_frac = 0.013

# ---------- Assembly ----------
water_holes = [
    (2, 5), (2, 8), (2, 11),
    (3, 3), (3, 13),
    (5, 2), (5, 5), (5, 8), (5, 11), (5, 14),
    (8, 2), (8, 5), (8, 8), (8, 11), (8, 14),
    (11, 2), (11, 5), (11, 8), (11, 11), (11, 14),
    (13, 3), (13, 13),
    (14, 5), (14, 8), (14, 11),
]
n_pins = 17 * 17 - len(water_holes)

def make_pin(fuel):
    fuel_surface = openmc.ZCylinder(r=fuel_radius)
    clad_surface = openmc.ZCylinder(r=clad_radius)
    fuel_cell = openmc.Cell(fill=fuel, region=-fuel_surface)
    clad_cell = openmc.Cell(fill=clad, region=+fuel_surface & -clad_surface)
    water_cell = openmc.Cell(fill=water, region=+clad_surface)
    return openmc.Universe(cells=[fuel_cell, clad_cell, water_cell])

def make_assembly(fuel):
    pin = make_pin(fuel)
    hole = openmc.Universe(cells=[openmc.Cell(fill=water)])
    n = 17

    grid = [[pin] * n for _ in range(n)]
    for row, col in water_holes:
        grid[row][col] = hole

    lattice = openmc.RectLattice()
    lattice.lower_left = (-half, -half)
    lattice.pitch = (pitch, pitch)
    lattice.universes = grid

    box = openmc.model.RectangularPrism(2 * half, 2 * half, boundary_type="reflective")
    root_cell = openmc.Cell(fill=lattice, region=-box)
    return openmc.Geometry([root_cell])

# ---------- Build model ----------
fuel = make_th_fuel(u233_frac)
fuel.depletable = True
fuel.volume = n_pins * math.pi * fuel_radius**2

materials = openmc.Materials([fuel, clad, water])
geometry = make_assembly(fuel)

settings = openmc.Settings()
settings.run_mode = "fixed source"
settings.particles = 2000
settings.batches = 10
settings.source = openmc.IndependentSource(
    space=openmc.stats.Point((0, 0, 0)),
    angle=openmc.stats.Isotropic(),
    energy=openmc.stats.Discrete([14.1e6], [1.0]),
)

total = openmc.Tally(name="total")
total.scores = ["nu-fission", "fission"]
tallies = openmc.Tallies([total])

model = openmc.Model(geometry=geometry, materials=materials, settings=settings, tallies=tallies)

# ---------- Depletion (driver held constant) ----------
source_rate = 1.5e14  # driver neutrons per second (per cm of height)
timesteps = [1, 4, 5, 20, 30, 30, 90, 185]  # days
results_path = os.path.join(OUT, "depletion_results.h5")

if not os.path.exists(results_path):
    os.makedirs(OUT, exist_ok=True)
    os.chdir(OUT)
    operator = openmc.deplete.CoupledOperator(model, normalization_mode="source-rate")
    integrator = openmc.deplete.PredictorIntegrator(
        operator, timesteps,
        source_rates=[source_rate] * len(timesteps),
        timestep_units="d",
    )
    integrator.integrate()
    os.chdir(ROOT)
else:
    print("Loaded saved driver depletion results")

# ---------- Read results ----------
results = openmc.deplete.Results(results_path)
mat_id = list(results[0].index_mat.keys())[0]

days, u233 = results.get_atoms(mat_id, "U233", time_units="d")
_, pa233 = results.get_atoms(mat_id, "Pa233", time_units="d")
_, th232 = results.get_atoms(mat_id, "Th232", time_units="d")

made, fissions = [], []
for i in range(len(days)):
    with openmc.StatePoint(os.path.join(OUT, f"openmc_simulation_n{i}.h5")) as sp:
        t = sp.get_tally(name="total")
        made.append(t.get_values(scores=["nu-fission"]).flatten()[0])
        fissions.append(t.get_values(scores=["fission"]).flatten()[0])

made = np.array(made)
fissions = np.array(fissions)

M = 1 + made
k_guess = 1 - 1 / M
joules_per_fission = 200e6 * 1.602e-19
power_per_pin = fissions * source_rate * joules_per_fission / n_pins  # W per cm

# ---------- Plots ----------
fig, axes = plt.subplots(2, 2, figsize=(13, 10))

ax = axes[0, 0]
ax.plot(days, M, "o-", color="purple")
ax.set_title("Multiplication M vs time")
ax.set_xlabel("Time (days)")
ax.set_ylabel("Neutrons per driver neutron")

ax = axes[0, 1]
ax.plot(days, k_guess, "o-", color="red")
ax.set_title("k guess (1 - 1/M) vs time")
ax.set_xlabel("Time (days)")
ax.set_ylabel("k")

ax = axes[1, 0]
ax.plot(days, power_per_pin, "o-", color="black")
ax.set_title("Power per pin (driver held constant)")
ax.set_xlabel("Time (days)")
ax.set_ylabel("W per cm")

ax = axes[1, 1]
ax.plot(days, u233, "o-", color="blue", label="U-233 (fuel)")
ax.plot(days, pa233, "o-", color="orange", label="Pa-233 (waiting room)")
ax.set_title("Breeding")
ax.set_xlabel("Time (days)")
ax.set_ylabel("Atoms (per cm, whole assembly)")
ax.legend()

plt.tight_layout()
os.makedirs(FIG, exist_ok=True)
plt.savefig(os.path.join(FIG, "depletion_driver.png"), dpi=150)
print("Saved to results/figures/depletion_driver.png")

# ---------- Quick numbers ----------
print(f"M: {M[0]:.1f}x -> {M[-1]:.1f}x")
print(f"k guess: {k_guess[0]:.4f} -> {k_guess[-1]:.4f}")
print(f"Power per pin: {power_per_pin[0]:.0f} -> {power_per_pin[-1]:.0f} W/cm")
print(f"U-233: {u233[0]:.3e} -> {u233[-1]:.3e} atoms")
print(f"Pa-233 at end: {pa233[-1]:.3e} atoms")