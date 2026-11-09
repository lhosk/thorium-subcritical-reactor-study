import os
import numpy as np
import matplotlib.pyplot as plt
import openmc

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
true_k = 0.9553

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

# ---------- Detectors (in water holes) ----------
detectors = {
    "D1 near":   (8, 11),
    "D2 diag":   (5, 5),
    "D3 far":    (2, 8),
    "D4 corner": (3, 3),
}
det_colors = ["magenta", "orange", "yellow", "white"]
line_colors = ["magenta", "orange", "gold", "black"]

def hole_center(row, col):
    x = -half + (col + 0.5) * pitch
    y = half - (row + 0.5) * pitch
    return x, y

# ---------- Picture of detector spots ----------
fuel = make_th_fuel(u233_frac)
geometry = make_assembly(fuel)

colors = {fuel: "green", clad: "gray", water: "skyblue"}
ax = geometry.root_universe.plot(
    width=(2 * half, 2 * half), pixels=(600, 600),
    color_by="material", colors=colors,
)

ax.plot(0, 0, "r*", markersize=18, markeredgecolor="black", label="Driver")
for (name, (row, col)), c in zip(detectors.items(), det_colors):
    x, y = hole_center(row, col)
    ax.plot(x, y, "o", markersize=12, color=c, markeredgecolor="black", label=name)

ax.legend(loc="upper right", fontsize=8)
ax.set_title("Driver and detector locations")
os.makedirs("results/figures", exist_ok=True)
plt.savefig("results/figures/detector_locations.png", dpi=150)
plt.close()
print("Saved to results/figures/detector_locations.png")

# ---------- Time-dependent driver run ----------
def run_time(label, time_dist, time_bins, fission_on):
    fuel = make_th_fuel(u233_frac)
    materials = openmc.Materials([fuel, clad, water])
    geometry = make_assembly(fuel)

    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.particles = 5000
    settings.batches = 10
    settings.create_fission_neutrons = fission_on
    settings.cutoff = {"time_neutron": time_bins[-1]}
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((0, 0, 0)),
        angle=openmc.stats.Isotropic(),
        energy=openmc.stats.Discrete([14.1e6], [1.0]),
        time=time_dist,
    )

    rod_mesh = openmc.RegularMesh()
    rod_mesh.dimension = (17, 17)
    rod_mesh.lower_left = (-half, -half)
    rod_mesh.upper_right = (half, half)

    tally = openmc.Tally(name="det")
    tally.filters = [openmc.TimeFilter(time_bins), openmc.MeshFilter(rod_mesh)]
    tally.scores = ["flux"]
    tallies = openmc.Tallies([tally])

    model = openmc.Model(geometry=geometry, materials=materials, settings=settings, tallies=tallies)
    out_dir = f"04_subcritical/output/time_{label}_{'on' if fission_on else 'off'}"
    os.makedirs(out_dir, exist_ok=True)
    sp_path = model.run(threads=8, cwd=out_dir, output=False)

    with openmc.StatePoint(sp_path) as sp:
        vals = sp.get_tally(name="det").mean.reshape(len(time_bins) - 1, 17, 17)

    return np.array([vals[:, 16 - row, col] for row, col in detectors.values()])

def get_data(label, time_dist, time_bins):
    path = f"04_subcritical/output/{label}_data.npz"
    if os.path.exists(path):
        d = np.load(path)
        print(f"Loaded saved {label} data")
        return d["bins"], d["on"], d["off"]
    on = run_time(label, time_dist, time_bins, fission_on=True)
    off = run_time(label, time_dist, time_bins, fission_on=False)
    np.savez(path, bins=time_bins, on=on, off=off)
    return time_bins, on, off

# ---------- Pulse: all driver neutrons at t = 0 ----------
pulse_bins, pulse_on, pulse_off = get_data(
    "pulse",
    openmc.stats.Discrete([0.0], [1.0]),
    np.concatenate(([0.0], np.logspace(-9, -1, 41))),
)

# ---------- Continuous: driver neutrons spread evenly over 20 ms ----------
T = 0.02
cont_bins, cont_on, cont_off = get_data(
    "continuous",
    openmc.stats.Uniform(0.0, T),
    np.linspace(0.0, T, 41),
)

# ---------- Turn readings into flux and k guesses ----------
def flux_rate(arr, bins):
    return arr / np.diff(bins)

def k_guess(on, off):
    with np.errstate(divide="ignore", invalid="ignore"):
        M = on / np.where(off > 0, off, np.nan)
        return 1 - 1 / M

# Pulse: k from everything counted up to time t
pulse_k = k_guess(np.cumsum(pulse_on, axis=1), np.cumsum(pulse_off, axis=1))
pulse_flux = flux_rate(pulse_on, pulse_bins)
pulse_t = pulse_bins[1:]

# Continuous: k from each time bucket on its own
cont_k = k_guess(cont_on, cont_off)
cont_flux = flux_rate(cont_on, cont_bins)
cont_t = 0.5 * (cont_bins[:-1] + cont_bins[1:]) * 1000  # ms

# ---------- Plots ----------
fig, axes = plt.subplots(2, 2, figsize=(13, 10))

ax = axes[0, 0]
for i, name in enumerate(detectors):
    ax.loglog(pulse_t, pulse_flux[i], color=line_colors[i], label=name)
ax.set_title("Pulse: neutron flux vs time")
ax.set_xlabel("Time (s)")
ax.set_ylabel("Flux (per driver neutron per s)")
ax.legend()

ax = axes[0, 1]
for i, name in enumerate(detectors):
    ax.semilogx(pulse_t, pulse_k[i], color=line_colors[i], label=name)
ax.axhline(true_k, color="red", linestyle="--", label=f"true k = {true_k}")
ax.set_ylim(0.5, 1.0)
ax.set_title("Pulse: k guess vs time")
ax.set_xlabel("Time (s)")
ax.set_ylabel("k guess")
ax.legend()

ax = axes[1, 0]
for i, name in enumerate(detectors):
    ax.plot(cont_t, cont_flux[i], color=line_colors[i], label=name)
ax.set_title("Continuous: neutron flux vs time")
ax.set_xlabel("Time (ms)")
ax.set_ylabel("Flux (per driver neutron per s)")
ax.legend()

ax = axes[1, 1]
for i, name in enumerate(detectors):
    ax.plot(cont_t, cont_k[i], color=line_colors[i], label=name)
ax.axhline(true_k, color="red", linestyle="--", label=f"true k = {true_k}")
ax.set_ylim(0.5, 1.0)
ax.set_title("Continuous: k guess vs time")
ax.set_xlabel("Time (ms)")
ax.set_ylabel("k guess")
ax.legend()

plt.tight_layout()
plt.savefig("results/figures/detector_time_plots.png", dpi=150)
print("Saved to results/figures/detector_time_plots.png")