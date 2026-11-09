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






    # ---------- Find U-233 amount for k = 0.95 ----------
def run_k(u233_frac):
    fuel = make_th_fuel(u233_frac)
    materials = openmc.Materials([fuel, clad, water])
    geometry = make_assembly(fuel)

    settings = openmc.Settings()
    settings.run_mode = "eigenvalue"
    settings.particles = 10000
    settings.batches = 100
    settings.inactive = 20
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Box((-half, -half, -1), (half, half, 1)),
        constraints={"fissionable": True},
    )

    model = openmc.Model(geometry=geometry, materials=materials, settings=settings)
    out_dir = f"04_subcritical/output/sweep_{u233_frac:.3f}"
    os.makedirs(out_dir, exist_ok=True)
    sp_path = model.run(threads=8, cwd=out_dir, output=False)
    with openmc.StatePoint(sp_path) as sp:
        return sp.keff.nominal_value

# Uncomment to rerun the U-233 sweep (slow)
# fractions = [0.010, 0.015, 0.020, 0.025, 0.030]
# for f in fractions:
#     k = run_k(f)
#     print(f"U-233 = {f * 100:.1f}%  ->  k = {k:.4f}")
u233_frac = 0.013
k_sub = run_k(u233_frac)
print(f"U-233 = {u233_frac * 100:.1f}%  ->  k = {k_sub:.4f}")
print(f"Expected multiplication 1/(1-k) = {1 / (1 - k_sub):.1f}x")









# ---------- Fixed-source run (driver ON) ----------
def run_driver(u233_frac):
    fuel = make_th_fuel(u233_frac)
    materials = openmc.Materials([fuel, clad, water])
    geometry = make_assembly(fuel)

    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.particles = 10000
    settings.batches = 20
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((0, 0, 0)),
        angle=openmc.stats.Isotropic(),
        energy=openmc.stats.Discrete([14.1e6], [1.0]),
    )

    flux_mesh = openmc.RegularMesh()
    flux_mesh.dimension = (85, 85)
    flux_mesh.lower_left = (-half, -half)
    flux_mesh.upper_right = (half, half)

    rod_mesh = openmc.RegularMesh()
    rod_mesh.dimension = (17, 17)
    rod_mesh.lower_left = (-half, -half)
    rod_mesh.upper_right = (half, half)

    total = openmc.Tally(name="total")
    total.scores = ["nu-fission"]

    flux_tally = openmc.Tally(name="flux_map")
    flux_tally.filters = [openmc.MeshFilter(flux_mesh)]
    flux_tally.scores = ["flux"]

    fission_tally = openmc.Tally(name="fission_map")
    fission_tally.filters = [openmc.MeshFilter(rod_mesh)]
    fission_tally.scores = ["fission"]

    tallies = openmc.Tallies([total, flux_tally, fission_tally])

    model = openmc.Model(geometry=geometry, materials=materials, settings=settings, tallies=tallies)
    out_dir = "04_subcritical/output/driver"
    os.makedirs(out_dir, exist_ok=True)
    sp_path = model.run(threads=8, cwd=out_dir, output=False)
    with openmc.StatePoint(sp_path) as sp:
        made = sp.get_tally(name="total").mean.flatten()[0]
    return made, sp_path

made, drv_sp = run_driver(u233_frac)
M_sim = 1 + made
print(f"Fission neutrons made per driver neutron = {made:.1f}")
print(f"Simulated multiplication = {M_sim:.1f}x  (formula said {1 / (1 - k_sub):.1f}x)")

# ---------- Driver heatmaps ----------
with openmc.StatePoint(drv_sp) as sp:
    flux = sp.get_tally(name="flux_map").get_values(scores=["flux"]).reshape(85, 85)
    power = sp.get_tally(name="fission_map").get_values(scores=["fission"]).reshape(17, 17)

flux = flux / flux.max()
power = np.ma.masked_equal(power, 0)
power = power / power.mean()

extent = (-half, half, -half, half)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6))

im1 = ax1.imshow(flux, origin="lower", extent=extent, cmap="inferno", vmax=np.percentile(flux, 99))
fig.colorbar(im1, ax=ax1, shrink=0.8)
ax1.set_title("Driver ON: flux")

cmap = plt.get_cmap("RdYlGn_r").copy()
cmap.set_bad("black")
im2 = ax2.imshow(power, origin="lower", extent=extent, cmap=cmap)
fig.colorbar(im2, ax=ax2, shrink=0.8)
ax2.set_title("Driver ON: power per rod")

os.makedirs("results/figures", exist_ok=True)
plt.savefig("results/figures/subcritical_heatmaps.png", dpi=150)
print("Saved to results/figures/subcritical_heatmaps.png")



# ---------- Particle tracks ----------
def run_tracks(u233_frac, n_tracks=3):
    fuel = make_th_fuel(u233_frac)
    materials = openmc.Materials([fuel, clad, water])
    geometry = make_assembly(fuel)

    settings = openmc.Settings()
    settings.run_mode = "fixed source"
    settings.particles = n_tracks
    settings.batches = 1
    settings.source = openmc.IndependentSource(
        space=openmc.stats.Point((0, 0, 0)),
        angle=openmc.stats.Isotropic(),
        energy=openmc.stats.Discrete([14.1e6], [1.0]),
    )
    settings.track = [(1, 1, i) for i in range(1, n_tracks + 1)]

    model = openmc.Model(geometry=geometry, materials=materials, settings=settings)
    out_dir = "04_subcritical/output/tracks"
    os.makedirs(out_dir, exist_ok=True)
    model.run(threads=1, cwd=out_dir, output=False)
    return openmc.Tracks(os.path.join(out_dir, "tracks.h5")), geometry, fuel

tracks, geom, fuel = run_tracks(0.013)

colors = {fuel: "green", clad: "gray", water: "skyblue"}
ax = geom.root_universe.plot(
    width=(2 * half, 2 * half), pixels=(600, 600),
    color_by="material", colors=colors,
)

for track in tracks:
    for i, pt in enumerate(track.particle_tracks):
        r = pt.states["r"]
        if i == 0:
            ax.plot(r["x"], r["y"], color="red", linewidth=2)
        else:
            ax.plot(r["x"], r["y"], color="yellow", linewidth=0.7)

ax.set_title("3 driver neutrons (red) and their fission babies (yellow)")
os.makedirs("results/figures", exist_ok=True)
plt.savefig("results/figures/particle_tracks.png", dpi=150)
print("Saved to results/figures/particle_tracks.png")
