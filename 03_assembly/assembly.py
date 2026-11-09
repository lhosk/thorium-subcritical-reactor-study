import os
import numpy as np
import matplotlib.pyplot as plt
import openmc

# ---------- Materials ----------
u_fuel = openmc.Material(name="UO2 fuel")
u_fuel.add_element("U", 1.0, enrichment=3.0)
u_fuel.add_element("O", 2.0)
u_fuel.set_density("g/cm3", 10.4)

th_fuel = openmc.Material(name="ThO2 + U-233 fuel")
th_fuel.add_element("Th", 0.95)
th_fuel.add_nuclide("U233", 0.05)
th_fuel.add_element("O", 2.0)
th_fuel.set_density("g/cm3", 10.0)

clad = openmc.Material(name="Zirconium cladding")
clad.add_element("Zr", 1.0)
clad.set_density("g/cm3", 6.55)

water = openmc.Material(name="Water")
water.add_element("H", 2.0)
water.add_element("O", 1.0)
water.set_density("g/cm3", 1.0)
water.add_s_alpha_beta("c_H_in_H2O")

materials = openmc.Materials([u_fuel, th_fuel, clad, water])

# ---------- Sizes (cm) ----------
fuel_radius = 0.39
clad_radius = 0.46
pitch = 1.26

# ---------- One rod (a reusable stamp) ----------
def make_pin(fuel):
    fuel_surface = openmc.ZCylinder(r=fuel_radius)
    clad_surface = openmc.ZCylinder(r=clad_radius)
    fuel_cell = openmc.Cell(fill=fuel, region=-fuel_surface)
    clad_cell = openmc.Cell(fill=clad, region=+fuel_surface & -clad_surface)
    water_cell = openmc.Cell(fill=water, region=+clad_surface)
    return openmc.Universe(cells=[fuel_cell, clad_cell, water_cell])


# # ---------- Grid of rods (17 x 17) ----------
# def make_assembly(fuel):
#     pin = make_pin(fuel)
#     n = 17
#     lattice = openmc.RectLattice()
#     lattice.lower_left = (-n * pitch / 2, -n * pitch / 2)
#     lattice.pitch = (pitch, pitch)
#     lattice.universes = [[pin] * n for _ in range(n)]

#     box = openmc.model.RectangularPrism(n * pitch, n * pitch, boundary_type="reflective")
#     root_cell = openmc.Cell(fill=lattice, region=-box)
#     return openmc.Geometry([root_cell])

# ---------- Grid of rods (17 x 17) ----------
water_holes = [
    (2, 5), (2, 8), (2, 11),
    (3, 3), (3, 13),
    (5, 2), (5, 5), (5, 8), (5, 11), (5, 14),
    (8, 2), (8, 5), (8, 8), (8, 11), (8, 14),
    (11, 2), (11, 5), (11, 8), (11, 11), (11, 14),
    (13, 3), (13, 13),
    (14, 5), (14, 8), (14, 11),
]

def make_assembly(fuel):
    pin = make_pin(fuel)
    hole = openmc.Universe(cells=[openmc.Cell(fill=water)])
    n = 17

    grid = [[pin] * n for _ in range(n)]
    for row, col in water_holes:
        grid[row][col] = hole

    lattice = openmc.RectLattice()
    lattice.lower_left = (-n * pitch / 2, -n * pitch / 2)
    lattice.pitch = (pitch, pitch)
    lattice.universes = grid

    box = openmc.model.RectangularPrism(n * pitch, n * pitch, boundary_type="reflective")
    root_cell = openmc.Cell(fill=lattice, region=-box)
    return openmc.Geometry([root_cell])






# ---------- Build both assemblies ----------
u_geom = make_assembly(u_fuel)
th_geom = make_assembly(th_fuel)

# ---------- Picture (side by side) ----------
width = 17 * pitch
colors = {u_fuel: "orange", th_fuel: "green", clad: "gray", water: "skyblue"}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))
u_geom.root_universe.plot(
    width=(width, width), pixels=(600, 600),
    color_by="material", colors=colors, axes=ax1,
)
th_geom.root_universe.plot(
    width=(width, width), pixels=(600, 600),
    color_by="material", colors=colors, axes=ax2,
)
ax1.set_title("Uranium assembly")
ax2.set_title("Thorium assembly")

os.makedirs("results/figures", exist_ok=True)
plt.savefig("results/figures/assembly_geometry.png", dpi=150)
print("Saved to results/figures/assembly_geometry.png")





# ---------- Simulation settings ----------
half = 17 * pitch / 2

settings = openmc.Settings()
settings.run_mode = "eigenvalue"
settings.particles = 10000
settings.batches = 100
settings.inactive = 20
settings.source = openmc.IndependentSource(
    space=openmc.stats.Box((-half, -half, -1), (half, half, 1)),
    constraints={"fissionable": True},
)

# # ---------- Heatmap tally ----------
# mesh = openmc.RegularMesh()
# mesh.dimension = (170, 170)
# mesh.lower_left = (-half, -half)
# mesh.upper_right = (half, half)

# mesh_tally = openmc.Tally(name="heatmap")
# mesh_tally.filters = [openmc.MeshFilter(mesh)]
# mesh_tally.scores = ["flux", "fission"]
# tallies = openmc.Tallies([mesh_tally])



# ---------- Heatmap tallies ----------
flux_mesh = openmc.RegularMesh()
flux_mesh.dimension = (85, 85)
flux_mesh.lower_left = (-half, -half)
flux_mesh.upper_right = (half, half)

rod_mesh = openmc.RegularMesh()
rod_mesh.dimension = (17, 17)
rod_mesh.lower_left = (-half, -half)
rod_mesh.upper_right = (half, half)

flux_tally = openmc.Tally(name="flux_map")
flux_tally.filters = [openmc.MeshFilter(flux_mesh)]
flux_tally.scores = ["flux"]

fission_tally = openmc.Tally(name="fission_map")
fission_tally.filters = [openmc.MeshFilter(rod_mesh)]
fission_tally.scores = ["fission"]

tallies = openmc.Tallies([flux_tally, fission_tally])




# ---------- Run both ----------
def run_case(geom, name):
    model = openmc.Model(geometry=geom, materials=materials, settings=settings, tallies=tallies)
    out_dir = f"03_assembly/output/{name}"
    os.makedirs(out_dir, exist_ok=True)
    sp_path = model.run(threads=8, cwd=out_dir)
    with openmc.StatePoint(sp_path) as sp:
        k = sp.keff
    print(f"{name}: k-eff = {k.nominal_value:.5f} +/- {k.std_dev:.5f}")
    return sp_path

u_sp = run_case(u_geom, "uranium")
th_sp = run_case(th_geom, "thorium")




# # ---------- Heatmaps (2 x 2) ----------
# def get_map(sp_path, score):
#     with openmc.StatePoint(sp_path) as sp:
#         t = sp.get_tally(name="heatmap")
#         values = t.get_values(scores=[score]).reshape(170, 170)
#     return values / values.max()

# extent = (-half, half, -half, half)
# fig, axes = plt.subplots(2, 2, figsize=(12, 12))

# cases = [(u_sp, "Uranium"), (th_sp, "Thorium")]
# scores = ["flux", "fission"]

# # this just makes it so teh fissile colors are dramaticized and look better
# for col, (sp_path, fuel_name) in enumerate(cases):
#     for row, score in enumerate(scores):
#         ax = axes[row][col]
#         values = get_map(sp_path, score)
#         if score == "fission":
#             values = np.ma.masked_equal(values, 0)
#             vmin = np.percentile(values.compressed(), 25)
#             cmap = plt.get_cmap("RdYlGn_r").copy()
#             cmap.set_bad("black")
#         else:
#             vmin = values.min()
#             cmap = "inferno"
#         im = ax.imshow(values, origin="lower", extent=extent, cmap=cmap, vmin=vmin, vmax=values.max())
#         fig.colorbar(im, ax=ax, shrink=0.8)
#         ax.set_title(f"{fuel_name}: {score}")

# plt.savefig("results/figures/assembly_heatmaps.png", dpi=150)
# print("Saved to results/figures/assembly_heatmaps.png")


# ---------- Heatmaps (2 x 2) ----------
def get_flux_map(sp_path):
    with openmc.StatePoint(sp_path) as sp:
        t = sp.get_tally(name="flux_map")
        values = t.get_values(scores=["flux"]).reshape(85, 85)
    return values / values.max()

def get_rod_power_map(sp_path):
    with openmc.StatePoint(sp_path) as sp:
        t = sp.get_tally(name="fission_map")
        values = t.get_values(scores=["fission"]).reshape(17, 17)
    values = np.ma.masked_equal(values, 0)
    return values / values.mean()

extent = (-half, half, -half, half)
fig, axes = plt.subplots(2, 2, figsize=(12, 12))
cases = [(u_sp, "Uranium"), (th_sp, "Thorium")]

for col, (sp_path, fuel_name) in enumerate(cases):
    flux = get_flux_map(sp_path)
    im = axes[0][col].imshow(flux, origin="lower", extent=extent, cmap="inferno")
    fig.colorbar(im, ax=axes[0][col], shrink=0.8)
    axes[0][col].set_title(f"{fuel_name}: flux")

    power = get_rod_power_map(sp_path)
    cmap = plt.get_cmap("RdYlGn_r").copy()
    cmap.set_bad("black")
    im = axes[1][col].imshow(power, origin="lower", extent=extent, cmap=cmap)
    fig.colorbar(im, ax=axes[1][col], shrink=0.8)
    axes[1][col].set_title(f"{fuel_name}: power per rod")

plt.savefig("results/figures/assembly_heatmaps.png", dpi=150)
print("Saved to results/figures/assembly_heatmaps.png")