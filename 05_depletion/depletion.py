import os
import math
import openmc
import openmc.deplete

openmc.deplete.pool.USE_MULTIPROCESSING = False

# ---------- Chain file ----------
openmc.config["chain_file"] = os.path.expanduser("~/nuclear_data/chain_endfb80_pwr.xml")

# ---------- Materials ----------
fuel = openmc.Material(name="ThO2 + U233 fuel")
fuel.add_nuclide("Th232", 0.987)
fuel.add_nuclide("U233", 0.013)
fuel.add_element("O", 2.0)
fuel.set_density("g/cm3", 10.0)
fuel.depletable = True

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

fuel.volume = math.pi * fuel_radius**2

materials = openmc.Materials([fuel, clad, water])

# ---------- Geometry ----------
fuel_surface = openmc.ZCylinder(r=fuel_radius)
clad_surface = openmc.ZCylinder(r=clad_radius)
box = openmc.model.RectangularPrism(pitch, pitch, boundary_type="reflective")

fuel_cell = openmc.Cell(fill=fuel, region=-fuel_surface)
clad_cell = openmc.Cell(fill=clad, region=+fuel_surface & -clad_surface)
water_cell = openmc.Cell(fill=water, region=+clad_surface & -box)

geometry = openmc.Geometry([fuel_cell, clad_cell, water_cell])

print("Fuel volume (cm^2 per cm of height):", fuel.volume)

# ---------- Simulation settings ----------
settings = openmc.Settings()
settings.run_mode = "eigenvalue"
settings.particles = 1000
settings.batches = 30
settings.inactive = 10
settings.source = openmc.IndependentSource(space=openmc.stats.Point())

model = openmc.Model(geometry=geometry, materials=materials, settings=settings)

# ---------- Depletion ----------
os.makedirs("05_depletion/output", exist_ok=True)
os.chdir("05_depletion/output")

operator = openmc.deplete.CoupledOperator(model)

timesteps = [1, 4, 5, 20, 30, 30, 90, 185]  # days
power = 200.0  # watts per cm of pin

integrator = openmc.deplete.PredictorIntegrator(
    operator, timesteps, power=power, timestep_units="d"
)
integrator.integrate()
print("Done. Results in 05_depletion/output/depletion_results.h5")