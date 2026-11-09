import os
import numpy as np
import matplotlib.pyplot as plt
import openmc

# ---------- Fixed sizes (cm) ----------
fuel_radius = 0.39
clad_radius = 0.46
u233_frac = 0.013

# ---------- One pin run ----------
def run_pin(fuel_temp=294, water_density=1.0, pitch=1.26, tag="base"):
    fuel = openmc.Material(name="ThO2 + U233 fuel")
    fuel.add_element("Th", 1.0 - u233_frac)
    fuel.add_nuclide("U233", u233_frac)
    fuel.add_element("O", 2.0)
    fuel.set_density("g/cm3", 10.0)
    fuel.temperature = fuel_temp

    clad = openmc.Material(name="Zirconium cladding")
    clad.add_element("Zr", 1.0)
    clad.set_density("g/cm3", 6.55)

    water = openmc.Material(name="Water")
    water.add_element("H", 2.0)
    water.add_element("O", 1.0)
    water.set_density("g/cm3", water_density)
    water.add_s_alpha_beta("c_H_in_H2O")

    materials = openmc.Materials([fuel, clad, water])

    fuel_surface = openmc.ZCylinder(r=fuel_radius)
    clad_surface = openmc.ZCylinder(r=clad_radius)
    box = openmc.model.RectangularPrism(pitch, pitch, boundary_type="reflective")

    fuel_cell = openmc.Cell(fill=fuel, region=-fuel_surface)
    clad_cell = openmc.Cell(fill=clad, region=+fuel_surface & -clad_surface)
    water_cell = openmc.Cell(fill=water, region=+clad_surface & -box)
    geometry = openmc.Geometry([fuel_cell, clad_cell, water_cell])

    settings = openmc.Settings()
    settings.run_mode = "eigenvalue"
    settings.particles = 5000
    settings.batches = 60
    settings.inactive = 10
    settings.source = openmc.IndependentSource(space=openmc.stats.Point())
    settings.temperature = {"method": "interpolation"}

    model = openmc.Model(geometry=geometry, materials=materials, settings=settings)
    out_dir = f"06_sensitivity/output/{tag}"
    os.makedirs(out_dir, exist_ok=True)
    sp_path = model.run(threads=8, cwd=out_dir, output=False)
    with openmc.StatePoint(sp_path) as sp:
        return sp.keff.nominal_value, sp.keff.std_dev

def sweep(name, values, **fixed):
    ks, errs = [], []
    for v in values:
        k, e = run_pin(**{name: v}, **fixed, tag=f"{name}_{v}")
        print(f"{name} = {v}  ->  k = {k:.4f} +/- {e:.4f}")
        ks.append(k)
        errs.append(e)
    return np.array(ks), np.array(errs)

# ---------- Sweeps ----------
temps = [294, 600, 900, 1200]
densities = [1.0, 0.9, 0.8, 0.7, 0.6]
pitches = [1.10, 1.20, 1.26, 1.40, 1.60]

k_T, e_T = sweep("fuel_temp", temps)
k_rho, e_rho = sweep("water_density", densities)
k_p, e_p = sweep("pitch", pitches)

# ---------- Coefficients ----------
def pcm(k1, k2):
    return (k2 - k1) / (k1 * k2) * 1e5

doppler = pcm(k_T[0], k_T[-1]) / (temps[-1] - temps[0])
print(f"\nDoppler coefficient: {doppler:.2f} pcm per K")
print(f"Water 1.0 -> 0.6 g/cm3: {pcm(k_rho[0], k_rho[-1]):.0f} pcm")
print(f"Pitch 1.10 -> 1.60 cm: {pcm(k_p[0], k_p[-1]):.0f} pcm")

# ---------- Plots ----------
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 5))

ax1.errorbar(temps, k_T, yerr=e_T, fmt="o-", color="red", capsize=3)
ax1.set_title("Fuel temperature")
ax1.set_xlabel("Fuel temperature (K)")
ax1.set_ylabel("k")

ax2.errorbar(densities, k_rho, yerr=e_rho, fmt="o-", color="blue", capsize=3)
ax2.set_title("Water density")
ax2.set_xlabel("Water density (g/cm3)")
ax2.set_ylabel("k")
ax2.invert_xaxis()

ax3.errorbar(pitches, k_p, yerr=e_p, fmt="o-", color="green", capsize=3)
ax3.axvline(1.26, color="gray", linestyle="--", label="base pitch")
ax3.set_title("Pin pitch (more water ->)")
ax3.set_xlabel("Pitch (cm)")
ax3.set_ylabel("k")
ax3.legend()

plt.tight_layout()
os.makedirs("results/figures", exist_ok=True)
plt.savefig("results/figures/sensitivity.png", dpi=150)
print("Saved to results/figures/sensitivity.png")