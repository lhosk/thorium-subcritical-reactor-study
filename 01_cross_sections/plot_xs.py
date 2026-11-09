import os
import matplotlib.pyplot as plt
import openmc.data

# Find the nuclear data using the label we set earlier
library = openmc.data.DataLibrary.from_xml(os.environ["OPENMC_CROSS_SECTIONS"])

def load(nuclide):
    path = library.get_by_material(nuclide)["path"]
    return openmc.data.IncidentNeutron.from_hdf5(path)

nuclides = ["U235", "U233", "Th232"]
reactions = {18: "Fission", 102: "Capture"}  # MT numbers = reaction IDs

fig, axes = plt.subplots(1, 2, figsize=(13, 5))

for ax, (mt, name) in zip(axes, reactions.items()):
    for nuc_name in nuclides:
        nuc = load(nuc_name)
        energy = nuc.energy["294K"]
        xs = nuc[mt].xs["294K"](energy)
        ax.loglog(energy, xs, label=nuc_name)
    ax.set_title(f"{name} cross section (room temp)")
    ax.set_xlabel("Neutron energy (eV)")
    ax.set_ylabel("Cross section (barns)")
    ax.set_xlim(1e-3, 2e7)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()

plt.tight_layout()
os.makedirs("results/figures", exist_ok=True)
plt.savefig("results/figures/cross_sections.png", dpi=150)
print("Saved to results/figures/cross_sections.png")