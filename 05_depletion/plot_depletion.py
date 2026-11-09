import os
import matplotlib.pyplot as plt
import openmc.deplete

# ---------- Load results ----------
results = openmc.deplete.Results("05_depletion/output/depletion_results.h5")
mat_id = list(results[0].index_mat.keys())[0]

days, k = results.get_keff(time_units="d")
_, th232 = results.get_atoms(mat_id, "Th232", time_units="d")
_, pa233 = results.get_atoms(mat_id, "Pa233", time_units="d")
_, u233 = results.get_atoms(mat_id, "U233", time_units="d")

# ---------- Plots ----------
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 5))

ax1.plot(days, th232, "o-", color="green")
ax1.set_title("Th-232 (used up)")
ax1.set_xlabel("Time (days)")
ax1.set_ylabel("Atoms (per cm of pin)")

ax2.plot(days, u233, "o-", color="blue", label="U-233 (fuel)")
ax2.plot(days, pa233, "o-", color="orange", label="Pa-233 (waiting room)")
ax2.set_title("Breeding")
ax2.set_xlabel("Time (days)")
ax2.set_ylabel("Atoms (per cm of pin)")
ax2.legend()

ax3.errorbar(days, k[:, 0], yerr=k[:, 1], fmt="o-", color="red", capsize=3)
ax3.set_title("k vs time")
ax3.set_xlabel("Time (days)")
ax3.set_ylabel("k")

plt.tight_layout()
os.makedirs("results/figures", exist_ok=True)
plt.savefig("results/figures/depletion.png", dpi=150)
print("Saved to results/figures/depletion.png")

# ---------- Quick numbers ----------
print(f"k: {k[0, 0]:.4f} -> {k[-1, 0]:.4f}")
print(f"U-233: {u233[0]:.3e} -> {u233[-1]:.3e} atoms")
print(f"Pa-233 at end: {pa233[-1]:.3e} atoms")