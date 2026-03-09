import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

data = np.loadtxt("thermal_history.txt", comments="#")
if data.ndim == 1: data = data[np.newaxis, :]
t = data[:, 0] * 1e3  # ms
labels = [
    "P1: track1 start surface",
    "P2: track1 centre surface",
    "P3: track1 end surface",
    "P4: track1 centre 40um",
    "P5: track1 centre 100um",
    "P6: 100um offset surface",
    "P7: midtrack surface",
    "P8: track2 centre surface",
    "P9: far substrate surface",
    "P10: deep substrate",
]
fig, ax = plt.subplots(figsize=(12, 6))
for i in range(10):
    ax.plot(t, data[:, i+1], label=labels[i])
ax.axhline(1563, color="gray", linestyle="--", linewidth=0.8, label="T_solid")
ax.axhline(2650, color="red",  linestyle="--", linewidth=0.8, label="T_boiling")
ax.set_xlabel("Time (ms)")
ax.set_ylabel("Temperature (K)")
ax.set_title("Thermal History - 10 Monitoring Points")
ax.legend(fontsize=7, loc="upper right")
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("thermal_history.png", dpi=150)
print("Saved thermal_history.png")
