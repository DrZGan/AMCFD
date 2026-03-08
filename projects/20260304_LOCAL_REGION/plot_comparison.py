#!/usr/bin/env python3
"""Plot thermal history comparison: baseline vs localnum=3."""
import numpy as np
import os

PROJ = os.path.dirname(os.path.abspath(__file__))

POINT_NAMES = [
    "P1 (1.0,0.50,0.695)", "P2 (2.0,0.50,0.695)", "P3 (3.0,0.50,0.695)",
    "P4 (1.0,0.50,0.655)", "P5 (2.0,0.50,0.655)", "P6 (3.0,0.50,0.655)",
    "P7 (2.0,0.25,0.695)", "P8 (2.0,0.75,0.695)",
    "P9 (2.0,1.50,0.695)", "P10 (2.0,0.50,0.200)",
]

def load(fname):
    data = np.loadtxt(fname, comments="#")
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return data

def main():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available, skipping plot generation")
        return

    base_file = os.path.join(PROJ, "baseline_thermal_history.txt")
    loc3_file = os.path.join(PROJ, "local3_thermal_history.txt")
    if not os.path.exists(base_file) or not os.path.exists(loc3_file):
        print("Missing thermal history files")
        return

    base = load(base_file)
    loc3 = load(loc3_file)
    n = min(len(base), len(loc3))
    base = base[:n]
    loc3 = loc3[:n]
    t = base[:, 0] * 1000  # convert to ms

    # Plot 1: Temperature histories (2x5 grid)
    fig, axes = plt.subplots(2, 5, figsize=(20, 8))
    fig.suptitle("Thermal History: Baseline (localnum=0) vs Local (localnum=3)", fontsize=14)
    for p in range(10):
        ax = axes[p // 5, p % 5]
        ax.plot(t, base[:, p+1], 'b-', linewidth=0.8, label='Baseline')
        ax.plot(t, loc3[:, p+1], 'r--', linewidth=0.8, label='Local3')
        ax.set_title(POINT_NAMES[p], fontsize=8)
        ax.set_xlabel("Time (ms)", fontsize=7)
        ax.set_ylabel("T (K)", fontsize=7)
        ax.tick_params(labelsize=6)
        if p == 0:
            ax.legend(fontsize=6)
    plt.tight_layout()
    out1 = os.path.join(PROJ, "thermal_history_comparison.png")
    plt.savefig(out1, dpi=150)
    print(f"Saved {out1}")
    plt.close()

    # Plot 2: Temperature differences
    fig, axes = plt.subplots(2, 5, figsize=(20, 8))
    fig.suptitle("Temperature Difference (Local3 - Baseline)", fontsize=14)
    for p in range(10):
        ax = axes[p // 5, p % 5]
        diff = loc3[:, p+1] - base[:, p+1]
        ax.plot(t, diff, 'k-', linewidth=0.8)
        ax.axhline(0, color='gray', linewidth=0.5)
        ax.set_title(POINT_NAMES[p], fontsize=8)
        ax.set_xlabel("Time (ms)", fontsize=7)
        ax.set_ylabel("dT (K)", fontsize=7)
        ax.tick_params(labelsize=6)
    plt.tight_layout()
    out2 = os.path.join(PROJ, "thermal_history_diff.png")
    plt.savefig(out2, dpi=150)
    print(f"Saved {out2}")
    plt.close()

if __name__ == "__main__":
    main()
