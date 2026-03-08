#!/usr/bin/env python3
"""compare_results.py — Task 6: Compare baseline vs localnum=3 thermal histories and timing."""
import numpy as np
import os, re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PROJ = os.path.dirname(os.path.abspath(__file__))

POINT_NAMES = [
    "P1  (1.0mm,0.50mm,0.695mm) track1 start surface",
    "P2  (2.0mm,0.50mm,0.695mm) track1 centre surface",
    "P3  (3.0mm,0.50mm,0.695mm) track1 end surface",
    "P4  (1.0mm,0.50mm,0.655mm) track1 start -40um",
    "P5  (2.0mm,0.50mm,0.655mm) track1 centre -40um",
    "P6  (3.0mm,0.50mm,0.655mm) track1 end -40um",
    "P7  (2.0mm,0.25mm,0.695mm) offset surface",
    "P8  (2.0mm,0.75mm,0.695mm) offset surface",
    "P9  (2.0mm,1.50mm,0.695mm) far-field surface",
    "P10 (2.0mm,0.50mm,0.200mm) deep substrate",
]

def load_thermal(fname):
    data = np.loadtxt(fname, comments="#")
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return data  # columns: time, T1..T10

def extract_timing_total(fname):
    with open(fname) as f:
        text = f.read()
    m = re.search(r'Total CPU time:\s+([\d.]+)\s+s', text)
    cpu = float(m.group(1)) if m else None
    m2 = re.search(r'Total wall time:\s+([\d.]+)\s+s', text)
    wall = float(m2.group(1)) if m2 else None
    return cpu, wall

SHORT_NAMES = ["P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9", "P10"]

COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

def plot_comparison(base, loc3):
    """Plot thermal histories: solid=baseline, dashed=local."""
    time_base = base[:, 0] * 1000  # convert to ms
    time_loc = loc3[:, 0] * 1000

    fig, axes = plt.subplots(2, 5, figsize=(20, 8))
    axes = axes.flatten()

    for p in range(10):
        ax = axes[p]
        T_base = base[:, p + 1]
        T_loc = loc3[:, p + 1]
        ax.plot(time_base, T_base, '-', color=COLORS[p], linewidth=1.2, label='Baseline')
        ax.plot(time_loc, T_loc, '--', color=COLORS[p], linewidth=1.2, label='Local (N=3)')
        ax.set_title(SHORT_NAMES[p], fontsize=10)
        ax.set_xlabel('Time (ms)', fontsize=8)
        ax.set_ylabel('Temperature (K)', fontsize=8)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=7, loc='best')
        ax.grid(True, alpha=0.3)

    fig.suptitle('Thermal History Comparison: Baseline (solid) vs Local (dashed)', fontsize=13)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out_png = os.path.join(PROJ, "thermal_history_comparison.png")
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"Comparison plot saved to {out_png}")

def main():
    base = load_thermal(os.path.join(PROJ, "baseline_thermal_history.txt"))
    loc3 = load_thermal(os.path.join(PROJ, "local3_thermal_history.txt"))

    # Align to same number of timesteps
    n = min(len(base), len(loc3))
    base = base[:n]
    loc3 = loc3[:n]

    lines = []
    lines.append("=" * 90)
    lines.append("  LOCAL_REGION Verification Report")
    lines.append("  Comparing localnum=0 (baseline) vs localnum=3")
    lines.append("=" * 90)
    lines.append("")

    # Header
    lines.append(f"{'Point':<6} | {'MaxAbsDiff(K)':>13} | {'MaxRelDiff(%)':>13} | {'RelRMSD(%)':>10} | {'PeakTDiff(K)':>12}")
    lines.append("-" * 72)

    max_rel_diff_all = 0.0
    max_rel_rmsd_all = 0.0

    for p in range(10):
        T_base = base[:, p + 1]
        T_loc = loc3[:, p + 1]
        diff = T_loc - T_base

        max_abs = np.max(np.abs(diff))
        T_peak_base = np.max(T_base)
        max_rel = max_abs / T_peak_base * 100.0 if T_peak_base > 0 else 0.0
        # Relative RMSD: sqrt(mean(((T_local - T_global) / T_global)^2))
        rel_diff = diff / T_base
        rel_rmsd = np.sqrt(np.mean(rel_diff ** 2)) * 100.0
        peak_diff = abs(np.max(T_loc) - T_peak_base)

        max_rel_diff_all = max(max_rel_diff_all, max_rel)
        max_rel_rmsd_all = max(max_rel_rmsd_all, rel_rmsd)

        lines.append(f"P{p+1:<5} | {max_abs:13.3f} | {max_rel:13.4f} | {rel_rmsd:10.4f} | {peak_diff:12.3f}")

    lines.append("-" * 72)
    lines.append("")

    # Acceptance
    pass_rel_diff = max_rel_diff_all < 5.0
    pass_rel_rmsd = max_rel_rmsd_all < 3.0
    all_pass = pass_rel_diff and pass_rel_rmsd

    lines.append("Acceptance criteria:")
    lines.append(f"  Max relative difference (max|dT|/T_peak): {max_rel_diff_all:.4f}% (threshold < 5%)  => {'PASS' if pass_rel_diff else 'FAIL'}")
    lines.append(f"  Max relative RMSD (rms(dT/T_base)):       {max_rel_rmsd_all:.4f}% (threshold < 3%)  => {'PASS' if pass_rel_rmsd else 'FAIL'}")
    lines.append(f"  Overall result: {'PASS' if all_pass else 'FAIL'}")
    lines.append("")

    # Timing comparison
    lines.append("=" * 90)
    lines.append("  Timing Comparison")
    lines.append("=" * 90)
    cpu_base, wall_base = extract_timing_total(os.path.join(PROJ, "baseline_timing_report.txt"))
    cpu_loc3, wall_loc3 = extract_timing_total(os.path.join(PROJ, "local3_timing_report.txt"))
    if cpu_base and cpu_loc3:
        speedup_cpu = cpu_base / cpu_loc3 if cpu_loc3 > 0 else float('inf')
        lines.append(f"  Baseline (localnum=0) total CPU time:  {cpu_base:.3f} s")
        lines.append(f"  Local    (localnum=3) total CPU time:  {cpu_loc3:.3f} s")
        lines.append(f"  CPU speedup ratio: {speedup_cpu:.2f}x")
        if wall_base and wall_loc3:
            speedup_wall = wall_base / wall_loc3 if wall_loc3 > 0 else float('inf')
            lines.append(f"  Baseline (localnum=0) total wall time: {wall_base:.3f} s")
            lines.append(f"  Local    (localnum=3) total wall time: {wall_loc3:.3f} s")
            lines.append(f"  Wall speedup ratio: {speedup_wall:.2f}x")
        if speedup_cpu <= 1.0:
            lines.append("  WARNING: No CPU speedup achieved!")
    else:
        lines.append("  Could not parse timing reports.")

    lines.append("")
    lines.append("=" * 90)

    # Timing report details side by side
    for label, fname in [("Baseline (localnum=0)", "baseline_timing_report.txt"),
                          ("Local (localnum=3)", "local3_timing_report.txt")]:
        fpath = os.path.join(PROJ, fname)
        if os.path.exists(fpath):
            lines.append(f"\n--- {label} timing_report.txt ---")
            with open(fpath) as f:
                lines.append(f.read())

    report = "\n".join(lines)
    out_path = os.path.join(PROJ, "local_region_verification.txt")
    with open(out_path, "w") as f:
        f.write(report)
    print(report)
    print(f"\nReport written to {out_path}")

    # Generate comparison PNG
    plot_comparison(base, loc3)

if __name__ == "__main__":
    main()
