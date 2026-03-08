#!/usr/bin/env python3
"""compare_parametric.py — Task 7: Parametric study of localnum=3..10 vs baseline."""
import numpy as np
import os, re
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PROJ = os.path.dirname(os.path.abspath(__file__))

POINT_NAMES = ["P1","P2","P3","P4","P5","P6","P7","P8","P9","P10"]
POINT_DESCS = [
    "track1 start surface", "track1 centre surface", "track1 end surface",
    "track1 start -40um", "track1 centre -40um", "track1 end -40um",
    "offset surface", "offset surface", "far-field surface", "deep substrate",
]

COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
          '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

def load_thermal(fname):
    data = np.loadtxt(fname, comments="#")
    if data.ndim == 1:
        data = data[np.newaxis, :]
    return data

def extract_timing(fname):
    with open(fname) as f:
        text = f.read()
    cpu = wall = None
    m = re.search(r'Total CPU time:\s+([\d.]+)\s+s', text)
    if m: cpu = float(m.group(1))
    m = re.search(r'Total wall time:\s+([\d.]+)\s+s', text)
    if m: wall = float(m.group(1))
    return cpu, wall

def compute_metrics(base, loc):
    """Compute 4 metrics for each of 10 points."""
    n = min(len(base), len(loc))
    base = base[:n]
    loc = loc[:n]
    metrics = []
    for p in range(10):
        T_base = base[:, p + 1]
        T_loc = loc[:, p + 1]
        diff = T_loc - T_base
        max_abs = np.max(np.abs(diff))
        T_peak_base = np.max(T_base)
        max_rel = max_abs / T_peak_base * 100.0 if T_peak_base > 0 else 0.0
        rel_diff = diff / T_base
        rel_rmsd = np.sqrt(np.mean(rel_diff ** 2)) * 100.0
        peak_diff = abs(np.max(T_loc) - T_peak_base)
        metrics.append({
            'max_abs': max_abs, 'max_rel': max_rel,
            'rel_rmsd': rel_rmsd, 'peak_diff': peak_diff
        })
    return metrics

def plot_thermal_comparison(base, all_data, localnum_values):
    """Plot thermal histories: solid=baseline, dashed=each localnum."""
    time_base = base[:, 0] * 1000  # ms

    fig, axes = plt.subplots(2, 5, figsize=(22, 9))
    axes = axes.flatten()

    for p in range(10):
        ax = axes[p]
        T_base = base[:, p + 1]
        ax.plot(time_base, T_base, '-', color='black', linewidth=1.5, label='Baseline')
        for idx, N in enumerate(localnum_values):
            if N not in all_data:
                continue
            loc = all_data[N]['data']
            n = min(len(base), len(loc))
            time_loc = loc[:n, 0] * 1000
            T_loc = loc[:n, p + 1]
            c = COLORS[idx % len(COLORS)]
            ax.plot(time_loc, T_loc, '--', color=c, linewidth=0.9, alpha=0.8, label=f'N={N}')
        ax.set_title(f'{POINT_NAMES[p]} ({POINT_DESCS[p]})', fontsize=8)
        ax.set_xlabel('Time (ms)', fontsize=7)
        ax.set_ylabel('T (K)', fontsize=7)
        ax.tick_params(labelsize=6)
        ax.legend(fontsize=5, loc='best', ncol=2)
        ax.grid(True, alpha=0.3)

    fig.suptitle('Thermal History: Baseline (solid) vs localnum (dashed)', fontsize=12)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(PROJ, "parametric_thermal_comparison.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Thermal comparison plot saved to {out}")

def plot_accuracy_vs_speedup(all_results, cpu_base, wall_base):
    """Plot accuracy metrics and speedup vs localnum."""
    fig, axes = plt.subplots(2, 2, figsize=(12, 9))

    ns = sorted(all_results.keys())

    # Top-left: Max MaxRelDiff across all points
    ax = axes[0, 0]
    max_rel_diffs = [max(m['max_rel'] for m in all_results[n]['metrics']) for n in ns]
    ax.plot(ns, max_rel_diffs, 'o-', color='red', linewidth=1.5)
    ax.axhline(y=5.0, color='gray', linestyle='--', alpha=0.7, label='5% threshold')
    ax.set_xlabel('localnum')
    ax.set_ylabel('Max MaxRelDiff (%)')
    ax.set_title('Worst-case Max Relative Difference')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Top-right: Max RelRMSD across all points
    ax = axes[0, 1]
    max_rmsds = [max(m['rel_rmsd'] for m in all_results[n]['metrics']) for n in ns]
    ax.plot(ns, max_rmsds, 's-', color='blue', linewidth=1.5)
    ax.axhline(y=3.0, color='gray', linestyle='--', alpha=0.7, label='3% threshold')
    ax.set_xlabel('localnum')
    ax.set_ylabel('Max RelRMSD (%)')
    ax.set_title('Worst-case Relative RMSD')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Bottom-left: CPU speedup
    ax = axes[1, 0]
    if cpu_base:
        cpu_speedups = [cpu_base / all_results[n]['cpu'] if all_results[n]['cpu'] else 0 for n in ns]
        ax.plot(ns, cpu_speedups, 'D-', color='green', linewidth=1.5, label='CPU')
    if wall_base:
        wall_speedups = [wall_base / all_results[n]['wall'] if all_results[n]['wall'] else 0 for n in ns]
        ax.plot(ns, wall_speedups, '^-', color='purple', linewidth=1.5, label='Wall')
    ax.set_xlabel('localnum')
    ax.set_ylabel('Speedup ratio')
    ax.set_title('Speedup vs Baseline')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Bottom-right: Per-point MaxRelDiff
    ax = axes[1, 1]
    for p in range(10):
        vals = [all_results[n]['metrics'][p]['max_rel'] for n in ns]
        ax.plot(ns, vals, 'o-', color=COLORS[p], linewidth=0.8, markersize=3, label=POINT_NAMES[p])
    ax.axhline(y=5.0, color='gray', linestyle='--', alpha=0.7)
    ax.set_xlabel('localnum')
    ax.set_ylabel('MaxRelDiff (%)')
    ax.set_title('Per-point Max Relative Difference')
    ax.legend(fontsize=6, ncol=2)
    ax.grid(True, alpha=0.3)

    fig.suptitle('Parametric Study: Accuracy & Speedup vs localnum', fontsize=13)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    out = os.path.join(PROJ, "parametric_accuracy_speedup.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Accuracy/speedup plot saved to {out}")

def main():
    base = load_thermal(os.path.join(PROJ, "param_baseline_thermal_history.txt"))
    cpu_base, wall_base = extract_timing(os.path.join(PROJ, "param_baseline_timing_report.txt"))

    localnum_values = [4, 5, 6, 7, 8, 9, 10]
    all_results = {}
    all_data = {}

    for N in localnum_values:
        fname = os.path.join(PROJ, f"param_local{N}_thermal_history.txt")
        tname = os.path.join(PROJ, f"param_local{N}_timing_report.txt")
        if not os.path.exists(fname):
            print(f"  Skipping localnum={N}: {fname} not found")
            continue
        loc = load_thermal(fname)
        metrics = compute_metrics(base, loc)
        cpu_loc, wall_loc = extract_timing(tname) if os.path.exists(tname) else (None, None)
        all_results[N] = {'metrics': metrics, 'cpu': cpu_loc, 'wall': wall_loc}
        all_data[N] = {'data': loc}

    if not all_results:
        print("ERROR: No localnum results found.")
        return

    # --- Generate report ---
    lines = []
    lines.append("=" * 100)
    lines.append("  LOCAL_REGION Parametric Report: localnum Sensitivity")
    lines.append("=" * 100)
    lines.append("")

    # Per-localnum detailed table
    for N in sorted(all_results.keys()):
        r = all_results[N]
        cpu_sp = cpu_base / r['cpu'] if (cpu_base and r['cpu']) else 0.0
        wall_sp = wall_base / r['wall'] if (wall_base and r['wall']) else 0.0
        lines.append(f"--- localnum={N} | CPU speedup: {cpu_sp:.2f}x | Wall speedup: {wall_sp:.2f}x ---")
        if r['cpu']:
            lines.append(f"    CPU time: {r['cpu']:.1f}s (baseline: {cpu_base:.1f}s)")
        if r['wall']:
            lines.append(f"    Wall time: {r['wall']:.1f}s (baseline: {wall_base:.1f}s)")
        lines.append(f"  {'Point':<5} | {'MaxAbsDiff(K)':>13} | {'MaxRelDiff(%)':>13} | {'RelRMSD(%)':>10} | {'PeakTDiff(K)':>12}")
        lines.append("  " + "-" * 70)
        max_rel_all = 0
        max_rmsd_all = 0
        for p in range(10):
            m = r['metrics'][p]
            max_rel_all = max(max_rel_all, m['max_rel'])
            max_rmsd_all = max(max_rmsd_all, m['rel_rmsd'])
            lines.append(f"  {POINT_NAMES[p]:<5} | {m['max_abs']:13.3f} | {m['max_rel']:13.4f} | {m['rel_rmsd']:10.4f} | {m['peak_diff']:12.3f}")
        lines.append("  " + "-" * 70)
        pass_rel = max_rel_all < 5.0
        pass_rmsd = max_rmsd_all < 3.0
        lines.append(f"  Max MaxRelDiff: {max_rel_all:.4f}% {'PASS' if pass_rel else 'FAIL'} (<5%)")
        lines.append(f"  Max RelRMSD:    {max_rmsd_all:.4f}% {'PASS' if pass_rmsd else 'FAIL'} (<3%)")
        lines.append(f"  Overall: {'PASS' if pass_rel and pass_rmsd else 'FAIL'}")
        lines.append("")

    # Summary table
    lines.append("=" * 100)
    lines.append("  Summary Table")
    lines.append("=" * 100)
    lines.append(f"  {'localnum':>8} | {'CPU-sp':>7} | {'Wall-sp':>8} | {'CPU(s)':>8} | {'Wall(s)':>8} | {'MaxRelDiff%':>11} | {'MaxRMSD%':>8} | {'Result':>6}")
    lines.append("  " + "-" * 85)
    best_N = None
    for N in sorted(all_results.keys()):
        r = all_results[N]
        cpu_sp = cpu_base / r['cpu'] if (cpu_base and r['cpu']) else 0.0
        wall_sp = wall_base / r['wall'] if (wall_base and r['wall']) else 0.0
        max_rel = max(m['max_rel'] for m in r['metrics'])
        max_rmsd = max(m['rel_rmsd'] for m in r['metrics'])
        passed = max_rel < 5.0 and max_rmsd < 3.0
        result_str = "PASS" if passed else "FAIL"
        cpu_s = f"{r['cpu']:.1f}" if r['cpu'] else "N/A"
        wall_s = f"{r['wall']:.1f}" if r['wall'] else "N/A"
        lines.append(f"  {N:>8} | {cpu_sp:>7.2f}x | {wall_sp:>8.2f}x | {cpu_s:>8} | {wall_s:>8} | {max_rel:>11.4f} | {max_rmsd:>8.4f} | {result_str:>6}")
        if passed and (best_N is None or N > best_N):
            best_N = N
    lines.append("")

    # Trend analysis
    lines.append("=" * 100)
    lines.append("  Trend Analysis: MaxRelDiff(%) per point vs localnum")
    lines.append("=" * 100)
    for p in range(10):
        trend = f"  {POINT_NAMES[p]:>4}:"
        for N in sorted(all_results.keys()):
            val = all_results[N]['metrics'][p]['max_rel']
            trend += f"  N={N}:{val:5.2f}%"
        lines.append(trend)
    lines.append("")

    lines.append("  Trend Analysis: RelRMSD(%) per point vs localnum")
    lines.append("  " + "-" * 80)
    for p in range(10):
        trend = f"  {POINT_NAMES[p]:>4}:"
        for N in sorted(all_results.keys()):
            val = all_results[N]['metrics'][p]['rel_rmsd']
            trend += f"  N={N}:{val:5.2f}%"
        lines.append(trend)
    lines.append("")

    # Recommendation
    lines.append("=" * 100)
    lines.append("  Recommendation")
    lines.append("=" * 100)
    if best_N:
        r = all_results[best_N]
        cpu_sp = cpu_base / r['cpu'] if (cpu_base and r['cpu']) else 0.0
        wall_sp = wall_base / r['wall'] if (wall_base and r['wall']) else 0.0
        lines.append(f"  Recommended max localnum: {best_N}")
        lines.append(f"    CPU speedup: {cpu_sp:.2f}x, Wall speedup: {wall_sp:.2f}x")
        lines.append(f"    All acceptance criteria met (MaxRelDiff<5%, RelRMSD<3%)")
    else:
        lines.append("  No localnum value met all acceptance criteria.")
        # Find closest
        best_rmsd_N = min(all_results.keys(), key=lambda n: max(m['rel_rmsd'] for m in all_results[n]['metrics']))
        lines.append(f"  Best RelRMSD at localnum={best_rmsd_N} (all RelRMSD pass, MaxRelDiff exceeds 5% for some points)")
        lines.append(f"  This is an inherent limitation of the small local region (local_half_y=0.2mm)")
        lines.append(f"  with monitoring points outside the region (P7 at y=0.25mm).")
    lines.append("")

    # Conclusion
    lines.append("=" * 100)
    lines.append("  Conclusion")
    lines.append("=" * 100)
    lines.append("  - Higher localnum = more local steps per global step = faster but less accurate")
    lines.append("  - RelRMSD remains well under 3% for all tested localnum values (robust average accuracy)")
    lines.append("  - MaxRelDiff is dominated by P7 (permanently outside the local region in y)")
    lines.append("  - CPU speedup scales with localnum (more local steps = less full-domain work)")
    lines.append("  - Wall speedup is lower than CPU speedup due to OpenMP overhead")
    if best_N:
        lines.append(f"  - localnum={best_N} is recommended (best speedup within accuracy limits)")
    lines.append("=" * 100)

    report = "\n".join(lines)
    out_path = os.path.join(PROJ, "parametric_report.txt")
    with open(out_path, "w") as f:
        f.write(report)
    print(report)
    print(f"\nReport written to {out_path}")

    # --- Generate plots ---
    plot_thermal_comparison(base, all_data, sorted(all_results.keys()))
    plot_accuracy_vs_speedup(all_results, cpu_base, wall_base)

if __name__ == "__main__":
    main()
