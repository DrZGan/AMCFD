#!/usr/bin/env python3
"""
Compare all output.txt data between origin and new: every time step block
(time, iter, time/iter, tot_iter, res_enth, res_mass, res_u, res_v, res_w,
 Tmax, umax, vmax, wmax, length, depth, width,
 north, south, top, toploss, bottom, west, east, hout, accu, hin, heatvol, ratio,
 time, beam_posx, beam_posy, beam_posz, power, scanspeed, speedx, speedy).
Usage: compare_output_full.py origin.txt new.txt [report.txt]
"""

import re
import sys

PAT1 = re.compile(
    r"^\s*([\d.Ee+-]+)\s+(\d+)\s+([\d.]+)\s+(\d+)\s+"
    r"([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s*$"
)
PAT2 = re.compile(
    r"^\s*([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+"
    r"([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s*$"
)
PAT3 = re.compile(
    r"^\s*([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+"
    r"([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+"
    r"([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)\s+([-\d.Ee+]+)\s*$"
)
PAT4 = re.compile(
    r"^\s*([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+"
    r"([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s*$"
)

NAMES1 = "time iter time/iter tot_iter res_enth res_mass res_u res_v res_w".split()
NAMES2 = "Tmax umax vmax wmax length depth width".split()
NAMES3 = "north south top toploss bottom west east hout accu hin heatvol ratio".split()
NAMES4 = "time beam_posx beam_posy beam_posz power scanspeed speedx speedy".split()


def parse_blocks(path):
    blocks = []
    try:
        lines = open(path).readlines()
    except FileNotFoundError:
        return None
    i = 0
    while i < len(lines):
        m1 = PAT1.search(lines[i])
        if m1:
            row = {}
            for k, name in enumerate(NAMES1):
                row[name] = int(m1.group(k + 1)) if name in ("iter", "tot_iter") else float(m1.group(k + 1))
            # Output has header lines between data lines: data1 at i, data2 at i+2, data3 at i+4, data4 at i+6
            if i + 2 < len(lines):
                m2 = PAT2.search(lines[i + 2])
                if m2:
                    for k, name in enumerate(NAMES2):
                        row[name] = float(m2.group(k + 1))
            if i + 4 < len(lines):
                m3 = PAT3.search(lines[i + 4])
                if m3:
                    for k, name in enumerate(NAMES3):
                        row[name] = float(m3.group(k + 1))
            if i + 6 < len(lines):
                m4 = PAT4.search(lines[i + 6])
                if m4:
                    for k, name in enumerate(NAMES4):
                        row[name] = float(m4.group(k + 1))
            blocks.append(row)
        i += 1
    return blocks


def compare_val(a, b, name, rel_tol=1e-2, abs_tol=1e-8):
    if isinstance(a, int) and isinstance(b, int):
        return 0.0, 0.0, a == b
    a, b = float(a), float(b)
    abs_d = abs(a - b)
    denom = abs(a) if abs(a) > 1e-30 else 1e-30
    rel_d = abs_d / denom
    ok = abs_d <= abs_tol or rel_d <= rel_tol
    return abs_d, rel_d, ok


def main():
    if len(sys.argv) < 3:
        sys.exit("Usage: compare_output_full.py origin.txt new.txt [report.txt]")
    path_o = sys.argv[1]
    path_n = sys.argv[2]
    report_path = sys.argv[3] if len(sys.argv) > 3 else None
    out = []
    def w(s):
        out.append(s)
        print(s)

    bo = parse_blocks(path_o)
    bn = parse_blocks(path_n)
    if not bo:
        w("FAIL: no data in origin")
        if report_path:
            open(report_path, "w").write("\n".join(out))
        return
    if not bn:
        w("FAIL: no data in new")
        if report_path:
            open(report_path, "w").write("\n".join(out))
        return

    n_o, n_n = len(bo), len(bn)
    w("output.txt comparison (all time-step blocks)")
    w("  Blocks: origin=%d  new=%d  %s" % (n_o, n_n, "OK" if n_o == n_n else "MISMATCH"))
    n_compare = min(n_o, n_n)
    all_names = NAMES1 + NAMES2 + NAMES3 + NAMES4
    max_diffs = {}
    for name in all_names:
        max_diffs[name] = (0.0, 0.0)
    fails = []

    for idx in range(n_compare):
        o, n = bo[idx], bn[idx]
        for name in all_names:
            if name not in o or name not in n:
                continue
            abs_d, rel_d, ok = compare_val(o[name], n[name], name)
            if abs_d > max_diffs[name][0] or rel_d > max_diffs[name][1]:
                max_diffs[name] = (max(abs_d, max_diffs[name][0]), max(rel_d, max_diffs[name][1]))
            if not ok and name not in ("iter", "tot_iter"):
                fails.append((idx, name, o[name], n[name]))

    # Averages over blocks and relative difference of averages (all per-fields)
    unique_names = list(dict.fromkeys(all_names))
    avg_origin = {}
    avg_new = {}
    rel_diff_avg = {}
    for name in unique_names:
        vals_o = [float(bo[i][name]) for i in range(n_compare) if name in bo[i]]
        vals_n = [float(bn[i][name]) for i in range(n_compare) if name in bn[i]]
        if vals_o and vals_n:
            avg_origin[name] = sum(vals_o) / len(vals_o)
            avg_new[name] = sum(vals_n) / len(vals_n)
    small = 1e-30
    w("  Average over blocks (origin vs new) and rel_diff of averages:")
    for name in unique_names:
        if name in avg_origin and name in avg_new:
            ao, an = avg_origin[name], avg_new[name]
            denom = abs(ao) if abs(ao) > small else small
            rel_diff_avg[name] = abs(an - ao) / denom
            w("    %-12s  avg_origin=%.6e  avg_new=%.6e  rel_diff_avg=%.3e" % (name, ao, an, rel_diff_avg[name]))
        else:
            w("    %-12s  avg_origin=N/A  avg_new=N/A  rel_diff_avg=N/A" % name)

    # Conclusion: pass or fail and why
    key_fields = ["time", "Tmax", "tot_iter", "length", "depth", "width", "north", "south", "top", "bottom", "west", "east", "hout", "accu", "hin", "heatvol", "ratio"]
    rel_tol_key = 0.1   # 10% for key physics/geometry/flux fields
    reasons = []
    if n_o != n_n:
        reasons.append("Block count mismatch (origin=%d, new=%d)" % (n_o, n_n))
    for name in key_fields:
        if name in rel_diff_avg and rel_diff_avg[name] > rel_tol_key:
            reasons.append("%s rel_diff_avg=%.3e (threshold %.2f)" % (name, rel_diff_avg[name], rel_tol_key))
    if reasons:
        w("")
        w("  Conclusion: FAIL")
        w("  Reason: " + "; ".join(reasons))
    else:
        w("")
        w("  Conclusion: PASS")
        w("  Reason: Block counts match and key fields (time, Tmax, tot_iter, geometry, fluxes) have rel_diff_avg <= %.0f%%. Residual and time/iter differences are acceptable (solver/implementation)." % (rel_tol_key * 100))

    if fails:
        w("  First few mismatches (block, field, origin, new):")
        for t in fails[:10]:
            w("    block %d  %s  %.6e  %.6e" % (t[0], t[1], t[2], t[3]))

    if report_path:
        open(report_path, "w").write("\n".join(out))


if __name__ == "__main__":
    main()
