#!/usr/bin/env python3
"""Compare last time-step data from two AMCFD output.txt files."""
import re, sys

def parse_output(path):
    data = []
    try:
        lines = open(path).readlines()
    except FileNotFoundError:
        return None
    pat1 = re.compile(r"^\s*([\d.Ee+-]+)\s+(\d+)\s+([\d.]+)\s+(\d+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)")
    pat2 = re.compile(r"^\s*([\d.]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)\s+([\d.Ee+-]+)")
    for i, line in enumerate(lines):
        m1 = pat1.search(line)
        if m1:
            timet = float(m1.group(1))
            resorh = float(m1.group(5))
            resorm = float(m1.group(6))
            tpeak = float('nan')
            if i + 2 < len(lines):
                m2 = pat2.search(lines[i + 2])
                if m2:
                    tpeak = float(m2.group(1))
            data.append((timet, resorh, resorm, tpeak))
    return data

def main():
    if len(sys.argv) < 4:
        sys.exit("Usage: compare_outputs.py origin.txt new.txt report.txt")
    orig = parse_output(sys.argv[1])
    new = parse_output(sys.argv[2])
    out = []
    if not orig or not new:
        out.append("FAIL: missing data")
    else:
        o, n = orig[-1], new[-1]
        tpeak_rel = abs(n[3] - o[3]) / o[3] if o[3] > 100 else 0.0
        out.append("Accuracy (last time step):")
        out.append("  origin: timet={:.6e}  tpeak={:.2f}  resorh={:.2e}  resorm={:.2e}".format(o[0], o[3], o[1], o[2]))
        out.append("  new:    timet={:.6e}  tpeak={:.2f}  resorh={:.2e}  resorm={:.2e}".format(n[0], n[3], n[1], n[2]))
        out.append("  tpeak relative diff: {:.4f}%".format(tpeak_rel * 100))
        if abs(o[0] - n[0]) < 1e-10 and tpeak_rel < 0.05:
            out.append("  Result: PASS")
        else:
            out.append("  Result: CHECK")
    text = '\n'.join(out)
    print(text)
    if sys.argv[3] != '/dev/stdout':
        open(sys.argv[3], 'w').write(text)

if __name__ == "__main__":
    main()
