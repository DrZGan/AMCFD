#!/usr/bin/env python3
"""Compare BINARY legacy VTK files from two result directories.
For each file, report sum of T, Velocity, den, diff, solidID, vis for origin and new;
then conclusion (pass/fail) and reason."""
import struct, re, sys, os

FIELDS = ['T', 'Velocity', 'den', 'diff', 'solidID', 'vis']
REL_TOL = 0.1   # 10% relative difference threshold for pass

def read_vtk(path):
    with open(path, 'rb') as f:
        raw = f.read()
    idx = raw.find(b'POINTS ')
    if idx < 0:
        return None
    end = raw.find(b'\n', idx)
    m = re.match(r'POINTS\s+(\d+)\s+float', raw[idx:end].decode('ascii'))
    if not m:
        return None
    npts = int(m.group(1))
    ps = end + 1
    points = list(struct.unpack_from('>' + 'f' * (3 * npts), raw, ps))
    datasets = {}
    pos = raw.find(b'POINT_DATA ')
    if pos < 0:
        return dict(npts=npts, points=points, datasets=datasets)
    pos = raw.find(b'\n', pos) + 1
    while pos < len(raw):
        if raw[pos:pos+1] != b'\n':
            pos += 1
            continue
        pos += 1
        if raw[pos:pos+8] == b'VECTORS ':
            el = raw.find(b'\n', pos)
            name = raw[pos:el].decode('ascii').split()[1]
            pos = el + 1
            if pos + 3 * npts * 4 <= len(raw):
                datasets[name] = list(struct.unpack_from('>' + 'f' * (3 * npts), raw, pos))
            pos += 3 * npts * 4
        elif raw[pos:pos+8] == b'SCALARS ':
            el = raw.find(b'\n', pos)
            name = raw[pos:el].decode('ascii').split()[1]
            pos = el + 1
            if raw[pos:pos+20] == b'LOOKUP_TABLE default':
                pos = raw.find(b'\n', pos) + 1
            if pos + npts * 4 <= len(raw):
                datasets[name] = list(struct.unpack_from('>' + 'f' * npts, raw, pos))
            pos += npts * 4
        else:
            pos += 1
    return dict(npts=npts, points=points, datasets=datasets)

def sum_field(data):
    return sum(data)

def sum_velocity_magnitude(data):
    """Velocity is 3*npts: vx,vy,vz per point. Return sum of magnitudes."""
    n = len(data) // 3
    total = 0.0
    for i in range(n):
        vx, vy, vz = data[3*i], data[3*i+1], data[3*i+2]
        total += (vx*vx + vy*vy + vz*vz) ** 0.5
    return total

def main():
    d1, d2 = sys.argv[1], sys.argv[2]
    common = sorted(set(f for f in os.listdir(d1) if f.endswith('.vtk')) & set(f for f in os.listdir(d2) if f.endswith('.vtk')))
    if not common:
        print("No common VTK files found.")
        return
    print("VTK comparison: sums of T, Velocity, den, diff, solidID, vis (origin vs new)")
    print("")
    small = 1e-30
    failures = []
    for f in common:
        a = read_vtk(os.path.join(d1, f))
        b = read_vtk(os.path.join(d2, f))
        if not a or not b:
            print("  %s: FAIL (parse)" % f)
            failures.append("%s: parse error" % f)
            continue
        if a['npts'] != b['npts']:
            print("  %s: npts mismatch %d vs %d" % (f, a['npts'], b['npts']))
            failures.append("%s: npts mismatch" % f)
            continue
        print("  %s  npts=%d" % (f, a['npts']))
        for name in FIELDS:
            if name not in a['datasets'] or name not in b['datasets']:
                print("      %s: missing in one or both" % name)
                failures.append("%s %s: missing" % (f, name))
                continue
            va, vb = a['datasets'][name], b['datasets'][name]
            if name == 'Velocity':
                sa, sb = sum_velocity_magnitude(va), sum_velocity_magnitude(vb)
            else:
                sa, sb = sum_field(va), sum_field(vb)
            denom = abs(sa) if abs(sa) > small else small
            rel = abs(sb - sa) / denom
            print("      %-8s  sum_origin=%.6e  sum_new=%.6e  rel_diff=%.3e" % (name, sa, sb, rel))
            if rel > REL_TOL:
                failures.append("%s %s rel_diff=%.3e (threshold %.2f)" % (f, name, rel, REL_TOL))
        print("")
    # Conclusion
    if failures:
        print("Conclusion: FAIL")
        print("Reason: " + "; ".join(failures))
    else:
        print("Conclusion: PASS")
        print("Reason: All VTK files have matching npts and sums of T, Velocity, den, diff, solidID, vis within %.0f%% relative difference." % (REL_TOL * 100))

if __name__ == '__main__':
    main()
