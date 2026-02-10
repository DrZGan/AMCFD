"""
Debug test to compare Fortran vs Taichi TDMA algorithm
"""
import os
os.environ['TI_LOG_LEVEL'] = 'error'
if os.path.exists('/usr/lib/wsl/lib'):
    os.environ['LD_LIBRARY_PATH'] = '/usr/lib/wsl/lib:' + os.environ.get('LD_LIBRARY_PATH', '')

import taichi as ti
import numpy as np

ti.init(arch=ti.cpu, cpu_max_num_threads=1)

# Small grid for debugging
NI, NJ, NK = 40, 40, 10
NIM1, NJM1, NKM1 = NI - 1, NJ - 1, NK - 1

# Allocate Taichi fields
enthalpy = ti.field(dtype=ti.f64, shape=(NI, NJ, NK))
ap = ti.field(dtype=ti.f64, shape=(NI, NJ, NK))
ae = ti.field(dtype=ti.f64, shape=(NI, NJ, NK))
aw = ti.field(dtype=ti.f64, shape=(NI, NJ, NK))
an = ti.field(dtype=ti.f64, shape=(NI, NJ, NK))
as_ = ti.field(dtype=ti.f64, shape=(NI, NJ, NK))
at = ti.field(dtype=ti.f64, shape=(NI, NJ, NK))
ab = ti.field(dtype=ti.f64, shape=(NI, NJ, NK))
su = ti.field(dtype=ti.f64, shape=(NI, NJ, NK))

# TDMA temp arrays - use 1D per line
pr = ti.field(dtype=ti.f64, shape=(NI,))
qr = ti.field(dtype=ti.f64, shape=(NI,))

@ti.kernel
def initialize():
    for i, j, k in ti.ndrange(NI, NJ, NK):
        enthalpy[i, j, k] = 1.0
        ap[i, j, k] = 6.0
        ae[i, j, k] = 1.0
        aw[i, j, k] = 1.0
        an[i, j, k] = 1.0
        as_[i, j, k] = 1.0
        at[i, j, k] = 1.0
        ab[i, j, k] = 1.0
        su[i, j, k] = 0.1

def tdma_solve_line_python(k, j):
    """Pure Python TDMA solve for one line - matching Fortran exactly"""
    # Get numpy arrays
    enth = enthalpy.to_numpy()
    ap_np = ap.to_numpy()
    ae_np = ae.to_numpy()
    aw_np = aw.to_numpy()
    an_np = an.to_numpy()
    as_np = as_.to_numpy()
    at_np = at.to_numpy()
    ab_np = ab.to_numpy()
    su_np = su.to_numpy()
    
    pr_local = np.zeros(NI)
    qr_local = np.zeros(NI)
    
    # Fortran uses 1-based indexing: pr(1), qr(1)
    # In 0-based: pr[0], qr[0]
    pr_local[0] = 0.0
    qr_local[0] = enth[0, j, k]
    
    # Forward elimination: Fortran i = 2 to nim1
    # In 0-based: i = 1 to NI-2
    for i in range(1, NI - 1):
        d = (at_np[i, j, k] * enth[i, j, k + 1] +
             ab_np[i, j, k] * enth[i, j, k - 1] +
             an_np[i, j, k] * enth[i, j + 1, k] +
             as_np[i, j, k] * enth[i, j - 1, k] +
             su_np[i, j, k])
        
        denom = ap_np[i, j, k] - aw_np[i, j, k] * pr_local[i - 1]
        
        if denom <= 1e-12 and denom >= 0.0:
            denom = denom + 1e-13
        if denom >= -1e-12 and denom < 0.0:
            denom = denom - 1e-13
        
        pr_local[i] = ae_np[i, j, k] / denom
        qr_local[i] = (d + aw_np[i, j, k] * qr_local[i - 1]) / denom
    
    # Back substitution: Fortran i = nim1 to 2, -1
    # In 0-based: i = NI-2 to 1, -1
    for i in range(NI - 2, 0, -1):
        enth[i, j, k] = pr_local[i] * enth[i + 1, j, k] + qr_local[i]
    
    # Copy back to Taichi field
    enthalpy.from_numpy(enth)

def tdma_solve_python():
    """Complete TDMA solve matching Fortran structure"""
    for ksweep in range(2):
        # Fortran: k = nkm1 to 2, -1
        # nkm1 = nk - 1 = 10 - 1 = 9 (1-based)
        # In 0-based: k = 8 to 1
        for k in range(NK - 2, 0, -1):
            for jsweep in range(2):
                # Fortran: j = 2 to njm1
                # njm1 = nj - 1 = 40 - 1 = 39 (1-based)
                # In 0-based: j = 1 to 38
                for j in range(1, NJ - 1):
                    tdma_solve_line_python(k, j)

# Initialize
initialize()
ti.sync()

# Run TDMA
tdma_solve_python()

# Compute stats
enth = enthalpy.to_numpy()
checksum = np.sum(enth)
max_val = np.max(enth)
min_val = np.min(enth)
avg_val = checksum / (NI * NJ * NK)

print(f"Taichi (Python TDMA): checksum={checksum:.12e}, max={max_val:.12e}, min={min_val:.12e}, avg={avg_val:.12e}")
