"""
Pure NumPy TDMA - matching Fortran's 10 iterations
"""
import numpy as np

# Grid dimensions
ni, nj, nk = 40, 40, 10
nim1, njm1, nkm1 = ni - 1, nj - 1, nk - 1
n_iterations = 10

# Initialize arrays
enthalpy = np.ones((ni, nj, nk), dtype=np.float64)
ap = np.full((ni, nj, nk), 6.0, dtype=np.float64)
ae = np.ones((ni, nj, nk), dtype=np.float64)
aw = np.ones((ni, nj, nk), dtype=np.float64)
an = np.ones((ni, nj, nk), dtype=np.float64)
as_ = np.ones((ni, nj, nk), dtype=np.float64)
at = np.ones((ni, nj, nk), dtype=np.float64)
ab = np.ones((ni, nj, nk), dtype=np.float64)
su = np.full((ni, nj, nk), 0.1, dtype=np.float64)

def tdma_solve():
    """TDMA solve matching Fortran structure exactly"""
    pr = np.zeros(ni, dtype=np.float64)
    qr = np.zeros(ni, dtype=np.float64)
    
    for ksweep in range(1, 3):  # ksweep = 1, 2
        for k in range(nkm1 - 1, 0, -1):  # k = 8 to 1 (0-based)
            for jsweep in range(1, 3):  # jsweep = 1, 2
                for j in range(1, njm1):  # j = 1 to 38 (0-based)
                    pr[0] = 0.0
                    qr[0] = enthalpy[0, j, k]
                    
                    for i in range(1, nim1):  # i = 1 to 38 (0-based)
                        d = (at[i, j, k] * enthalpy[i, j, k + 1] +
                             ab[i, j, k] * enthalpy[i, j, k - 1] +
                             an[i, j, k] * enthalpy[i, j + 1, k] +
                             as_[i, j, k] * enthalpy[i, j - 1, k] +
                             su[i, j, k])
                        
                        denom = ap[i, j, k] - aw[i, j, k] * pr[i - 1]
                        
                        if denom <= 1e-12 and denom >= 0.0:
                            denom = denom + 1e-13
                        if denom >= -1e-12 and denom < 0.0:
                            denom = denom - 1e-13
                        
                        pr[i] = ae[i, j, k] / denom
                        qr[i] = (d + aw[i, j, k] * qr[i - 1]) / denom
                    
                    for i in range(nim1 - 1, 0, -1):  # i = 38 to 1 (0-based)
                        enthalpy[i, j, k] = pr[i] * enthalpy[i + 1, j, k] + qr[i]

# Warmup run (like Fortran)
tdma_solve()

# Re-initialize (like Fortran)
enthalpy = np.ones((ni, nj, nk), dtype=np.float64)

# Run 10 iterations (like Fortran benchmark)
for _ in range(n_iterations):
    tdma_solve()

# Compute stats
checksum = np.sum(enthalpy)
max_val = np.max(enthalpy)
min_val = np.min(enthalpy)
avg_val = checksum / (ni * nj * nk)

print(f"NumPy (10 iter): checksum={checksum:.12e}, max={max_val:.12e}, min={min_val:.12e}, avg={avg_val:.12e}")
