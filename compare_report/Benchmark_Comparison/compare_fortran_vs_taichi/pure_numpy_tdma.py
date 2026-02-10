"""
Pure NumPy TDMA to verify algorithm matches Fortran
"""
import numpy as np

# Grid dimensions - same as Fortran (1-based will be handled in loops)
ni, nj, nk = 40, 40, 10
nim1, njm1, nkm1 = ni - 1, nj - 1, nk - 1

# Arrays - use 1-based indexing by making arrays 1 element larger
# Actually, let's keep 0-based but adjust loop ranges carefully
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
        # Fortran: k = nkm1 to 2, -1 (1-based: 9 to 2)
        # 0-based: k = 8 to 1
        for k in range(nkm1 - 1, 0, -1):  # nkm1-1 = 8, goes to 1
            for jsweep in range(1, 3):  # jsweep = 1, 2
                # Fortran: j = 2 to njm1 (1-based: 2 to 39)
                # 0-based: j = 1 to 38
                for j in range(1, njm1):  # j = 1 to 38
                    # Fortran: pr(1) = 0, qr(1) = enthalpy(1,j,k)
                    # 0-based: pr[0] = 0, qr[0] = enthalpy[0,j,k]
                    pr[0] = 0.0
                    qr[0] = enthalpy[0, j, k]
                    
                    # Forward elimination
                    # Fortran: i = 2 to nim1 (1-based: 2 to 39)
                    # 0-based: i = 1 to 38
                    for i in range(1, nim1):  # i = 1 to 38
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
                    
                    # Back substitution
                    # Fortran: i = nim1 to 2, -1 (1-based: 39 to 2)
                    # 0-based: i = 38 to 1
                    for i in range(nim1 - 1, 0, -1):  # i = 38 to 1
                        enthalpy[i, j, k] = pr[i] * enthalpy[i + 1, j, k] + qr[i]

# Run TDMA
tdma_solve()

# Compute stats
checksum = np.sum(enthalpy)
max_val = np.max(enthalpy)
min_val = np.min(enthalpy)
avg_val = checksum / (ni * nj * nk)

print(f"NumPy TDMA: checksum={checksum:.12e}, max={max_val:.12e}, min={min_val:.12e}, avg={avg_val:.12e}")
