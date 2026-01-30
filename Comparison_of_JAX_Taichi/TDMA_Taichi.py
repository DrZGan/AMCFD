"""
Taichi TDMA - Tridiagonal Matrix Algorithm on GPU.
Follows the exact loop logic from the provided Fortran code.
"""

import taichi as ti
import time, os

ti.init(arch=ti.cpu)

# os.environ["CUDA_VISIBLE_DEVICES"] = "1"
# ti.init(arch=ti.cuda)

# Grid dimensions (nim1, njm1, nkm1 from original Fortran)
NI, NJ, NK = 1000, 1000, 1000

# Use ghost cells to match Fortran-style 1-based indexing and neighbor access
GI, GJ, GK = NI + 2, NJ + 2, NK + 2

# Coefficients and fields (3D)
aw = ti.field(dtype=ti.f32, shape=(GI, GJ, GK))
ae = ti.field(dtype=ti.f32, shape=(GI, GJ, GK))
ap = ti.field(dtype=ti.f32, shape=(GI, GJ, GK))
an = ti.field(dtype=ti.f32, shape=(GI, GJ, GK))
as_ = ti.field(dtype=ti.f32, shape=(GI, GJ, GK))
at = ti.field(dtype=ti.f32, shape=(GI, GJ, GK))
ab = ti.field(dtype=ti.f32, shape=(GI, GJ, GK))
su = ti.field(dtype=ti.f32, shape=(GI, GJ, GK))
enthalpy = ti.field(dtype=ti.f32, shape=(GI, GJ, GK))

# TDMA temp arrays
pr = ti.field(dtype=ti.f32, shape=(GI, GJ, GK))
qr = ti.field(dtype=ti.f32, shape=(GI, GJ, GK))


@ti.kernel
def init_random():
    """Initialize coefficients and enthalpy with random values."""
    for i, j, k in ti.ndrange(GI, GJ, GK):
        aw[i, j, k] = 0.1 * ti.random()
        ae[i, j, k] = 0.1 * ti.random()
        an[i, j, k] = 0.1 * ti.random()
        as_[i, j, k] = 0.1 * ti.random()
        at[i, j, k] = 0.1 * ti.random()
        ab[i, j, k] = 0.1 * ti.random()
        su[i, j, k] = 0.1 * ti.random()
        ap[i, j, k] = 1.0 + 0.2 * ti.random()  # diagonal dominant
        enthalpy[i, j, k] = ti.random()


@ti.kernel
def tdma_solve():
    # do k=nkm1,2,-1 (1-based)
    for kk in range(NK - 1):
        k = NK - kk
        if k <= 1:
            continue
        # do j=2,njm1
        for j in range(2, NJ + 1):
            # pr(1)=0.0, qr(1)=enthalpy(1,j,k)
            pr[1, j, k] = 0.0
            qr[1, j, k] = enthalpy[1, j, k]
            # do i=2,nim1
            for i in range(2, NI + 1):
                d = (
                    at[i, j, k] * enthalpy[i, j, k + 1]
                    + ab[i, j, k] * enthalpy[i, j, k - 1]
                    + an[i, j, k] * enthalpy[i, j + 1, k]
                    + as_[i, j, k] * enthalpy[i, j - 1, k]
                    + su[i, j, k]
                )
                denom = ap[i, j, k] - aw[i, j, k] * pr[i - 1, j, k]
                if ti.abs(denom) < 1e-12:
                    denom = denom + 1e-13
                pr[i, j, k] = ae[i, j, k] / denom
                qr[i, j, k] = (d + aw[i, j, k] * qr[i - 1, j, k]) / denom
            # -----back
            for idx in range(NI - 1):
                i = NI - idx
                enthalpy[i, j, k] = pr[i, j, k] * enthalpy[i + 1, j, k] + qr[i, j, k]


if __name__ == "__main__":
    init_random()
    ti.sync()

    t0 = time.perf_counter()
    tdma_solve()
    ti.sync()
    t1 = time.perf_counter()

    elapsed_ms = (t1 - t0) * 1000
    print(f"TDMA solve time: {elapsed_ms:.4f} ms")

    result = enthalpy.to_numpy()
    print(f"\nEnthalpy shape: {result.shape}")
    print("Enthalpy sample (k=NK, j=2, i=1..10):")
    print(result[1:11, 2, NK])
    print("\nEnthalpy statistics (all):")
    print(f"  min: {result.min():.6f}, max: {result.max():.6f}")
    print(f"  mean: {result.mean():.6f}, std: {result.std():.6f}")
