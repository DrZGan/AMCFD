"""
TDMA Solver using Taichi - CPU and GPU versions.

Run with:
    python tdma_taichi.py --cpu      # CPU backend
    python tdma_taichi.py --gpu      # GPU backend (CUDA/Vulkan)
    python tdma_taichi.py            # Default: CPU
"""

import taichi as ti
import numpy as np
import time
import argparse

# Parse command line arguments
parser = argparse.ArgumentParser(description='TDMA Solver Benchmark')
parser.add_argument('--cpu', action='store_true', help='Use CPU backend')
parser.add_argument('--gpu', action='store_true', help='Use GPU backend')
parser.add_argument('--threads', type=int, default=8, help='Number of CPU threads')
args = parser.parse_args()

# Initialize Taichi backend
if args.gpu:
    try:
        ti.init(arch=ti.cuda)
        backend_name = "CUDA GPU"
    except:
        try:
            ti.init(arch=ti.vulkan)
            backend_name = "Vulkan GPU"
        except:
            ti.init(arch=ti.cpu, cpu_max_num_threads=args.threads)
            backend_name = f"CPU (fallback, {args.threads} threads)"
else:
    ti.init(arch=ti.cpu, cpu_max_num_threads=args.threads)
    backend_name = f"CPU ({args.threads} threads)"

print(f"Backend: {backend_name}")

# Grid dimensions
NI, NJ, NK = 128, 128, 128
NIM1, NJM1, NKM1 = NI - 1, NJ - 1, NK - 1
N_ITERATIONS = 10

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

# Temporary arrays for TDMA
pr = ti.field(dtype=ti.f64, shape=(NK, NJ, NI))
qr = ti.field(dtype=ti.f64, shape=(NK, NJ, NI))


@ti.kernel
def initialize_arrays():
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


@ti.kernel
def tdma_sweep_k(k: ti.i32):
    """
    TDMA sweep for a single k-plane.
    Parallel over j, sequential over i (TDMA dependency).
    """
    for j in range(1, NJM1):
        # Initialize
        pr[k, j, 0] = 0.0
        qr[k, j, 0] = enthalpy[0, j, k]

        # Forward elimination
        for i in range(1, NIM1):
            # Load coefficients
            at_ijk = at[i, j, k]
            ab_ijk = ab[i, j, k]
            an_ijk = an[i, j, k]
            as_ijk = as_[i, j, k]
            ae_ijk = ae[i, j, k]
            aw_ijk = aw[i, j, k]
            ap_ijk = ap[i, j, k]
            su_ijk = su[i, j, k]

            d = (at_ijk * enthalpy[i, j, k + 1] +
                 ab_ijk * enthalpy[i, j, k - 1] +
                 an_ijk * enthalpy[i, j + 1, k] +
                 as_ijk * enthalpy[i, j - 1, k] +
                 su_ijk)

            denom = ap_ijk - aw_ijk * pr[k, j, i - 1]

            # Avoid divide by zero
            if ti.abs(denom) < 1e-12:
                if denom >= 0:
                    denom = denom + 1e-13
                else:
                    denom = denom - 1e-13

            inv_denom = 1.0 / denom
            pr[k, j, i] = ae_ijk * inv_denom
            qr[k, j, i] = (d + aw_ijk * qr[k, j, i - 1]) * inv_denom

        # Back substitution
        for i in range(NIM1 - 1, 0, -1):
            enthalpy[i, j, k] = pr[k, j, i] * enthalpy[i + 1, j, k] + qr[k, j, i]


def tdma_solve():
    """Complete TDMA solve with k and j sweeps."""
    for ksweep in range(2):
        for k in range(NKM1 - 1, 0, -1):
            for jsweep in range(2):
                tdma_sweep_k(k)


@ti.kernel
def compute_checksum() -> ti.f64:
    total = 0.0
    for i, j, k in ti.ndrange(NI, NJ, NK):
        total += enthalpy[i, j, k]
    return total


def main():
    print("Initializing arrays...")
    initialize_arrays()
    ti.sync()

    # Warmup
    print("Warmup run...")
    tdma_solve()
    ti.sync()

    # Benchmark
    print("\nBenchmarking...")
    times = []

    # Re-initialize
    initialize_arrays()
    ti.sync()

    for iter in range(N_ITERATIONS):
        start = time.perf_counter()
        tdma_solve()
        ti.sync()
        end = time.perf_counter()
        elapsed = (end - start) * 1000
        times.append(elapsed)
        print(f"Iteration {iter + 1}: {elapsed:.6f} ms")

    print("-" * 40)
    print(f"Backend: {backend_name}")
    print(f"Grid size: {NI} x {NJ} x {NK}")
    print(f"Number of iterations: {N_ITERATIONS}")
    print(f"Average time per iteration: {np.mean(times):.6f} ms")
    print(f"Total time: {np.sum(times):.6f} ms")
    print(f"Min time: {np.min(times):.6f} ms")
    print(f"Max time: {np.max(times):.6f} ms")

    checksum = compute_checksum()
    print(f"Checksum (sum of enthalpy): {checksum:.8e}")


if __name__ == "__main__":
    main()
