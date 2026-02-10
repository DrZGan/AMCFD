"""
TDMA Solver Benchmark using Taichi - CPU and GPU version

Parallel j-loop (Jacobi style) for GPU performance.

Usage:
    CPU: python tdma_taichi_benchmark.py <num_threads>
    GPU: python tdma_taichi_benchmark.py --gpu
"""

import sys
import time
import os

os.environ['TI_LOG_LEVEL'] = 'error'

if os.path.exists('/usr/lib/wsl/lib'):
    current_ld = os.environ.get('LD_LIBRARY_PATH', '')
    os.environ['LD_LIBRARY_PATH'] = '/usr/lib/wsl/lib:' + current_ld

use_gpu = '--gpu' in sys.argv
num_threads = 1
arch_name = "GPU"

# Parse --size N for grid size (default 100)
grid_size = 100
for idx, arg in enumerate(sys.argv):
    if arg == '--size' and idx + 1 < len(sys.argv):
        grid_size = int(sys.argv[idx + 1])

if use_gpu:
    sys.argv = [arg for arg in sys.argv if arg not in ('--gpu', '--size') and not arg.isdigit()]
else:
    # Find positional arg for threads (not --size value)
    args_clean = [a for i, a in enumerate(sys.argv[1:], 1)
                  if a not in ('--gpu', '--size') and
                  not (i > 1 and sys.argv[i-1] == '--size')]
    num_threads = int(args_clean[0]) if args_clean else 1
    arch_name = f"CPU({num_threads}t)"

import taichi as ti
import numpy as np

if use_gpu:
    ti.init(arch=ti.cuda, default_fp=ti.f64)
    arch_name = "Taichi-GPU"
else:
    ti.init(arch=ti.cpu, cpu_max_num_threads=num_threads)

# Grid dimensions - parameterized for benchmark
NI, NJ, NK = grid_size, grid_size, grid_size
NIM1, NJM1, NKM1 = NI - 1, NJ - 1, NK - 1
N_ITERATIONS = 50

# Allocate fields
enthalpy = ti.field(dtype=ti.f64, shape=(NI, NJ, NK))
ap = ti.field(dtype=ti.f64, shape=(NI, NJ, NK))
ae = ti.field(dtype=ti.f64, shape=(NI, NJ, NK))
aw = ti.field(dtype=ti.f64, shape=(NI, NJ, NK))
an = ti.field(dtype=ti.f64, shape=(NI, NJ, NK))
as_ = ti.field(dtype=ti.f64, shape=(NI, NJ, NK))
at = ti.field(dtype=ti.f64, shape=(NI, NJ, NK))
ab = ti.field(dtype=ti.f64, shape=(NI, NJ, NK))
su = ti.field(dtype=ti.f64, shape=(NI, NJ, NK))

# Per-thread pr/qr for parallel j-loop (Jacobi style)
pr = ti.field(dtype=ti.f64, shape=(NJ, NI))
qr = ti.field(dtype=ti.f64, shape=(NJ, NI))


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
def tdma_solve_k_jsweep(k: ti.i32):
    """TDMA solve for all j at given k - parallel j (Jacobi style)"""
    for j in range(1, NJM1):
        pr[j, 0] = 0.0
        qr[j, 0] = enthalpy[0, j, k]

        for i in range(1, NIM1):
            d = (at[i, j, k] * enthalpy[i, j, k + 1] +
                 ab[i, j, k] * enthalpy[i, j, k - 1] +
                 an[i, j, k] * enthalpy[i, j + 1, k] +
                 as_[i, j, k] * enthalpy[i, j - 1, k] +
                 su[i, j, k])

            denom = ap[i, j, k] - aw[i, j, k] * pr[j, i - 1]
            if denom <= 1e-12 and denom >= 0.0:
                denom = denom + 1e-13
            if denom >= -1e-12 and denom < 0.0:
                denom = denom - 1e-13

            pr[j, i] = ae[i, j, k] / denom
            qr[j, i] = (d + aw[i, j, k] * qr[j, i - 1]) / denom

        for ii in range(NIM1 - 1):
            i = NIM1 - 1 - ii
            enthalpy[i, j, k] = pr[j, i] * enthalpy[i + 1, j, k] + qr[j, i]


def tdma_solve():
    """Complete TDMA solve with parallel j-loop (Jacobi)"""
    for ksweep in range(2):
        for k in range(NKM1 - 1, 0, -1):
            for jsweep in range(2):
                tdma_solve_k_jsweep(k)
                ti.sync()


@ti.kernel
def compute_checksum() -> ti.f64:
    total = 0.0
    for i, j, k in ti.ndrange(NI, NJ, NK):
        total += enthalpy[i, j, k]
    return total


@ti.kernel
def compute_max() -> ti.f64:
    max_val = enthalpy[0, 0, 0]
    for i, j, k in ti.ndrange(NI, NJ, NK):
        ti.atomic_max(max_val, enthalpy[i, j, k])
    return max_val


@ti.kernel
def compute_min() -> ti.f64:
    min_val = enthalpy[0, 0, 0]
    for i, j, k in ti.ndrange(NI, NJ, NK):
        ti.atomic_min(min_val, enthalpy[i, j, k])
    return min_val


def main():
    initialize_arrays()
    ti.sync()

    tdma_solve()
    ti.sync()

    initialize_arrays()
    ti.sync()

    total_time = 0.0
    for _ in range(N_ITERATIONS):
        start = time.perf_counter()
        tdma_solve()
        ti.sync()
        end = time.perf_counter()
        total_time += (end - start)

    avg_time = (total_time / N_ITERATIONS) * 1000
    checksum = compute_checksum()
    max_val = compute_max()
    min_val = compute_min()
    avg_val = checksum / (NI * NJ * NK)

    print(f"{arch_name},{avg_time:12.6f},{checksum:20.12e},{max_val:20.12e},{min_val:20.12e},{avg_val:20.12e}")


if __name__ == "__main__":
    main()
