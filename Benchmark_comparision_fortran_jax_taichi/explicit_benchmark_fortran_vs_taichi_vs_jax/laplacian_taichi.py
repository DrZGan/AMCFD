"""
Explicit Euler with Laplacian - Taichi GPU Implementation
Usage: python laplacian_taichi.py <grid_size>
"""
import os
if os.path.exists('/usr/lib/wsl/lib'):
    os.environ['LD_LIBRARY_PATH'] = '/usr/lib/wsl/lib:' + os.environ.get('LD_LIBRARY_PATH', '')

import taichi as ti
import numpy as np
import time
import sys

ti.init(arch=ti.cuda)  # 或 ti.vulkan
print("Taichi arch:", ti.cfg.arch)
print("Taichi debug:", ti.cfg.debug)

# Configuration
NSTEPS = 100
DX = 0.1
DT = 0.0001
ALPHA = 0.1

# Get grid size from command line
N = int(sys.argv[1]) if len(sys.argv) > 1 else 100

# Initialize Taichi with CUDA backend
ti.init(arch=ti.cuda, default_fp=ti.f64)

# Allocate fields
u = ti.field(dtype=ti.f64, shape=(N, N, N))
u_new = ti.field(dtype=ti.f64, shape=(N, N, N))

@ti.kernel
def initialize(n: ti.i32):
    for i, j, k in u:
        u[i, j, k] = ti.sin((i + 1) * 0.1) * ti.sin((j + 1) * 0.1) * ti.sin((k + 1) * 0.1)
        u_new[i, j, k] = u[i, j, k]

@ti.kernel
def laplacian_step(coeff: ti.f64, n: ti.i32):
    for i, j, k in u:
        if 1 <= i < n - 1 and 1 <= j < n - 1 and 1 <= k < n - 1:
            laplacian = (
                u[i + 1, j, k] + u[i - 1, j, k] +
                u[i, j + 1, k] + u[i, j - 1, k] +
                u[i, j, k + 1] + u[i, j, k - 1] -
                6.0 * u[i, j, k]
            )
            u_new[i, j, k] = u[i, j, k] + coeff * laplacian
        else:
            u_new[i, j, k] = u[i, j, k]

@ti.kernel
def swap_arrays():
    for i, j, k in u:
        u[i, j, k] = u_new[i, j, k]

@ti.kernel
def compute_sum() -> ti.f64:
    total = 0.0
    for i, j, k in u:
        total += u[i, j, k]
    return total

@ti.kernel
def compute_min() -> ti.f64:
    min_val = u[0, 0, 0]
    for i, j, k in u:
        ti.atomic_min(min_val, u[i, j, k])
    return min_val

@ti.kernel
def compute_max() -> ti.f64:
    max_val = u[0, 0, 0]
    for i, j, k in u:
        ti.atomic_max(max_val, u[i, j, k])
    return max_val

def main():
    coeff = ALPHA * DT / (DX * DX)

    # Initialize
    initialize(N)
    ti.sync()

    # Warm-up run
    for _ in range(2):
        laplacian_step(coeff, N)
        swap_arrays()
    ti.sync()

    # Re-initialize for actual benchmark
    initialize(N)
    ti.sync()

    # Benchmark
    start_time = time.perf_counter()
    for step in range(NSTEPS):
        laplacian_step(coeff, N)
        swap_arrays()
    ti.sync()
    end_time = time.perf_counter()

    elapsed_time = end_time - start_time

    result_sum = compute_sum()
    result_min = compute_min()
    result_max = compute_max()
    result_avg = result_sum / (N * N * N)

    print(f"CSV: Taichi,{result_sum:.14e},{result_min:.14e},{result_max:.14e},{result_avg:.14e},{elapsed_time:.6f}")

if __name__ == "__main__":
    main()
