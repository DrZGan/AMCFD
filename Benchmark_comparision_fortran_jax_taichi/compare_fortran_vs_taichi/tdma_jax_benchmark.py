"""
TDMA Solver Benchmark using JAX with jax.lax.linalg.tridiagonal_solve

Parallel j-loop (Jacobi style) using batched tridiagonal_solve.

Usage:
    CPU: python tdma_jax_benchmark.py
    GPU: python tdma_jax_benchmark.py --gpu
"""

import sys
import time
import os

os.environ['XLA_PYTHON_CLIENT_PREALLOCATE'] = 'false'
os.environ['JAX_PLATFORMS'] = 'cpu'

use_gpu = '--gpu' in sys.argv
if use_gpu:
    os.environ.pop('JAX_PLATFORMS', None)

# Parse --size N for grid size (default 100)
grid_size = 100
for idx, arg in enumerate(sys.argv):
    if arg == '--size' and idx + 1 < len(sys.argv):
        grid_size = int(sys.argv[idx + 1])

import jax
import jax.numpy as jnp
from functools import partial

# Grid dimensions - parameterized for benchmark
NI, NJ, NK = grid_size, grid_size, grid_size
NIM1, NJM1, NKM1 = NI - 1, NJ - 1, NK - 1
N_ITERATIONS = 50

arch_name = "JAX-GPU" if use_gpu else "JAX-CPU"

jax.config.update("jax_enable_x64", True)


def initialize_arrays():
    """Initialize all arrays"""
    enthalpy = jnp.ones((NI, NJ, NK), dtype=jnp.float64)
    ap = jnp.full((NI, NJ, NK), 6.0, dtype=jnp.float64)
    ae = jnp.ones((NI, NJ, NK), dtype=jnp.float64)
    aw = jnp.ones((NI, NJ, NK), dtype=jnp.float64)
    an = jnp.ones((NI, NJ, NK), dtype=jnp.float64)
    as_ = jnp.ones((NI, NJ, NK), dtype=jnp.float64)
    at = jnp.ones((NI, NJ, NK), dtype=jnp.float64)
    ab = jnp.ones((NI, NJ, NK), dtype=jnp.float64)
    su = jnp.full((NI, NJ, NK), 0.1, dtype=jnp.float64)
    return enthalpy, ap, ae, aw, an, as_, at, ab, su


@partial(jax.jit, static_argnums=(9,))
def tdma_solve_all_j_lines(enthalpy, ap, ae, aw, an, as_, at, ab, su, k):
    """
    Solve TDMA for all j-lines at a given k using jax.lax.linalg.tridiagonal_solve.
    Parallel j-loop (Jacobi style) - all j-lines solved simultaneously.
    """
    n_lines = NJM1 - 1  # j = 1 to NJM1-1
    n_unknowns = NIM1 - 1  # i = 1 to NIM1-1

    # JAX tridiagonal_solve expects all diagonals with same shape (n_lines, n_unknowns)
    # with dl[0]=0 and du[-1]=0

    # Diagonal: ap[i,j,k] for i=1..NIM1-1, j=1..NJM1-1
    diag = ap[1:NIM1, 1:NJM1, k].T  # (n_lines, n_unknowns)

    # Lower diagonal: -aw[i,j,k] for i=1..NIM1-1, with dl[:,0]=0
    # aw[i,j,k] represents coefficient for x[i-1], so dl[i] = -aw[i+1,j,k]
    dl = jnp.zeros((n_lines, n_unknowns), dtype=jnp.float64)
    dl = dl.at[:, 1:].set(-aw[2:NIM1, 1:NJM1, k].T)  # dl[:,0]=0, dl[:,1:] = -aw[2:,...]

    # Upper diagonal: -ae[i,j,k] for i=1..NIM1-1, with du[:,-1]=0
    du = jnp.zeros((n_lines, n_unknowns), dtype=jnp.float64)
    du = du.at[:, :-1].set(-ae[1:NIM1-1, 1:NJM1, k].T)  # du[:,:-1] = -ae[1:-1,...], du[:,-1]=0

    # RHS - uses current enthalpy values (Jacobi: same iteration values for all j)
    rhs = (at[1:NIM1, 1:NJM1, k] * enthalpy[1:NIM1, 1:NJM1, k + 1] +
           ab[1:NIM1, 1:NJM1, k] * enthalpy[1:NIM1, 1:NJM1, k - 1] +
           an[1:NIM1, 1:NJM1, k] * enthalpy[1:NIM1, 2:NJM1 + 1, k] +
           as_[1:NIM1, 1:NJM1, k] * enthalpy[1:NIM1, 0:NJM1 - 1, k] +
           su[1:NIM1, 1:NJM1, k]).T  # (n_lines, n_unknowns)

    # Add boundary contributions
    rhs = rhs.at[:, 0].add(aw[1, 1:NJM1, k] * enthalpy[0, 1:NJM1, k])
    rhs = rhs.at[:, -1].add(ae[NIM1 - 1, 1:NJM1, k] * enthalpy[NIM1, 1:NJM1, k])

    # Reshape rhs to (n_lines, n_unknowns, 1) as required by tridiagonal_solve
    rhs_3d = rhs[:, :, jnp.newaxis]

    # Solve
    solution = jax.lax.linalg.tridiagonal_solve(dl, diag, du, rhs_3d)

    # Update enthalpy - solution has shape (n_lines, n_unknowns, 1)
    enthalpy = enthalpy.at[1:NIM1, 1:NJM1, k].set(solution[:, :, 0].T)

    return enthalpy


def tdma_solve(enthalpy, ap, ae, aw, an, as_, at, ab, su):
    """Complete TDMA solve with k and j sweeps - parallel j (Jacobi)"""
    for ksweep in range(2):
        for k in range(NKM1 - 1, 0, -1):
            for jsweep in range(2):
                enthalpy = tdma_solve_all_j_lines(enthalpy, ap, ae, aw, an, as_, at, ab, su, k)
    return enthalpy


def main():
    print(f"JAX devices: {jax.devices()}", file=sys.stderr)

    enthalpy, ap, ae, aw, an, as_, at, ab, su = initialize_arrays()

    # Warmup
    enthalpy = tdma_solve(enthalpy, ap, ae, aw, an, as_, at, ab, su)
    enthalpy.block_until_ready()

    # Re-initialize
    enthalpy, ap, ae, aw, an, as_, at, ab, su = initialize_arrays()

    # Benchmark
    total_time = 0.0
    for _ in range(N_ITERATIONS):
        start = time.perf_counter()
        enthalpy = tdma_solve(enthalpy, ap, ae, aw, an, as_, at, ab, su)
        enthalpy.block_until_ready()
        end = time.perf_counter()
        total_time += (end - start)

    avg_time = (total_time / N_ITERATIONS) * 1000
    checksum = float(jnp.sum(enthalpy))
    max_val = float(jnp.max(enthalpy))
    min_val = float(jnp.min(enthalpy))
    avg_val = checksum / (NI * NJ * NK)

    print(f"{arch_name},{avg_time:12.6f},{checksum:20.12e},{max_val:20.12e},{min_val:20.12e},{avg_val:20.12e}")


if __name__ == "__main__":
    main()
