"""
TDMA Solver Benchmark using JAX with jax.lax.linalg.tridiagonal_solve

Scheme B:
- Batch over BOTH k and j (Jacobi in k and j):
  Each global sweep uses enthalpy_old to build RHS for all interior k-planes,
  then solves all (k,j) i-lines in one batched tridiagonal_solve.

Usage:
    CPU: python tdma_jax_benchmark.py
    GPU: python tdma_jax_benchmark.py --gpu
    Size: python tdma_jax_benchmark.py --size 128
"""

import sys
import time
import os

# Set env BEFORE importing jax
os.environ["XLA_PYTHON_CLIENT_PREALLOCATE"] = "false"
os.environ["JAX_PLATFORMS"] = "cpu"

use_gpu = "--gpu" in sys.argv
if use_gpu:
    os.environ.pop("JAX_PLATFORMS", None)

# Parse --size N for grid size (default 100)
grid_size = 100
for idx, arg in enumerate(sys.argv):
    if arg == "--size" and idx + 1 < len(sys.argv):
        grid_size = int(sys.argv[idx + 1])

import jax
import jax.numpy as jnp

# Grid dimensions - parameterized for benchmark
NI, NJ, NK = grid_size, grid_size, grid_size
NIM1, NJM1, NKM1 = NI - 1, NJ - 1, NK - 1
N_ITERATIONS = 10

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


@jax.jit
def tdma_solve(enthalpy, ap, ae, aw, an, as_, at, ab, su):
    """
    Scheme B: batch over k and j (Jacobi in k and j).
    One "global sweep" updates all interior (i,j,k) using enthalpy_old everywhere.
    """
    # interior slices
    islc = slice(1, NIM1)   # i = 1..NIM1-1  (NI-2 unknowns)
    jslc = slice(1, NJM1)   # j = 1..NJM1-1  (NJ-2 lines)
    kslc = slice(1, NKM1)   # k = 1..NKM1-1  (NK-2 planes)

    # neighbor slices
    j_p1 = slice(2, NJM1 + 1)
    j_m1 = slice(0, NJM1 - 1)
    k_p1 = slice(2, NKM1 + 1)   # k+1
    k_m1 = slice(0, NKM1 - 1)   # k-1

    # Build diagonals for ALL k at once.
    # Shape convention for tridiagonal_solve:
    #   dl, d, du, b  all have shape (..., n)
    # We'll use (...)= (k, j) and n = i-unknowns.
    diag = ap[islc, jslc, kslc].transpose(2, 1, 0)  # (k, j, i)

    # Lower diag dl[...,0]=0, dl[...,1:]=-aw[i+1]
    aw_i = aw[2:NIM1, jslc, kslc].transpose(2, 1, 0)  # (k, j, i-1)
    zeros_1 = jnp.zeros(diag.shape[:-1] + (1,), dtype=diag.dtype)
    dl = jnp.concatenate([zeros_1, -aw_i], axis=-1)

    # Upper diag du[...,-1]=0, du[...,:-1]=-ae[i]
    ae_i = ae[1:NIM1 - 1, jslc, kslc].transpose(2, 1, 0)  # (k, j, i-1)
    du = jnp.concatenate([-ae_i, zeros_1], axis=-1)

    def one_global_sweep(h):
        h0 = h  # Jacobi: RHS uses old field everywhere

        rhs = (
            at[islc, jslc, kslc] * h0[islc, jslc, k_p1] +
            ab[islc, jslc, kslc] * h0[islc, jslc, k_m1] +
            an[islc, jslc, kslc] * h0[islc, j_p1, kslc] +
            as_[islc, jslc, kslc] * h0[islc, j_m1, kslc] +
            su[islc, jslc, kslc]
        ).transpose(2, 1, 0)  # (k, j, i)

        # i-boundary contributions (same as original logic)
        bc_lo = (aw[1, jslc, kslc] * h0[0, jslc, kslc]).transpose(1, 0)            # (k, j)
        bc_hi = (ae[NIM1 - 1, jslc, kslc] * h0[NIM1, jslc, kslc]).transpose(1, 0)  # (k, j)
        rhs = rhs.at[..., 0].add(bc_lo)
        rhs = rhs.at[..., -1].add(bc_hi)

        sol = jax.lax.linalg.tridiagonal_solve(dl, diag, du, rhs[..., None])  # (k,j,i,1)
        update = sol[..., 0].transpose(2, 1, 0)  # back to (i,j,k)

        return h0.at[islc, jslc, kslc].set(update)

    # Original did 2*(ksweep) * 2*(jsweep) passes, each pass covered all k planes sequentially.
    # Here we approximate that with 4 global Jacobi sweeps over all k planes.
    enthalpy = jax.lax.fori_loop(0, 4, lambda _, h: one_global_sweep(h), enthalpy)
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

    print(
        f"{arch_name},{avg_time:12.6f},{checksum:20.12e},"
        f"{max_val:20.12e},{min_val:20.12e},{avg_val:20.12e}"
    )


if __name__ == "__main__":
    main()
