"""
JAX TDMA - Tridiagonal Matrix Algorithm.
Follows the exact loop logic from the provided Fortran code.
"""

# import time
# # os.environ["CUDA_VISIBLE_DEVICES"] = "2"

# import os
# os.environ["JAX_PLATFORM_NAME"] = "cpu"
# os.environ["CUDA_VISIBLE_DEVICES"] = ""   # 禁用 GPU（可选但推荐）

# import jax
# import jax.numpy as jnp
# from jax import lax
# print("backend:", jax.default_backend())
# print("devices:", jax.devices())

import os, time
os.environ["JAX_PLATFORMS"] = "cpu"      
# os.environ["CUDA_VISIBLE_DEVICES"] = ""

import jax
import jax.numpy as jnp
from jax import lax

print("backend:", jax.default_backend())
print("devices:", jax.devices())



# Grid dimensions (nim1, njm1, nkm1 from original Fortran)
NI, NJ, NK = 500, 500, 500

# Use ghost cells to match Fortran-style 1-based indexing and neighbor access
GI, GJ, GK = NI + 2, NJ + 2, NK + 2


def init_random(key):
    """Initialize coefficients and enthalpy with random values."""
    keys = jax.random.split(key, 9)
    aw = 0.1 * jax.random.uniform(keys[0], (GI, GJ, GK), dtype=jnp.float32)
    ae = 0.1 * jax.random.uniform(keys[1], (GI, GJ, GK), dtype=jnp.float32)
    an = 0.1 * jax.random.uniform(keys[2], (GI, GJ, GK), dtype=jnp.float32)
    as_ = 0.1 * jax.random.uniform(keys[3], (GI, GJ, GK), dtype=jnp.float32)
    at = 0.1 * jax.random.uniform(keys[4], (GI, GJ, GK), dtype=jnp.float32)
    ab = 0.1 * jax.random.uniform(keys[5], (GI, GJ, GK), dtype=jnp.float32)
    su = 0.1 * jax.random.uniform(keys[6], (GI, GJ, GK), dtype=jnp.float32)
    ap = 1.0 + 0.2 * jax.random.uniform(keys[7], (GI, GJ, GK), dtype=jnp.float32)
    enthalpy = jax.random.uniform(keys[8], (GI, GJ, GK), dtype=jnp.float32)
    return aw, ae, ap, an, as_, at, ab, su, enthalpy


def tdma_solve(aw, ae, ap, an, as_, at, ab, su, enthalpy):
    """Thomas algorithm: forward sweep then backward substitution."""
    pr = jnp.zeros_like(enthalpy)
    qr = jnp.zeros_like(enthalpy)

    def k_body(kk, state):
        enthalpy, pr, qr = state
        k = NK - kk

        def j_body(j, state):
            enthalpy, pr, qr = state
            pr = pr.at[1, j, k].set(0.0)
            qr = qr.at[1, j, k].set(enthalpy[1, j, k])

            def i_body(i, state):
                enthalpy, pr, qr = state
                d = (
                    at[i, j, k] * enthalpy[i, j, k + 1]
                    + ab[i, j, k] * enthalpy[i, j, k - 1]
                    + an[i, j, k] * enthalpy[i, j + 1, k]
                    + as_[i, j, k] * enthalpy[i, j - 1, k]
                    + su[i, j, k]
                )
                denom = ap[i, j, k] - aw[i, j, k] * pr[i - 1, j, k]
                denom = jnp.where(jnp.abs(denom) < 1e-12, denom + 1e-13, denom)
                pr = pr.at[i, j, k].set(ae[i, j, k] / denom)
                qr = qr.at[i, j, k].set((d + aw[i, j, k] * qr[i - 1, j, k]) / denom)
                return enthalpy, pr, qr

            enthalpy, pr, qr = lax.fori_loop(2, NI + 1, i_body, (enthalpy, pr, qr))

            def back_body(idx, enthalpy):
                i = NI - idx
                return enthalpy.at[i, j, k].set(
                    pr[i, j, k] * enthalpy[i + 1, j, k] + qr[i, j, k]
                )

            enthalpy = lax.fori_loop(0, NI - 1, back_body, enthalpy)
            return enthalpy, pr, qr

        enthalpy, pr, qr = lax.fori_loop(2, NJ + 1, j_body, (enthalpy, pr, qr))
        return enthalpy, pr, qr

    enthalpy, pr, qr = lax.fori_loop(0, NK - 1, k_body, (enthalpy, pr, qr))
    return enthalpy


tdma_solve_jit = jax.jit(tdma_solve)


if __name__ == "__main__":
    key = jax.random.PRNGKey(0)
    aw, ae, ap, an, as_, at, ab, su, enthalpy = init_random(key)

    t0 = time.perf_counter()
    enthalpy = tdma_solve_jit(aw, ae, ap, an, as_, at, ab, su, enthalpy)
    enthalpy.block_until_ready()
    t1 = time.perf_counter()

    elapsed_ms = (t1 - t0) * 1000
    print(f"Backend: {jax.default_backend()}")
    print(f"TDMA solve time: {elapsed_ms:.4f} ms")

    result = jnp.asarray(enthalpy)
    print(f"\nEnthalpy shape: {result.shape}")
    print("Enthalpy sample (k=NK, j=2, i=1..10):")
    print(result[1:11, 2, NK])
    print("\nEnthalpy statistics (all):")
    print(f"  min: {result.min():.6f}, max: {result.max():.6f}")
    print(f"  mean: {result.mean():.6f}, std: {result.std():.6f}")
