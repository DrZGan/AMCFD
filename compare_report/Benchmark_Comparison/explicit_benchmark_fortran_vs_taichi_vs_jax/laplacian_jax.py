"""
Explicit Euler with Laplacian - JAX GPU Implementation
Usage: python laplacian_jax.py <grid_size>
"""
import jax
import jax.numpy as jnp
from jax import jit
import time
import numpy as np
import sys

print("JAX devices:", jax.devices())
print("JAX default backend:", jax.default_backend())

# Configuration
NSTEPS = 100
DX = 0.1
DT = 0.0001
ALPHA = 0.1

def main():
    # Enable 64-bit precision
    jax.config.update("jax_enable_x64", True)

    # Get grid size from command line
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 100

    # Print device info
    devices = jax.devices()
    print(f"JAX devices: {devices}", flush=True)
    print(f"Grid size: {N}x{N}x{N}", flush=True)

    coeff = ALPHA * DT / (DX * DX)

    # Initialize the field
    i = jnp.arange(1, N + 1, dtype=jnp.float64)
    j = jnp.arange(1, N + 1, dtype=jnp.float64)
    k = jnp.arange(1, N + 1, dtype=jnp.float64)
    I, J, K = jnp.meshgrid(i, j, k, indexing='ij')
    u_init = jnp.sin(I * 0.1) * jnp.sin(J * 0.1) * jnp.sin(K * 0.1)

    @jit
    def laplacian_step(u, coeff):
        laplacian = (
            u[2:, 1:-1, 1:-1] + u[:-2, 1:-1, 1:-1] +
            u[1:-1, 2:, 1:-1] + u[1:-1, :-2, 1:-1] +
            u[1:-1, 1:-1, 2:] + u[1:-1, 1:-1, :-2] -
            6.0 * u[1:-1, 1:-1, 1:-1]
        )
        u_new = u.at[1:-1, 1:-1, 1:-1].set(
            u[1:-1, 1:-1, 1:-1] + coeff * laplacian
        )
        return u_new

    @jit
    def run_simulation(u, coeff):
        def body_fn(i, u):
            return laplacian_step(u, coeff)
        return jax.lax.fori_loop(0, NSTEPS, body_fn, u)

    # Warm-up run (compile the JIT functions)
    _ = run_simulation(u_init, coeff)
    _ = _.block_until_ready()

    # Re-initialize for actual benchmark
    u = u_init

    # Benchmark
    start_time = time.perf_counter()
    u_final = run_simulation(u, coeff)
    u_final = u_final.block_until_ready()
    end_time = time.perf_counter()

    elapsed_time = end_time - start_time

    # Convert to numpy for verification
    u_np = np.array(u_final)

    result_sum = float(np.sum(u_np))
    result_min = float(np.min(u_np))
    result_max = float(np.max(u_np))
    result_avg = float(np.mean(u_np))

    print(f"CSV: JAX,{result_sum:.14e},{result_min:.14e},{result_max:.14e},{result_avg:.14e},{elapsed_time:.6f}")

if __name__ == "__main__":
    main()
