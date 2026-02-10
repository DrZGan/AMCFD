"""
Benchmark Runner - Runs all grid sizes and writes benchmark_results.txt
"""
import subprocess
import os
import sys

WORK_DIR = os.path.dirname(os.path.abspath(__file__))
GRID_SIZES = [50, 100, 200, 300, 500]
THREAD_COUNTS = list(range(1, 32, 2))  # 1, 3, 5, ..., 31
NSTEPS = 100
RESULT_FILE = os.path.join(WORK_DIR, "benchmark_results.txt")

JAX_GPU_ID = 1
TAICHI_GPU_ID = 2

def compile_fortran():
    print("Compiling Fortran code...")
    result = subprocess.run(
        ["gfortran", "-O3", "-fopenmp", "-o",
         os.path.join(WORK_DIR, "laplacian_fortran"),
         os.path.join(WORK_DIR, "laplacian_fortran.f90")],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        print(f"Compilation failed: {result.stderr}")
        return False
    print("Compilation successful!")
    return True

def parse_fortran_csv(line):
    parts = line.strip().split(",")
    if len(parts) >= 6:
        return {
            'threads': int(parts[0]),
            'sum': float(parts[1]),
            'min': float(parts[2]),
            'max': float(parts[3]),
            'avg': float(parts[4]),
            'time': float(parts[5])
        }
    return None

def parse_python_csv(stdout):
    for line in stdout.split('\n'):
        if line.startswith("CSV:"):
            parts = line.split(",")
            if len(parts) >= 6:
                return {
                    'sum': float(parts[1]),
                    'min': float(parts[2]),
                    'max': float(parts[3]),
                    'avg': float(parts[4]),
                    'time': float(parts[5])
                }
    return None

def run_fortran(grid_size):
    results = []
    for threads in THREAD_COUNTS:
        print(f"  Fortran {threads:>2} threads ...", end=" ", flush=True)
        r = subprocess.run(
            [os.path.join(WORK_DIR, "laplacian_fortran"), str(threads), str(grid_size)],
            capture_output=True, text=True, timeout=600
        )
        if r.returncode != 0:
            print(f"FAILED")
            continue
        parsed = parse_fortran_csv(r.stdout)
        if parsed:
            print(f"{parsed['time']:.6f}s")
            results.append(parsed)
    return results

def run_jax(grid_size):
    print(f"  JAX GPU ...", end=" ", flush=True)
    r = subprocess.run(
        [sys.executable, os.path.join(WORK_DIR, "laplacian_jax.py"), str(grid_size)],
        capture_output=True, text=True, timeout=600
    )
    
    env = os.environ.copy()
    env["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"   # 让编号更稳定（可选）
    env["CUDA_VISIBLE_DEVICES"] = str(JAX_GPU_ID)
    env["JAX_PLATFORMS"] = "gpu"              # 可选：强制只初始化 gpu backend
    
    if r.returncode != 0:
        print(f"FAILED: {r.stderr[-200:]}")
        return None
    parsed = parse_python_csv(r.stdout)
    if parsed:
        print(f"{parsed['time']:.6f}s")
    return parsed

def run_taichi(grid_size):
    print(f"  Taichi GPU ...", end=" ", flush=True)
    env = os.environ.copy()
    
    env = os.environ.copy()
    env["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"   # 可选
    env["CUDA_VISIBLE_DEVICES"] = str(TAICHI_GPU_ID)
    
    if os.path.exists('/usr/lib/wsl/lib'):
        env['LD_LIBRARY_PATH'] = '/usr/lib/wsl/lib:' + env.get('LD_LIBRARY_PATH', '')
    r = subprocess.run(
        [sys.executable, os.path.join(WORK_DIR, "laplacian_taichi.py"), str(grid_size)],
        capture_output=True, text=True, timeout=600, env=env
    )
    if r.returncode != 0:
        print(f"FAILED: {r.stderr[-200:]}")
        return None
    parsed = parse_python_csv(r.stdout)
    if parsed:
        print(f"{parsed['time']:.6f}s")
    return parsed

def write_section(f, section_num, grid_size, fortran_results, jax_result, taichi_result):
    N = grid_size
    sep = "=" * 90
    dash = "-" * 90

    f.write(f"\n{sep}\n")
    f.write(f"  SECTION {section_num}: Grid Size {N}x{N}x{N}  ({N**3:,} cells, {NSTEPS} time steps)\n")
    f.write(f"{sep}\n\n")

    # --- Verification ---
    f.write(f"  VERIFICATION\n")
    f.write(f"  {dash}\n")
    f.write(f"  {'Implementation':<28} {'Sum':>22} {'Min':>22} {'Max':>22} {'Average':>22}\n")
    f.write(f"  {dash}\n")

    if fortran_results:
        ref = fortran_results[0]
        f.write(f"  {'Fortran (reference)':<28} {ref['sum']:>22.14e} {ref['min']:>22.14e} {ref['max']:>22.14e} {ref['avg']:>22.14e}\n")
    if jax_result:
        f.write(f"  {'JAX GPU':<28} {jax_result['sum']:>22.14e} {jax_result['min']:>22.14e} {jax_result['max']:>22.14e} {jax_result['avg']:>22.14e}\n")
    if taichi_result:
        f.write(f"  {'Taichi GPU':<28} {taichi_result['sum']:>22.14e} {taichi_result['min']:>22.14e} {taichi_result['max']:>22.14e} {taichi_result['avg']:>22.14e}\n")

    f.write(f"\n")

    # --- Differences from reference ---
    if fortran_results:
        ref = fortran_results[0]
        f.write(f"  DIFFERENCES FROM FORTRAN REFERENCE\n")
        f.write(f"  {'-'*70}\n")
        f.write(f"  {'Implementation':<28} {'Sum diff':>18} {'Min diff':>18} {'Max diff':>18} {'Avg diff':>18}\n")
        f.write(f"  {'-'*70}\n")
        if jax_result:
            f.write(f"  {'JAX GPU':<28} {abs(jax_result['sum']-ref['sum']):>18.6e} {abs(jax_result['min']-ref['min']):>18.6e} {abs(jax_result['max']-ref['max']):>18.6e} {abs(jax_result['avg']-ref['avg']):>18.6e}\n")
        if taichi_result:
            f.write(f"  {'Taichi GPU':<28} {abs(taichi_result['sum']-ref['sum']):>18.6e} {abs(taichi_result['min']-ref['min']):>18.6e} {abs(taichi_result['max']-ref['max']):>18.6e} {abs(taichi_result['avg']-ref['avg']):>18.6e}\n")
        f.write(f"\n")

    # --- Computational Time ---
    f.write(f"  COMPUTATIONAL TIME\n")
    f.write(f"  {dash}\n")

    ref_time = fortran_results[0]['time'] if fortran_results else None

    f.write(f"  {'Implementation':<28} {'Time (s)':>12} {'Speedup vs 1-thread':>22}\n")
    f.write(f"  {dash}\n")

    if fortran_results:
        for r in fortran_results:
            speedup = ref_time / r['time'] if r['time'] > 0 else 0
            f.write(f"  Fortran ({r['threads']:>2} threads)        {r['time']:>12.6f} {speedup:>22.2f}x\n")

    if jax_result:
        if ref_time:
            speedup = ref_time / jax_result['time'] if jax_result['time'] > 0 else 0
            f.write(f"  {'JAX GPU':<28} {jax_result['time']:>12.6f} {speedup:>22.2f}x\n")
        else:
            f.write(f"  {'JAX GPU':<28} {jax_result['time']:>12.6f}\n")

    if taichi_result:
        if ref_time:
            speedup = ref_time / taichi_result['time'] if taichi_result['time'] > 0 else 0
            f.write(f"  {'Taichi GPU':<28} {taichi_result['time']:>12.6f} {speedup:>22.2f}x\n")
        else:
            f.write(f"  {'Taichi GPU':<28} {taichi_result['time']:>12.6f}\n")

    # Summary line
    if fortran_results:
        best_f = min(fortran_results, key=lambda x: x['time'])
        f.write(f"\n  Best Fortran: {best_f['threads']} threads @ {best_f['time']:.6f}s\n")
    if jax_result:
        f.write(f"  JAX GPU:     {jax_result['time']:.6f}s\n")
    if taichi_result:
        f.write(f"  Taichi GPU:  {taichi_result['time']:.6f}s\n")

    if fortran_results and (jax_result or taichi_result):
        best_f = min(fortran_results, key=lambda x: x['time'])
        gpu_entries = []
        if jax_result:
            gpu_entries.append(('JAX GPU', jax_result['time']))
        if taichi_result:
            gpu_entries.append(('Taichi GPU', taichi_result['time']))
        best_gpu = min(gpu_entries, key=lambda x: x[1])
        ratio = best_f['time'] / best_gpu[1] if best_gpu[1] > 0 else 0
        f.write(f"  Best GPU ({best_gpu[0]}) vs Best Fortran speedup: {ratio:.2f}x\n")

    f.write(f"\n")

def main():
    if not compile_fortran():
        print("Cannot proceed without Fortran compilation.")
        return

    all_results = {}

    for i, grid_size in enumerate(GRID_SIZES):
        section_num = i + 1
        print(f"\n{'='*60}")
        print(f"SECTION {section_num}: Grid {grid_size}x{grid_size}x{grid_size}")
        print(f"{'='*60}")

        fortran_results = run_fortran(grid_size)
        jax_result = run_jax(grid_size)
        taichi_result = run_taichi(grid_size)

        all_results[grid_size] = (fortran_results, jax_result, taichi_result)

    # Write results file
    print(f"\nWriting results to {RESULT_FILE} ...")

    with open(RESULT_FILE, 'w') as f:
        f.write("=" * 90 + "\n")
        f.write("  EXPLICIT EULER LAPLACIAN BENCHMARK RESULTS\n")
        f.write("  Fortran (OpenMP) vs JAX (GPU/CUDA) vs Taichi (GPU/CUDA)\n")
        f.write("=" * 90 + "\n")
        f.write(f"  Equation:    du/dt = alpha * Laplacian(u)\n")
        f.write(f"  Method:      Explicit Euler, dt={0.0001}, dx={0.1}, alpha={0.1}\n")
        f.write(f"  Time steps:  {NSTEPS}\n")
        f.write(f"  Fortran threads tested: {THREAD_COUNTS}\n")
        f.write(f"  Grid sizes tested: {['{}^3'.format(n) for n in GRID_SIZES]}\n")

        for i, grid_size in enumerate(GRID_SIZES):
            section_num = i + 1
            fortran_results, jax_result, taichi_result = all_results[grid_size]
            write_section(f, section_num, grid_size, fortran_results, jax_result, taichi_result)

        # Final cross-section summary
        f.write("=" * 90 + "\n")
        f.write("  CROSS-SECTION SUMMARY\n")
        f.write("=" * 90 + "\n\n")

        f.write(f"  {'Grid':<12} {'Fortran 1T (s)':>16} {'Best Fortran (s)':>18} {'(threads)':>10} {'JAX GPU (s)':>14} {'Taichi GPU (s)':>16}\n")
        f.write(f"  {'-'*86}\n")

        for grid_size in GRID_SIZES:
            fortran_results, jax_result, taichi_result = all_results[grid_size]
            f1 = fortran_results[0]['time'] if fortran_results else float('nan')
            if fortran_results:
                best_f = min(fortran_results, key=lambda x: x['time'])
                bf_time = best_f['time']
                bf_threads = best_f['threads']
            else:
                bf_time = float('nan')
                bf_threads = 0
            jt = jax_result['time'] if jax_result else float('nan')
            tt = taichi_result['time'] if taichi_result else float('nan')
            f.write(f"  {grid_size:>3}^3        {f1:>16.6f} {bf_time:>18.6f} {bf_threads:>10d} {jt:>14.6f} {tt:>16.6f}\n")

        f.write(f"\n{'='*90}\n")
        f.write(f"  END OF BENCHMARK RESULTS\n")
        f.write(f"{'='*90}\n")

    print(f"Done! Results written to benchmark_results.txt")

if __name__ == "__main__":
    main()
