#!/bin/bash
# Run TDMA benchmarks for multiple grid sizes
# Usage: bash run_all_benchmarks.sh

WORKDIR="$(cd "$(dirname "$0")" && pwd)"
cd "$WORKDIR"

SIZES=(50 100)  #(50 100 200 400)
FORTRAN_THREADS=(1 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30)
N_ITER=10

echo "=== TDMA Multi-Grid Benchmark ==="
echo "Hardware:"
lscpu | grep "Model name" | head -1
nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null || echo "No GPU info"
echo ""

for SIZE in "${SIZES[@]}"; do
    echo "========================================="
    echo "Running grid size: ${SIZE}x${SIZE}x${SIZE}"
    echo "========================================="

    OUTFILE="${WORKDIR}/bench_${SIZE}.txt"
    > "$OUTFILE"

    # --- Fortran ---
    echo "  Compiling Fortran for ${SIZE}^3..."
    # Generate Fortran source with this grid size
    sed "s/integer, parameter :: ni = [0-9]*, nj = [0-9]*, nk = [0-9]*/integer, parameter :: ni = ${SIZE}, nj = ${SIZE}, nk = ${SIZE}/" \
        tdma_fortran_benchmark.f90 > "/tmp/tdma_fortran_bench_${SIZE}.f90"

    gfortran -O3 -fopenmp -march=native "/tmp/tdma_fortran_bench_${SIZE}.f90" -o "/tmp/tdma_fortran_bench_${SIZE}"

    if [ $? -eq 0 ]; then
        echo "  Running Fortran benchmarks..."
        echo "FORTRAN_RESULTS" >> "$OUTFILE"
        for T in "${FORTRAN_THREADS[@]}"; do
            echo "    threads=$T ..."
            RESULT=$("/tmp/tdma_fortran_bench_${SIZE}" "$T" 2>/dev/null)
            echo "$RESULT" >> "$OUTFILE"
            echo "    -> $RESULT"
        done
    else
        echo "  Fortran compilation failed for size $SIZE!"
        echo "FORTRAN_RESULTS" >> "$OUTFILE"
        echo "COMPILATION_FAILED" >> "$OUTFILE"
    fi

    # --- Taichi GPU ---
    echo "  Running Taichi GPU benchmark for ${SIZE}^3..."
    echo "TAICHI_RESULTS" >> "$OUTFILE"
    # RESULT=$(LD_LIBRARY_PATH=/usr/lib/wsl/lib:$LD_LIBRARY_PATH python3 tdma_taichi_benchmark.py --gpu --size "$SIZE" 2>/dev/null)
    RESULT=$(CUDA_VISIBLE_DEVICES=1 python3 tdma_taichi_benchmark.py --gpu --size "$SIZE")
    echo "$RESULT" >> "$OUTFILE"
    echo "    -> $RESULT"

    # --- JAX GPU ---
    echo "  Running JAX GPU benchmark for ${SIZE}^3..."
    echo "JAX_RESULTS" >> "$OUTFILE"
    # RESULT=$(python3 tdma_jax_benchmark.py --gpu --size "$SIZE" 2>/dev/null)
    # RESULT=$(CUDA_VISIBLE_DEVICES=2 python3 tdma_jax_benchmark.py --gpu --size "$SIZE")
    RESULT=$(CUDA_VISIBLE_DEVICES=2 python3 tdma_jax_benchmark.py --gpu --size "$SIZE" 2>&1)
    echo "$RESULT" >> "$OUTFILE"
    echo "    -> $RESULT"

    # --- JAX GPU 2 vectorize---
    echo "  Running JAX GPU benchmark for ${SIZE}^3..."
    echo "JAX_Vectorized_RESULTS" >> "$OUTFILE"
    RESULT=$(CUDA_VISIBLE_DEVICES=2 python3 tdma_jax_vectorize_benchmark.py --gpu --size "$SIZE" 2>&1)
    echo "$RESULT" >> "$OUTFILE"
    echo "    -> $RESULT"

    echo "  Done with ${SIZE}^3"
    echo ""
done

echo "All benchmarks complete. Results in bench_*.txt files."
