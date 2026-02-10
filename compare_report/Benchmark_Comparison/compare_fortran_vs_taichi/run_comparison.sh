#!/bin/bash

# Set CUDA library path for WSL2
export LD_LIBRARY_PATH=/usr/lib/wsl/lib:$LD_LIBRARY_PATH

echo "=========================================="
echo "TDMA Solver Performance Comparison"
echo "Grid: 400 x 400 x 100"
echo "Iterations: 10 (averaged)"
echo "=========================================="
echo ""

# Compile Fortran
echo "Compiling Fortran..."
gfortran -O3 -fopenmp -march=native -ffast-math -o tdma_fortran_benchmark tdma_fortran_benchmark.f90

if [ $? -ne 0 ]; then
    echo "Fortran compilation failed!"
    exit 1
fi
echo "Compilation successful!"
echo ""

# Create results file
RESULTS_FILE="benchmark_results.txt"
echo "TDMA Solver Performance Comparison" > $RESULTS_FILE
echo "Grid: 400 x 400 x 100" >> $RESULTS_FILE
echo "Iterations: 10 (averaged)" >> $RESULTS_FILE
echo "Date: $(date)" >> $RESULTS_FILE
echo "" >> $RESULTS_FILE

# Hardware specifications
echo "=========================================="
echo "HARDWARE SPECIFICATIONS"
echo "=========================================="
echo ""

echo "HARDWARE SPECIFICATIONS" >> $RESULTS_FILE
echo "----------------------------" >> $RESULTS_FILE

# CPU info
echo "CPU Information:"
CPU_MODEL=$(lscpu 2>/dev/null | grep "Model name" | sed 's/Model name:[[:space:]]*//' || echo "Unknown")
CPU_CORES=$(nproc 2>/dev/null || echo "Unknown")
CPU_THREADS=$(lscpu 2>/dev/null | grep "^CPU(s):" | awk '{print $2}' || echo "Unknown")
echo "  Model: $CPU_MODEL"
echo "  Cores: $CPU_CORES"
echo "  Threads: $CPU_THREADS"

echo "CPU: $CPU_MODEL" >> $RESULTS_FILE
echo "CPU Cores: $CPU_CORES" >> $RESULTS_FILE
echo "CPU Threads: $CPU_THREADS" >> $RESULTS_FILE

# GPU info (try nvidia-smi first, then lspci)
echo ""
echo "GPU Information:"
if command -v nvidia-smi &> /dev/null; then
    GPU_MODEL=$(nvidia-smi --query-gpu=name --format=csv,noheader 2>/dev/null | head -1 || echo "Unknown")
    GPU_MEMORY=$(nvidia-smi --query-gpu=memory.total --format=csv,noheader 2>/dev/null | head -1 || echo "Unknown")
    GPU_DRIVER=$(nvidia-smi --query-gpu=driver_version --format=csv,noheader 2>/dev/null | head -1 || echo "Unknown")
    echo "  Model: $GPU_MODEL"
    echo "  Memory: $GPU_MEMORY"
    echo "  Driver: $GPU_DRIVER"
    echo "GPU: $GPU_MODEL" >> $RESULTS_FILE
    echo "GPU Memory: $GPU_MEMORY" >> $RESULTS_FILE
    echo "GPU Driver: $GPU_DRIVER" >> $RESULTS_FILE
else
    GPU_INFO=$(lspci 2>/dev/null | grep -i "vga\|3d\|display" | head -1 || echo "No GPU detected")
    echo "  $GPU_INFO"
    echo "GPU: $GPU_INFO" >> $RESULTS_FILE
fi
echo "" >> $RESULTS_FILE
echo ""

# Run Fortran benchmarks
echo "Running Fortran benchmarks..."
echo "Fortran Results (threads, avg_time_ms, checksum):" >> $RESULTS_FILE

declare -a FORTRAN_TIMES
declare -a FORTRAN_CHECKSUMS
declare -a FORTRAN_MAX
declare -a FORTRAN_MIN
declare -a FORTRAN_AVG

for threads in 1 2 3 4 5 6 7 8 12 16 20 24 28 32; do
    result=$(./tdma_fortran_benchmark $threads)
    echo "  Threads $threads: $result"
    echo "$result" >> $RESULTS_FILE

    # Extract values for later comparison
    time_val=$(echo $result | cut -d',' -f2 | tr -d ' ')
    checksum_val=$(echo $result | cut -d',' -f3 | tr -d ' ')
    max_val=$(echo $result | cut -d',' -f4 | tr -d ' ')
    min_val=$(echo $result | cut -d',' -f5 | tr -d ' ')
    avg_val=$(echo $result | cut -d',' -f6 | tr -d ' ')
    FORTRAN_TIMES[$threads]=$time_val
    FORTRAN_CHECKSUMS[$threads]=$checksum_val
    FORTRAN_MAX[$threads]=$max_val
    FORTRAN_MIN[$threads]=$min_val
    FORTRAN_AVG[$threads]=$avg_val
done

echo "" >> $RESULTS_FILE

# Run Taichi GPU benchmark (skip CPU tests)
echo ""
echo "Running Taichi GPU benchmark..."
echo "Taichi GPU Results:" >> $RESULTS_FILE

GPU_RESULT=$(python3 tdma_taichi_benchmark.py --gpu 2>&1 | grep -E "^GPU|^CUDA")
if [ ! -z "$GPU_RESULT" ]; then
    echo "  $GPU_RESULT"
    echo "$GPU_RESULT" >> $RESULTS_FILE

    GPU_TIME=$(echo $GPU_RESULT | cut -d',' -f2 | tr -d ' ')
    GPU_CHECKSUM=$(echo $GPU_RESULT | cut -d',' -f3 | tr -d ' ')
    GPU_MAX=$(echo $GPU_RESULT | cut -d',' -f4 | tr -d ' ')
    GPU_MIN=$(echo $GPU_RESULT | cut -d',' -f5 | tr -d ' ')
    GPU_AVG=$(echo $GPU_RESULT | cut -d',' -f6 | tr -d ' ')
else
    echo "  GPU not available or failed"
    echo "GPU not available" >> $RESULTS_FILE
    GPU_TIME=""
fi

echo "" >> $RESULTS_FILE

# Print comparison table
echo ""
echo "=========================================="
echo "COMPARISON TABLE"
echo "=========================================="
echo ""

# Header - all speedups relative to Fortran 1-thread baseline
f_baseline=${FORTRAN_TIMES[1]}
printf "| %-7s | %-15s | %-12s |\n" "Threads" "Fortran (ms)" "Speedup"
printf "|---------|-----------------|--------------|--|\n"

for threads in 1 2 3 4 5 6 7 8 12 16 20 24 28 32; do
    f_time=${FORTRAN_TIMES[$threads]}

    if [ ! -z "$f_time" ]; then
        # Calculate speedup relative to Fortran 1-thread baseline
        f_speedup=$(awk "BEGIN {printf \"%.2f\", $f_baseline / $f_time}" 2>/dev/null || echo "N/A")
        printf "| %-7d | %15.3f | %12s |\n" $threads $f_time "${f_speedup}x"
    fi
done

# Add GPU row if available
if [ ! -z "$GPU_TIME" ]; then
    gpu_speedup=$(awk "BEGIN {printf \"%.2f\", $f_baseline / $GPU_TIME}" 2>/dev/null || echo "N/A")
    printf "|---------|-----------------|--------------|--|\n"
    printf "| %-7s | %15.3f | %12s |\n" "GPU" $GPU_TIME "${gpu_speedup}x"
fi

echo ""
echo "All speedups relative to Fortran 1-thread baseline (${f_baseline} ms)"
echo "Speedup > 1.0 means faster than baseline, < 1.0 means slower"
echo ""

# Verify results match
echo "=========================================="
echo "DETAILED VERIFICATION"
echo "=========================================="
echo ""
echo "Using identical Gauss-Seidel algorithm (sequential j-loop) for both"
echo "Fortran and Taichi to ensure matching results."
echo ""

# Set Fortran reference values
f_check=${FORTRAN_CHECKSUMS[1]}
f_max=${FORTRAN_MAX[1]}
f_min=${FORTRAN_MIN[1]}
f_avg=${FORTRAN_AVG[1]}

# GPU verification
if [ ! -z "$GPU_TIME" ]; then
    echo "=== Fortran (1 thread) vs Taichi GPU ==="
    echo ""
    echo "                        Fortran                    Taichi GPU"
    echo "--------------------------------------------------------------------------------------------"
    printf "Sum (checksum): %22s %22s\n" "$f_check" "$GPU_CHECKSUM"
    printf "Max value:      %22s %22s\n" "$f_max" "$GPU_MAX"
    printf "Min value:      %22s %22s\n" "$f_min" "$GPU_MIN"
    printf "Average value:  %22s %22s\n" "$f_avg" "$GPU_AVG"
    echo ""

    echo "Numerical Differences (Fortran vs Taichi GPU):"
    echo "--------------------------------------------------------------------------------------------"

    # GPU Checksum difference
    gpu_check_diff=$(echo "$f_check $GPU_CHECKSUM" | awk '{printf "%.6e", $1 - $2}')
    gpu_check_rel=$(echo "$f_check $GPU_CHECKSUM" | awk '{if($1!=0) printf "%.6e", ($1 - $2)/$1; else print "N/A"}')
    echo "  Checksum: absolute diff = $gpu_check_diff, relative diff = $gpu_check_rel"

    # GPU Max difference
    gpu_max_diff=$(echo "$f_max $GPU_MAX" | awk '{printf "%.6e", $1 - $2}')
    gpu_max_rel=$(echo "$f_max $GPU_MAX" | awk '{if($1!=0) printf "%.6e", ($1 - $2)/$1; else print "N/A"}')
    echo "  Max:      absolute diff = $gpu_max_diff, relative diff = $gpu_max_rel"

    # GPU Min difference
    gpu_min_diff=$(echo "$f_min $GPU_MIN" | awk '{printf "%.6e", $1 - $2}')
    gpu_min_rel=$(echo "$f_min $GPU_MIN" | awk '{if($1!=0) printf "%.6e", ($1 - $2)/$1; else print "N/A"}')
    echo "  Min:      absolute diff = $gpu_min_diff, relative diff = $gpu_min_rel"

    # GPU Average difference
    gpu_avg_diff=$(echo "$f_avg $GPU_AVG" | awk '{printf "%.6e", $1 - $2}')
    gpu_avg_rel=$(echo "$f_avg $GPU_AVG" | awk '{if($1!=0) printf "%.6e", ($1 - $2)/$1; else print "N/A"}')
    echo "  Average:  absolute diff = $gpu_avg_diff, relative diff = $gpu_avg_rel"
    echo ""

    # Check if results match (within floating-point tolerance)
    match_status="MISMATCH"
    rel_diff=$(echo "$f_check $GPU_CHECKSUM" | awk '{diff=($1-$2)/$1; if(diff<0) diff=-diff; print diff}')
    is_match=$(echo "$rel_diff" | awk '{if($1 < 0.001) print "yes"; else print "no"}')
    if [ "$is_match" == "yes" ]; then
        match_status="MATCH (within 0.1%)"
    fi
    echo "  Verification Status: $match_status"
    echo ""
else
    echo ""
    echo "GPU not available - no verification performed"
    echo ""
fi

# Also save table to file
echo "" >> $RESULTS_FILE
echo "COMPARISON TABLE (all speedups relative to Fortran 1-thread: ${f_baseline} ms)" >> $RESULTS_FILE
printf "| %-7s | %-15s | %-12s |\n" "Config" "Time (ms)" "Speedup" >> $RESULTS_FILE
printf "|---------|-----------------|--------------|--|\n" >> $RESULTS_FILE

for threads in 1 2 3 4 5 6 7 8 12 16 20 24 28 32; do
    f_time=${FORTRAN_TIMES[$threads]}

    if [ ! -z "$f_time" ]; then
        f_speedup=$(awk "BEGIN {printf \"%.2f\", $f_baseline / $f_time}" 2>/dev/null || echo "N/A")
        printf "| F-%02dt   | %15.3f | %12s |\n" $threads $f_time "${f_speedup}x" >> $RESULTS_FILE
    fi
done

# Add GPU row to file if available
if [ ! -z "$GPU_TIME" ]; then
    gpu_speedup=$(awk "BEGIN {printf \"%.2f\", $f_baseline / $GPU_TIME}" 2>/dev/null || echo "N/A")
    printf "|---------|-----------------|--------------|--|\n" >> $RESULTS_FILE
    printf "| T-GPU   | %15.3f | %12s |\n" $GPU_TIME "${gpu_speedup}x" >> $RESULTS_FILE
fi

# Save detailed verification to file
echo "" >> $RESULTS_FILE
echo "===========================================" >> $RESULTS_FILE
echo "DETAILED VERIFICATION" >> $RESULTS_FILE
echo "===========================================" >> $RESULTS_FILE
echo "" >> $RESULTS_FILE
echo "Algorithm: Gauss-Seidel (sequential j-loop) - identical for Fortran and Taichi" >> $RESULTS_FILE
echo "" >> $RESULTS_FILE

if [ ! -z "$GPU_TIME" ]; then
    echo "Fortran (1 thread) vs Taichi GPU:" >> $RESULTS_FILE
    echo "--------------------------------------------------------------------------------------------" >> $RESULTS_FILE
    printf "  Sum (checksum): Fortran = %s, Taichi = %s\n" "$f_check" "$GPU_CHECKSUM" >> $RESULTS_FILE
    printf "  Max value:      Fortran = %s, Taichi = %s\n" "$f_max" "$GPU_MAX" >> $RESULTS_FILE
    printf "  Min value:      Fortran = %s, Taichi = %s\n" "$f_min" "$GPU_MIN" >> $RESULTS_FILE
    printf "  Average value:  Fortran = %s, Taichi = %s\n" "$f_avg" "$GPU_AVG" >> $RESULTS_FILE
    echo "" >> $RESULTS_FILE
    echo "Numerical Differences:" >> $RESULTS_FILE
    echo "  Checksum: absolute diff = $gpu_check_diff, relative diff = $gpu_check_rel" >> $RESULTS_FILE
    echo "  Max:      absolute diff = $gpu_max_diff, relative diff = $gpu_max_rel" >> $RESULTS_FILE
    echo "  Min:      absolute diff = $gpu_min_diff, relative diff = $gpu_min_rel" >> $RESULTS_FILE
    echo "  Average:  absolute diff = $gpu_avg_diff, relative diff = $gpu_avg_rel" >> $RESULTS_FILE
    echo "" >> $RESULTS_FILE
    echo "Verification Status: $match_status" >> $RESULTS_FILE
fi

echo ""
echo "Results saved to: $RESULTS_FILE"
