#!/bin/bash

echo "=========================================="
echo "TDMA Solver Performance Comparison"
echo "Fortran (OpenMP) vs Taichi (CPU/GPU)"
echo "=========================================="
echo ""

# Set OpenMP threads
export OMP_NUM_THREADS=8

# Compile Fortran versions
echo "Compiling Fortran code..."
gfortran -O3 -fopenmp -march=native -o tdma_fortran tdma_fortran.f90
gfortran -O3 -fopenmp -march=native -funroll-loops -ffast-math -o tdma_fortran_optimized tdma_fortran_optimized.f90

if [ $? -eq 0 ]; then
    echo "Fortran compilation successful!"
else
    echo "Fortran compilation failed!"
    exit 1
fi

echo ""
echo "=========================================="
echo "1. Running Fortran (basic OpenMP)..."
echo "=========================================="
./tdma_fortran

echo ""
echo "=========================================="
echo "2. Running Fortran (optimized OpenMP)..."
echo "=========================================="
./tdma_fortran_optimized

echo ""
echo "=========================================="
echo "3. Running Taichi CPU version..."
echo "=========================================="
python tdma_taichi.py --cpu --threads 8

echo ""
echo "=========================================="
echo "4. Running Taichi GPU version..."
echo "=========================================="
python tdma_taichi.py --gpu

echo ""
echo "=========================================="
echo "Comparison complete!"
echo "=========================================="
