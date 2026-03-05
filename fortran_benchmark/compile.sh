#!/bin/bash
# AM-CFD Build Script
# Compiles the Fortran CFD code with OpenMP parallelization

set -e  # Exit on error

# Full clean (clean.sh)
bash clean.sh

echo "Compiling AM-CFD..."

# Compile all modules and main program
# Order matters: dependencies must be compiled first
gfortran -fopenmp -O3 -march=native -c \
    mod_precision.f90 \
    mod_const.f90 \
    mod_cfd_utils.f90 \
    mod_param.f90 \
    mod_geom.f90 \
    mod_field_data.f90 \
    mod_coeff_data.f90 \
    mod_sim_state.f90 \
    mod_local_enthalpy.f90 \
    mod_init.f90 \
    mod_laser.f90 \
    mod_dimen.f90 \
    mod_bound.f90 \
    mod_discret.f90 \
    mod_entot.f90 \
    mod_sour.f90 \
    mod_flux.f90 \
    mod_prop.f90 \
    mod_resid.f90 \
    mod_revise.f90 \
    mod_solve.f90 \
    mod_print.f90 \
    mod_converge.f90 \
    mod_toolpath.f90 \
    mod_timing.f90 \
    main.f90

# Link all object files
gfortran -fopenmp -O3 -march=native *.o -o cluster_main

echo "Build complete: cluster_main"
echo ""
echo "To run: export OMP_NUM_THREADS=12 && ./cluster_main"
