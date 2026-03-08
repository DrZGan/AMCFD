#!/bin/bash
# AM-CFD Run Script
# Usage: bash run.sh <case_name> [omp_threads]
set -e

CASE_NAME="${1:?Usage: bash run.sh <case_name> [omp_threads]}"
OMP_THREADS="${2:-4}"

# Update case_name in input file
sed -i "s/case_name='[^']*'/case_name='$CASE_NAME'/" ./inputfile/input_param.txt

export OMP_NUM_THREADS=$OMP_THREADS
./cluster_main
