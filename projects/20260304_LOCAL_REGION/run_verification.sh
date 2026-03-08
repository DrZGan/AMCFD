#!/bin/bash
# run_verification.sh — Task 6: Compare localnum=0 (baseline) vs localnum=3
set -e

PROJ_DIR="$(cd "$(dirname "$0")" && pwd)"
SRC_DIR="$(cd "$PROJ_DIR/../../fortran_new" && pwd)"
INPUT_FILE="$SRC_DIR/inputfile/input_param.txt"
OMP_THREADS=4

echo "=== LOCAL_REGION Verification ==="

# Backup original input
cp "$INPUT_FILE" "$PROJ_DIR/input_param_backup.txt"

# Compile once
cd "$SRC_DIR"
bash compile.sh

# --- Run 1: localnum=0 (baseline) ---
echo ""
echo "=== Run 1: localnum=0 (baseline) ==="
sed -i "s/localnum=[0-9]*/localnum=0/" "$INPUT_FILE"
cd "$SRC_DIR"
bash run.sh baseline $OMP_THREADS
cp result/baseline/baseline_thermal_history.txt "$PROJ_DIR/baseline_thermal_history.txt"
cp result/baseline/baseline_timing_report.txt "$PROJ_DIR/baseline_timing_report.txt"
[ -f result/baseline/baseline_thermal_history.png ] && cp result/baseline/baseline_thermal_history.png "$PROJ_DIR/baseline_thermal_history.png"

# --- Run 2: localnum=3 ---
echo ""
echo "=== Run 2: localnum=3 ==="
sed -i "s/localnum=[0-9]*/localnum=3/" "$INPUT_FILE"
cd "$SRC_DIR"
bash run.sh local3 $OMP_THREADS
cp result/local3/local3_thermal_history.txt "$PROJ_DIR/local3_thermal_history.txt"
cp result/local3/local3_timing_report.txt "$PROJ_DIR/local3_timing_report.txt"
[ -f result/local3/local3_thermal_history.png ] && cp result/local3/local3_thermal_history.png "$PROJ_DIR/local3_thermal_history.png"

# Restore original input
cp "$PROJ_DIR/input_param_backup.txt" "$INPUT_FILE"
echo ""
echo "=== Verification runs complete ==="
