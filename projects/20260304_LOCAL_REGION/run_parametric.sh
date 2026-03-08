#!/bin/bash
# run_parametric.sh — Task 7: Sweep localnum=0,4..10 with timax=0.01
# Runs cases sequentially (each case modifies the shared input file)
set -e

PROJ_DIR="$(cd "$(dirname "$0")" && pwd)"
SRC_DIR="$(cd "$PROJ_DIR/../../fortran_new" && pwd)"
INPUT_FILE="$SRC_DIR/inputfile/input_param.txt"
OMP_THREADS=4

echo "=== LOCAL_REGION Parametric Study (timax=0.01) ==="

# Backup original input
cp "$INPUT_FILE" "$PROJ_DIR/input_param_backup.txt"

# Compile once
cd "$SRC_DIR"
bash compile.sh

# Run one case
run_case() {
    local LOCALNUM=$1
    local CASE="param_local${LOCALNUM}"
    [ "$LOCALNUM" -eq 0 ] && CASE="param_baseline"

    # Set localnum and timax in input file
    sed -i "s/localnum=[0-9]*/localnum=$LOCALNUM/" "$INPUT_FILE"
    sed -i "s/timax=[0-9.eE+-]*/timax=0.01/" "$INPUT_FILE"

    echo ""
    echo "[${CASE}] Starting localnum=$LOCALNUM with $OMP_THREADS threads..."
    cd "$SRC_DIR"
    bash run.sh "$CASE" $OMP_THREADS

    # Copy results to project dir
    local RDIR="$SRC_DIR/result/$CASE"
    cp "$RDIR/${CASE}_thermal_history.txt" "$PROJ_DIR/${CASE}_thermal_history.txt"
    cp "$RDIR/${CASE}_timing_report.txt" "$PROJ_DIR/${CASE}_timing_report.txt"
    [ -f "$RDIR/${CASE}_thermal_history.png" ] && cp "$RDIR/${CASE}_thermal_history.png" "$PROJ_DIR/${CASE}_thermal_history.png"

    echo "[${CASE}] Done."
}

# Run all cases sequentially
for N in 0 4 5 6 7 8 9 10; do
    run_case $N
done

# Restore original input
cp "$PROJ_DIR/input_param_backup.txt" "$INPUT_FILE"
echo ""
echo "=== Parametric study complete ==="
