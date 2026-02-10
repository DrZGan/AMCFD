#!/bin/bash
# Verify localnum=0..5: run origin ONCE, then for each localnum run only fortran_new and compare.
# Run from AMCFD folder: bash verify_nl_report.sh
BASEDIR="$(cd "$(dirname "$0")" && pwd)"
REPORT="$BASEDIR/compare_report/nl_verification_report.txt"
REPORT_DIR="$BASEDIR/compare_report"
INPUT_PARAM="$BASEDIR/fortran_new/inputfile/input_param.txt"
TIMING_TXT="$BASEDIR/fortran_new/result/timing_report.txt"
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-12}

mkdir -p "$REPORT_DIR"

echo "=============================================="
echo "  verify_nl_report: build -> run origin once -> run new for each localnum"
echo "=============================================="
echo ""
echo "1. Killing any running cluster_main..."
killall -9 cluster_main 2>/dev/null || true
sleep 1

echo ""
echo "2. Building fortran_origin and fortran_new..."
cd "$BASEDIR/fortran_origin"
mkdir -p result
bash clean.sh 2>/dev/null || true
bash compile.sh 2>&1
[ $? -ne 0 ] && { echo "ERROR: fortran_origin build failed."; exit 1; }
cd "$BASEDIR/fortran_new"
mkdir -p result
bash compile.sh 2>&1
[ $? -ne 0 ] && { echo "ERROR: fortran_new build failed."; exit 1; }

echo ""
echo "3. Cleaning old results..."
rm -f "$BASEDIR/fortran_origin/result/output.txt" "$BASEDIR/fortran_origin/result/"*.vtk
rm -f "$BASEDIR/fortran_new/result/output.txt" "$BASEDIR/fortran_new/result/"*.vtk

echo ""
echo "4. Running fortran_origin once (reference for all localnum)..."
cd "$BASEDIR/fortran_origin"
t0_orig=$(date +%s.%N)
./cluster_main > /dev/null 2>&1
t1_orig=$(date +%s.%N)
orig_wall=$(awk "BEGIN {printf \"%.1f\", $t1_orig - $t0_orig}")
echo "   fortran_origin finished in ${orig_wall}s (used as reference for all localnum)"

{
  echo "=============================================="
  echo "  localnum=0..5 Verification Report"
  echo "  Generated: $(date -Iseconds 2>/dev/null || date)"
  echo "  Local solver: local_half_x=1e-3, local_half_y=2e-4, local_depth_z=2e-4"
  echo "  Origin run once; new run per localnum (comparison vs same origin result)"
  echo "=============================================="
  echo ""
  printf "  %8s  %10s  %10s  %12s  %-8s  %-8s\n" "localnum" "origin(s)" "new(s)" "new_CPU(s)" "output" "vtk"
  echo "  -------------------------------------------------------------------------"
} > "$REPORT"

echo ""
echo "5. Running fortran_new for each localnum and comparing..."
for localnum in 0 1 2 3 4 5; do
  # Update only the local_solver localnum entry (avoid touching denl in material properties)
  sed -E -i "s/^(&local_solver[[:space:]]+localnum=)[0-9]+/\\1$localnum/" "$INPUT_PARAM"

  cd "$BASEDIR/fortran_new"
  rm -f result/output.txt result/*.vtk result/timing_report.txt
  t0_new=$(date +%s.%N)
  ./cluster_main > /dev/null 2>&1
  t1_new=$(date +%s.%N)
  new_wall=$(awk "BEGIN {printf \"%.1f\", $t1_new - $t0_new}")

  python3 "$BASEDIR/compare_output_full.py" \
    "$BASEDIR/fortran_origin/result/output.txt" \
    "$BASEDIR/fortran_new/result/output.txt" \
    "$REPORT_DIR/output_compare.txt" 2>/dev/null || true
  out_conc=$(sed -n 's/.*Conclusion:[[:space:]]*\(PASS\|FAIL\).*/\1/p' "$REPORT_DIR/output_compare.txt" 2>/dev/null | head -1)

  python3 "$BASEDIR/compare_vtk.py" \
    "$BASEDIR/fortran_origin/result" \
    "$BASEDIR/fortran_new/result" 2>/dev/null | tee "$REPORT_DIR/vtk_compare.txt" > /dev/null
  vtk_conc=$(sed -n 's/^Conclusion:[[:space:]]*\(PASS\|FAIL\).*/\1/p' "$REPORT_DIR/vtk_compare.txt" 2>/dev/null | head -1)

  new_cpu=""
  [ -f "$TIMING_TXT" ] && new_cpu=$(sed -n 's/.*Total CPU time:[[:space:]]*\([0-9.]*\).*/\1/p' "$TIMING_TXT" | head -1)

  printf "  %8d  %10s  %10s  %12s  %-8s  %-8s\n" "$localnum" "$orig_wall" "${new_wall}" "${new_cpu:---}" "${out_conc:---}" "${vtk_conc:---}" >> "$REPORT"
  echo "   localnum=$localnum done: new=${new_wall}s  output=$out_conc  vtk=$vtk_conc"
done

{
  echo ""
  echo "=============================================="
  echo "  End of report. Last comparison details:"
  echo "    $REPORT_DIR/output_compare.txt"
  echo "    $REPORT_DIR/vtk_compare.txt"
  echo "=============================================="
} >> "$REPORT"

echo ""
echo "Report written to: $REPORT"
cat "$REPORT"
