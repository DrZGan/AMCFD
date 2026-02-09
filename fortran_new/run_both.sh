#!/bin/bash
# Kill, clean, compile, and run both fortran_new and fortran_origin
# Both run in parallel; script waits for both to finish, then summarizes.
SCRIPTDIR="$(cd "$(dirname "$0")" && pwd)"
BASEDIR="$(cd "$SCRIPTDIR/.." && pwd)"
# Comparison scripts live in AMCFD (BASEDIR)
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-12}
export GFORTRAN_UNBUFFERED_ALL=1
echo "=============================================="
echo "  run_both: kill -> clean -> compile -> run"
echo "=============================================="
echo ""
echo "1. Killing any running cluster_main..."
killall -9 cluster_main 2>/dev/null || true
sleep 1
echo ""
echo "2. Building fortran_origin..."
cd "$BASEDIR/fortran_origin"
mkdir -p result
bash clean.sh 2>/dev/null || true
bash compile.sh 2>&1
ORIG_BUILD=$?
echo ""
echo "3. Building fortran_new..."
cd "$BASEDIR/fortran_new"
mkdir -p result
bash compile.sh 2>&1
NEW_BUILD=$?
if [ $ORIG_BUILD -ne 0 ]; then echo "ERROR: fortran_origin build failed."; exit 1; fi
if [ $NEW_BUILD -ne 0 ]; then echo "ERROR: fortran_new build failed."; exit 1; fi
echo ""
echo "4. Cleaning old results..."
rm -f "$BASEDIR/fortran_new/result/output.txt" "$BASEDIR/fortran_new/result/"*.vtk
rm -f "$BASEDIR/fortran_origin/result/output.txt" "$BASEDIR/fortran_origin/result/"*.vtk
echo ""
echo "5. Running both simulations..."
cd "$BASEDIR/fortran_origin"
t0_orig=$(date +%s.%N)
./cluster_main > /dev/null 2>&1 &
ORIG_PID=$!
echo "   fortran_origin PID=$ORIG_PID"
cd "$BASEDIR/fortran_new"
t0_new=$(date +%s.%N)
./cluster_main > /dev/null 2>&1 &
NEW_PID=$!
echo "   fortran_new    PID=$NEW_PID"
echo "   Waiting for both to finish..."
wait $NEW_PID 2>/dev/null
t1_new=$(date +%s.%N)
new_wall=$(awk "BEGIN {printf \"%.1f\", $t1_new - $t0_new}")
echo "   fortran_new    finished in ${new_wall}s"
wait $ORIG_PID 2>/dev/null
t1_orig=$(date +%s.%N)
orig_wall=$(awk "BEGIN {printf \"%.1f\", $t1_orig - $t0_orig}")
echo "   fortran_origin finished in ${orig_wall}s"
echo ""
echo "=============================================="
echo "  Results"
echo "=============================================="
echo ""
echo "--- fortran_origin ---"
echo "  Wall time: ${orig_wall}s"
orig_lines=$(wc -l < "$BASEDIR/fortran_origin/result/output.txt" 2>/dev/null || echo 0)
orig_vtk=$(ls "$BASEDIR/fortran_origin/result/"*.vtk 2>/dev/null | wc -l)
echo "  output.txt: ${orig_lines} lines"
echo "  VTK files:  ${orig_vtk}"
echo ""
echo "--- fortran_new ---"
echo "  Wall time: ${new_wall}s"
new_lines=$(wc -l < "$BASEDIR/fortran_new/result/output.txt" 2>/dev/null || echo 0)
new_vtk=$(ls "$BASEDIR/fortran_new/result/"*.vtk 2>/dev/null | wc -l)
echo "  output.txt: ${new_lines} lines"
echo "  VTK files:  ${new_vtk}"
echo ""
ratio=$(awk "BEGIN {if($new_wall>0) printf \"%.2f\", $orig_wall/$new_wall; else print 0}")
echo "  Speedup: ${ratio}x"
if [ -f "$BASEDIR/fortran_origin/result/output.txt" ] && [ -f "$BASEDIR/fortran_new/result/output.txt" ]; then
  echo ""
  echo "--- output.txt comparison ---"
  python3 "$BASEDIR/compare_outputs.py" "$BASEDIR/fortran_origin/result/output.txt" "$BASEDIR/fortran_new/result/output.txt" /dev/stdout 2>/dev/null || true
fi
if [ "$orig_vtk" -gt 0 ] 2>/dev/null && [ "$new_vtk" -gt 0 ] 2>/dev/null; then
  echo ""
  echo "--- VTK comparison ---"
  python3 "$BASEDIR/compare_vtk.py" "$BASEDIR/fortran_origin/result" "$BASEDIR/fortran_new/result" 2>/dev/null || true
fi
echo ""
echo "DONE"
