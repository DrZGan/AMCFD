#!/bin/bash
# Run both fortran_new and fortran_origin for 3 minutes
# Use GFORTRAN_UNBUFFERED_ALL=1 to ensure output is flushed even when killed

BASEDIR="$(cd "$(dirname "$0")" && pwd)"
export GFORTRAN_UNBUFFERED_ALL=1

# Clean old results
rm -f "$BASEDIR/fortran_new/result/output.txt" "$BASEDIR/fortran_new/result/"*.vtk
rm -f "$BASEDIR/fortran_origin/result/output.txt" "$BASEDIR/fortran_origin/result/"*.vtk

echo "Starting both simulations with unbuffered I/O..."

# Run fortran_new in background
cd "$BASEDIR/fortran_new"
./cluster_main > /dev/null 2>&1 &
NEW_PID=$!
echo "fortran_new  PID=$NEW_PID"

# Run fortran_origin in background
cd "$BASEDIR/fortran_origin"
./cluster_main > /dev/null 2>&1 &
ORIG_PID=$!
echo "fortran_origin PID=$ORIG_PID"

echo "Waiting 180 seconds..."
sleep 180

echo "Killing processes..."
kill -INT $NEW_PID 2>/dev/null
kill -INT $ORIG_PID 2>/dev/null
sleep 3
kill -9 $NEW_PID 2>/dev/null
kill -9 $ORIG_PID 2>/dev/null
wait $NEW_PID 2>/dev/null
wait $ORIG_PID 2>/dev/null

echo ""
echo "=== fortran_new results ==="
wc -l "$BASEDIR/fortran_new/result/output.txt" 2>/dev/null || echo "NO output.txt"
ls "$BASEDIR/fortran_new/result/"*.vtk 2>/dev/null | wc -l
echo " VTK files"

echo "=== fortran_origin results ==="
wc -l "$BASEDIR/fortran_origin/result/output.txt" 2>/dev/null || echo "NO output.txt"
ls "$BASEDIR/fortran_origin/result/"*.vtk 2>/dev/null | wc -l
echo " VTK files"

echo ""
echo "DONE"
