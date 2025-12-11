#!/bin/bash

# Debug script to check actual threading behavior
# Usage: ./debug_threading.sh

echo "=== Environment Variables Check ==="
echo "OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo "OPENBLAS_NUM_THREADS: $OPENBLAS_NUM_THREADS"
echo "OMP_PROC_BIND: $OMP_PROC_BIND"
echo "OMP_PLACES: $OMP_PLACES"
echo ""

echo "=== All OpenMP/BLAS Environment Variables ==="
env | grep -E "(OMP_|OPENBLAS_|MKL_|BLAS_)" | sort
echo ""

echo "=== Starting program with monitoring ==="
echo "Setting strict thread limits..."

# Set ALL possible threading environment variables
export OMP_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export MKL_NUM_THREADS=2
export BLAS_NUM_THREADS=2
export LAPACK_NUM_THREADS=2
export OMP_PROC_BIND=true
export OMP_PLACES=cores
export GOMP_CPU_AFFINITY="0,1"

echo "Environment set to:"
env | grep -E "(OMP_|OPENBLAS_|MKL_|BLAS_|LAPACK_)" | sort
echo ""

# Start the program in background
echo "Starting ./catcharea 1.0 100.0 0.001 ..."
./catcharea 1.0 100.0 0.001 &
PID=$!

# Monitor thread usage
sleep 2
if ps -p $PID > /dev/null; then
    echo "=== Process Information ==="
    echo "PID: $PID"
    echo "Number of threads:"
    ps -o pid,ppid,lwp,nlwp,comm -p $PID
    echo ""
    
    echo "=== Thread Details ==="
    ps -eLf | grep -E "(PID|$PID)" | head -10
    echo ""
    
    echo "=== CPU Affinity ==="
    taskset -cp $PID 2>/dev/null || echo "taskset not available"
    echo ""
    
    echo "=== Monitoring CPU usage for 10 seconds ==="
    echo "Press Ctrl+C to stop monitoring"
    top -p $PID -d 1 -n 10 -b | grep -E "(PID|$PID|%CPU)"
    
    # Wait for process to complete
    wait $PID
else
    echo "Process completed too quickly to monitor"
fi

echo ""
echo "=== Diagnosis ==="
echo "If you see more than 2 threads, possible causes:"
echo "1. OpenBLAS internal threading (check OPENBLAS_NUM_THREADS)"
echo "2. Nested parallelism in your code"
echo "3. Default threading in LAPACK routines"
echo "4. System default overriding your settings"
