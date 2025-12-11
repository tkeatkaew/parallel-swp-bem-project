#!/bin/bash

# Script to build and run catcharea program with optimal parallel settings
# Usage: ./run_catcharea.sh [num_threads]

# Set default number of threads (or use command line argument)
NUM_THREADS=${1:-8}

echo "=== System Information ==="
echo "Available CPU cores: $(nproc)"
echo "Setting OMP_NUM_THREADS to: $NUM_THREADS"

# Set OpenMP environment variables for optimal performance
export OMP_NUM_THREADS=$NUM_THREADS
export OMP_PROC_BIND=true
export OMP_PLACES=cores
export OPENBLAS_NUM_THREADS=$NUM_THREADS

echo "OpenMP Configuration:"
echo "  OMP_NUM_THREADS=$OMP_NUM_THREADS"
echo "  OMP_PROC_BIND=$OMP_PROC_BIND"
echo "  OMP_PLACES=$OMP_PLACES"
echo "  OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS"
echo ""

echo "=== Building catcharea program ==="
echo "Step 1: First make..."
make
if [ $? -ne 0 ]; then
    echo "Error: First make failed!"
    exit 1
fi

echo ""
echo "Step 2: Second make..."
make
if [ $? -ne 0 ]; then
    echo "Error: Second make failed!"
    exit 1
fi

echo ""
echo "=== Running catcharea program ==="
echo "Executing: ./catcharea 1.0 100.0 0.001"
echo "Start time: $(date)"
time ./catcharea 1.0 100.0 0.001
if [ $? -ne 0 ]; then
    echo "Error: catcharea execution failed!"
    exit 1
fi

echo ""
echo "End time: $(date)"
echo "=== Process completed successfully ==="
