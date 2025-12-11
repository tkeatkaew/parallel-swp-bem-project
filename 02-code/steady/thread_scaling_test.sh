#!/bin/bash
# thread_scaling_test.sh

for threads in 1 2 4 6; do
    export OMP_NUM_THREADS=$threads
    export OPENBLAS_NUM_THREADS=$threads
    
    echo "===== Testing $threads threads ====="
    
    # OpenBLAS
    echo "--- OpenBLAS ---"
    ./catcharea 1.0 100.0 0.001 0 3 32 1 | tee "openblas_${threads}t.txt"
    
    # Hybrid
    echo "--- Hybrid ---"
    ./catcharea 1.0 100.0 0.001 0 3 32 0 | tee "hybrid_${threads}t.txt"
done
