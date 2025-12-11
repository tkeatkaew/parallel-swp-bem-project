#!/bin/bash
# block_size_test.sh

export OMP_NUM_THREADS=6
export OPENBLAS_NUM_THREADS=6
export OMP_PROC_BIND=spread
export OMP_PLACES=cores

for bs in 16 32 64 128; do
    echo "===== Testing Block Size = $bs ====="
    ./catcharea 1.0 100.0 0.001 0 3 $bs 0 | tee "hybrid_b${bs}.txt"
    echo ""
done

echo "===== Done! ====="
