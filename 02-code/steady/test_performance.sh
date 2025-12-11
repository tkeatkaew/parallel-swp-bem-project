#!/bin/bash

# Performance testing script for different thread counts
# Usage: ./test_performance.sh

echo "=== Performance Testing for Catcharea ==="
echo "System has $(nproc) CPU cores"
echo ""

# Test different thread counts
for threads in 1 2 4 8 16; do
    if [ $threads -le $(nproc) ]; then
        echo "=== Testing with $threads threads ==="
        export OMP_NUM_THREADS=$threads
        export OPENBLAS_NUM_THREADS=$threads
        export OMP_PROC_BIND=true
        export OMP_PLACES=cores
        
        echo "Running with OMP_NUM_THREADS=$threads"
        echo "Start: $(date)"
        
        # Run and time the execution
        /usr/bin/time -f "Real time: %e seconds, CPU: %P, Memory: %M KB" \
            ./catcharea 1.0 100.0 0.001 2>&1 | grep -E "(Real time|Error|Cannot)"
        
        echo "End: $(date)"
        echo "----------------------------------------"
        echo ""
    fi
done

echo "=== Performance test completed ==="
