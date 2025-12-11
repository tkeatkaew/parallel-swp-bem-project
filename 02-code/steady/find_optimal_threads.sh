#!/bin/bash

# Find optimal thread configuration for your 6-core system
# Usage: ./find_optimal_threads.sh

echo "=== CPU Architecture Analysis ==="
echo "Physical cores: 6"
echo "Logical cores: 12 (with hyperthreading)"
echo "Current OMP_NUM_THREADS: $OMP_NUM_THREADS"
echo ""

# Ensure the program is compiled
echo "=== Building program ==="
make clean 2>/dev/null
make

if [ $? -ne 0 ]; then
    echo "Build failed!"
    exit 1
fi

echo ""
echo "=== Testing Different Thread Configurations ==="
echo "Format: Threads | Real Time | CPU Usage | Memory"
echo "=================================================="

# Test different thread counts
for threads in 1 2 4 6 8 12; do
    echo "Testing with $threads threads..."
    
    # Set environment variables
    export OMP_NUM_THREADS=$threads
    export OPENBLAS_NUM_THREADS=$threads
    export OMP_PROC_BIND=true
    export OMP_PLACES=cores
    
    # Run 3 times and take average
    total_time=0
    for run in 1 2 3; do
        echo "  Run $run/3..."
        start_time=$(date +%s.%N)
        ./catcharea 1.0 100.0 0.001 > /dev/null 2>&1
        end_time=$(date +%s.%N)
        run_time=$(echo "$end_time - $start_time" | bc -l)
        total_time=$(echo "$total_time + $run_time" | bc -l)
    done
    
    avg_time=$(echo "scale=3; $total_time / 3" | bc -l)
    speedup=$(echo "scale=2; $baseline_time / $avg_time" | bc -l 2>/dev/null || echo "1.00")
    
    if [ $threads -eq 1 ]; then
        baseline_time=$avg_time
        speedup="1.00"
    fi
    
    printf "%2d threads | %8.3f s | Speedup: %5.2fx\n" $threads $avg_time $speedup
done

echo ""
echo "=== Memory Bandwidth Test ==="
echo "Testing memory-intensive operations..."

for threads in 1 6 12; do
    export OMP_NUM_THREADS=$threads
    export OPENBLAS_NUM_THREADS=$threads
    
    echo "Testing $threads threads with memory monitoring..."
    /usr/bin/time -f "Memory: %M KB, CPU: %P" ./catcharea 1.0 100.0 0.001 2>&1 | grep "Memory:"
done

echo ""
echo "=== Recommendations ==="
echo "Based on your 6-core system:"
echo "1. Try OMP_NUM_THREADS=6 (one thread per physical core)"
echo "2. If memory-bound: consider OMP_NUM_THREADS=4"
echo "3. For CPU-intensive tasks: test OMP_NUM_THREADS=8"
echo "4. Avoid OMP_NUM_THREADS=12 unless specifically beneficial"
echo ""
echo "Optimal settings to add to /etc/environment:"
echo "OMP_NUM_THREADS=6"
echo "OMP_PROC_BIND=true" 
echo "OMP_PLACES=cores"
echo "OPENBLAS_NUM_THREADS=6"
