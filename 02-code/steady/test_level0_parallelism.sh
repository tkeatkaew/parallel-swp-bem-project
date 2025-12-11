#!/usr/bin/env bash
#
# test_level0_parallelism.sh
# Test Level 0 (Zone-Level) Parallelism with different thread configurations
#
# Usage: ./test_level0_parallelism.sh

set -u

echo ""
echo "╔══════════════════════════════════════════════════════════════╗"
echo "║  LEVEL 0 PARALLELISM TEST SUITE                             ║"
echo "╚══════════════════════════════════════════════════════════════╝"
echo ""

# Check if executable exists
if [ ! -f "./catcharea_level0" ]; then
    echo "ERROR: catcharea_level0 not found!"
    echo "Please compile first:"
    echo "  gcc -O3 -fopenmp -mavx2 -mfma -march=native \\"
    echo "      catcharea_level0.c matrix.c matrix_inv.c matrix_multiply_optimized.c \\"
    echo "      performance_summary.c bsolve.c area.c catchment.c co_matrix.c \\"
    echo "      -o catcharea_level0 -lm -lopenblas -llapacke"
    exit 1
fi

# Create output directory
OUTPUT_DIR="level0_test_results_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTPUT_DIR"

echo "Output directory: $OUTPUT_DIR"
echo ""

# Test configurations
declare -a CONFIGS=(
    "1:1:Sequential baseline"
    "2:1:2 zones parallel"
    "3:1:3 zones parallel"
    "6:1:6 zones parallel (Strategy 1)"
    "3:2:3 zones × 2 threads (Strategy 2 - Hybrid)"
)

echo "Test Configurations:"
echo "────────────────────────────────────────────────────────────"
for config in "${CONFIGS[@]}"; do
    IFS=':' read -r zone_threads lapack_threads desc <<< "$config"
    echo "  $desc"
    echo "    Zone threads: $zone_threads, LAPACK threads: $lapack_threads"
done
echo ""
echo "════════════════════════════════════════════════════════════"
echo ""

# Run tests
test_num=1
for config in "${CONFIGS[@]}"; do
    IFS=':' read -r zone_threads lapack_threads desc <<< "$config"
    
    echo ""
    echo "────────────────────────────────────────────────────────────"
    echo "Test $test_num: $desc"
    echo "────────────────────────────────────────────────────────────"
    echo ""
    
    # Set environment
    export OMP_NUM_THREADS=$zone_threads
    export OPENBLAS_NUM_THREADS=$lapack_threads
    export OMP_PROC_BIND=close
    export OMP_PLACES=cores
    
    echo "Environment:"
    echo "  OMP_NUM_THREADS=$OMP_NUM_THREADS"
    echo "  OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS"
    echo ""
    
    # Output file
    output_file="$OUTPUT_DIR/test${test_num}_${zone_threads}zones_${lapack_threads}threads.txt"
    timing_file="$OUTPUT_DIR/test${test_num}_timing.txt"
    
    # Run test
    echo "Running test..."
    { time ./catcharea_level0 1.0 99.0 0.001 0 3 32 > "$output_file" 2>&1; } 2> "$timing_file"
    
    # Extract results
    if [ -f "$output_file" ]; then
        catchment_area=$(grep "Catchment area:" "$output_file" | tail -1 | awk '{print $NF}')
        bem_time=$(grep "BEM computation time:" "$output_file" | tail -1 | awk '{print $4}')
        total_time=$(grep "Total execution time:" "$output_file" | tail -1 | awk '{print $4}')
        
        echo ""
        echo "Results:"
        echo "  Catchment area: $catchment_area"
        echo "  BEM time: ${bem_time}s"
        echo "  Total time: ${total_time}s"
        
        # Save to summary CSV
        if [ $test_num -eq 1 ]; then
            echo "TestNum,ZoneThreads,LapackThreads,Description,CatchmentArea,BEMTime,TotalTime" > "$OUTPUT_DIR/summary.csv"
        fi
        echo "$test_num,$zone_threads,$lapack_threads,\"$desc\",$catchment_area,$bem_time,$total_time" >> "$OUTPUT_DIR/summary.csv"
    else
        echo "ERROR: Output file not created"
    fi
    
    test_num=$((test_num + 1))
    
    # Wait between tests
    sleep 2
done

echo ""
echo "════════════════════════════════════════════════════════════"
echo "  ALL TESTS COMPLETED"
echo "════════════════════════════════════════════════════════════"
echo ""

# Generate summary report
if [ -f "$OUTPUT_DIR/summary.csv" ]; then
    echo "Performance Summary:"
    echo "────────────────────────────────────────────────────────────"
    echo ""
    
    # Read baseline
    baseline=$(sed -n '2p' "$OUTPUT_DIR/summary.csv" | cut -d',' -f7)
    
    # Print table
    printf "%-8s %-10s %-12s %-10s %-10s %-8s\n" \
        "Test" "Zones" "LAPACK" "BEM(s)" "Total(s)" "Speedup"
    printf "%-8s %-10s %-12s %-10s %-10s %-8s\n" \
        "────" "──────" "────────────" "────────" "────────" "───────"
    
    tail -n +2 "$OUTPUT_DIR/summary.csv" | while IFS=',' read -r num zones lapack desc area bem total; do
        if [ -n "$baseline" ] && [ -n "$total" ]; then
            speedup=$(awk -v b="$baseline" -v t="$total" 'BEGIN {printf "%.2fx", b/t}')
        else
            speedup="N/A"
        fi
        
        printf "%-8s %-10s %-12s %-10.2f %-10.2f %-8s\n" \
            "$num" "$zones" "$lapack" "$bem" "$total" "$speedup"
    done
    
    echo ""
    echo "────────────────────────────────────────────────────────────"
    echo ""
    echo "Detailed results saved to: $OUTPUT_DIR/"
    echo "  - summary.csv: Performance data"
    echo "  - test*.txt: Full output for each test"
    echo "  - test*_timing.txt: Wall-clock timing"
    echo ""
    
    # Best configuration
    best_test=$(tail -n +2 "$OUTPUT_DIR/summary.csv" | \
                awk -F',' 'NR==1 {min=$7; best=$1} $7<min {min=$7; best=$1} END {print best}')
    
    if [ -n "$best_test" ]; then
        best_config=$(sed -n "${best_test}p" "$OUTPUT_DIR/summary.csv")
        IFS=',' read -r num zones lapack desc area bem total <<< "$best_config"
        
        echo "╔══════════════════════════════════════════════════════════╗"
        echo "║  BEST CONFIGURATION                                      ║"
        echo "╚══════════════════════════════════════════════════════════╝"
        echo ""
        echo "  Test $num: $desc"
        echo "  Zone threads: $zones"
        echo "  LAPACK threads: $lapack"
        echo "  Total time: ${total}s"
        echo ""
        
        if [ -n "$baseline" ]; then
            speedup=$(awk -v b="$baseline" -v t="$total" 'BEGIN {printf "%.2fx", b/t}')
            echo "  Speedup vs baseline: $speedup"
            
            time_saved=$(awk -v b="$baseline" -v t="$total" 'BEGIN {printf "%.0f", b-t}')
            echo "  Time saved: ${time_saved}s"
            echo ""
        fi
    fi
fi

echo "════════════════════════════════════════════════════════════"
echo ""

exit 0
