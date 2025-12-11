#!/usr/bin/env bash
################################################################################
# benchmark_openblas.sh - Benchmark OpenBLAS-based Matrix Multiplication
#
# This script tests the OpenBLAS implementation of matrix operations.
# Methods 1-3 all use the same OpenBLAS DGEMM implementation.
#
# Usage: ./benchmark_openblas.sh [NUM_THREADS] [ARG1 ARG2 ARG3]
#
# Example:
#   ./benchmark_openblas.sh 6 1.0 100.0 0.001
#
# Author: Optimized for BEM Flow Path Analysis
# Date: 2025-11-21
################################################################################

set -euo pipefail

# Configuration
NUM_THREADS="${1:-6}"
ARG1="${2:-1.0}"
ARG2="${3:-100.0}"
ARG3="${4:-0.001}"
INVERSION_METHOD=0  # Always use LAPACK (parallel)

# Output directory
OUTPUT_DIR="benchmark_openblas_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTPUT_DIR"

# Logging
LOG_FILE="$OUTPUT_DIR/benchmark.log"
exec > >(tee -a "$LOG_FILE") 2>&1

################################################################################
# Helper functions
################################################################################

ts() { printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*"; }

print_header() {
    echo ""
    echo "════════════════════════════════════════════════════════════════════════════"
    echo "  $1"
    echo "════════════════════════════════════════════════════════════════════════════"
    echo ""
}

print_subheader() {
    echo ""
    echo "────────────────────────────────────────────────────────────────────────────"
    echo "  $1"
    echo "────────────────────────────────────────────────────────────────────────────"
    echo ""
}

################################################################################
# Check prerequisites
################################################################################

check_prerequisites() {
    print_subheader "Checking Prerequisites"
    
    local missing=0
    
    # Check executable
    if [ ! -f "./catcharea" ]; then
        ts "❌ Error: ./catcharea not found!"
        ts "   Please build first: make clean && make"
        missing=1
    else
        ts "✅ Executable found: ./catcharea"
    fi
    
   
    echo ""
}

################################################################################
# Setup environment
################################################################################

setup_environment() {
    print_subheader "Environment Configuration"
    
    # OpenBLAS threading (critical for performance!)
    export OPENBLAS_NUM_THREADS="$NUM_THREADS"
    
    # OpenMP threading (for our custom code sections)
    export OMP_NUM_THREADS="$NUM_THREADS"
    export OMP_PROC_BIND=close
    export OMP_PLACES=cores
    
    # Memory optimization
    export MALLOC_ARENA_MAX=4
    
    # Disable CPU frequency scaling for consistent results
    export CPU_FREQ_GOV_SAVED=$(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor 2>/dev/null || echo "unknown")
    
    ts "Configuration:"
    ts "  OPENBLAS_NUM_THREADS: $OPENBLAS_NUM_THREADS"
    ts "  OMP_NUM_THREADS:      $OMP_NUM_THREADS"
    ts "  OMP_PROC_BIND:        $OMP_PROC_BIND"
    ts "  OMP_PLACES:           $OMP_PLACES"
    ts "  MALLOC_ARENA_MAX:     $MALLOC_ARENA_MAX"
    
    echo ""
}

################################################################################
# Run single test
################################################################################

run_test() {
    local method=$1
    local method_name=$2
    
    local output_file="$OUTPUT_DIR/method${method}.txt"
    local timing_file="$OUTPUT_DIR/timing_method${method}.txt"
    
    print_subheader "Test: $method_name"
    
    ts "Command: ./catcharea $ARG1 $ARG2 $ARG3 $INVERSION_METHOD $method 0"
    ts "Output: $output_file"
    echo ""
    
    # Run with timing
    local start_time=$(date +%s.%N)
    
    { time ./catcharea $ARG1 $ARG2 $ARG3 $INVERSION_METHOD $method 0 > "$output_file" 2>&1; } \
        2> "$timing_file"
    
    local exit_code=$?
    local end_time=$(date +%s.%N)
    
    if [ $exit_code -ne 0 ]; then
        ts "❌ Test FAILED with exit code $exit_code"
        return 1
    fi
    
    # Parse results from CSV section
    local csv_section=$(awk '/^Parameter,Value$/,/^Performance data exported/' "$output_file")
    
    local area=$(echo "$csv_section" | awk -F, '/^Catchment_Area,/{print $2}')
    local dgemm_time=$(echo "$csv_section" | awk -F, '/^Multiply_Time_sec,/{print $2}')
    local inv_time=$(echo "$csv_section" | awk -F, '/^Inversion_Time_sec,/{print $2}')
    local comp_time=$(echo "$csv_section" | awk -F, '/^Computation_Time_sec,/{print $2}')
    local total_time=$(echo "$csv_section" | awk -F, '/^Total_Time_sec,/{print $2}')
    
    # Parse wall-clock time
    local real_time=$(grep "^real" "$timing_file" | awk '{print $2}')
    local real_time_sec=$(echo "$real_time" | awk -F'[ms]' '{
        if (NF==3) print $1*60 + $2;
        else if ($0 ~ /m/) print $1*60 + $2;
        else print $1;
    }')
    
    # Calculate execution time from timestamps
    local exec_time=$(awk -v start="$start_time" -v end="$end_time" \
        'BEGIN {printf "%.6f", end - start}')
    
    # Parse GFLOPS if available
    local gflops=$(tail -100 "$output_file" | grep "DGEMM GFLOPS:" | tail -1 | \
        awk '{for(i=1;i<=NF;i++) if($i ~ /^[0-9]+\.?[0-9]*$/) {val=$i}} END {print val}')
    
    # Display results
    ts "✅ Test completed successfully"
    echo ""
    ts "Results:"
    ts "  Catchment Area:        $area"
    ts "  Matrix Multiply Time:  ${dgemm_time}s"
    ts "  Matrix Inversion Time: ${inv_time}s"
    ts "  Matrix Ops Total:      ${comp_time}s"
    ts "  Total Runtime:         ${total_time}s"
    ts "  Wall-Clock Time:       ${real_time_sec}s"
    
    if [ -n "$gflops" ]; then
        ts "  DGEMM Performance:     ${gflops} GFLOPS"
    fi
    
    echo ""
    
    # Save to results CSV
    echo "$method|$method_name|$area|$dgemm_time|$inv_time|$comp_time|$real_time_sec|$gflops" \
        >> "$OUTPUT_DIR/results.csv"
    
    # Copy performance CSV if exists
    if [ -f "performance_results.csv" ]; then
        cp -f "performance_results.csv" "$OUTPUT_DIR/method${method}_perf.csv" 2>/dev/null || true
    fi
    
    return 0
}

################################################################################
# Main benchmark routine
################################################################################

main() {
    print_header "OpenBLAS Matrix Multiplication Benchmark Suite"
    
    ts "Benchmark Configuration:"
    ts "  Number of threads: $NUM_THREADS"
    ts "  Parameters: $ARG1 $ARG2 $ARG3"
    ts "  Output directory: $OUTPUT_DIR"
    echo ""
    
    # Check prerequisites
    check_prerequisites
    
    # Setup environment
    setup_environment
    
    # Create results CSV header
    echo "Method|Name|Area|DGEMM_Time|Inv_Time|Comp_Time|Wall_Time|GFLOPS" \
        > "$OUTPUT_DIR/results.csv"
    
    # Test configurations
    declare -A methods=(
        [0]="Sequential (baseline)"
        [1]="OpenBLAS DGEMM"
        [2]="OpenBLAS DGEMM (same as Method 1)"
        [3]="OpenBLAS DGEMM (same as Method 1)"
    )
    
    local all_success=true
    
    # Run tests
    for method in 0 1 2 3; do
        if ! run_test $method "${methods[$method]}"; then
            all_success=false
        fi
        
        # Small delay between tests
        if [ $method -lt 3 ]; then
            sleep 2
        fi
    done
    
    # Generate comparison report
    print_header "Performance Comparison"
    
    # Read baseline
    local baseline=$(grep "^0|" "$OUTPUT_DIR/results.csv")
    local baseline_wall=$(echo "$baseline" | cut -d'|' -f7)
    local baseline_dgemm=$(echo "$baseline" | cut -d'|' -f4)
    local baseline_area=$(echo "$baseline" | cut -d'|' -f3)
    
    ts "Baseline (Method 0 - Sequential):"
    ts "  Wall-clock time: ${baseline_wall}s"
    ts "  DGEMM time: ${baseline_dgemm}s"
    ts "  Catchment area: ${baseline_area}"
    echo ""
    
    # Results table
    printf "%-40s %12s %12s %12s %12s %12s\n" \
        "Method" "Wall(s)" "DGEMM(s)" "Inv(s)" "Speedup" "GFLOPS"
    printf "%-40s %12s %12s %12s %12s %12s\n" \
        "------" "-------" "--------" "-----" "-------" "-------"
    
    tail -n +2 "$OUTPUT_DIR/results.csv" | while IFS='|' read -r method name area dgemm inv comp wall gflops; do
        local speedup="N/A"
        if [ -n "$wall" ] && [ -n "$baseline_wall" ] && \
           awk -v w="$wall" 'BEGIN {exit !(w > 0)}'; then
            speedup=$(awk -v b="$baseline_wall" -v w="$wall" \
                'BEGIN {printf "%.2fx", b/w}')
        fi
        
        local gflops_fmt="${gflops:-N/A}"
        
        printf "%-40s %12.2f %12.2f %12.2f %12s %12s\n" \
            "$name" "$wall" "$dgemm" "$inv" "$speedup" "$gflops_fmt"
    done
    
    echo ""
    
    # Correctness verification
    print_subheader "Correctness Verification"
    
    local areas=$(tail -n +2 "$OUTPUT_DIR/results.csv" | cut -d'|' -f3)
    local first_area=$(echo "$areas" | head -1)
    local all_match=true
    
    while read -r area; do
        local diff=$(awk -v a="$area" -v b="$first_area" \
            'BEGIN {d=a-b; printf "%.10f", (d<0)?-d:d}')
        local is_close=$(awk -v d="$diff" 'BEGIN {print (d < 0.001) ? "yes" : "no"}')
        
        if [ "$is_close" = "no" ]; then
            all_match=false
            ts "⚠️  Warning: Area mismatch: $area vs $first_area (diff: $diff)"
        fi
    done <<< "$areas"
    
    if $all_match; then
        ts "✅ All methods produce identical results (within tolerance)"
    else
        ts "❌ Some methods produce different results - check output files"
    fi
    
    echo ""
    
    # Recommendations
    print_header "Analysis & Recommendations"
    
    ts "OpenBLAS Implementation Notes:"
    ts "  • Methods 1-3 use identical OpenBLAS cblas_dgemm() implementation"
    ts "  • All optimizations (threading, cache, SIMD) are built into OpenBLAS"
    ts "  • Expected speedup: $(get_expected_speedup 1 $NUM_THREADS)× with $NUM_THREADS threads"
    echo ""
    
    ts "Performance Tuning:"
    ts "  • Adjust threads: export OPENBLAS_NUM_THREADS=<N>"
    ts "  • Thread binding: export OMP_PROC_BIND=close"
    ts "  • Thread placement: export OMP_PLACES=cores"
    echo ""
    
    # Find best performer
    local best_line=$(tail -n +2 "$OUTPUT_DIR/results.csv" | \
        sort -t'|' -k7 -n | head -1)
    local best_method=$(echo "$best_line" | cut -d'|' -f1)
    local best_name=$(echo "$best_line" | cut -d'|' -f2)
    local best_time=$(echo "$best_line" | cut -d'|' -f7)
    local best_speedup=$(awk -v b="$baseline_wall" -v t="$best_time" \
        'BEGIN {printf "%.2f", b/t}')
    
    ts "Best Performer:"
    ts "  Method $best_method: $best_name"
    ts "  Wall-clock time: ${best_time}s"
    ts "  Speedup: ${best_speedup}× vs baseline"
    echo ""
    
    # Summary
    print_header "Benchmark Complete"
    
    ts "Output files saved to: $OUTPUT_DIR/"
    ts "  • benchmark.log        - Full execution log"
    ts "  • results.csv          - Summary results table"
    ts "  • method*.txt          - Detailed output for each method"
    ts "  • timing_method*.txt   - Timing information"
    ts "  • method*_perf.csv     - Performance CSV files"
    echo ""
    
    if $all_success; then
        ts "✅ All tests completed successfully!"
    else
        ts "⚠️  Some tests failed - check log files"
    fi
    
    echo ""
}

################################################################################
# Helper function for expected speedup calculation
################################################################################

get_expected_speedup() {
    local method=$1
    local threads=$2
    
    if [ $method -eq 0 ]; then
        echo "1.00"
    else
        # OpenBLAS expected speedup
        if [ $threads -le 4 ]; then
            awk -v t=$threads 'BEGIN {printf "%.2f", t * 0.95}'
        elif [ $threads -le 8 ]; then
            awk -v t=$threads 'BEGIN {printf "%.2f", t * 0.90}'
        else
            awk -v t=$threads 'BEGIN {printf "%.2f", t * 0.80}'
        fi
    fi
}

################################################################################
# Run main
################################################################################

main

exit 0
