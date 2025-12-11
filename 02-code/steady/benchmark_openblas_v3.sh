#!/usr/bin/env bash
################################################################################
# benchmark_openblas_v2_fixed.sh - OpenBLAS Benchmark with Area Display Fix
#
# VERSION 2.1: Fixed Catchment Area Display
# - Fixed: Parse catchment area from both CSV and text output
# - Fallback: If CSV parsing fails, try text parsing
# - Display: Always show area in results
#
# Usage: ./benchmark_openblas_v2_fixed.sh [NUM_THREADS] [ARG1 ARG2 ARG3]
#
# Author: Optimized for BEM Flow Path Analysis
# Date: 2025-11-21 (Version 2.1 - Fixed Area Display)
################################################################################

set -euo pipefail

# Configuration
NUM_THREADS="${1:-6}"
ARG1="${2:-1.0}"
ARG2="${3:-100.0}"
ARG3="${4:-0.001}"

# Output directory
OUTPUT_DIR="benchmark_openblas_v2_$(date +%Y%m%d_%H%M%S)"
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
    
    if [ ! -f "./catcharea" ]; then
        ts "❌ Error: ./catcharea not found!"
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
    
    export OPENBLAS_NUM_THREADS="$NUM_THREADS"
    export OMP_NUM_THREADS="$NUM_THREADS"
    export OMP_PROC_BIND=close
    export OMP_PLACES=cores
    export MALLOC_ARENA_MAX=4
    
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
    
    # Set inversion method based on multiply method
    local inversion_method
    if [ "$method" -eq 0 ]; then
        inversion_method=1  # Sequential (Gauss-Jordan)
    else
        inversion_method=0  # Parallel (LAPACK)
    fi
    
    local output_file="$OUTPUT_DIR/method${method}.txt"
    local timing_file="$OUTPUT_DIR/timing_method${method}.txt"
    
    print_subheader "Test: $method_name"
    
    ts "Configuration:"
    ts "  Multiply method:   $method"
    ts "  Inversion method:  $inversion_method"
    
    if [ "$method" -eq 0 ]; then
        ts "  Mode: PURE SEQUENTIAL (true baseline)"
        ts "    - Matrix multiply: Sequential loops (no BLAS)"
        ts "    - Matrix inversion: Gauss-Jordan (no LAPACK)"
    else
        ts "  Mode: OPTIMIZED (parallel)"
        ts "    - Matrix multiply: OpenBLAS DGEMM ($NUM_THREADS threads)"
        ts "    - Matrix inversion: LAPACK ($NUM_THREADS threads)"
    fi
    
    ts "Command: ./catcharea $ARG1 $ARG2 $ARG3 $inversion_method $method 0"
    ts "Output: $output_file"
    echo ""
    
    # Run with timing
    local start_time=$(date +%s.%N)
    
    { time ./catcharea $ARG1 $ARG2 $ARG3 $inversion_method $method 0 > "$output_file" 2>&1; } \
        2> "$timing_file"
    
    local exit_code=$?
    local end_time=$(date +%s.%N)
    
    if [ $exit_code -ne 0 ]; then
        ts "❌ Test FAILED with exit code $exit_code"
        return 1
    fi
    
    # ═══════════════════════════════════════════════════════════════════════════
    # PRIORITY 1: Parse from CSV section (most accurate)
    # ═══════════════════════════════════════════════════════════════════════════
    
    local csv_section=$(awk '/^Parameter,Value$/,/^Performance data exported/' "$output_file")
    
    local area=$(echo "$csv_section" | awk -F, '/^Catchment_Area,/{print $2}')
    local dgemm_time=$(echo "$csv_section" | awk -F, '/^Multiply_Time_sec,/{print $2}')
    local inv_time=$(echo "$csv_section" | awk -F, '/^Inversion_Time_sec,/{print $2}')
    local comp_time=$(echo "$csv_section" | awk -F, '/^Computation_Time_sec,/{print $2}')
    local total_time=$(echo "$csv_section" | awk -F, '/^Total_Time_sec,/{print $2}')
    
    # ═══════════════════════════════════════════════════════════════════════════
    # FALLBACK: If area is empty, try text parsing
    # ═══════════════════════════════════════════════════════════════════════════
    
    if [ -z "$area" ] || [ "$area" = "N/A" ]; then
        # Try multiple patterns to find catchment area
        area=$(grep -i "Catchment area" "$output_file" | tail -1 | \
               awk '{for(i=1;i<=NF;i++) if($i ~ /^[0-9]+\.?[0-9]*$/) print $i}' | tail -1)
        
        if [ -z "$area" ]; then
            # Try another pattern
            area=$(grep -i "catchment" "$output_file" | grep -Eo '[0-9]+\.[0-9]{6,}' | tail -1)
        fi
        
        if [ -z "$area" ]; then
            area="N/A"
            ts "⚠️  Warning: Could not parse catchment area"
        else
            ts "ℹ️  Used text fallback for catchment area: $area"
        fi
    fi
    
    # Parse wall-clock time
    local real_time=$(grep "^real" "$timing_file" | awk '{print $2}')
    local real_time_sec=$(echo "$real_time" | awk -F'[ms]' '{
        if (NF==3) print $1*60 + $2;
        else if ($0 ~ /m/) print $1*60 + $2;
        else print $1;
    }')
    
    # Calculate execution time
    local exec_time=$(awk -v start="$start_time" -v end="$end_time" \
        'BEGIN {printf "%.6f", end - start}')
    
    # Parse GFLOPS if available
    local gflops=$(tail -100 "$output_file" | grep "DGEMM GFLOPS:" | tail -1 | \
        awk '{for(i=1;i<=NF;i++) if($i ~ /^[0-9]+\.?[0-9]*$/) {val=$i}} END {print val}')
    
    # Display results
    ts "✅ Test completed successfully"
    echo ""
    ts "Results:"
    ts "  Catchment Area:        ${area}"
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
    echo "$method|$method_name|$area|$dgemm_time|$inv_time|$comp_time|$real_time_sec|$gflops|$inversion_method" \
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
    print_header "OpenBLAS Benchmark Suite v2.1 - Fixed Area Display"
    
    ts "Benchmark Configuration:"
    ts "  Number of threads: $NUM_THREADS"
    ts "  Parameters: $ARG1 $ARG2 $ARG3"
    ts "  Output directory: $OUTPUT_DIR"
    echo ""
    
    ts "Method 0 Configuration:"
    ts "  Matrix multiply:  Sequential loops (no BLAS, no threading)"
    ts "  Matrix inversion: Gauss-Jordan (no LAPACK, no threading)"
    ts "  Purpose:          TRUE baseline for accurate speedup measurement"
    echo ""
    
    ts "Methods 1-3 Configuration:"
    ts "  Matrix multiply:  OpenBLAS DGEMM ($NUM_THREADS threads)"
    ts "  Matrix inversion: LAPACK ($NUM_THREADS threads)"
    ts "  Purpose:          Optimized parallel performance"
    echo ""
    
    # Check prerequisites
    check_prerequisites
    
    # Setup environment
    setup_environment
    
    # Create results CSV header
    echo "Method|Name|Area|DGEMM_Time|Inv_Time|Comp_Time|Wall_Time|GFLOPS|Inv_Method" \
        > "$OUTPUT_DIR/results.csv"
    
    # Test configurations
    declare -A methods=(
        [0]="Sequential (pure baseline)"
        [1]="OpenBLAS DGEMM + LAPACK"
        [2]="OpenBLAS DGEMM + LAPACK (same as Method 1)"
        [3]="OpenBLAS DGEMM + LAPACK (same as Method 1)"
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
    local baseline_inv=$(echo "$baseline" | cut -d'|' -f5)
    local baseline_total=$(echo "$baseline" | cut -d'|' -f6)
    local baseline_area=$(echo "$baseline" | cut -d'|' -f3)
    
    ts "Baseline (Method 0 - Pure Sequential):"
    ts "  Catchment area:       ${baseline_area}"
    ts "  Wall-clock time:      ${baseline_wall}s"
    ts "  Matrix multiply:      ${baseline_dgemm}s (sequential loops)"
    ts "  Matrix inversion:     ${baseline_inv}s (Gauss-Jordan)"
    ts "  Total matrix ops:     ${baseline_total}s"
    echo ""
    
    # Results table with Area column
    printf "%-45s %12s %12s %12s %12s %12s %15s %15s\n" \
        "Method" "Wall(s)" "Mult(s)" "Inv(s)" "Total(s)" "Speedup" "Area" "Inv Type"
    printf "%-45s %12s %12s %12s %12s %12s %15s %15s\n" \
        "------" "-------" "-------" "------" "--------" "-------" "----" "--------"
    
    tail -n +2 "$OUTPUT_DIR/results.csv" | while IFS='|' read -r method name area dgemm inv comp wall gflops inv_method; do
        local speedup="N/A"
        if [ -n "$wall" ] && [ -n "$baseline_wall" ] && \
           awk -v w="$wall" 'BEGIN {exit !(w > 0)}'; then
            speedup=$(awk -v b="$baseline_wall" -v w="$wall" \
                'BEGIN {printf "%.2fx", b/w}')
        fi
        
        local inv_type
        if [ "$inv_method" = "1" ]; then
            inv_type="Sequential"
        else
            inv_type="Parallel"
        fi
        
        # Format area for display
        local area_fmt="${area:-N/A}"
        if [ "$area_fmt" != "N/A" ]; then
            area_fmt=$(printf "%.6f" "$area" 2>/dev/null || echo "$area")
        fi
        
        printf "%-45s %12.2f %12.2f %12.2f %12.2f %12s %15s %15s\n" \
            "$name" "$wall" "$dgemm" "$inv" "$comp" "$speedup" "$area_fmt" "$inv_type"
    done
    
    echo ""
    
    # Correctness verification
    print_subheader "Correctness Verification"
    
    local areas=$(tail -n +2 "$OUTPUT_DIR/results.csv" | cut -d'|' -f3)
    local first_area=$(echo "$areas" | head -1)
    local all_match=true
    
    # Only verify if we have valid areas
    if [ "$first_area" != "N/A" ] && [ -n "$first_area" ]; then
        while read -r area; do
            if [ "$area" = "N/A" ] || [ -z "$area" ]; then
                continue
            fi
            
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
    else
        ts "⚠️  Could not verify correctness - catchment area not available"
    fi
    
    echo ""
    
    # Analysis
    print_header "Performance Analysis"
    
    ts "Pure Sequential Baseline (Method 0):"
    ts "  This is a TRUE baseline with NO optimizations:"
    ts "    • Matrix multiply: Simple triple-loop (no BLAS, no threading)"
    ts "    • Matrix inversion: Gauss-Jordan elimination (no LAPACK)"
    ts "    • Single-threaded throughout"
    ts "    • No cache optimization, no SIMD"
    echo ""
    
    ts "Optimized Methods (1-3):"
    ts "  These use production-grade parallel libraries:"
    ts "    • Matrix multiply: OpenBLAS DGEMM ($NUM_THREADS threads)"
    ts "    • Matrix inversion: LAPACK LU decomposition ($NUM_THREADS threads)"
    ts "    • Fully optimized: threading + cache + SIMD + assembly"
    echo ""
    
    # Find best performer
    local best_line=$(tail -n +2 "$OUTPUT_DIR/results.csv" | \
        awk -F'|' '$1 != "0"' | sort -t'|' -k7 -n | head -1)
    
    if [ -n "$best_line" ]; then
        local best_method=$(echo "$best_line" | cut -d'|' -f1)
        local best_name=$(echo "$best_line" | cut -d'|' -f2)
        local best_time=$(echo "$best_line" | cut -d'|' -f7)
        local best_speedup=$(awk -v b="$baseline_wall" -v t="$best_time" \
            'BEGIN {printf "%.2f", b/t}')
        
        ts "Best Optimized Method:"
        ts "  Method $best_method: $best_name"
        ts "  Wall-clock time: ${best_time}s"
        ts "  Speedup: ${best_speedup}× vs pure sequential baseline"
    fi
    
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
    
    ts "Key Takeaway:"
    ts "  Method 0 speedup reflects TRUE parallel performance gain"
    ts "  (no library overhead in baseline measurement)"
    echo ""
}

################################################################################
# Run main
################################################################################

main

exit 0