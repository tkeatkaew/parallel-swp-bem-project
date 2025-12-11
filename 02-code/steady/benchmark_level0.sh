#!/usr/bin/env bash
#
# benchmark_level0.sh
#
# Comprehensive benchmark script for Level 0 Contour-Layer Parallelism
# Tests zone-level parallelism with different thread counts
#
# Usage: ./benchmark_level0.sh [NUM_THREADS] [ARG1 ARG2 ARG3]
#
# Author: Generated for BEM SWP Research
# Date: 2025-11-14
# Version: 1.0

set -u
set -e

#══════════════════════════════════════════════════════════════════════════
# Configuration
#══════════════════════════════════════════════════════════════════════════

NUM_THREADS="${1:-6}"
ARG1="${2:-1.0}"
ARG2="${3:-99.0}"
ARG3="${4:-0.001}"
INVERSION_METHOD=1  # 1=Sequential LAPACK (recommended for Level 0)
MULTIPLY_METHOD=3   # 3=Full optimization (OpenMP+Cache+SIMD)
BLOCK_SIZE=32

# Output directory
OUTPUT_DIR="benchmark_level0_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTPUT_DIR"

# Logging
LOG_FILE="$OUTPUT_DIR/benchmark.log"
exec > >(tee -a "$LOG_FILE") 2>&1

#══════════════════════════════════════════════════════════════════════════
# Helper functions
#══════════════════════════════════════════════════════════════════════════

ts() { printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*"; }

print_header() {
    echo ""
    echo "═══════════════════════════════════════════════════════════════════"
    echo "  $1"
    echo "═══════════════════════════════════════════════════════════════════"
    echo ""
}

#══════════════════════════════════════════════════════════════════════════
# Setup environment
#══════════════════════════════════════════════════════════════════════════

setup_environment() {
    # CRITICAL: Set OPENBLAS to use 1 thread
    # This prevents oversubscription when we parallelize zones
    export OPENBLAS_NUM_THREADS=1
    export MKL_NUM_THREADS=1
    
    # OpenMP configuration for zone-level parallelism
    export OMP_NUM_THREADS="$NUM_THREADS"
    export OMP_PROC_BIND=close
    export OMP_PLACES=cores
    export OMP_NESTED=FALSE  # Disable nested parallelism
    
    # Memory optimization
    export MALLOC_ARENA_MAX=4
    
    ts "Environment configured:"
    ts "  OMP_NUM_THREADS:      $OMP_NUM_THREADS"
    ts "  OPENBLAS_NUM_THREADS: 1 (to prevent oversubscription)"
    ts "  OMP_NESTED:           FALSE"
}

#══════════════════════════════════════════════════════════════════════════
# Run single test
#══════════════════════════════════════════════════════════════════════════

run_test() {
    local num_threads=$1
    local test_name=$2
    
    export OMP_NUM_THREADS="$num_threads"
    
    local output_file="$OUTPUT_DIR/${test_name}_t${num_threads}.txt"
    local timing_file="$OUTPUT_DIR/timing_${test_name}_t${num_threads}.txt"
    local csv_file="$OUTPUT_DIR/${test_name}_t${num_threads}_zones.csv"
    
    ts "Testing: $test_name with $num_threads threads"
    ts "  Command: ./catcharea_level0 $ARG1 $ARG2 $ARG3 $INVERSION_METHOD $MULTIPLY_METHOD $BLOCK_SIZE"
    
    # Run test
    { time ./catcharea_level0 $ARG1 $ARG2 $ARG3 $INVERSION_METHOD $MULTIPLY_METHOD $BLOCK_SIZE \
        > "$output_file" 2>&1; } 2> "$timing_file"
    
    local exit_code=$?
    
    if [ $exit_code -eq 0 ]; then
        # Extract results
        local area=$(awk '/Total catchment area:/{print $NF}' "$output_file")
        local parallel_time=$(awk '/Level 0.*parallel processing time:/{print $(NF-1)}' "$output_file")
        
        # Get wall-clock time from timing file
        local real_time=$(grep "^real" "$timing_file" | awk '{print $2}')
        local real_time_sec=$(echo "$real_time" | awk -F'[ms]' '{
          if (NF==3) print $1*60 + $2;
          else if ($0 ~ /m/) print $1*60 + $2;
          else print $1;
        }')
        
        ts "  Status: SUCCESS ✅"
        ts "  Catchment area:    $area"
        ts "  Parallel time:     ${parallel_time}s"
        ts "  Wall-clock time:   ${real_time_sec}s"
        
        # Save to results CSV
        echo "$num_threads|$test_name|$area|$parallel_time|$real_time_sec" \
            >> "$OUTPUT_DIR/results.csv"
        
        # Copy per-zone CSV if it exists
        if [ -f "zone_results.csv" ]; then
            mv "zone_results.csv" "$csv_file"
        fi
        
        return 0
    else
        ts "  Status: FAILED ❌"
        return 1
    fi
}

#══════════════════════════════════════════════════════════════════════════
# Run baseline comparison
#══════════════════════════════════════════════════════════════════════════

run_baseline() {
    print_header "BASELINE: Sequential Zone Processing"
    
    ts "Running original catcharea for comparison..."
    
    if [ ! -f "./catcharea" ]; then
        ts "Warning: ./catcharea not found, skipping baseline"
        return 1
    fi
    
    export OMP_NUM_THREADS=6
    export OPENBLAS_NUM_THREADS=6
    
    local output_file="$OUTPUT_DIR/baseline_sequential.txt"
    local timing_file="$OUTPUT_DIR/timing_baseline.txt"
    
    { time ./catcharea $ARG1 $ARG2 $ARG3 $INVERSION_METHOD $MULTIPLY_METHOD $BLOCK_SIZE \
        > "$output_file" 2>&1; } 2> "$timing_file"
    
    local area=$(awk '/Catchment area:/{val=$NF} END{print val}' "$output_file")
    local real_time=$(grep "^real" "$timing_file" | awk '{print $2}')
    
    ts "Baseline results:"
    ts "  Catchment area:  $area"
    ts "  Time:            $real_time"
    
    # Save baseline for comparison
    echo "baseline|sequential_zones|$area|0|$real_time" \
        >> "$OUTPUT_DIR/results.csv"
    
    return 0
}

#══════════════════════════════════════════════════════════════════════════
# Main benchmark routine
#══════════════════════════════════════════════════════════════════════════

main() {
    print_header "LEVEL 0 CONTOUR-LAYER PARALLELISM BENCHMARK"
    
    ts "Configuration:"
    ts "  Max threads:       $NUM_THREADS"
    ts "  Parameters:        $ARG1 $ARG2 $ARG3"
    ts "  Inversion method:  $INVERSION_METHOD (Sequential LAPACK)"
    ts "  Multiply method:   $MULTIPLY_METHOD (Full optimization)"
    ts "  Output directory:  $OUTPUT_DIR"
    echo ""
    
    # Check executable
    if [ ! -f "./catcharea_level0" ]; then
        ts "Error: ./catcharea_level0 not found!"
        ts "Please build the project first:"
        ts "  make -f Makefile.level0 clean"
        ts "  make -f Makefile.level0 catcharea_level0"
        exit 1
    fi
    
    # Setup environment
    setup_environment
    echo ""
    
    # Run baseline if available
    # run_baseline || ts "Baseline comparison skipped"
    # echo ""
    
    # Initialize results CSV
    echo "Threads|Name|Area|Parallel_Time|Wall_Clock" > "$OUTPUT_DIR/results.csv"
    
    # Test different thread counts
    print_header "SCALABILITY TESTS"
    
    local thread_counts=(1 2 3 4 6)
    local all_success=true
    
    for nt in "${thread_counts[@]}"; do
        if [ $nt -le $NUM_THREADS ]; then
            echo ""
            if ! run_test $nt "level0_${nt}threads"; then
                all_success=false
            fi
            echo ""
            sleep 2
        fi
    done
    
    # Analyze results
    print_header "BENCHMARK RESULTS"
    
    if [ ! -f "$OUTPUT_DIR/results.csv" ]; then
        ts "Error: Results file not created"
        exit 1
    fi
    
    ts "Reading results..."
    echo ""
    
    # Display results table
    printf "%-15s %12s %18s %18s %12s\n" \
        "Configuration" "Threads" "Area" "Parallel_Time(s)" "Speedup"
    printf "%-15s %12s %18s %18s %12s\n" \
        "---------------" "-------" "----" "-----------------" "-------"
    
    # Get baseline time (1 thread)
    local baseline_line=$(grep "^1|" "$OUTPUT_DIR/results.csv" | tail -1)
    local baseline_time=$(echo "$baseline_line" | cut -d'|' -f5)
    
    if [ -z "$baseline_time" ]; then
        ts "Warning: Could not extract baseline time"
        baseline_time=1.0
    fi
    
    ts "Baseline (1 thread): ${baseline_time}s"
    echo ""
    
    while IFS='|' read -r threads name area parallel_time wall_time; do
        if [ "$threads" = "Threads" ]; then continue; fi
        
        # Calculate speedup
        if [ -n "$wall_time" ] && [ -n "$baseline_time" ] && [ "$wall_time" != "N/A" ]; then
            speedup=$(awk -v b="$baseline_time" -v t="$wall_time" 'BEGIN {
              if (t > 0) printf "%.2f", b/t
              else print "N/A"
            }')
        else
            speedup="N/A"
        fi
        
        printf "%-15s %12s %18s %18s %11sx\n" \
            "$name" "$threads" "$area" "$wall_time" "$speedup"
    done < "$OUTPUT_DIR/results.csv"
    
    echo ""
    
    # Performance analysis
    print_header "PERFORMANCE ANALYSIS"
    
    ts "Speedup Analysis:"
    ts "  Expected (Amdahl's Law):"
    ts "    2 threads → 2.0× speedup"
    ts "    4 threads → 4.0× speedup"
    ts "    6 threads → 6.0× speedup"
    echo ""
    
    ts "Efficiency Analysis:"
    for nt in "${thread_counts[@]}"; do
        if [ $nt -le 1 ]; then continue; fi
        
        local line=$(grep "^${nt}|" "$OUTPUT_DIR/results.csv" | tail -1)
        if [ -z "$line" ]; then continue; fi
        
        local time=$(echo "$line" | cut -d'|' -f5)
        if [ -z "$time" ] || [ "$time" = "N/A" ]; then continue; fi
        
        local speedup=$(awk -v b="$baseline_time" -v t="$time" 'BEGIN {
            if (t > 0) printf "%.2f", b/t
            else print "0"
        }')
        
        local efficiency=$(awk -v s="$speedup" -v n="$nt" 'BEGIN {
            printf "%.1f", (s/n)*100
        }')
        
        ts "  $nt threads: Speedup=${speedup}×, Efficiency=${efficiency}%"
    done
    echo ""
    
    # Memory usage
    ts "Memory Usage:"
    for nt in "${thread_counts[@]}"; do
        local output_file="$OUTPUT_DIR/level0_${nt}threads_t${nt}.txt"
        if [ -f "$output_file" ]; then
            local peak_mem=$(grep "Peak memory" "$output_file" | tail -1 | awk '{print $(NF-1), $NF}')
            if [ -n "$peak_mem" ]; then
                ts "  $nt threads: $peak_mem"
            fi
        fi
    done
    echo ""
    
    # Correctness verification
    print_header "CORRECTNESS VERIFICATION"
    
    ts "Checking consistency of catchment areas..."
    
    local areas=$(cut -d'|' -f3 "$OUTPUT_DIR/results.csv" | tail -n +2)
    local first_area=$(echo "$areas" | head -1)
    local all_match=true
    
    while read -r area; do
        local diff=$(awk -v a="$area" -v b="$first_area" 'BEGIN {
            printf "%.10f", (a-b)>0?(a-b):(b-a)
        }')
        local is_close=$(awk -v d="$diff" 'BEGIN {
            print (d < 0.001) ? "yes" : "no"
        }')
        
        if [ "$is_close" = "no" ]; then
            all_match=false
            ts "  ⚠️  Warning: Area mismatch: $area vs $first_area (diff: $diff)"
        fi
    done <<< "$areas"
    
    if $all_match; then
        ts "  ✅ All configurations produce identical results!"
    else
        ts "  ❌ Some configurations produce different results - investigate!"
    fi
    echo ""
    
    # Final summary
    print_header "BENCHMARK COMPLETE"
    
    ts "Results saved to: $OUTPUT_DIR/"
    ts "  - benchmark.log:           Full log"
    ts "  - results.csv:             Summary table"
    ts "  - level0_*threads_t*.txt:  Detailed output"
    ts "  - timing_*.txt:            Timing information"
    ts "  - *_zones.csv:             Per-zone statistics"
    echo ""
    
    if $all_success; then
        ts "All tests completed successfully! ✅"
    else
        ts "Some tests failed - check log files ⚠️"
    fi
    echo ""
    
    # Recommendations
    print_header "RECOMMENDATIONS"
    
    local best_line=$(grep -v "^Threads" "$OUTPUT_DIR/results.csv" | \
                      sort -t'|' -k5 -n | head -1)
    local best_threads=$(echo "$best_line" | cut -d'|' -f1)
    local best_time=$(echo "$best_line" | cut -d'|' -f5)
    
    ts "Optimal configuration: $best_threads threads"
    ts "  Best time: ${best_time}s"
    
    local final_speedup=$(awk -v b="$baseline_time" -v t="$best_time" 'BEGIN {
        printf "%.2f", b/t
    }')
    ts "  Speedup: ${final_speedup}× vs single thread"
    echo ""
    
    ts "💡 For production use:"
    ts "   export OMP_NUM_THREADS=$best_threads"
    ts "   export OPENBLAS_NUM_THREADS=1"
    ts "   ./catcharea_level0 $ARG1 $ARG2 $ARG3 $INVERSION_METHOD $MULTIPLY_METHOD $BLOCK_SIZE"
    echo ""
}

#══════════════════════════════════════════════════════════════════════════
# Run main
#══════════════════════════════════════════════════════════════════════════

main

exit 0
