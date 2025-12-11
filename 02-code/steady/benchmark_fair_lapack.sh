#!/usr/bin/env bash
# CORRECTED: Fair LAPACK Benchmark Script
# Properly controls OPENBLAS_NUM_THREADS for fair comparison
# Usage: ./benchmark_fair_lapack.sh [NUM_THREADS] [ARG1 ARG2 ARG3]

set -u

# Configuration
NUM_THREADS="${1:-6}"
ARG1="${2:-1.0}"
ARG2="${3:-100.0}"
ARG3="${4:-0.001}"
MULTIPLY_METHOD=3  # Use best matrix multiply method
BLOCK_SIZE=32

# Output directory
OUTPUT_DIR="lapack_fair_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTPUT_DIR"

# Logging
LOG_FILE="$OUTPUT_DIR/benchmark.log"
exec > >(tee -a "$LOG_FILE") 2>&1

# Helper functions
ts() { printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*"; }

print_header() {
    echo ""
    echo "================================================================================"
    echo "  $1"
    echo "================================================================================"
    echo ""
}

# =============================================================================
# Fair test function - Sets OPENBLAS_NUM_THREADS correctly for each test
# =============================================================================

run_test() {
    local test_name=$1
    local openblas_threads=$2
    local output_file="$OUTPUT_DIR/${test_name}.txt"
    local timing_file="$OUTPUT_DIR/${test_name}_timing.txt"
    
    print_header "Test: $test_name"
    
    # === CRITICAL: Set thread count for THIS specific test ===
    export OMP_NUM_THREADS=$openblas_threads
    export OPENBLAS_NUM_THREADS=$openblas_threads
    export OMP_PROC_BIND=close
    export OMP_PLACES=cores
    export MALLOC_ARENA_MAX=4
    
    ts "Configuration:"
    ts "  OMP_NUM_THREADS:      $OMP_NUM_THREADS"
    ts "  OPENBLAS_NUM_THREADS: $OPENBLAS_NUM_THREADS"
    ts "  Matrix multiply:      Method $MULTIPLY_METHOD"
    ts "  Block size:           $BLOCK_SIZE"
    echo ""
    
    # INVERSION_METHOD=0 means use LAPACK (parallel-capable)
    # Parallelism controlled by OPENBLAS_NUM_THREADS
    ts "Running: ./catcharea $ARG1 $ARG2 $ARG3 0 $MULTIPLY_METHOD $BLOCK_SIZE"
    
    { time ./catcharea $ARG1 $ARG2 $ARG3 0 $MULTIPLY_METHOD $BLOCK_SIZE \
      > "$output_file" 2>&1; } 2> "$timing_file"
    
    local exit_code=$?
    
    if [ $exit_code -eq 0 ]; then
        local real_time=$(grep "^real" "$timing_file" | awk '{print $2}')
        
        # Extract inversion time
        #local inv_time=$(sed -n 's/.*Matrix inversion completed in \([0-9.]\+\) seconds.*/\1/p' "$output_file" | head -1)
        
        # Extract GFLOPS
        #local gflops=$(sed -n 's/.*Matrix inversion.*(\([0-9.]\+\) GFLOPS).*/\1/p' "$output_file" | head -1)
        
        # Extract actual OpenBLAS threads used
        #local actual_threads=$(sed -n 's/.*OpenBLAS actual threads: \([0-9]\+\).*/\1/p' "$output_file" | head -1)
        

        # Extract inversion time  (from the "MATRIX INVERSION COMPLETE" summary)
        local inv_time=$(sed -n 's/^[[:space:]]*Time:[[:space:]]*\([0-9.]\+\) seconds.*/\1/p' \
                        "$output_file" | head -1)

        # Extract GFLOPS          (same summary block)
        local gflops=$(sed -n 's/^[[:space:]]*GFLOPS:[[:space:]]*\([0-9.]\+\).*/\1/p' \
                    "$output_file" | head -1)

        # Extract actual OpenBLAS threads used (from diag in mat_inv)
        local actual_threads=$(sed -n 's/.*Actual threads in use: \([0-9]\+\).*/\1/p' \
                            "$output_file" | head -1)

        ts "Status: SUCCESS ✅"
        ts "  Total time:          $real_time"
        ts "  Inversion time:      ${inv_time}s"
        ts "  GFLOPS:              $gflops"
        [ -n "$actual_threads" ] && ts "  OpenBLAS threads:    $actual_threads"
        if [ -n "$gflops" ] && [ -n "$openblas_threads" ]; then
            local gflops_per_thread=$(awk -v g="$gflops" -v t="$openblas_threads" \
                                      'BEGIN {printf "%.2f", g/t}')
            ts "  GFLOPS/thread:       $gflops_per_thread"
        fi
        echo ""
        
        echo "$test_name|$openblas_threads|$real_time|$inv_time|$gflops|$actual_threads" \
             >> "$OUTPUT_DIR/results.csv"
        
        return 0
    else
        ts "Status: FAILED ❌"
        return 1
    fi
}

# =============================================================================
# Main benchmark routine
# =============================================================================

main() {
    print_header "FAIR LAPACK BENCHMARK: CONTROLLED THREAD COMPARISON"
    
    ts "Configuration:"
    ts "  Comparing:    1 thread vs $NUM_THREADS threads"
    ts "  Algorithm:    LAPACK (same for both)"
    ts "  Parameters:   $ARG1 $ARG2 $ARG3"
    ts "  Output dir:   $OUTPUT_DIR"
    echo ""
    
    # Check binary
    if [ ! -f "./catcharea" ]; then
        ts "Error: ./catcharea not found!"
        ts "Please build: make clean && make all"
        exit 1
    fi
    
    # Check CPU info
    if grep -q avx2 /proc/cpuinfo 2>/dev/null; then
        ts "AVX2 support: YES ✅"
    else
        ts "AVX2 support: NO ⚠️"
    fi
    
    local total_cores=$(nproc)
    ts "Total CPU cores: $total_cores"
    echo ""
    
    # Create results CSV
    echo "Test|Threads|Total_Time|Inversion_Time|GFLOPS|Actual_Threads" > "$OUTPUT_DIR/results.csv"
    
    print_header "RUNNING BENCHMARKS"
    
    # =============================================================================
    # Test 1: Single-threaded LAPACK (baseline)
    # =============================================================================
    
    ts "Test 1 of 2: Single-threaded baseline"
    run_test "sequential_1thread" 1
    
    sleep 3  # Let system stabilize
    
    # =============================================================================
    # Test 2: Multi-threaded LAPACK
    # =============================================================================
    
    ts "Test 2 of 2: Multi-threaded parallel"
    run_test "parallel_${NUM_THREADS}threads" $NUM_THREADS
    
    # =============================================================================
    # Results Analysis
    # =============================================================================
    
    print_header "RESULTS ANALYSIS"
    
    if [ ! -f "$OUTPUT_DIR/results.csv" ]; then
        ts "Error: Results file not created"
        exit 1
    fi
    
    # Read results
    local baseline=$(grep "^sequential_1thread" "$OUTPUT_DIR/results.csv")
    local parallel=$(grep "^parallel_${NUM_THREADS}threads" "$OUTPUT_DIR/results.csv")
    
    if [ -z "$baseline" ] || [ -z "$parallel" ]; then
        ts "Error: Could not read results"
        cat "$OUTPUT_DIR/results.csv"
        exit 1
    fi
    
    # Extract data
    local base_inv_time=$(echo "$baseline" | cut -d'|' -f4)
    local par_inv_time=$(echo "$parallel" | cut -d'|' -f4)
    
    local base_gflops=$(echo "$baseline" | cut -d'|' -f5)
    local par_gflops=$(echo "$parallel" | cut -d'|' -f5)
    
    local base_threads=$(echo "$baseline" | cut -d'|' -f6)
    local par_threads=$(echo "$parallel" | cut -d'|' -f6)
    
    # Verify thread counts
    ts "Thread Verification:"
    ts "  Sequential requested: 1, actual: ${base_threads:-N/A}"
    ts "  Parallel requested:   $NUM_THREADS, actual: ${par_threads:-N/A}"
    echo ""
    
    if [ "$base_threads" != "1" ] 2>/dev/null; then
        ts "⚠️  WARNING: Sequential test did not use 1 thread!"
        ts "⚠️  Results may not be fair comparison"
        echo ""
    fi
    
    if [ "$par_threads" != "$NUM_THREADS" ] 2>/dev/null; then
        ts "⚠️  WARNING: Parallel test did not use $NUM_THREADS threads!"
        ts "⚠️  OpenBLAS may not be properly configured"
        echo ""
    fi
    
    # Calculate speedup
    if [ -n "$base_inv_time" ] && [ -n "$par_inv_time" ]; then
        local speedup=$(awk -v b="$base_inv_time" -v p="$par_inv_time" \
                        'BEGIN {printf "%.2f", b/p}')
        
        local efficiency=$(awk -v s="$speedup" -v t="$NUM_THREADS" \
                           'BEGIN {printf "%.1f", s/t*100}')
    else
        local speedup="N/A"
        local efficiency="N/A"
    fi
    
    # Display results table
    printf "\n%-30s %12s %12s %10s %10s\n" "Test" "Inv Time(s)" "GFLOPS" "Threads" "Speedup"
    printf "%-30s %12s %12s %10s %10s\n" "----" "-----------" "-------" "-------" "-------"
    printf "%-30s %12s %12s %10s %10s\n" \
        "1 thread (baseline)" "$base_inv_time" "$base_gflops" "${base_threads:-1}" "1.00×"
    printf "%-30s %12s %12s %10s %10s\n" \
        "$NUM_THREADS threads" "$par_inv_time" "$par_gflops" "${par_threads:-$NUM_THREADS}" "${speedup}×"
    echo ""
    
    # Summary
    ts "===== SUMMARY ====="
    ts ""
    ts "Speedup:    ${speedup}×"
    ts "Efficiency: ${efficiency}%"
    ts ""
    
    # Expected values
    local expected_min=$(awk -v t="$NUM_THREADS" 'BEGIN {printf "%.1f", t*0.70}')
    local expected_max=$(awk -v t="$NUM_THREADS" 'BEGIN {printf "%.1f", t*1.0}')
    
    ts "Expected range for $NUM_THREADS threads: ${expected_min}× - ${expected_max}×"
    ts ""
    
    # Interpretation
    ts "Interpretation:"
    
    if [ "$speedup" = "N/A" ]; then
        ts "  ❌ Could not calculate speedup (missing data)"
    elif (( $(echo "$speedup < 1.5" | bc -l) )); then
        ts "  ❌ PROBLEM: Speedup < 1.5× indicates:"
        ts "     - LAPACK not using multiple threads"
        ts "     - OpenBLAS misconfigured"
        ts "     - Check OPENBLAS_NUM_THREADS setting"
    elif (( $(echo "$speedup < $expected_min" | bc -l) )); then
        ts "  ⚠️  Speedup is LOWER than expected"
        ts "     Possible causes:"
        ts "     - Memory bandwidth saturation"
        ts "     - NUMA effects"
        ts "     - Cache contention"
    elif (( $(echo "$speedup > $expected_max" | bc -l) )); then
        ts "  ⚠️  Speedup is HIGHER than expected"
        if (( $(echo "$speedup > 15" | bc -l) )); then
            ts "  🚨 WARNING: Speedup > 15× is SUSPICIOUS!"
            ts "     This suggests:"
            ts "     - Unfair comparison (different algorithms)"
            ts "     - Measurement error"
            ts "     - Caching effects"
            ts "     - Invalid timing"
        else
            ts "     This is unusual but possible if:"
            ts "     - Super-linear speedup from cache effects"
            ts "     - NUMA optimization"
        fi
    else
        ts "  ✅ Speedup is REASONABLE for $NUM_THREADS threads"
        ts "  ✅ LAPACK is properly using multiple threads"
        ts "  ✅ Efficiency is ${efficiency}%"
        
        if (( $(echo "$efficiency > 80" | bc -l) )); then
            ts "  ✅ Excellent scaling!"
        elif (( $(echo "$efficiency > 60" | bc -l) )); then
            ts "  ✅ Good scaling"
        else
            ts "  ⚠️  Moderate scaling (typical for dense linear algebra)"
        fi
    fi
    echo ""
    
    # Recommendations
    ts "Recommendations:"
    
    if [ "$speedup" != "N/A" ] && (( $(echo "$speedup < $expected_min" | bc -l) )); then
        ts "  1. Verify OpenBLAS configuration:"
        ts "     export OPENBLAS_VERBOSE=1"
        ts "     ./catcharea ... (check output for thread usage)"
        ts ""
        ts "  2. Check NUMA settings:"
        ts "     numactl --hardware"
        ts "     numactl --cpubind=0 ./catcharea ..."
        ts ""
        ts "  3. Profile with perf:"
        ts "     perf stat -e cache-misses,cycles ./catcharea ..."
    elif [ "$speedup" != "N/A" ] && (( $(echo "$speedup > 15" | bc -l) )); then
        ts "  1. ⚠️  RE-VERIFY MEASUREMENT:"
        ts "     - Check both tests use same algorithm"
        ts "     - Verify timing excludes I/O"
        ts "     - Run tests multiple times"
        ts ""
        ts "  2. Check for caching:"
        ts "     - Clear caches between runs"
        ts "     - Randomize input data"
    else
        ts "  ✅ Results look good!"
        ts "  1. Run multiple times to confirm consistency"
        ts "  2. Test with different problem sizes"
        ts "  3. Document results in paper"
    fi
    echo ""
    
    print_header "BENCHMARK COMPLETE"
    
    ts "Results saved to: $OUTPUT_DIR/"
    ts "  - benchmark.log:              Full log"
    ts "  - results.csv:                Summary table"
    ts "  - sequential_1thread.txt:     1-thread detailed output"
    ts "  - parallel_${NUM_THREADS}threads.txt:    $NUM_THREADS-thread detailed output"
    echo ""
    
    ts "Next steps:"
    ts "  1. Review detailed outputs for verification"
    ts "  2. Run with perf for hardware counters"
    ts "  3. Update Section III with measured speedup"
    echo ""
}

# =============================================================================
# Run main
# =============================================================================

main

exit 0
