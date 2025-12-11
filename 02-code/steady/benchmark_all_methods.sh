#!/usr/bin/env bash
# Benchmark script - Test all 4 matrix multiplication methods
# Compares performance and validates correctness
# Usage: ./benchmark_all_methods.sh [NUM_THREADS] [ARG1 ARG2 ARG3]
#
# ═══════════════════════════════════════════════════════════════════════════
# CORRECTED VERSION - Uses CSV data to avoid double-counting issue
# ═══════════════════════════════════════════════════════════════════════════

set -u

# Configuration
NUM_THREADS="${1:-6}"
ARG1="${2:-1.0}"
ARG2="${3:-100.0}"
ARG3="${4:-0.001}"
INVERSION_METHOD=0 

# Output directory
OUTPUT_DIR="benchmark_results_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTPUT_DIR"

# Logging
LOG_FILE="$OUTPUT_DIR/benchmark.log"
exec > >(tee -a "$LOG_FILE") 2>&1

# ===============================================================================
# Helper functions
# ===============================================================================

ts() { printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*"; }

print_header() {
    echo ""
    echo "═══════════════════════════════════════════════════════════════════════════"
    echo "  $1"
    echo "═══════════════════════════════════════════════════════════════════════════"
    echo ""
}

# ===============================================================================
# Setup environment
# ===============================================================================

setup_environment() {
    export OMP_NUM_THREADS="$NUM_THREADS"
    export OMP_PROC_BIND=close
    export OMP_PLACES=cores
    export OPENBLAS_NUM_THREADS="$NUM_THREADS"
    export MALLOC_ARENA_MAX=4
}

# ===============================================================================
# Run single test
# ===============================================================================

run_test() {
    local method=$1
    local block_size=$2
    local method_name=$3
    
    local current_inv_method=0
    if [ "$method" -eq 0 ]; then
        current_inv_method=1
    fi
    
    local output_file="$OUTPUT_DIR/method${method}.txt"
    local timing_file="$OUTPUT_DIR/timing_method${method}.txt"
    
    ts "Testing: $method_name"
    ts "  Command: ./catcharea $ARG1 $ARG2 $ARG3 $current_inv_method $method $block_size"

    { time ./catcharea $ARG1 $ARG2 $ARG3 $current_inv_method $method $block_size > "$output_file" 2>&1; } 2> "$timing_file"
    
    local exit_code=$?
    
    if [ $exit_code -eq 0 ]; then
        # Extract catchment area
        local area=$(awk '/Catchment area/{val=$NF} END{if(val!="") printf "%s", val}' "$output_file")
        
        # ═══════════════════════════════════════════════════════════════════════════
        # PRIORITY 1: Parse from CSV embedded in output (MOST ACCURATE!)
        # ═══════════════════════════════════════════════════════════════════════════
        
        # Extract CSV section from output file
        local csv_section=$(awk '/^Parameter,Value$/,/^Performance data exported/' "$output_file")
        
        local dgemm_time=$(echo "$csv_section" | awk -F, '/^Multiply_Time_sec,/{print $2}')
        local inversion_time=$(echo "$csv_section" | awk -F, '/^Inversion_Time_sec,/{print $2}')
        local total_comp_time=$(echo "$csv_section" | awk -F, '/^Computation_Time_sec,/{print $2}')
        local csv_total_time=$(echo "$csv_section" | awk -F, '/^Total_Time_sec,/{print $2}')
        
        # ═══════════════════════════════════════════════════════════════════════════
        # FALLBACK: If CSV parsing failed, try text (LAST occurrence to avoid duplicates)
        # ═══════════════════════════════════════════════════════════════════════════
        
        local used_fallback=false
        
        if [ -z "$dgemm_time" ] || [ -z "$inversion_time" ] || [ -z "$total_comp_time" ]; then
          ts "  ⚠️  CSV parsing failed, trying text fallback"
          used_fallback=true
          
          # Use tail to get LAST occurrence (avoid double-counting from first section)
          dgemm_time=$(tail -200 "$output_file" | awk '
            /Total DGEMM.*time:.*seconds/ {
              for(i=1; i<=NF; i++) {
                if($i ~ /^[0-9]+\.?[0-9]*$/ && $(i+1) == "seconds") {
                  value=$i
                }
              }
            }
            END {if(value) print value}
          ')
          
          inversion_time=$(tail -200 "$output_file" | awk '
            /Total Matrix Inversion time:.*seconds/ {
              for(i=1; i<=NF; i++) {
                if($i ~ /^[0-9]+\.?[0-9]*$/ && $(i+1) == "seconds") {
                  value=$i
                }
              }
            }
            END {if(value) print value}
          ')
          
          total_comp_time=$(tail -200 "$output_file" | awk '
            /Total computation time:.*seconds/ {
              for(i=1; i<=NF; i++) {
                if($i ~ /^[0-9]+\.?[0-9]*$/ && $(i+1) == "seconds") {
                  value=$i
                }
              }
            }
            END {if(value) print value}
          ')
        fi
        
        # ═══════════════════════════════════════════════════════════════════════════
        # Calculate total if we have parts but not sum
        # ═══════════════════════════════════════════════════════════════════════════
        
        if [ -n "$dgemm_time" ] && [ -n "$inversion_time" ] && [ -z "$total_comp_time" ]; then
          total_comp_time=$(awk -v d="$dgemm_time" -v i="$inversion_time" 'BEGIN {printf "%.6f", d+i}')
        fi
        
        # ═══════════════════════════════════════════════════════════════════════════
        # Extract wall-clock time
        # ═══════════════════════════════════════════════════════════════════════════
        
        local real_time=$(grep "^real" "$timing_file" | awk '{print $2}')
        
        # Convert m:s format to just seconds
        local real_time_sec=$(echo "$real_time" | awk -F'[ms]' '{
          if (NF==3) print $1*60 + $2;
          else if ($0 ~ /m/) print $1*60 + $2;
          else print $1;
        }')
        
        # ═══════════════════════════════════════════════════════════════════════════
        # Display results
        # ═══════════════════════════════════════════════════════════════════════════
        
        ts "  Status: SUCCESS ✅"
        ts "  Catchment area: $area"
        
        if [ -n "$dgemm_time" ] && [ -n "$inversion_time" ] && [ -n "$total_comp_time" ]; then
          ts "  Matrix Multiplication: ${dgemm_time}s"
          ts "  Matrix Inversion: ${inversion_time}s"
          ts "  Matrix Ops Total: ${total_comp_time}s"
          
          if $used_fallback; then
            ts "  (⚠️  Used text fallback - verify accuracy)"
          fi
          
          # Verify that total <= wall-clock
          local is_valid=$(awk -v t="$total_comp_time" -v w="$real_time_sec" 'BEGIN {
            print (t <= w * 1.01) ? "yes" : "no"
          }')
          
          if [ "$is_valid" = "no" ]; then
            ts "  ⚠️  WARNING: Matrix Ops (${total_comp_time}s) > Wall-Clock (${real_time_sec}s)!"
            ts "             This indicates timing measurement error (double-counting)"
          fi
        else
          ts "  ⚠️  Could not extract all timing values"
        fi
        
        ts "  Total wall-clock time: ${real_time_sec}s"
        
        # Copy performance_results.csv if it exists
        if [ -f "performance_results.csv" ]; then
            cp -f "performance_results.csv" "$OUTPUT_DIR/method${method}_perf.csv" 2>/dev/null || true
        fi
        
        # ═══════════════════════════════════════════════════════════════════════════
        # Save to CSV
        # ═══════════════════════════════════════════════════════════════════════════
        
        local dgemm_csv="${dgemm_time:-N/A}"
        local inv_csv="${inversion_time:-N/A}"
        local total_csv="${total_comp_time:-N/A}"
        
        echo "$method|$method_name|$area|$dgemm_csv|$inv_csv|$total_csv|$real_time_sec" >> "$OUTPUT_DIR/results.csv"
        
        return 0
    else
        ts "  Status: FAILED ❌"
        return 1
    fi
}

# ===============================================================================
# Main benchmark routine
# ===============================================================================

main() {
    print_header "CATCHAREA MATRIX OPTIMIZATION BENCHMARK SUITE"
    
    ts "Configuration:"
    ts "  Threads: $NUM_THREADS"
    ts "  Parameters: $ARG1 $ARG2 $ARG3"
    ts "  Output directory: $OUTPUT_DIR"
    echo ""
    
    if [ ! -f "./catcharea" ]; then
        ts "Error: ./catcharea not found!"
        ts "Please build the project first: make clean && make"
        exit 1
    fi
    
    if grep -q avx2 /proc/cpuinfo 2>/dev/null; then
        ts "AVX2 support: YES ✅"
    else
        ts "AVX2 support: NO ⚠️  (method 3 will fall back to method 2)"
    fi
    echo ""
    
    setup_environment
    
    # Create results CSV header
    echo "Method|Name|Catchment_Area|DGEMM_Time|Inversion_Time|Total_Matrix_Ops|Wall_Clock" > "$OUTPUT_DIR/results.csv"
    
    print_header "RUNNING BENCHMARKS"
    
    declare -A methods=(
        [0]="Sequential (baseline)"
        [1]="OpenMP"
        [2]="OpenMP + Cache Blocking"
        [3]="OpenMP + Cache + SIMD"
    )
    
    declare -A block_sizes=(
        [0]=32
        [1]=32
        [2]=32
        [3]=32
    )
    
    local all_success=true
    
    for method in 0 1 2 3; do
        echo ""
        if ! run_test $method ${block_sizes[$method]} "${methods[$method]}"; then
            all_success=false
        fi
        echo ""
        sleep 2
    done
    
    print_header "BENCHMARK RESULTS"
    
    if [ ! -f "$OUTPUT_DIR/results.csv" ]; then
        ts "Error: Results file not created"
        exit 1
    fi
    
    ts "Reading results..."
    echo ""
    
    # Read baseline
    local baseline_line=$(grep "^0|" "$OUTPUT_DIR/results.csv")
    local baseline_wall=$(echo "$baseline_line" | cut -d'|' -f7)
    local baseline_dgemm=$(echo "$baseline_line" | cut -d'|' -f4)
    local baseline_inv=$(echo "$baseline_line" | cut -d'|' -f5)
    local baseline_total=$(echo "$baseline_line" | cut -d'|' -f6)
    local baseline_area=$(echo "$baseline_line" | cut -d'|' -f3)
    
    if [ -z "$baseline_wall" ]; then
        ts "Warning: Could not extract baseline time"
        baseline_wall=1.0
    fi
    
    ts "Baseline (Method 0 - Sequential):"
    ts "  Wall-clock time: ${baseline_wall}s"
    if [ "$baseline_dgemm" != "N/A" ]; then
      ts "  Matrix Multiply: ${baseline_dgemm}s"
    fi
    if [ "$baseline_inv" != "N/A" ]; then
      ts "  Matrix Inversion: ${baseline_inv}s"
    fi
    if [ "$baseline_total" != "N/A" ]; then
      ts "  Matrix Ops Total: ${baseline_total}s"
    fi
    ts "  Catchment area: $baseline_area"
    echo ""
    
    # ═══════════════════════════════════════════════════════════════════════════
    # Display results table
    # ═══════════════════════════════════════════════════════════════════════════
    
    printf "%-35s %15s %18s %18s %18s %12s %15s\n" \
        "Method" "Wall-Clock(s)" "Matrix Mult(s)" "Matrix Inv(s)" "Matrix Ops(s)" "Speedup" "Area"
    printf "%-35s %15s %18s %18s %18s %12s %15s\n" \
        "------" "-------------" "---------------" "--------------" "--------------" "-------" "----"
    
    while IFS='|' read -r method name area dgemm_time inv_time total_comp wall_time; do
        if [ "$method" = "Method" ]; then continue; fi
        
        # Calculate speedup
        if [ -n "$wall_time" ] && [ -n "$baseline_wall" ] && [ "$wall_time" != "N/A" ]; then
            speedup=$(awk -v b="$baseline_wall" -v t="$wall_time" 'BEGIN {
              if (t > 0) printf "%.2f", b/t
              else print "N/A"
            }')
        else
            speedup="N/A"
        fi
        
        # Format times
        if [ "$dgemm_time" != "N/A" ]; then
          dgemm_fmt=$(printf "%.2f" "$dgemm_time" 2>/dev/null || echo "$dgemm_time")
        else
          dgemm_fmt="N/A"
        fi
        
        if [ "$inv_time" != "N/A" ]; then
          inv_fmt=$(printf "%.2f" "$inv_time" 2>/dev/null || echo "$inv_time")
        else
          inv_fmt="N/A"
        fi
        
        if [ "$total_comp" != "N/A" ]; then
          total_fmt=$(printf "%.2f" "$total_comp" 2>/dev/null || echo "$total_comp")
        else
          total_fmt="N/A"
        fi
        
        wall_fmt=$(printf "%.2f" "$wall_time" 2>/dev/null || echo "$wall_time")
        
        printf "%-35s %15s %18s %18s %18s %11sx %15s\n" \
            "$name" \
            "$wall_fmt" \
            "$dgemm_fmt" \
            "$inv_fmt" \
            "$total_fmt" \
            "$speedup" \
            "$area"
    done < "$OUTPUT_DIR/results.csv"
    
    echo ""
    
    # Verify correctness
    ts "Correctness Verification:"
    
    local areas=$(cut -d'|' -f3 "$OUTPUT_DIR/results.csv" | tail -n +2)
    local first_area=$(echo "$areas" | head -1)
    local all_match=true
    
    while read -r area; do
        local diff=$(awk -v a="$area" -v b="$first_area" 'BEGIN {printf "%.10f", (a-b)>0?(a-b):(b-a)}')
        local is_close=$(awk -v d="$diff" 'BEGIN {print (d < 0.001) ? "yes" : "no"}')
        
        if [ "$is_close" = "no" ]; then
            all_match=false
            ts "  ⚠️  Warning: Area mismatch detected: $area vs $first_area (diff: $diff)"
        fi
    done <<< "$areas"
    
    if $all_match; then
        ts "  ✅ All methods produce identical results (within tolerance)"
    else
        ts "  ❌ Some methods produce different results - check output files"
    fi
    echo ""
    
    # Performance recommendations
    print_header "RECOMMENDATIONS"
    
    local best_method=0
    local best_time=999999
    
    for method in 0 1 2 3; do
        local line=$(grep "^${method}|" "$OUTPUT_DIR/results.csv")
        local time=$(echo "$line" | cut -d'|' -f7)
        
        if [ -n "$time" ] && [ "$time" != "N/A" ]; then
            local is_better=$(awk -v t="$time" -v b="$best_time" 'BEGIN {print (t < b) ? "yes" : "no"}')
            if [ "$is_better" = "yes" ]; then
                best_method=$method
                best_time=$time
            fi
        fi
    done
    
    ts "Best performing method: Method $best_method (${methods[$best_method]})"
    ts "  Wall-clock time: ${best_time}s"
    ts "  Speedup: $(awk -v b="$baseline_wall" -v t="$best_time" 'BEGIN {printf "%.2fx", b/t}') vs baseline"
    echo ""
    
    if [ $best_method -eq 3 ]; then
        ts "✅ You're already using the fastest method!"
    else
        ts "💡 Try method 3 for better performance:"
        ts "   ./catcharea $ARG1 $ARG2 $ARG3 0 3 32"
    fi
    echo ""
    
    if [ $best_method -ge 2 ]; then
        ts "💡 Fine-tune performance by testing different block sizes:"
        ts "   for bs in 16 24 32 48 64; do"
        ts "     ./catcharea $ARG1 $ARG2 $ARG3 0 $best_method \$bs"
        ts "   done"
        echo ""
    fi
    
    print_header "BENCHMARK COMPLETE"
    
    ts "Results saved to: $OUTPUT_DIR/"
    ts "  - benchmark.log: Full log"
    ts "  - results.csv: Summary table"
    ts "  - method*.txt: Detailed output for each method"
    ts "  - timing_method*.txt: Timing information"
    ts "  - method*_perf.csv: Performance CSV files (if available)"
    echo ""
    
    if $all_success; then
        ts "All tests completed successfully! ✅"
    else
        ts "Some tests failed - check log files ⚠️"
    fi
    echo ""
    
    # Print note
    echo ""
    ts "📊 NOTE: Matrix Operations breakdown:"
    ts "   • Matrix Mult(s):  DGEMM (dense matrix-matrix multiplication)"
    ts "   • Matrix Inv(s):   LU factorization + solve"
    ts "   • Matrix Ops(s):   Sum of both operations"
    ts "   Data source: CSV embedded in output (most accurate!)"
    echo ""
    
    # Check for timing issues
    local any_warnings=$(grep "WARNING: Matrix Ops" "$LOG_FILE" 2>/dev/null)
    if [ -n "$any_warnings" ]; then
        echo ""
        ts "⚠️  TIMING WARNINGS DETECTED:"
        ts "   Some methods show Matrix Ops > Wall-Clock time"
        ts "   This indicates double-counting in timing measurement"
        ts "   CSV values are correct; text parsing may be wrong"
        echo ""
    fi
}

# ===============================================================================
# Run main
# ===============================================================================

main

exit 0