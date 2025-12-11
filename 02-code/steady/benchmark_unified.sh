#!/usr/bin/env bash
################################################################################
# benchmark_unified.sh - Unified Benchmark Script for Section IV
#
# VERSION 3.0: Supports both Hybrid and OpenBLAS testing in single script
#
# Features:
# - Tier 1: Synthetic benchmarks (controlled experiments)
# - Tier 2: Real watershed validation (Dontako)
# - Hybrid vs OpenBLAS comparison
# - All optimization levels (0-3)
# - Thread scaling analysis
# - CSV output for paper tables
#
# Usage:
#   ./benchmark_unified.sh [MODE] [THREADS]
#
# Modes:
#   all      - Run all experiments (default)
#   quick    - Quick test (methods 0 and 3 only)
#   hybrid   - Test Hybrid implementation only
#   openblas - Test OpenBLAS implementation only
#   compare  - Side-by-side Hybrid vs OpenBLAS
#   scaling  - Thread scaling analysis
#
# Author: BEM Flow Path Analysis - Parallel Optimization Research
# Date: 2025-12-05
################################################################################

set -euo pipefail

#═══════════════════════════════════════════════════════════════════════════════
# Configuration
#═══════════════════════════════════════════════════════════════════════════════

MODE="${1:-all}"
NUM_THREADS="${2:-6}"

# Watershed parameters (Dontako default)
ARG_STEP="${3:-1.0}"
ARG_RM="${4:-100.0}"
ARG_DR="${5:-0.001}"

# Output directory
TIMESTAMP=$(date +%Y%m%d_%H%M%S)
OUTPUT_DIR="benchmark_unified_${TIMESTAMP}"
mkdir -p "$OUTPUT_DIR"

# Logging
LOG_FILE="$OUTPUT_DIR/benchmark.log"
exec > >(tee -a "$LOG_FILE") 2>&1

#═══════════════════════════════════════════════════════════════════════════════
# Helper Functions
#═══════════════════════════════════════════════════════════════════════════════

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

#═══════════════════════════════════════════════════════════════════════════════
# Environment Setup
#═══════════════════════════════════════════════════════════════════════════════

setup_environment() {
    export OMP_NUM_THREADS="$NUM_THREADS"
    export OPENBLAS_NUM_THREADS="$NUM_THREADS"
    export OMP_PROC_BIND=close
    export OMP_PLACES=cores
    export MALLOC_ARENA_MAX=4
    
    ts "Environment configured:"
    ts "  OMP_NUM_THREADS:      $OMP_NUM_THREADS"
    ts "  OPENBLAS_NUM_THREADS: $OPENBLAS_NUM_THREADS"
    ts "  OMP_PROC_BIND:        $OMP_PROC_BIND"
    ts "  OMP_PLACES:           $OMP_PLACES"
}

#═══════════════════════════════════════════════════════════════════════════════
# Core Test Function
#═══════════════════════════════════════════════════════════════════════════════

run_test() {
    local dgemm_type=$1      # 0=Hybrid, 1=OpenBLAS
    local mult_method=$2     # 0=Seq, 1=OMP, 2=Cache, 3=SIMD
    local inv_method=$3      # 0=Parallel, 1=Sequential
    local block_size=${4:-64}
    local test_name=$5
    
    local output_file="$OUTPUT_DIR/${test_name}.txt"
    local timing_file="$OUTPUT_DIR/${test_name}_timing.txt"
    
    ts "Running: $test_name"
    ts "  Config: dgemm=$dgemm_type, method=$mult_method, inv=$inv_method, block=$block_size"
    ts "  Command: ./catcharea $ARG_STEP $ARG_RM $ARG_DR $inv_method $mult_method $block_size $dgemm_type"
    
    # Run with timing
    local start_time=$(date +%s.%N)
    
    { time ./catcharea $ARG_STEP $ARG_RM $ARG_DR $inv_method $mult_method $block_size $dgemm_type \
        > "$output_file" 2>&1; } 2> "$timing_file"
    
    local exit_code=$?
    local end_time=$(date +%s.%N)
    
    if [ $exit_code -ne 0 ]; then
        ts "  ❌ FAILED (exit code: $exit_code)"
        return 1
    fi
    
    # Parse results from CSV section
    local csv_section=$(awk '/^Parameter,Value$/,/^Performance data exported/' "$output_file")
    
    local area=$(echo "$csv_section" | awk -F, '/^Catchment_Area,/{print $2}')
    local dgemm_time=$(echo "$csv_section" | awk -F, '/^Multiply_Time_sec,/{print $2}')
    local inv_time=$(echo "$csv_section" | awk -F, '/^Inversion_Time_sec,/{print $2}')
    local comp_time=$(echo "$csv_section" | awk -F, '/^Computation_Time_sec,/{print $2}')
    local total_time=$(echo "$csv_section" | awk -F, '/^Total_Time_sec,/{print $2}')
    
    # Fallback for area if CSV parsing fails
    if [ -z "$area" ] || [ "$area" = "N/A" ]; then
        area=$(grep -i "Catchment area" "$output_file" | tail -1 | \
               awk '{for(i=1;i<=NF;i++) if($i ~ /^[0-9]+\.?[0-9]*$/) print $i}' | tail -1)
    fi
    
    # Get wall-clock time
    local real_time=$(grep "^real" "$timing_file" | awk '{print $2}')
    local real_time_sec=$(echo "$real_time" | awk -F'[ms]' '{
        if (NF==3) print $1*60 + $2;
        else if ($0 ~ /m/) print $1*60 + $2;
        else print $1;
    }')
    
    ts "  ✅ SUCCESS"
    ts "     Area:        ${area:-N/A}"
    ts "     DGEMM:       ${dgemm_time:-N/A}s"
    ts "     Inversion:   ${inv_time:-N/A}s"
    ts "     Computation: ${comp_time:-N/A}s"
    ts "     Wall-clock:  ${real_time_sec}s"
    
    # Copy performance CSV
    if [ -f "performance_results.csv" ]; then
        cp -f "performance_results.csv" "$OUTPUT_DIR/${test_name}_perf.csv" 2>/dev/null || true
    fi
    
    # Append to results CSV
    echo "$test_name|$dgemm_type|$mult_method|$inv_method|$block_size|${area:-N/A}|${dgemm_time:-N/A}|${inv_time:-N/A}|${comp_time:-N/A}|${real_time_sec}" \
        >> "$OUTPUT_DIR/results.csv"
    
    return 0
}

#═══════════════════════════════════════════════════════════════════════════════
# Test Modes
#═══════════════════════════════════════════════════════════════════════════════

#───────────────────────────────────────────────────────────────────────────────
# Mode: All Experiments (Comprehensive)
#───────────────────────────────────────────────────────────────────────────────
run_all_experiments() {
    print_header "COMPREHENSIVE BENCHMARK SUITE"
    
    # Initialize results CSV
    echo "Test|DGEMM_Type|Mult_Method|Inv_Method|Block_Size|Area|DGEMM_Time|Inv_Time|Comp_Time|Wall_Time" \
        > "$OUTPUT_DIR/results.csv"
    
    #───────────────────────────────────────────────────────────────────────────
    # Part 1: Sequential Baseline (True baseline for both)
    #───────────────────────────────────────────────────────────────────────────
    print_subheader "Part 1: Sequential Baseline"
    
    # Method 0 with sequential inversion - same for both Hybrid and OpenBLAS
    run_test 0 0 1 32 "baseline_sequential"
    
    #───────────────────────────────────────────────────────────────────────────
    # Part 2: Hybrid Implementation (All optimization levels)
    #───────────────────────────────────────────────────────────────────────────
    print_subheader "Part 2: Hybrid Implementation (dgemm_type=0)"
    
    run_test 0 1 0 64 "hybrid_method1_openmp"
    run_test 0 2 0 32 "hybrid_method2_cache_b32"
    run_test 0 2 0 64 "hybrid_method2_cache_b64"
    run_test 0 3 0 32 "hybrid_method3_simd_b32"
    run_test 0 3 0 64 "hybrid_method3_simd_b64"
    
    #───────────────────────────────────────────────────────────────────────────
    # Part 3: OpenBLAS Implementation
    #───────────────────────────────────────────────────────────────────────────
    print_subheader "Part 3: OpenBLAS Implementation (dgemm_type=1)"
    
    run_test 1 1 0 64 "openblas_method1"
    run_test 1 2 0 64 "openblas_method2"
    run_test 1 3 0 64 "openblas_method3"
    
    #───────────────────────────────────────────────────────────────────────────
    # Part 4: Block Size Tuning (Hybrid only)
    #───────────────────────────────────────────────────────────────────────────
    print_subheader "Part 4: Block Size Tuning (Hybrid Method 3)"
    
    for BS in 16 24 32 48 64 96 128; do
        run_test 0 3 0 $BS "hybrid_block_${BS}"
    done
}

#───────────────────────────────────────────────────────────────────────────────
# Mode: Quick Test
#───────────────────────────────────────────────────────────────────────────────
run_quick_test() {
    print_header "QUICK BENCHMARK (Baseline + Best Methods)"
    
    echo "Test|DGEMM_Type|Mult_Method|Inv_Method|Block_Size|Area|DGEMM_Time|Inv_Time|Comp_Time|Wall_Time" \
        > "$OUTPUT_DIR/results.csv"
    
    # Sequential baseline
    run_test 0 0 1 32 "baseline_sequential"
    
    # Best Hybrid (Method 3)
    run_test 0 3 0 64 "hybrid_best"
    
    # Best OpenBLAS
    run_test 1 3 0 64 "openblas_best"
}

#───────────────────────────────────────────────────────────────────────────────
# Mode: Hybrid vs OpenBLAS Comparison
#───────────────────────────────────────────────────────────────────────────────
run_comparison() {
    print_header "HYBRID vs OPENBLAS COMPARISON"
    
    echo "Test|DGEMM_Type|Mult_Method|Inv_Method|Block_Size|Area|DGEMM_Time|Inv_Time|Comp_Time|Wall_Time" \
        > "$OUTPUT_DIR/results.csv"
    
    # Baseline
    run_test 0 0 1 32 "baseline_sequential"
    
    print_subheader "Method 1: OpenMP Only"
    run_test 0 1 0 64 "hybrid_method1"
    run_test 1 1 0 64 "openblas_method1"
    
    print_subheader "Method 2: OpenMP + Cache"
    run_test 0 2 0 64 "hybrid_method2"
    run_test 1 2 0 64 "openblas_method2"
    
    print_subheader "Method 3: Full Optimization"
    run_test 0 3 0 64 "hybrid_method3"
    run_test 1 3 0 64 "openblas_method3"
}

#───────────────────────────────────────────────────────────────────────────────
# Mode: Thread Scaling Analysis
#───────────────────────────────────────────────────────────────────────────────
run_scaling() {
    print_header "THREAD SCALING ANALYSIS"
    
    echo "Test|DGEMM_Type|Mult_Method|Inv_Method|Block_Size|Threads|Area|DGEMM_Time|Inv_Time|Comp_Time|Wall_Time" \
        > "$OUTPUT_DIR/scaling_results.csv"
    
    # Save original thread count
    local orig_threads=$NUM_THREADS
    
    for T in 1 2 3 4 5 6; do
        print_subheader "Threads: $T"
        
        export OMP_NUM_THREADS=$T
        export OPENBLAS_NUM_THREADS=$T
        
        # Hybrid
        run_test 0 3 0 64 "hybrid_t${T}"
        
        # OpenBLAS
        run_test 1 3 0 64 "openblas_t${T}"
    done
    
    # Restore original
    export OMP_NUM_THREADS=$orig_threads
    export OPENBLAS_NUM_THREADS=$orig_threads
}

#═══════════════════════════════════════════════════════════════════════════════
# Results Analysis
#═══════════════════════════════════════════════════════════════════════════════

analyze_results() {
    print_header "RESULTS ANALYSIS"
    
    if [ ! -f "$OUTPUT_DIR/results.csv" ]; then
        ts "No results file found"
        return
    fi
    
    # Read baseline
    local baseline_line=$(grep "baseline_sequential" "$OUTPUT_DIR/results.csv" | head -1)
    local baseline_wall=$(echo "$baseline_line" | cut -d'|' -f10)
    local baseline_area=$(echo "$baseline_line" | cut -d'|' -f6)
    
    ts "Baseline (Sequential):"
    ts "  Wall-clock time: ${baseline_wall}s"
    ts "  Catchment area:  ${baseline_area}"
    echo ""
    
    # Results table
    printf "%-30s %8s %8s %8s %10s %10s\n" \
        "Test" "DGEMM" "Method" "Wall(s)" "Speedup" "Area"
    printf "%-30s %8s %8s %8s %10s %10s\n" \
        "----" "----" "------" "-------" "-------" "----"
    
    tail -n +2 "$OUTPUT_DIR/results.csv" | while IFS='|' read -r name dgemm method inv block area dgemm_t inv_t comp_t wall_t; do
        local speedup="N/A"
        if [ -n "$wall_t" ] && [ -n "$baseline_wall" ] && \
           awk -v w="$wall_t" 'BEGIN {exit !(w > 0)}' 2>/dev/null; then
            speedup=$(awk -v b="$baseline_wall" -v w="$wall_t" 'BEGIN {printf "%.2fx", b/w}')
        fi
        
        local dgemm_name
        [ "$dgemm" = "0" ] && dgemm_name="Hybrid" || dgemm_name="OpenBLAS"
        
        local area_short
        if [ -n "$area" ] && [ "$area" != "N/A" ]; then
            area_short=$(printf "%.4f" "$area" 2>/dev/null || echo "$area")
        else
            area_short="N/A"
        fi
        
        printf "%-30s %8s %8s %8.2f %10s %10s\n" \
            "$name" "$dgemm_name" "$method" "${wall_t:-0}" "$speedup" "$area_short"
    done
    
    echo ""
    
    # Correctness check
    ts "Correctness Verification:"
    local areas=$(tail -n +2 "$OUTPUT_DIR/results.csv" | cut -d'|' -f6 | grep -v "N/A")
    local first_area=$(echo "$areas" | head -1)
    local all_match=true
    
    while read -r area; do
        if [ -z "$area" ] || [ "$area" = "N/A" ]; then continue; fi
        local diff=$(awk -v a="$area" -v b="$first_area" 'BEGIN {d=a-b; printf "%.10f", (d<0)?-d:d}')
        local is_close=$(awk -v d="$diff" 'BEGIN {print (d < 0.001) ? "yes" : "no"}')
        if [ "$is_close" = "no" ]; then
            all_match=false
            ts "  ⚠️  Area mismatch: $area vs $first_area (diff: $diff)"
        fi
    done <<< "$areas"
    
    if $all_match; then
        ts "  ✅ All methods produce identical results (within tolerance)"
    fi
}

#═══════════════════════════════════════════════════════════════════════════════
# Generate LaTeX Table (for Section IV)
#═══════════════════════════════════════════════════════════════════════════════

generate_latex_table() {
    print_subheader "Generating LaTeX Table"
    
    local tex_file="$OUTPUT_DIR/table_results.tex"
    
    cat > "$tex_file" << 'EOF'
\begin{table}[htbp]
\centering
\caption{Performance Comparison: Hybrid vs OpenBLAS Implementation}
\label{tab:dgemm_comparison}
\begin{tabular}{lcccccc}
\hline
\textbf{Configuration} & \textbf{DGEMM} & \textbf{Method} & \textbf{DGEMM (s)} & \textbf{Inv (s)} & \textbf{Total (s)} & \textbf{Speedup} \\
\hline
EOF
    
    # Read baseline for speedup calculation
    local baseline_wall=$(grep "baseline_sequential" "$OUTPUT_DIR/results.csv" | cut -d'|' -f10)
    
    tail -n +2 "$OUTPUT_DIR/results.csv" | while IFS='|' read -r name dgemm method inv block area dgemm_t inv_t comp_t wall_t; do
        local speedup="1.00"
        if [ -n "$wall_t" ] && [ -n "$baseline_wall" ]; then
            speedup=$(awk -v b="$baseline_wall" -v w="$wall_t" 'BEGIN {printf "%.2f", b/w}')
        fi
        
        local dgemm_name
        [ "$dgemm" = "0" ] && dgemm_name="Hybrid" || dgemm_name="OpenBLAS"
        
        # Clean name for LaTeX
        local clean_name=$(echo "$name" | sed 's/_/\\_/g')
        
        printf "%s & %s & %s & %.3f & %.3f & %.2f & %.2f\\\\times \\\\\\\\\n" \
            "$clean_name" "$dgemm_name" "$method" \
            "${dgemm_t:-0}" "${inv_t:-0}" "${wall_t:-0}" "$speedup" >> "$tex_file"
    done
    
    cat >> "$tex_file" << 'EOF'
\hline
\end{tabular}
\end{table}
EOF
    
    ts "LaTeX table saved: $tex_file"
}

#═══════════════════════════════════════════════════════════════════════════════
# Main
#═══════════════════════════════════════════════════════════════════════════════

main() {
    print_header "BEM UNIFIED BENCHMARK SUITE v3.0"
    
    ts "Configuration:"
    ts "  Mode:       $MODE"
    ts "  Threads:    $NUM_THREADS"
    ts "  Parameters: step=$ARG_STEP, rm=$ARG_RM, dr=$ARG_DR"
    ts "  Output:     $OUTPUT_DIR"
    echo ""
    
    # Check executable
    if [ ! -f "./catcharea" ]; then
        ts "❌ Error: ./catcharea not found!"
        ts "   Please build first: make clean && make"
        exit 1
    fi
    
    # Check help
    if [ "$MODE" = "-h" ] || [ "$MODE" = "--help" ]; then
        echo ""
        echo "Usage: $0 [MODE] [THREADS] [STEP] [RM] [DR]"
        echo ""
        echo "Modes:"
        echo "  all      - Run all experiments (default)"
        echo "  quick    - Quick test (baseline + best methods)"
        echo "  hybrid   - Test Hybrid implementation only"
        echo "  openblas - Test OpenBLAS implementation only"
        echo "  compare  - Side-by-side Hybrid vs OpenBLAS"
        echo "  scaling  - Thread scaling analysis"
        echo ""
        echo "Examples:"
        echo "  $0 quick 6                # Quick test with 6 threads"
        echo "  $0 compare 6              # Compare Hybrid vs OpenBLAS"
        echo "  $0 scaling 6              # Thread scaling (1-6)"
        echo "  $0 all 6 1.0 100.0 0.001  # Full test with custom params"
        echo ""
        exit 0
    fi
    
    # Setup
    setup_environment
    
    # Run selected mode
    case "$MODE" in
        all)
            run_all_experiments
            ;;
        quick)
            run_quick_test
            ;;
        hybrid)
            print_header "HYBRID IMPLEMENTATION TEST"
            echo "Test|DGEMM_Type|Mult_Method|Inv_Method|Block_Size|Area|DGEMM_Time|Inv_Time|Comp_Time|Wall_Time" \
                > "$OUTPUT_DIR/results.csv"
            run_test 0 0 1 32 "baseline_sequential"
            for M in 1 2 3; do
                run_test 0 $M 0 64 "hybrid_method${M}"
            done
            ;;
        openblas)
            print_header "OPENBLAS IMPLEMENTATION TEST"
            echo "Test|DGEMM_Type|Mult_Method|Inv_Method|Block_Size|Area|DGEMM_Time|Inv_Time|Comp_Time|Wall_Time" \
                > "$OUTPUT_DIR/results.csv"
            run_test 0 0 1 32 "baseline_sequential"
            for M in 1 2 3; do
                run_test 1 $M 0 64 "openblas_method${M}"
            done
            ;;
        compare)
            run_comparison
            ;;
        scaling)
            run_scaling
            ;;
        *)
            ts "Unknown mode: $MODE"
            ts "Use --help for usage information"
            exit 1
            ;;
    esac
    
    # Analyze results
    analyze_results
    
    # Generate LaTeX
    generate_latex_table
    
    # Summary
    print_header "BENCHMARK COMPLETE"
    
    ts "Output files:"
    ts "  $OUTPUT_DIR/"
    ls -la "$OUTPUT_DIR"/*.csv 2>/dev/null || true
    ls -la "$OUTPUT_DIR"/*.tex 2>/dev/null || true
    echo ""
    
    ts "Next steps for Section IV:"
    ts "  1. Review results.csv for data tables"
    ts "  2. Copy table_results.tex to paper"
    ts "  3. Run: python plot_results.py $OUTPUT_DIR"
    echo ""
}

main
