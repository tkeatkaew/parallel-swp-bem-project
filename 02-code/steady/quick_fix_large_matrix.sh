#!/usr/bin/env bash
# QUICK FIX: Test with larger matrix immediately
# This will give you proper LAPACK speedup measurement

set -u

ts() { printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*"; }

echo ""
echo "================================================================================"
echo "  QUICK FIX: Testing with LARGER Matrix (N~2000)"
echo "================================================================================"
echo ""

ts "Problem identified:"
ts "  Previous test: N=420 (too small)"
ts "  OpenBLAS used: 2 threads instead of 6"
ts "  Result: SLOWDOWN instead of speedup"
echo ""

ts "Solution:"
ts "  Use Dr=0.0002 to generate N~2000"
ts "  This will activate full OpenBLAS threading"
echo ""

# Check binary
if [ ! -f "./catcharea" ]; then
    ts "Error: ./catcharea not found!"
    exit 1
fi

# Create output directory
OUTPUT_DIR="quick_fix_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTPUT_DIR"

LOG_FILE="$OUTPUT_DIR/quick_fix.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "================================================================================"
echo "  TEST 1: Sequential (1 thread) with LARGE matrix"
echo "================================================================================"
echo ""

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OMP_PROC_BIND=close
export OMP_PLACES=cores

ts "Configuration:"
ts "  OPENBLAS_NUM_THREADS: 1"
ts "  Dr: 0.0002 (generates N~2000)"
echo ""

ts "Running sequential test..."
{ time ./catcharea 1.0 100.0 0.0002 0 3 32 \
  > "$OUTPUT_DIR/sequential_large.txt" 2>&1; } 2> "$OUTPUT_DIR/sequential_large_timing.txt"

# Extract results
SEQ_N=$(grep "Matrix size:" "$OUTPUT_DIR/sequential_large.txt" | head -1 | awk '{print $3}')
SEQ_TIME=$(sed -n 's/^[[:space:]]*Time:[[:space:]]*\([0-9.]\+\) seconds.*/\1/p' \
           "$OUTPUT_DIR/sequential_large.txt" | head -1)
SEQ_GFLOPS=$(sed -n 's/^[[:space:]]*GFLOPS:[[:space:]]*\([0-9.]\+\).*/\1/p' \
             "$OUTPUT_DIR/sequential_large.txt" | head -1)
SEQ_THREADS=$(sed -n 's/.*Actual threads in use: \([0-9]\+\).*/\1/p' \
              "$OUTPUT_DIR/sequential_large.txt" | head -1)

if [ -n "$SEQ_N" ]; then
    ts "✅ Sequential test completed"
    ts "   Matrix size: $SEQ_N"
    ts "   Time: ${SEQ_TIME}s"
    ts "   GFLOPS: $SEQ_GFLOPS"
    ts "   Threads: ${SEQ_THREADS:-1}"
    
    if [ "$SEQ_N" -lt 1500 ]; then
        echo ""
        ts "⚠️  WARNING: Matrix still too small (N=$SEQ_N < 1500)"
        ts "   Try even smaller Dr (e.g., 0.0001)"
        echo ""
    else
        echo ""
        ts "✅ Matrix size is good (N=$SEQ_N ≥ 1500)"
        echo ""
    fi
else
    ts "❌ Could not extract results"
    exit 1
fi

echo ""
ts "Waiting 3 seconds for system to stabilize..."
sleep 3
echo ""

echo "================================================================================"
echo "  TEST 2: Parallel (6 threads) with LARGE matrix"
echo "================================================================================"
echo ""

export OPENBLAS_NUM_THREADS=6
export OMP_NUM_THREADS=6

ts "Configuration:"
ts "  OPENBLAS_NUM_THREADS: 6"
ts "  Dr: 0.0002 (same as sequential)"
echo ""

ts "Running parallel test..."
{ time ./catcharea 1.0 100.0 0.0002 0 3 32 \
  > "$OUTPUT_DIR/parallel_large.txt" 2>&1; } 2> "$OUTPUT_DIR/parallel_large_timing.txt"

# Extract results
PAR_N=$(grep "Matrix size:" "$OUTPUT_DIR/parallel_large.txt" | head -1 | awk '{print $3}')
PAR_TIME=$(sed -n 's/^[[:space:]]*Time:[[:space:]]*\([0-9.]\+\) seconds.*/\1/p' \
           "$OUTPUT_DIR/parallel_large.txt" | head -1)
PAR_GFLOPS=$(sed -n 's/^[[:space:]]*GFLOPS:[[:space:]]*\([0-9.]\+\).*/\1/p' \
             "$OUTPUT_DIR/parallel_large.txt" | head -1)
PAR_THREADS=$(sed -n 's/.*Actual threads in use: \([0-9]\+\).*/\1/p' \
              "$OUTPUT_DIR/parallel_large.txt" | head -1)

if [ -n "$PAR_N" ]; then
    ts "✅ Parallel test completed"
    ts "   Matrix size: $PAR_N"
    ts "   Time: ${PAR_TIME}s"
    ts "   GFLOPS: $PAR_GFLOPS"
    ts "   Threads: ${PAR_THREADS:-6}"
    echo ""
else
    ts "❌ Could not extract results"
    exit 1
fi

echo "================================================================================"
echo "  RESULTS COMPARISON"
echo "================================================================================"
echo ""

# Verify thread counts
ts "Thread Verification:"
ts "  Sequential requested: 1, actual: ${SEQ_THREADS:-N/A}"
ts "  Parallel requested:   6, actual: ${PAR_THREADS:-N/A}"
echo ""

if [ "${SEQ_THREADS:-1}" != "1" ]; then
    ts "⚠️  WARNING: Sequential didn't use 1 thread!"
fi

if [ "${PAR_THREADS:-0}" != "6" ]; then
    ts "⚠️  WARNING: Parallel didn't use 6 threads!"
    if [ "${PAR_THREADS:-0}" = "2" ]; then
        ts "    Matrix still too small (OpenBLAS limited to 2 threads)"
        ts "    Need N > 1500 for full threading"
    fi
fi
echo ""

# Calculate speedup
if [ -n "$SEQ_TIME" ] && [ -n "$PAR_TIME" ]; then
    SPEEDUP=$(awk -v s="$SEQ_TIME" -v p="$PAR_TIME" 'BEGIN {printf "%.2f", s/p}')
    EFFICIENCY=$(awk -v sp="$SPEEDUP" 'BEGIN {printf "%.1f", sp/6*100}')
    
    printf "%-25s %12s %12s %10s\n" "Configuration" "Time(s)" "GFLOPS" "Threads"
    printf "%-25s %12s %12s %10s\n" "-------------------------" "------------" "------------" "----------"
    printf "%-25s %12s %12s %10s\n" "Sequential (1 thread)" "$SEQ_TIME" "$SEQ_GFLOPS" "${SEQ_THREADS:-1}"
    printf "%-25s %12s %12s %10s\n" "Parallel (6 threads)" "$PAR_TIME" "$PAR_GFLOPS" "${PAR_THREADS:-6}"
    echo ""
    
    ts "===== PERFORMANCE SUMMARY ====="
    ts ""
    ts "Matrix size:     N=$SEQ_N"
    ts "Speedup:         ${SPEEDUP}×"
    ts "Efficiency:      ${EFFICIENCY}%"
    ts ""
    
    # Interpret results
    if (( $(echo "$SPEEDUP < 1.0" | bc -l) )); then
        ts "❌ RESULT: SLOWDOWN (parallel is slower!)"
        ts "   Cause: Matrix still too small or other issues"
        echo ""
    elif (( $(echo "$SPEEDUP < 2.0" | bc -l) )); then
        ts "⚠️  RESULT: Poor speedup"
        ts "   Matrix size N=$SEQ_N may still be too small"
        ts "   Try Dr=0.0001 for larger matrix"
        echo ""
    elif (( $(echo "$SPEEDUP < 3.5" | bc -l) )); then
        ts "✅ RESULT: Moderate speedup"
        ts "   Acceptable but not optimal"
        ts "   Consider larger problem for better scaling"
        echo ""
    elif (( $(echo "$SPEEDUP < 5.0" | bc -l) )); then
        ts "✅ RESULT: Good speedup!"
        ts "   This is acceptable for paper"
        echo ""
    else
        ts "✅ RESULT: Excellent speedup!"
        ts "   Great scaling for dense linear algebra"
        echo ""
    fi
    
    # Recommendations
    ts "Recommendations:"
    if (( $(echo "$SPEEDUP < 3.0" | bc -l) )); then
        ts "  1. Try smaller Dr to get larger matrix:"
        ts "     ./catcharea 1.0 100.0 0.0001 0 3 32"
        ts ""
        ts "  2. Verify OpenBLAS configuration:"
        ts "     ldd ./catcharea | grep openblas"
        ts "     export OPENBLAS_VERBOSE=1"
        ts ""
        ts "  3. Check actual matrix sizes in output logs"
    else
        ts "  ✅ Results are good!"
        ts "  1. Use these results in your paper"
        ts "  2. Report speedup: ${SPEEDUP}×"
        ts "  3. Report efficiency: ${EFFICIENCY}%"
        ts ""
        ts "  4. Create table for Section III:"
        ts "     | Threads | Time(s) | GFLOPS | Speedup |"
        ts "     |---------|---------|---------|---------|"
        ts "     | 1       | $SEQ_TIME | $SEQ_GFLOPS | 1.00×   |"
        ts "     | 6       | $PAR_TIME | $PAR_GFLOPS | ${SPEEDUP}× |"
    fi
    echo ""
else
    ts "❌ Could not calculate speedup (missing timing data)"
    exit 1
fi

echo "================================================================================"
echo "  FILES CREATED"
echo "================================================================================"
echo ""

ts "Output directory: $OUTPUT_DIR/"
ts "  - quick_fix.log:               This log"
ts "  - sequential_large.txt:        Sequential detailed output"
ts "  - parallel_large.txt:          Parallel detailed output"
ts "  - sequential_large_timing.txt: Sequential timing"
ts "  - parallel_large_timing.txt:   Parallel timing"
echo ""

echo "================================================================================"
echo "  NEXT STEPS"
echo "================================================================================"
echo ""

if (( $(echo "$SPEEDUP < 3.0" | bc -l) )); then
    ts "❌ Speedup still not satisfactory (${SPEEDUP}×)"
    echo ""
    ts "Run the full scaling study to find optimal Dr:"
    ts "  chmod +x benchmark_scaling_study.sh"
    ts "  ./benchmark_scaling_study.sh"
    echo ""
else
    ts "✅ Good speedup achieved (${SPEEDUP}×)!"
    echo ""
    ts "You can now:"
    ts "  1. Update benchmark_fair_lapack.sh to use Dr=0.0002"
    ts "  2. Run comprehensive benchmarks with this parameter"
    ts "  3. Update Section III with measured speedup: ${SPEEDUP}×"
    echo ""
fi

ts "To examine detailed output:"
ts "  less $OUTPUT_DIR/sequential_large.txt"
ts "  less $OUTPUT_DIR/parallel_large.txt"
echo ""

exit 0
