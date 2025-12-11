#!/usr/bin/env bash
# Fixed LAPACK Benchmark - Tests with LARGER matrices
# This script tests with different problem sizes to find optimal threading point

set -u

# Helper functions
ts() { printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*"; }

print_header() {
    echo ""
    echo "================================================================================"
    echo "  $1"
    echo "================================================================================"
    echo ""
}

print_header "LAPACK BENCHMARK - MATRIX SIZE SCALING STUDY"

ts "Problem: Previous benchmark used N=420, too small for threading"
ts "Solution: Test with multiple problem sizes to find crossover point"
echo ""

# Check binary exists
if [ ! -f "./catcharea" ]; then
    ts "Error: ./catcharea not found!"
    ts "Please build: make clean && make all"
    exit 1
fi

# Create output directory
OUTPUT_DIR="scaling_study_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$OUTPUT_DIR"

LOG_FILE="$OUTPUT_DIR/scaling_study.log"
exec > >(tee -a "$LOG_FILE") 2>&1

ts "Output directory: $OUTPUT_DIR"
echo ""

# =============================================================================
# Test different problem sizes
# =============================================================================

# Dr parameter controls matrix size:
# Dr = 0.001  → N ≈ 420  (too small)
# Dr = 0.0005 → N ≈ 800  (borderline)
# Dr = 0.0002 → N ≈ 2000 (good for threading)
# Dr = 0.0001 → N ≈ 4000 (excellent for threading)

declare -a DR_VALUES=("0.001" "0.0005" "0.0003" "0.0002")
declare -a EXPECTED_N=("~420" "~800" "~1300" "~2000")
declare -a NUM_THREADS_VALUES=("1" "2" "4" "6")

print_header "PHASE 1: MATRIX SIZE DISCOVERY"

ts "Testing different Dr values to find matrix sizes..."
echo ""

# Quick test to find actual matrix sizes
for i in "${!DR_VALUES[@]}"; do
    DR="${DR_VALUES[$i]}"
    EXPECTED="${EXPECTED_N[$i]}"
    
    ts "Testing Dr=$DR (expected N$EXPECTED)..."
    
    export OPENBLAS_NUM_THREADS=1
    export OMP_NUM_THREADS=1
    
    # Run briefly and extract matrix size
    timeout 60s ./catcharea 1.0 100.0 "$DR" 0 3 32 > "$OUTPUT_DIR/size_test_$DR.txt" 2>&1 || true
    
    # Extract actual matrix size
    ACTUAL_N=$(grep "Matrix size:" "$OUTPUT_DIR/size_test_$DR.txt" | head -1 | awk '{print $3}')
    
    if [ -n "$ACTUAL_N" ]; then
        ts "  Dr=$DR → N=$ACTUAL_N ✅"
        
        # Determine threading recommendation
        if [ "$ACTUAL_N" -lt 400 ]; then
            ts "    → Too small for threading (use 1 thread)"
        elif [ "$ACTUAL_N" -lt 800 ]; then
            ts "    → OpenBLAS will use 2 threads max"
        elif [ "$ACTUAL_N" -lt 1500 ]; then
            ts "    → OpenBLAS will use 2-4 threads"
        else
            ts "    → OpenBLAS will use all threads (good!)"
        fi
    else
        ts "  Dr=$DR → Could not determine size ❌"
    fi
    echo ""
done

print_header "PHASE 2: THREADING SCALING STUDY"

ts "Testing thread scaling for each problem size..."
echo ""

# Create CSV header
echo "Dr,Matrix_Size,Threads,Inv_Time,GFLOPS,Actual_Threads,Speedup_vs_1T" > "$OUTPUT_DIR/scaling_results.csv"

# For each problem size
for i in "${!DR_VALUES[@]}"; do
    DR="${DR_VALUES[$i]}"
    EXPECTED="${EXPECTED_N[$i]}"
    
    print_header "Problem Size: Dr=$DR (N$EXPECTED)"
    
    # Get actual matrix size from previous test
    ACTUAL_N=$(grep "Matrix size:" "$OUTPUT_DIR/size_test_$DR.txt" | head -1 | awk '{print $3}')
    
    if [ -z "$ACTUAL_N" ]; then
        ts "Skipping Dr=$DR (size unknown)"
        continue
    fi
    
    # Baseline time (1 thread)
    BASELINE_TIME=""
    
    # Test each thread count
    for THREADS in "${NUM_THREADS_VALUES[@]}"; do
        ts "Running with $THREADS thread(s)..."
        
        export OPENBLAS_NUM_THREADS=$THREADS
        export OMP_NUM_THREADS=$THREADS
        export OMP_PROC_BIND=close
        export OMP_PLACES=cores
        
        OUTPUT_FILE="$OUTPUT_DIR/dr${DR}_t${THREADS}.txt"
        TIMING_FILE="$OUTPUT_DIR/dr${DR}_t${THREADS}_timing.txt"
        
        # Run benchmark
        { time ./catcharea 1.0 100.0 "$DR" 0 3 32 > "$OUTPUT_FILE" 2>&1; } 2> "$TIMING_FILE"
        
        if [ $? -eq 0 ]; then
            # Extract metrics
            INV_TIME=$(sed -n 's/^[[:space:]]*Time:[[:space:]]*\([0-9.]\+\) seconds.*/\1/p' "$OUTPUT_FILE" | head -1)
            GFLOPS=$(sed -n 's/^[[:space:]]*GFLOPS:[[:space:]]*\([0-9.]\+\).*/\1/p' "$OUTPUT_FILE" | head -1)
            ACTUAL_THREADS=$(sed -n 's/.*Actual threads in use: \([0-9]\+\).*/\1/p' "$OUTPUT_FILE" | head -1)
            
            # Save baseline
            if [ "$THREADS" = "1" ]; then
                BASELINE_TIME="$INV_TIME"
            fi
            
            # Calculate speedup
            if [ -n "$BASELINE_TIME" ] && [ -n "$INV_TIME" ]; then
                SPEEDUP=$(awk -v b="$BASELINE_TIME" -v t="$INV_TIME" 'BEGIN {printf "%.2f", b/t}')
            else
                SPEEDUP="N/A"
            fi
            
            ts "  Threads: $THREADS → Actual: ${ACTUAL_THREADS:-N/A}, Time: ${INV_TIME}s, GFLOPS: ${GFLOPS}, Speedup: ${SPEEDUP}×"
            
            # Save to CSV
            echo "$DR,$ACTUAL_N,$THREADS,$INV_TIME,$GFLOPS,$ACTUAL_THREADS,$SPEEDUP" >> "$OUTPUT_DIR/scaling_results.csv"
        else
            ts "  FAILED ❌"
        fi
    done
    echo ""
done

print_header "PHASE 3: RESULTS ANALYSIS"

ts "Analyzing scaling behavior..."
echo ""

# Read results and analyze
if [ -f "$OUTPUT_DIR/scaling_results.csv" ]; then
    ts "Scaling Results Summary:"
    echo ""
    
    printf "%-8s %-12s %-8s %-12s %-10s %-10s %-10s\n" \
        "Dr" "Matrix_N" "Threads" "Inv_Time(s)" "GFLOPS" "Actual_T" "Speedup"
    printf "%-8s %-12s %-8s %-12s %-10s %-10s %-10s\n" \
        "--------" "------------" "--------" "------------" "----------" "----------" "----------"
    
    tail -n +2 "$OUTPUT_DIR/scaling_results.csv" | while IFS=, read -r dr n threads time gflops actual speedup; do
        printf "%-8s %-12s %-8s %-12s %-10s %-10s %-10s\n" \
            "$dr" "$n" "$threads" "$time" "$gflops" "$actual" "${speedup}×"
    done
    echo ""
    
    # Find best configuration
    ts "Recommendations:"
    echo ""
    
    # Find largest N with good speedup
    BEST_LINE=$(tail -n +2 "$OUTPUT_DIR/scaling_results.csv" | \
                awk -F, '$3=="6" && $7!="N/A" {print $1,$2,$7}' | \
                sort -k3 -nr | head -1)
    
    if [ -n "$BEST_LINE" ]; then
        read DR_BEST N_BEST SPEEDUP_BEST <<< "$BEST_LINE"
        
        ts "✅ Best configuration found:"
        ts "   Dr = $DR_BEST"
        ts "   Matrix size = $N_BEST"
        ts "   Speedup with 6 threads = ${SPEEDUP_BEST}×"
        echo ""
        
        if (( $(echo "$SPEEDUP_BEST > 3.0" | bc -l) )); then
            ts "✅ Excellent scaling! Use this for paper."
        elif (( $(echo "$SPEEDUP_BEST > 2.0" | bc -l) )); then
            ts "✅ Good scaling. Acceptable for paper."
        elif (( $(echo "$SPEEDUP_BEST > 1.5" | bc -l) )); then
            ts "⚠️  Moderate scaling. Consider larger problem."
        else
            ts "❌ Poor scaling. Need larger matrix size."
        fi
        echo ""
        
        ts "To use this configuration:"
        ts "  ./catcharea 1.0 100.0 $DR_BEST 0 3 32"
        echo ""
    else
        ts "❌ No good configuration found. All matrices too small."
        echo ""
    fi
    
    # Matrix size recommendations
    ts "Matrix Size Guidelines:"
    ts "  N < 400:    Use 1 thread (no benefit from threading)"
    ts "  400-800:    OpenBLAS uses 2 threads"
    ts "  800-1500:   OpenBLAS uses 2-4 threads"
    ts "  N > 1500:   OpenBLAS uses all threads (best)"
    echo ""
    
    # Find which Dr gives N > 1500
    GOOD_DR=$(tail -n +2 "$OUTPUT_DIR/scaling_results.csv" | \
              awk -F, '$2 > 1500 {print $1; exit}')
    
    if [ -n "$GOOD_DR" ]; then
        ts "✅ For optimal threading, use Dr ≤ $GOOD_DR"
    else
        ts "⚠️  None of the tested Dr values gave N > 1500"
        ts "   Try smaller Dr values (e.g., 0.0001)"
    fi
    echo ""
fi

print_header "PHASE 4: CREATING VISUALIZATION DATA"

# Create data for plotting
cat > "$OUTPUT_DIR/plot_data.py" << 'PYTHON_EOF'
#!/usr/bin/env python3
"""
Generate plots from scaling study results
"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read data
df = pd.read_csv('scaling_results.csv')

# Convert types
df['Threads'] = pd.to_numeric(df['Threads'])
df['Speedup_vs_1T'] = pd.to_numeric(df['Speedup_vs_1T'], errors='coerce')
df['GFLOPS'] = pd.to_numeric(df['GFLOPS'], errors='coerce')

# Create figure with subplots
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
fig.suptitle('LAPACK Threading Scaling Study', fontsize=16)

# Plot 1: Speedup vs Threads (for each Dr)
ax1 = axes[0, 0]
for dr in df['Dr'].unique():
    data = df[df['Dr'] == dr]
    n = data['Matrix_Size'].iloc[0]
    ax1.plot(data['Threads'], data['Speedup_vs_1T'], 
             marker='o', label=f'Dr={dr} (N={n})')
ax1.set_xlabel('Number of Threads')
ax1.set_ylabel('Speedup vs 1 Thread')
ax1.set_title('Thread Scaling by Problem Size')
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.axhline(y=1, color='r', linestyle='--', alpha=0.5, label='No speedup')

# Plot 2: GFLOPS vs Threads
ax2 = axes[0, 1]
for dr in df['Dr'].unique():
    data = df[df['Dr'] == dr]
    n = data['Matrix_Size'].iloc[0]
    ax2.plot(data['Threads'], data['GFLOPS'], 
             marker='s', label=f'Dr={dr} (N={n})')
ax2.set_xlabel('Number of Threads')
ax2.set_ylabel('GFLOPS')
ax2.set_title('Performance vs Thread Count')
ax2.legend()
ax2.grid(True, alpha=0.3)

# Plot 3: Efficiency (Speedup/Threads)
ax3 = axes[1, 0]
for dr in df['Dr'].unique():
    data = df[df['Dr'] == dr]
    n = data['Matrix_Size'].iloc[0]
    efficiency = data['Speedup_vs_1T'] / data['Threads'] * 100
    ax3.plot(data['Threads'], efficiency, 
             marker='^', label=f'Dr={dr} (N={n})')
ax3.set_xlabel('Number of Threads')
ax3.set_ylabel('Parallel Efficiency (%)')
ax3.set_title('Threading Efficiency')
ax3.legend()
ax3.grid(True, alpha=0.3)
ax3.axhline(y=100, color='g', linestyle='--', alpha=0.5, label='Perfect')
ax3.set_ylim(0, 120)

# Plot 4: Speedup vs Matrix Size (at 6 threads)
ax4 = axes[1, 1]
data_6t = df[df['Threads'] == 6]
ax4.scatter(data_6t['Matrix_Size'], data_6t['Speedup_vs_1T'], 
            s=100, c='red', marker='o')
for idx, row in data_6t.iterrows():
    ax4.annotate(f"Dr={row['Dr']}", 
                 (row['Matrix_Size'], row['Speedup_vs_1T']),
                 fontsize=8, ha='right')
ax4.set_xlabel('Matrix Size (N)')
ax4.set_ylabel('Speedup (6 threads vs 1 thread)')
ax4.set_title('Impact of Matrix Size on Parallel Speedup')
ax4.grid(True, alpha=0.3)
ax4.axhline(y=4, color='g', linestyle='--', alpha=0.5, label='Target (4×)')
ax4.axvline(x=1500, color='b', linestyle='--', alpha=0.5, label='OpenBLAS threshold')

plt.tight_layout()
plt.savefig('scaling_study_plots.png', dpi=150)
print("Plots saved to: scaling_study_plots.png")
PYTHON_EOF

chmod +x "$OUTPUT_DIR/plot_data.py"

ts "Visualization script created: $OUTPUT_DIR/plot_data.py"
ts "To generate plots (requires matplotlib):"
ts "  cd $OUTPUT_DIR && python3 plot_data.py"
echo ""

print_header "STUDY COMPLETE"

ts "All results saved to: $OUTPUT_DIR/"
ts ""
ts "Files:"
ts "  - scaling_study.log:       Full log"
ts "  - scaling_results.csv:     Summary data"
ts "  - size_test_*.txt:         Matrix size tests"
ts "  - dr*_t*.txt:              Detailed runs"
ts "  - plot_data.py:            Plotting script"
echo ""

ts "Next steps:"
ts "  1. Review scaling_results.csv to find best Dr value"
ts "  2. Run plot_data.py to visualize results"
ts "  3. Use recommended Dr in benchmark_fair_lapack.sh"
ts "  4. Update Section III with actual measured speedup"
echo ""

exit 0
