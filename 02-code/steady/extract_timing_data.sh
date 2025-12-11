#!/bin/bash
# Extract timing data from method files for comparison

echo "==================================================================="
echo "TIMING DATA EXTRACTION FROM METHOD OUTPUT FILES"
echo "==================================================================="
echo ""

for i in 0 1 2 3; do
    file="/mnt/user-data/uploads/method${i}.txt"
    
    if [ ! -f "$file" ]; then
        echo "Method $i: File not found"
        continue
    fi
    
    echo "Method $i:"
    echo "----------"
    
    # Extract DGEMM time
    dgemm=$(grep "Total DGEMM (matrix multiply) time:" "$file" | awk '{print $6}')
    
    # Extract Inversion time
    inversion=$(grep "Total Matrix Inversion time:" "$file" | awk '{print $5}')
    
    # Extract Total computation time
    computation=$(grep "Total computation time:" "$file" | awk '{print $4}')
    
    # Extract from CSV
    csv_multiply=$(grep "^Multiply_Time_sec," "$file" | cut -d',' -f2)
    csv_inversion=$(grep "^Inversion_Time_sec," "$file" | cut -d',' -f2)
    
    echo "  DGEMM time:             ${dgemm}s"
    echo "  Inversion time:         ${inversion}s"
    echo "  Total computation:      ${computation}s"
    echo "  ---"
    echo "  CSV Multiply_Time_sec:  ${csv_multiply}s"
    echo "  CSV Inversion_Time_sec: ${csv_inversion}s"
    echo ""
    
    # Verify computation = dgemm + inversion
    if [ -n "$dgemm" ] && [ -n "$inversion" ] && [ -n "$computation" ]; then
        sum=$(awk -v d="$dgemm" -v i="$inversion" 'BEGIN {printf "%.6f", d+i}')
        match=$(awk -v s="$sum" -v c="$computation" 'BEGIN {
            diff = (s-c > 0) ? s-c : c-s
            if (diff < 0.01) print "✅ MATCH"
            else print "❌ MISMATCH"
        }')
        echo "  Verification: $dgemm + $inversion = $sum  $match ($computation)"
    fi
    echo ""
done

echo "==================================================================="
echo "COMPARISON TABLE"
echo "==================================================================="
echo ""
printf "%-8s %12s %12s %12s %12s\n" "Method" "DGEMM(s)" "Inversion(s)" "Total(s)" "CSV Issue?"
printf "%-8s %12s %12s %12s %12s\n" "------" "--------" "------------" "--------" "-----------"

for i in 0 1 2 3; do
    file="/mnt/user-data/uploads/method${i}.txt"
    
    if [ ! -f "$file" ]; then
        continue
    fi
    
    dgemm=$(grep "Total DGEMM (matrix multiply) time:" "$file" | awk '{print $6}')
    inversion=$(grep "Total Matrix Inversion time:" "$file" | awk '{print $5}')
    computation=$(grep "Total computation time:" "$file" | awk '{print $4}')
    csv_inversion=$(grep "^Inversion_Time_sec," "$file" | cut -d',' -f2)
    
    # Check if CSV inversion matches text inversion
    if [ -n "$csv_inversion" ] && [ -n "$inversion" ]; then
        diff=$(awk -v c="$csv_inversion" -v t="$inversion" 'BEGIN {
            d = (c-t > 0) ? c-t : t-c
            if (d > 1.0) print "❌ BAD"
            else print "✅ OK"
        }')
    else
        diff="N/A"
    fi
    
    printf "%-8s %12.2f %12.2f %12.2f %12s\n" \
        "Method$i" \
        "${dgemm:-0}" \
        "${inversion:-0}" \
        "${computation:-0}" \
        "$diff"
done

echo ""
echo "==================================================================="
echo "KEY FINDINGS"
echo "==================================================================="
echo ""
echo "1. 'Total computation time' = DGEMM + Inversion (verified above)"
echo "2. CSV Inversion_Time_sec may be wrong (check '❌ BAD' entries)"
echo "3. Benchmark script should parse 'Total computation time' from text"
echo "4. This is the CORRECT value for Matrix Operations Time"
echo ""
