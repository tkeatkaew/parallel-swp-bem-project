#!/usr/bin/env python3
"""
compare_results.py

Compare results between sequential and Level 0 parallel implementations
to verify correctness and analyze performance.

Usage:
    python3 compare_results.py sequential.csv level0.csv
"""

import sys
import csv
from pathlib import Path

def load_results(filename):
    """Load results from CSV file."""
    results = []
    with open(filename, 'r') as f:
        reader = csv.DictReader(f, delimiter='|')
        for row in reader:
            results.append(row)
    return results

def compare_areas(seq_area, par_area, tolerance=1e-6):
    """Compare two area values within tolerance."""
    try:
        seq = float(seq_area)
        par = float(par_area)
        diff = abs(seq - par)
        rel_error = diff / abs(seq) if seq != 0 else diff
        
        return {
            'match': rel_error < tolerance,
            'diff': diff,
            'rel_error': rel_error,
            'tolerance': tolerance
        }
    except (ValueError, TypeError):
        return {'match': False, 'diff': None, 'rel_error': None}

def print_header(text):
    """Print formatted header."""
    print()
    print("═" * 70)
    print(f"  {text}")
    print("═" * 70)
    print()

def main():
    if len(sys.argv) != 3:
        print("Usage: python3 compare_results.py sequential.csv level0.csv")
        sys.exit(1)
    
    seq_file = Path(sys.argv[1])
    par_file = Path(sys.argv[2])
    
    # Check files exist
    if not seq_file.exists():
        print(f"Error: Sequential file not found: {seq_file}")
        sys.exit(1)
    if not par_file.exists():
        print(f"Error: Parallel file not found: {par_file}")
        sys.exit(1)
    
    # Load results
    print_header("LEVEL 0 CORRECTNESS VALIDATION")
    print(f"Sequential file: {seq_file}")
    print(f"Parallel file:   {par_file}")
    
    seq_results = load_results(seq_file)
    par_results = load_results(par_file)
    
    print(f"\nSequential results: {len(seq_results)} entries")
    print(f"Parallel results:   {len(par_results)} entries")
    
    # Compare catchment areas
    print_header("CATCHMENT AREA COMPARISON")
    
    print(f"{'Method':<30} {'Area':<15} {'Time (s)':<12} {'Speedup':<10}")
    print("-" * 70)
    
    all_match = True
    speedups = []
    
    # Get sequential baseline
    seq_baseline = None
    for r in seq_results:
        if 'Sequential' in r.get('Name', ''):
            seq_baseline = {
                'area': r.get('Catchment_Area', r.get('Area', 'N/A')),
                'time': float(r.get('Wall_Clock', r.get('Wall_Time', '0')))
            }
            print(f"{'Sequential (baseline)':<30} {seq_baseline['area']:<15} "
                  f"{seq_baseline['time']:<12.2f} {'1.00×':<10}")
            break
    
    # Compare parallel results
    for r in par_results:
        name = r.get('Name', 'Unknown')
        if 'threads' not in name.lower():
            continue
            
        area = r.get('Area', 'N/A')
        try:
            time = float(r.get('Wall_Clock', r.get('Wall_Time', '0')))
        except:
            time = 0.0
        
        # Calculate speedup
        if seq_baseline and time > 0:
            speedup = seq_baseline['time'] / time
            speedups.append((name, speedup))
            speedup_str = f"{speedup:.2f}×"
        else:
            speedup_str = "N/A"
        
        print(f"{name:<30} {area:<15} {time:<12.2f} {speedup_str:<10}")
        
        # Check area match
        if seq_baseline:
            cmp = compare_areas(seq_baseline['area'], area)
            if not cmp['match']:
                all_match = False
                print(f"  ⚠️  Area mismatch! Error: {cmp['rel_error']:.2e}")
    
    # Verification result
    print_header("VERIFICATION RESULT")
    
    if all_match:
        print("✅ PASSED: All parallel configurations produce correct results!")
        print(f"   Areas match sequential baseline within tolerance (1e-6)")
    else:
        print("❌ FAILED: Some configurations produce incorrect results!")
        print(f"   Please investigate discrepancies above.")
    
    # Performance analysis
    if speedups:
        print_header("PERFORMANCE ANALYSIS")
        
        print("Speedup Summary:")
        for name, speedup in speedups:
            threads = 1
            if 'threads' in name.lower():
                try:
                    threads = int(name.split('_')[1].replace('threads', ''))
                except:
                    pass
            
            efficiency = (speedup / threads) * 100 if threads > 0 else 0
            print(f"  {name:<30} Speedup: {speedup:>6.2f}×  "
                  f"Efficiency: {efficiency:>5.1f}%")
        
        # Find best configuration
        best = max(speedups, key=lambda x: x[1])
        print(f"\nBest configuration: {best[0]}")
        print(f"  Speedup: {best[1]:.2f}×")
        
        # Check if speedup is reasonable
        print("\nSpeedup Analysis:")
        for name, speedup in speedups:
            try:
                threads = int(name.split('_')[1].replace('threads', ''))
                expected_min = threads * 0.8  # 80% efficiency
                expected_max = threads * 1.0  # 100% efficiency
                
                if speedup >= expected_min:
                    status = "✅ Good"
                elif speedup >= threads * 0.6:
                    status = "⚠️  Acceptable"
                else:
                    status = "❌ Poor"
                
                print(f"  {threads} threads: Expected {expected_min:.1f}-{expected_max:.1f}×, "
                      f"Got {speedup:.2f}× {status}")
            except:
                pass
    
    print()

if __name__ == '__main__':
    main()
