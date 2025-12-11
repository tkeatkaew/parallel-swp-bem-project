#!/usr/bin/env bash
#
# quick_setup.sh
#
# Quick setup script for Level 0 Contour-Layer Parallelism
# This script helps you get started quickly by copying and setting up all files
#
# Usage: ./quick_setup.sh

set -e

print_header() {
    echo ""
    echo "═══════════════════════════════════════════════════════════════"
    echo "  $1"
    echo "═══════════════════════════════════════════════════════════════"
    echo ""
}

print_step() {
    echo ""
    echo "▶ $1"
    echo ""
}

print_header "LEVEL 0 CONTOUR-LAYER PARALLELISM - QUICK SETUP"

# Get current directory
CURRENT_DIR=$(pwd)
echo "Current directory: $CURRENT_DIR"

# Check if we're in the right place
if [ ! -f "../source/catcharea.c" ]; then
    echo "❌ Error: Cannot find ../source/catcharea.c"
    echo "   Please run this script from: ~/works/parallel-swp-bem/02-code/steady"
    exit 1
fi

print_step "Step 1: Checking source files location"

if [ -f "/home/thanit/works/parallel-swp-bem/02-code/level0_implement1/zone_processor.h" ]; then
    SRC_DIR="/home/thanit/works/parallel-swp-bem/02-code/level0_implement1"
    echo "✅ Found Level 0 files in /home/thanit/works/parallel-swp-bem/02-code/level0_implement1"
elif [ -f "/mnt/user-data/outputs/zone_processor.h" ]; then
    SRC_DIR="/mnt/user-data/outputs"
    echo "✅ Found Level 0 files in /mnt/user-data/outputs"
else
    echo "❌ Error: Cannot find Level 0 source files"
    echo "   Expected locations:"
    echo "   - /home/claude/zone_processor.h"
    echo "   - /mnt/user-data/outputs/zone_processor.h"
    exit 1
fi

print_step "Step 2: Copying Level 0 files"

# Copy all Level 0 files
echo "Copying zone_processor.h..."
cp "$SRC_DIR/zone_processor.h" .

echo "Copying zone_processor.c..."
cp "$SRC_DIR/zone_processor.c" .

echo "Copying Makefile.level0..."
cp "$SRC_DIR/Makefile.level0" .

echo "Copying benchmark_level0.sh..."
cp "$SRC_DIR/benchmark_level0.sh" .

echo "Copying compare_results.py..."
cp "$SRC_DIR/compare_results.py" .

echo "✅ All files copied successfully"

print_step "Step 3: Setting permissions"

chmod +x benchmark_level0.sh
chmod +x compare_results.py

echo "✅ Permissions set"

print_step "Step 4: Creating catcharea_level0.c"

if [ -f "catcharea_level0.c" ]; then
    echo "⚠️  catcharea_level0.c already exists"
    read -p "Overwrite? (y/n) " -n 1 -r
    echo
    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
        echo "Skipping catcharea_level0.c creation"
    else
        cp ../source/catcharea.c catcharea_level0.c
        echo "✅ Created catcharea_level0.c from catcharea.c"
        echo "⚠️  YOU MUST MANUALLY EDIT catcharea_level0.c:"
        echo "   1. Add: #include \"zone_processor.h\""
        echo "   2. Replace catchment_area() call with parallel zone loop"
        echo "   3. Add zone statistics output"
    fi
else
    cp ../source/catcharea.c catcharea_level0.c
    echo "✅ Created catcharea_level0.c from catcharea.c"
    echo ""
    echo "⚠️  IMPORTANT: YOU MUST MANUALLY EDIT catcharea_level0.c"
    echo ""
    echo "Required modifications:"
    echo "  1. Line ~10: Add #include \"zone_processor.h\""
    echo "  2. Line ~274: Replace catchment_area() with parallel loop"
    echo "  3. Add zone statistics output"
    echo ""
    echo "See LEVEL0_SUMMARY.md for detailed instructions"
fi

print_step "Step 5: Copying documentation"

for doc in VALIDATION_STRATEGY.md IMPLEMENTATION_GUIDE.md LEVEL0_SUMMARY.md; do
    if [ -f "$SRC_DIR/$doc" ]; then
        cp "$SRC_DIR/$doc" .
        echo "✅ Copied $doc"
    fi
done

print_step "Step 6: Verifying setup"

echo "Checking files..."
files=(
    "zone_processor.h"
    "zone_processor.c"
    "Makefile.level0"
    "benchmark_level0.sh"
    "compare_results.py"
    "catcharea_level0.c"
)

all_present=true
for f in "${files[@]}"; do
    if [ -f "$f" ]; then
        echo "  ✅ $f"
    else
        echo "  ❌ $f (missing)"
        all_present=false
    fi
done

if [ "$all_present" = false ]; then
    echo ""
    echo "⚠️  Some files are missing!"
    exit 1
fi

print_header "SETUP COMPLETE!"

echo "Next steps:"
echo ""
echo "1. Edit catcharea_level0.c:"
echo "   See LEVEL0_SUMMARY.md Section: 'Option A: Minimal Changes Approach'"
echo ""
echo "2. Build:"
echo "   make -f Makefile.level0 clean"
echo "   make -f Makefile.level0 catcharea_level0"
echo ""
echo "3. Test:"
echo "   export OMP_NUM_THREADS=2"
echo "   export OPENBLAS_NUM_THREADS=1"
echo "   ./catcharea_level0 1.0 99.0 0.001 1 3 32"
echo ""
echo "4. Benchmark:"
echo "   ./benchmark_level0.sh 6"
echo ""
echo "5. Validate:"
echo "   python3 compare_results.py \\"
echo "       benchmark_results_sequential/results.csv \\"
echo "       benchmark_results_level0_*/results.csv"
echo ""

print_header "DOCUMENTATION"

echo "Available documentation:"
echo "  • LEVEL0_SUMMARY.md          - Quick overview and status"
echo "  • IMPLEMENTATION_GUIDE.md    - Detailed implementation guide"
echo "  • VALIDATION_STRATEGY.md     - Correctness validation plan"
echo ""

echo "🚀 Ready to implement Level 0 Contour-Layer Parallelism!"
echo ""
