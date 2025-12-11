#!/usr/bin/env bash
# Quick Fix Script: แก้ไข LAPACK ipiv type mismatch bug
# จาก: int ipiv[n + 1]
# เป็น: lapack_int ipiv[n + 1]

set -e

echo "================================================================================"
echo "  LAPACK Bug Fix: ipiv Type Mismatch"
echo "================================================================================"
echo ""

# ตรวจสอบว่าอยู่ใน directory ที่ถูกต้อง
if [ ! -f "matrix_inv.c" ]; then
    echo "❌ Error: matrix_inv.c not found!"
    echo "   Please run this script in the directory containing matrix_inv.c"
    exit 1
fi

echo "✅ Found matrix_inv.c"
echo ""

# Backup file
BACKUP="matrix_inv.c.backup_$(date +%Y%m%d_%H%M%S)"
echo "📦 Creating backup: $BACKUP"
cp matrix_inv.c "$BACKUP"
echo "   ✅ Backup created successfully"
echo ""

# ตรวจสอบว่ามี bug หรือไม่
echo "🔍 Checking for bug..."
if grep -q "^[[:space:]]*int ipiv\[n + 1\];" matrix_inv.c; then
    echo "   ✅ Bug found at line 303: 'int ipiv[n + 1]'"
    echo ""
    
    # แสดง context
    echo "📄 Current code (before fix):"
    grep -n -B2 -A2 "int ipiv\[n + 1\]" matrix_inv.c | head -10
    echo ""
    
    # Apply fix
    echo "🔧 Applying fix..."
    sed -i.bak '303s/int ipiv\[n + 1\]/lapack_int ipiv[n + 1]/' matrix_inv.c
    
    # Verify fix
    echo "   ✅ Fix applied"
    echo ""
    echo "📄 Updated code (after fix):"
    grep -n -B2 -A2 "lapack_int ipiv\[n + 1\]" matrix_inv.c | head -10
    echo ""
    
    echo "✅ Bug fixed successfully!"
    echo ""
    
elif grep -q "^[[:space:]]*lapack_int ipiv\[n + 1\];" matrix_inv.c; then
    echo "   ℹ️  Code is already fixed (lapack_int ipiv[] exists)"
    echo "   No changes needed."
    echo ""
else
    echo "   ⚠️  Warning: Could not find expected pattern"
    echo "   Please check matrix_inv.c manually"
    echo ""
    exit 1
fi

echo "================================================================================"
echo "  Next Steps"
echo "================================================================================"
echo ""
echo "1. Compile the fixed code:"
echo "   cd /path/to/your/project"
echo "   make clean && make all"
echo ""
echo "2. Test the fix:"
echo "   ./catcharea 1.0 99.0 0.001 0 3 32 | grep -E '(ERROR|Catchment area)'"
echo ""
echo "3. Run full benchmark:"
echo "   ./benchmark_all_methods.sh 6 1.0 99.0 0.001"
echo ""
echo "Expected results after fix:"
echo "  - No 'ERROR: LAPACKE_dgetrf failed' messages"
echo "  - Catchment area is a normal number (not -nan)"
echo "  - All 4 methods produce consistent results"
echo ""
echo "If you need to restore the original:"
echo "  cp $BACKUP matrix_inv.c"
echo ""
echo "================================================================================"

exit 0
