#!/usr/bin/env bash
# Enhanced build and run script for catcharea with matrix optimization support
# Supports 4 matrix multiplication implementations + matrix inversion methods
# Usage: ./run_catcharea_optimized.sh [NUM_THREADS] [ARG1 ARG2 ARG3 [INVERSION_METHOD] [MATMUL_METHOD] [BLOCK_SIZE]]
#
# Arguments:
#   NUM_THREADS:       Number of OpenMP threads (default: 8)
#   ARG1:              Step size (default: 1.0)
#   ARG2:              Rm parameter (default: 100.0)
#   ARG3:              Dr parameter (default: 0.001)
#   INVERSION_METHOD:  0=Parallel/LAPACK (default), 1=Sequential
#   MATMUL_METHOD:     0=Sequential, 1=OpenMP, 2=OpenMP+Cache, 3=OpenMP+Cache+SIMD (default)
#   BLOCK_SIZE:        Cache block size 8-256 (default: 32)
#
# Examples:
#   ./run_catcharea_optimized.sh                    # Use all defaults
#   ./run_catcharea_optimized.sh 6                  # 6 threads, defaults for rest
#   ./run_catcharea_optimized.sh 6 1.0 99.0 0.001 0 3 32  # All parameters specified
#
set -u

# ===============================================================================
# CONFIGURATION AND ARGUMENT PARSING
# ===============================================================================

# Parse arguments with defaults
NUM_THREADS="${1:-8}"
ARG1="${2:-1.0}"
ARG2="${3:-100.0}"
ARG3="${4:-0.001}"
INVERSION_METHOD="${5:-0}"     # 0=Parallel (LAPACK), 1=Sequential
MATMUL_METHOD="${6:-3}"        # 0=Sequential, 1=OpenMP, 2=+Cache, 3=+Cache+SIMD
BLOCK_SIZE="${7:-32}"          # Cache block size

CMD="./catcharea $ARG1 $ARG2 $ARG3 $INVERSION_METHOD $MATMUL_METHOD $BLOCK_SIZE"

# Matrix method names for display
declare -A MATMUL_NAMES=(
    [0]="Sequential (baseline)"
    [1]="OpenMP (parallel)"
    [2]="OpenMP + Cache Blocking"
    [3]="OpenMP + Cache + SIMD (AVX2)"
)

# Validation
if [ "$MATMUL_METHOD" -lt 0 ] || [ "$MATMUL_METHOD" -gt 3 ]; then
    echo "Error: MATMUL_METHOD must be 0-3 (got: $MATMUL_METHOD)"
    exit 1
fi

if [ "$BLOCK_SIZE" -lt 8 ] || [ "$BLOCK_SIZE" -gt 256 ]; then
    echo "Warning: BLOCK_SIZE should be 8-256 (got: $BLOCK_SIZE), proceeding anyway..."
fi

# ===============================================================================
# HELPER FUNCTIONS
# ===============================================================================

ts() { printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*"; }
ts_block() {
  while IFS= read -r line; do ts "$line"; done
}

val_or_na() {
  if [ -n "${1:-}" ]; then echo "$1"; else echo "N/A"; fi
}

print_separator() {
    echo ""
    echo "================================================================================"
    echo "$1"
    echo "================================================================================"
    echo ""
}

# ===============================================================================
# SYSTEM INFORMATION GATHERING
# ===============================================================================

gather_system_info() {
    # OS information
    OS_KERNEL="$(uname -sr)"
    if command -v lsb_release >/dev/null 2>&1; then
        OS_NAME="$(lsb_release -ds 2>/dev/null | sed 's/^"//;s/"$//')"
    elif [ -r /etc/os-release ]; then
        # shellcheck disable=SC1091
        . /etc/os-release
        OS_NAME="${PRETTY_NAME:-Linux}"
    else
        OS_NAME="Linux"
    fi

    # CPU information
    CPU_MODEL="$(lscpu 2>/dev/null | awk -F: '/Model name/ {sub(/^[ \t]+/, "", $2); print $2; exit}')"
    if [ -z "$CPU_MODEL" ] && [ -r /proc/cpuinfo ]; then
        CPU_MODEL="$(awk -F: '/model name/ {print $2; exit}' /proc/cpuinfo | sed 's/^ //')"
    fi
    CPU_MODEL="$(val_or_na "$CPU_MODEL")"

    # Physical cores (unique Core,Socket pairs)
    if command -v lscpu >/dev/null 2>&1; then
        PHYS_CORES="$(lscpu -p=Core,Socket 2>/dev/null | awk -F, '!/^#/ {k=$1","$2; c[k]=1} END {print length(c)}')"
        LOGICAL_CORES="$(lscpu | awk '/^CPU\(s\):/ {print $2}')"
    else
        PHYS_CORES=""
        LOGICAL_CORES=""
    fi
    PHYS_CORES="$(val_or_na "$PHYS_CORES")"
    LOGICAL_CORES="$(val_or_na "$LOGICAL_CORES")"

    # Check for AVX2 support
    if grep -q avx2 /proc/cpuinfo 2>/dev/null; then
        AVX2_SUPPORT="Yes ✅"
    else
        AVX2_SUPPORT="No ❌ (method 3 will fall back to method 2)"
    fi

    # Check for FMA support
    if grep -q fma /proc/cpuinfo 2>/dev/null; then
        FMA_SUPPORT="Yes ✅"
    else
        FMA_SUPPORT="No"
    fi

    # Memory information
    if command -v free >/dev/null 2>&1; then
        MEM_TOTAL="$(free -h | awk '/^Mem:/ {print $2}')"
        MEM_AVAIL="$(free -h | awk '/^Mem:/ {print $7}')"
    elif [ -r /proc/meminfo ]; then
        kb=$(awk '/MemTotal:/ {print $2}' /proc/meminfo)
        if [ -n "$kb" ]; then
            MEM_TOTAL="$(awk -v kb="$kb" 'BEGIN { printf "%.1fGi", kb/1048576 }')"
        else
            MEM_TOTAL=""
        fi
        MEM_AVAIL=""
    else
        MEM_TOTAL=""
        MEM_AVAIL=""
    fi
    MEM_TOTAL="$(val_or_na "$MEM_TOTAL")"
    MEM_AVAIL="$(val_or_na "$MEM_AVAIL")"

    # Cache information
    CACHE_INFO="$(lscpu 2>/dev/null | awk -F: '
    /^(L1d cache|L1i cache|L2 cache|L3 cache)/ {
        gsub(/^[ \t]+/, "", $2);
        printf "%-35s %s\n", $1":", $2
    }')"

    # NUMA information
    NUMA_INFO="$(numactl --hardware 2>/dev/null | grep -E '(available|node.*cpus|node.*size)' || echo 'NUMA info not available')"
}

# ===============================================================================
# CHECK BUILD REQUIREMENTS
# ===============================================================================

check_build_requirements() {
    print_separator "CHECKING BUILD REQUIREMENTS"
    
    local all_ok=true
    
    # Check for gcc
    if command -v gcc >/dev/null 2>&1; then
        GCC_VERSION=$(gcc --version | head -n1)
        ts "GCC: $GCC_VERSION ✅"
    else
        ts "GCC: Not found ❌"
        all_ok=false
    fi
    
    # Check for OpenMP support
    if echo | cpp -fopenmp -dM 2>/dev/null | grep -q OPENMP; then
        OMP_VERSION=$(echo | cpp -fopenmp -dM 2>/dev/null | grep _OPENMP | awk '{print $3}')
        ts "OpenMP: Version $OMP_VERSION ✅"
    else
        ts "OpenMP: Not available ❌"
        all_ok=false
    fi
    
    # Check for OpenBLAS
    if ldconfig -p 2>/dev/null | grep -q libopenblas || [ -f /usr/lib/libopenblas.so ]; then
        ts "OpenBLAS: Found ✅"
    else
        ts "OpenBLAS: Not found ⚠️  (trying to continue)"
    fi
    
    # Check for LAPACK
    if ldconfig -p 2>/dev/null | grep -q liblapacke || [ -f /usr/lib/liblapacke.so ]; then
        ts "LAPACKE: Found ✅"
    else
        ts "LAPACKE: Not found ⚠️  (trying to continue)"
    fi
    
    # Check for make
    if command -v make >/dev/null 2>&1; then
        ts "Make: $(make --version | head -n1) ✅"
    else
        ts "Make: Not found ❌"
        all_ok=false
    fi
    
    echo ""
    
    if ! $all_ok; then
        ts "Error: Some required tools are missing!"
        ts "Please install: gcc, make, libopenblas-dev, liblapacke-dev"
        exit 1
    fi
}

# ===============================================================================
# SETUP ENVIRONMENT FOR OPTIMIZATION
# ===============================================================================

setup_environment() {
    print_separator "SETTING UP EXECUTION ENVIRONMENT"
    
    # OpenMP configuration
    export OMP_NUM_THREADS="$NUM_THREADS"
    #export OMP_PROC_BIND=close
    export OMP_PROC_BIND=spread
    export OMP_PLACES=cores
    export OMP_DISPLAY_ENV=VERBOSE
    
    # OpenBLAS configuration
    export OPENBLAS_NUM_THREADS=1  # Avoid oversubscription
    
    # Memory optimization
    export MALLOC_ARENA_MAX=4
    export MALLOC_MMAP_THRESHOLD_=131072
    export MALLOC_TRIM_THRESHOLD_=131072
    export MALLOC_TOP_PAD_=131072
    
    # Disable memory leak detection for performance
    export ASAN_OPTIONS=detect_leaks=0:malloc_context_size=0
    export LSAN_OPTIONS=exitcode=0
    
    ts "OpenMP Configuration:"
    ts "  OMP_NUM_THREADS=$OMP_NUM_THREADS"
    ts "  OMP_PROC_BIND=$OMP_PROC_BIND"
    ts "  OMP_PLACES=$OMP_PLACES"
    ts "  OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS"
    echo ""
    
    ts "Memory Optimization:"
    ts "  MALLOC_ARENA_MAX=$MALLOC_ARENA_MAX"
    ts "  MALLOC_MMAP_THRESHOLD_=$MALLOC_MMAP_THRESHOLD_"
    echo ""
    
    # CPU governor check
    if [ -r /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor ]; then
        GOVERNOR=$(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor)
        ts "CPU Governor: $GOVERNOR"
        if [ "$GOVERNOR" != "performance" ]; then
            ts "  ⚠️  Consider: sudo cpupower frequency-set -g performance"
        fi
        echo ""
    fi
    
    # Transparent Huge Pages (if possible)
    if [ -w /sys/kernel/mm/transparent_hugepage/enabled ]; then
        ts "Enabling Transparent Huge Pages..."
        echo always | sudo tee /sys/kernel/mm/transparent_hugepage/enabled > /dev/null 2>&1 || true
    fi
}

# ===============================================================================
# DISPLAY CONFIGURATION SUMMARY
# ===============================================================================

display_config_summary() {
    print_separator "EXECUTION CONFIGURATION"
    
    ts "Input Parameters:"
    ts "  ARG1 (step):              $ARG1"
    ts "  ARG2 (rm):                $ARG2"
    ts "  ARG3 (dr):                $ARG3"
    echo ""
    
    ts "Matrix Inversion Method:"
    if [ "$INVERSION_METHOD" = "0" ]; then
        ts "  Mode: PARALLEL (LAPACK) ⭐ [fastest]"
    elif [ "$INVERSION_METHOD" = "1" ]; then
        ts "  Mode: SEQUENTIAL (Gauss-Jordan) [for comparison]"
    else
        ts "  Mode: UNKNOWN (will default to PARALLEL)"
    fi
    echo ""
    
    ts "Matrix Multiplication Method:"
    ts "  Method $MATMUL_METHOD: ${MATMUL_NAMES[$MATMUL_METHOD]}"
    if [ "$MATMUL_METHOD" -ge 2 ]; then
        ts "  Block size: $BLOCK_SIZE"
    fi
    
    # Expected performance
    case $MATMUL_METHOD in
        0)
            ts "  Expected performance: Baseline (1.0×)"
            ;;
        1)
            ts "  Expected speedup: ~5.0× vs sequential"
            ;;
        2)
            ts "  Expected speedup: ~6.0× vs sequential"
            ;;
        3)
            ts "  Expected speedup: ~8-10× vs sequential ⭐"
            ;;
    esac
    
    if [ "$MATMUL_METHOD" = "3" ] && [ "$AVX2_SUPPORT" = "No ❌"* ]; then
        ts "  ⚠️  WARNING: AVX2 not available - will fall back to method 2"
    fi
    echo ""
    
    ts "Command to execute:"
    ts "  $CMD"
    echo ""
}

# ===============================================================================
# CLEAN BUILD ENVIRONMENT
# ===============================================================================

clean_build() {
    print_separator "CLEANING BUILD ENVIRONMENT"
    
    ts "Running: make clean"
    if make clean 2>&1 | ts_block; then
        ts "Clean completed successfully ✅"
    else
        ts "Warning: Clean failed (continuing anyway)"
    fi
    echo ""
}

# ===============================================================================
# BUILD PROJECT
# ===============================================================================

build_project() {
    print_separator "BUILDING PROJECT"
    
    # Check if Makefile exists
    if [ ! -f Makefile ]; then
        ts "Error: Makefile not found!"
        ts "Make sure you're in the correct directory with Makefile"
        exit 1
    fi
    
    # Show compilation flags
    ts "Compilation flags (from Makefile):"
    if grep -q "^CFLAGS" Makefile; then
        grep "^CFLAGS" Makefile | sed 's/^/  /' | ts_block
    fi
    echo ""
    
    # Build phase 1
    ts "Build Phase 1: Initial compilation..."
    if ! make all 2>&1 | ts_block; then
        ts "Error: Build phase 1 failed! ❌"
        exit 1
    fi
    ts "Build phase 1 completed ✅"
    echo ""
    
    # Build phase 2
    ts "Build Phase 2: Final compilation..."
    if ! make all 2>&1 | ts_block; then
        ts "Error: Build phase 2 failed! ❌"
        exit 1
    fi
    ts "Build phase 2 completed ✅"
    echo ""
    
    # Check binary
    if [ -f "./catcharea" ]; then
        BINARY_SIZE=$(du -h ./catcharea | cut -f1)
        ts "Binary created successfully:"
        ts "  File: ./catcharea"
        ts "  Size: $BINARY_SIZE"
        
        # Check for AVX2 instructions in binary
        if command -v objdump >/dev/null 2>&1; then
            if objdump -d ./catcharea 2>/dev/null | grep -q "vfmadd"; then
                ts "  AVX2 instructions: Detected ✅"
            else
                ts "  AVX2 instructions: Not found (method 3 will use scalar fallback)"
            fi
        fi
    else
        ts "Error: Binary ./catcharea not created! ❌"
        exit 1
    fi
    echo ""
}

# ===============================================================================
# RUN PROGRAM WITH TIMING
# ===============================================================================

run_program() {
    print_separator "EXECUTING CATCHAREA"
    
    ts "Command: $CMD"
    ts "Start time: $(date '+%Y-%m-%d %H:%M:%S')"
    ts "--------------------------------------------------------------------------------"
    echo ""
    
    # Save start state
    MEM_START=$(awk '/MemAvailable/ {print $2}' /proc/meminfo 2>/dev/null)
    
    # Run with detailed timing if available
    if command -v /usr/bin/time >/dev/null 2>&1; then
        ts "Running with detailed resource tracking..."
        echo ""
        if ! /usr/bin/time -v $CMD 2>&1 | ts_block; then
            ts "Error: Catcharea execution failed! ❌"
            exit 1
        fi
    else
        ts "Running (install 'time' package for detailed metrics)..."
        echo ""
        if ! { time $CMD; } 2>&1 | ts_block; then
            ts "Error: Catcharea execution failed! ❌"
            exit 1
        fi
    fi
    
    # Calculate memory usage
    MEM_END=$(awk '/MemAvailable/ {print $2}' /proc/meminfo 2>/dev/null)
    if [ -n "$MEM_START" ] && [ -n "$MEM_END" ]; then
        MEM_USED=$((MEM_START - MEM_END))
        MEM_USED_MB=$((MEM_USED / 1024))
        echo ""
        ts "Estimated memory used: ${MEM_USED_MB} MB"
    fi
    
    echo ""
    ts "--------------------------------------------------------------------------------"
    ts "End time: $(date '+%Y-%m-%d %H:%M:%S')"
    ts "Execution completed successfully ✅"
    echo ""
}

# ===============================================================================
# DISPLAY PERFORMANCE SUMMARY
# ===============================================================================

display_performance_summary() {
    print_separator "PERFORMANCE SUMMARY"
    
    # Extract timing from output files if they exist
    if [ -f "test.out" ] || [ -f "catchment.out" ]; then
        ts "Output files created successfully:"
        [ -f "test.out" ] && ts "  - test.out (streamlines)"
        [ -f "catchment.out" ] && ts "  - catchment.out (catchment data)"
        echo ""
    fi
    
    # Memory state after execution
    if command -v free >/dev/null 2>&1; then
        ts "Final memory state:"
        free -h | sed 's/^/  /' | ts_block
        echo ""
    fi
    
    # Check for memory pressure
    if [ -r /proc/pressure/memory ]; then
        ts "Memory pressure indicators:"
        cat /proc/pressure/memory | sed 's/^/  /' | ts_block
        echo ""
    fi
    
    # Performance recommendations
    ts "Performance Notes:"
    case $MATMUL_METHOD in
        0)
            ts "  ℹ️  You used sequential method (baseline)"
            ts "  💡 Try method 3 for ~8-10× speedup: add '0 3 32' to command"
            ;;
        1)
            ts "  ℹ️  You used OpenMP parallelization (~5× speedup)"
            ts "  💡 Try method 3 for additional ~70% improvement"
            ;;
        2)
            ts "  ℹ️  You used OpenMP + Cache blocking (~6× speedup)"
            ts "  💡 Try method 3 (add SIMD) for additional ~40% improvement"
            ;;
        3)
            ts "  ✅ You used the fastest method (~8-10× speedup)"
            if [ "$BLOCK_SIZE" != "32" ]; then
                ts "  💡 Consider testing block size 32 for optimal cache performance"
            fi
            ;;
    esac
    echo ""
}

# ===============================================================================
# RUN VALIDATION TEST (OPTIONAL)
# ===============================================================================

run_validation_test() {
    if [ "$RUN_VALIDATION" = "true" ]; then
        print_separator "VALIDATION TEST"
        
        ts "Running quick validation test..."
        ts "Testing all 4 matrix multiplication methods..."
        echo ""
        
        local test_cmd_base="./catcharea $ARG1 $ARG2 $ARG3 $INVERSION_METHOD"
        
        for method in 0 1 2 3; do
            ts "Testing method $method: ${MATMUL_NAMES[$method]}"
            
            # Run test
            local test_output="validation_method${method}.txt"
            if $test_cmd_base $method $BLOCK_SIZE > "$test_output" 2>&1; then
                # Extract catchment area
                local area=$(grep "Catchment area" "$test_output" | awk '{print $NF}')
                ts "  Result: Catchment area = $area"
            else
                ts "  Error: Test failed"
            fi
        done
        
        echo ""
        ts "Validation complete - check validation_method*.txt for details"
        echo ""
    fi
}

# ===============================================================================
# MAIN EXECUTION FLOW
# ===============================================================================

main() {
    # Print header
    echo ""
    echo "╔═══════════════════════════════════════════════════════════════════════════╗"
    echo "║                                                                           ║"
    echo "║              CATCHAREA OPTIMIZED BUILD AND RUN SCRIPT                    ║"
    echo "║              Enhanced with Matrix Optimization Support                   ║"
    echo "║                                                                           ║"
    echo "╚═══════════════════════════════════════════════════════════════════════════╝"
    echo ""
    
    # Gather system information
    gather_system_info
    
    # Display system information
    print_separator "SYSTEM INFORMATION"
    ts "Operating System:"
    ts "  $OS_NAME ($OS_KERNEL)"
    echo ""
    ts "CPU Information:"
    ts "  Model: $CPU_MODEL"
    ts "  Physical cores: $PHYS_CORES"
    ts "  Logical cores: $LOGICAL_CORES"
    ts "  AVX2 support: $AVX2_SUPPORT"
    ts "  FMA support: $FMA_SUPPORT"
    echo ""
    ts "Memory:"
    ts "  Total: $MEM_TOTAL"
    ts "  Available: $MEM_AVAIL"
    echo ""
    
    if [ -n "$CACHE_INFO" ]; then
        ts "Cache Hierarchy:"
        echo "$CACHE_INFO" | sed 's/^/  /' | ts_block
        echo ""
    fi
    
    ts "NUMA Configuration:"
    echo "$NUMA_INFO" | sed 's/^/  /' | ts_block
    
    # Check build requirements
    check_build_requirements
    
    # Setup environment
    setup_environment
    
    # Display configuration
    display_config_summary
    
    # Pause before build
    read -r -n1 -s -p "$(ts 'Press any key to start build, or Ctrl+C to cancel...')" _
    echo ""
    
    # Clean and build
    clean_build
    build_project
    
    # Run program
    run_program
    
    # Display summary
    display_performance_summary
    
    # Optional validation
    # Uncomment next line to enable validation test
    # RUN_VALIDATION=true run_validation_test
    
    # Final message
    print_separator "SCRIPT COMPLETED SUCCESSFULLY"
    ts "All operations completed successfully! ✅"
    ts "Check output files: test.out, catchment.out"
    echo ""
}

# ===============================================================================
# SCRIPT ENTRY POINT
# ===============================================================================

# Run main function
main

exit 0
