#!/usr/bin/env bash
# Enhanced build and run script for catcharea with memory optimization
# Usage: ./run_catcharea_optimized.sh [NUM_THREADS] [ARG1 ARG2 ARG3 [INVERSION_METHOD]]
#   INVERSION_METHOD: 0=Parallel (default), 1=Sequential
set -u

# ---- config / args ----
NUM_THREADS="${1:-8}"
ARG1="${2:-1.0}"
ARG2="${3:-100.0}"
ARG3="${4:-0.001}"
INVERSION_METHOD="${5:-0}"  # NEW: 0=Parallel (default), 1=Sequential
CMD="./catcharea $ARG1 $ARG2 $ARG3 $INVERSION_METHOD"

# ---- helpers ----
ts() { printf '[%(%Y-%m-%d %H:%M:%S)T] %s\n' -1 "$*"; }
ts_block() {
  while IFS= read -r line; do ts "$line"; done
}

val_or_na() {
  if [ -n "${1:-}" ]; then echo "$1"; else echo "N/A"; fi
}

# ---- gather system info (robust fallbacks) ----
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

CPU_MODEL="$(lscpu 2>/dev/null | awk -F: '/Model name/ {sub(/^[ \t]+/, "", $2); print $2; exit}')"
if [ -z "$CPU_MODEL" ] && [ -r /proc/cpuinfo ]; then
  CPU_MODEL="$(awk -F: '/model name/ {print $2; exit}' /proc/cpuinfo | sed 's/^ //')"
fi
CPU_MODEL="$(val_or_na "$CPU_MODEL")"

# physical cores (unique Core,Socket pairs)
if command -v lscpu >/dev/null 2>&1; then
  PHYS_CORES="$(lscpu -p=Core,Socket 2>/dev/null | awk -F, '!/^#/ {k=$1","$2; c[k]=1} END {print length(c)}')"
else
  PHYS_CORES=""
fi
PHYS_CORES="$(val_or_na "$PHYS_CORES")"

# memory total (human readable)
if command -v free >/dev/null 2>&1; then
  MEM_TOTAL="$(free -h | awk '/^Mem:/ {print $2}')"
elif [ -r /proc/meminfo ]; then
  kb=$(awk '/MemTotal:/ {print $2}' /proc/meminfo)
  if [ -n "$kb" ]; then
    awk -v kb="$kb" 'BEGIN { printf "%.1fGi", kb/1048576 }'
    echo
  else
    MEM_TOTAL=""
  fi
else
  MEM_TOTAL=""
fi
MEM_TOTAL="$(val_or_na "$MEM_TOTAL")"

# cache sizes via lscpu (L1d/L1i/L2/L3)
CACHE_INFO="$(lscpu 2>/dev/null | awk -F: '
/^(L1d cache|L1i cache|L2 cache|L3 cache)/ {
  gsub(/^[ \t]+/, "", $2);
  printf "%-35s %s\n", $1":", $2
}')"

# NUMA information
NUMA_INFO="$(numactl --hardware 2>/dev/null | grep -E '(available|node.*cpus|node.*size)' || echo 'NUMA info not available')"

# ---- print system/environment summary ----
echo ""
echo "================================================================================"
echo "                        SYSTEM INFORMATION SUMMARY"
echo "================================================================================"
echo ""
ts "Operating System:"
ts "  OS: $OS_NAME ($OS_KERNEL)"
echo ""
ts "CPU Information:"
ts "  Model: $CPU_MODEL"
ts "  Physical cores: $PHYS_CORES"
echo ""
ts "Memory:"
ts "  Total Memory: $MEM_TOTAL"
echo ""
if [ -n "$CACHE_INFO" ]; then 
  ts "Cache Hierarchy:"
  echo "$CACHE_INFO" | ts_block
  echo ""
fi

ts "NUMA Configuration:"
echo "$NUMA_INFO" | ts_block
echo ""

ts "Environment file: /etc/environment"
if [ -r /etc/environment ]; then
  sed 's/^/  /' /etc/environment | ts_block
else
  ts "  (not readable by current user)"
fi
echo ""

# ---- OpenMP / OpenBLAS env with OPTIMIZATIONS ----
export OMP_NUM_THREADS="$NUM_THREADS"
export OMP_PROC_BIND=spread
export OMP_PLACES=cores
export OMP_DISPLAY_ENV=VERBOSE
export OPENBLAS_NUM_THREADS="$NUM_THREADS"

# Memory optimization environment variables
# Disable memory sanitizers that can slow down execution
export ASAN_OPTIONS=detect_leaks=0:malloc_context_size=0
export LSAN_OPTIONS=exitcode=0

# Optimize malloc for multithreaded performance
export MALLOC_ARENA_MAX=2
export MALLOC_MMAP_THRESHOLD_=131072
export MALLOC_TRIM_THRESHOLD_=131072
export MALLOC_TOP_PAD_=131072

# For glibc malloc tuning - more arenas can help with parallel allocation
# but use more memory. Adjust based on your needs.
export MALLOC_ARENA_MAX=4

# Transparent Huge Pages can improve performance for large allocations
# Check if we can enable it
if [ -w /sys/kernel/mm/transparent_hugepage/enabled ]; then
  echo "Enabling Transparent Huge Pages for better memory performance..."
  echo always | sudo tee /sys/kernel/mm/transparent_hugepage/enabled > /dev/null 2>&1 || true
fi

echo "================================================================================"
echo "                    PERFORMANCE OPTIMIZATION SETTINGS"
echo "================================================================================"
echo ""
ts "OpenMP / OpenBLAS configuration:"
ts "  OMP_NUM_THREADS=$OMP_NUM_THREADS"
ts "  OMP_PROC_BIND=$OMP_PROC_BIND"
ts "  OMP_PLACES=$OMP_PLACES"
ts "  OPENBLAS_NUM_THREADS=$OPENBLAS_NUM_THREADS"
echo ""

ts "Matrix Inversion Method:"
if [ "$INVERSION_METHOD" = "0" ]; then
  ts "  Mode: PARALLEL (LAPACK - fastest)"
elif [ "$INVERSION_METHOD" = "1" ]; then
  ts "  Mode: SEQUENTIAL (manual Gauss-Jordan - for comparison)"
else
  ts "  Mode: UNKNOWN (will default to PARALLEL)"
fi
echo ""

ts "Memory optimization settings:"
ts "  MALLOC_ARENA_MAX=$MALLOC_ARENA_MAX"
ts "  MALLOC_MMAP_THRESHOLD_=$MALLOC_MMAP_THRESHOLD_"
ts "  MALLOC_TRIM_THRESHOLD_=$MALLOC_TRIM_THRESHOLD_"
ts "  MALLOC_TOP_PAD_=$MALLOC_TOP_PAD_"
echo ""

# CPU governor settings check
if [ -r /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor ]; then
  GOVERNOR=$(cat /sys/devices/system/cpu/cpu0/cpufreq/scaling_governor)
  ts "CPU Governor: $GOVERNOR"
  if [ "$GOVERNOR" != "performance" ]; then
    ts "  (Consider setting to 'performance' for benchmarking)"
  fi
  echo ""
fi

# Check available memory
AVAIL_MEM_KB=$(awk '/MemAvailable/ {print $2}' /proc/meminfo 2>/dev/null || echo "0")
AVAIL_MEM_MB=$((AVAIL_MEM_KB / 1024))
ts "Available memory: ${AVAIL_MEM_MB} MB"
echo ""

# ---- pause before build ----
ts "================================================================================"
read -r -n1 -s -p "$(ts 'Press any key to start build (make all)...')" _; echo
echo ""

# ---- build1 ----
ts "================================================================================"
ts "                              BUILD-1 PHASE"
ts "================================================================================"
ts "Running: make all"
if ! make all 2>&1 | ts_block; then
  ts "Error: build failed!"
  exit 1
fi
echo ""

# ---- build2 ----
ts "================================================================================"
ts "                              BUILD-2 PHASE"
ts "================================================================================"
ts "Running: make all"
if ! make all 2>&1 | ts_block; then
  ts "Error: build failed!"
  exit 1
fi
echo ""

# ---- check binary size ----
if [ -f "./catcharea" ]; then
  BINARY_SIZE=$(du -h ./catcharea | cut -f1)
  ts "Binary size: $BINARY_SIZE"
  echo ""
fi

# ---- run program ----
echo ""
ts "================================================================================"
ts "                           EXECUTION PHASE"
ts "================================================================================"
ts "Command: $CMD"
ts "Start time: $(date '+%Y-%m-%d %H:%M:%S')"
ts "================================================================================"
echo ""

# Save peak memory usage using /usr/bin/time if available
if command -v /usr/bin/time >/dev/null 2>&1; then
  ts "Running with detailed resource tracking..."
  echo ""
  if ! /usr/bin/time -v $CMD 2>&1 | ts_block; then
    ts "Error: catcharea execution failed!"
    exit 1
  fi
else
  ts "Running without /usr/bin/time (install 'time' package for more metrics)..."
  echo ""
  if ! { time $CMD; } 2>&1 | ts_block; then
    ts "Error: catcharea execution failed!"
    exit 1
  fi
fi

echo ""
ts "================================================================================"
ts "                          EXECUTION COMPLETED"
ts "================================================================================"
ts "End time: $(date '+%Y-%m-%d %H:%M:%S')"
ts "================================================================================"
echo ""

# ---- Final memory report ----
echo ""
ts "================================================================================"
ts "                        FINAL SYSTEM STATE"
ts "================================================================================"

# Memory info after run
if command -v free >/dev/null 2>&1; then
  echo ""
  ts "Memory state after execution:"
  free -h | ts_block
fi

# Check for any memory pressure indicators
if [ -r /proc/pressure/memory ]; then
  echo ""
  ts "Memory pressure indicators:"
  cat /proc/pressure/memory | ts_block
fi

echo ""
ts "================================================================================"
ts "Script completed successfully"
ts "================================================================================"
