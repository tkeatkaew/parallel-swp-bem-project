# ═══════════════════════════════════════════════════════════════════════════════
# Makefile for BEM Flow Path Analysis - Unified Version 3.0
# 
# Supports both Hybrid (custom OpenMP+Cache+SIMD) and OpenBLAS implementations
# selected at runtime via dgemm_type parameter
# ═══════════════════════════════════════════════════════════════════════════════

# OpenBLAS configuration
OPENBLAS_ROOT ?= /opt/OpenBLAS
OPENBLAS_INC  := -I$(OPENBLAS_ROOT)/include
OPENBLAS_LIBS := -L$(OPENBLAS_ROOT)/lib -lopenblas -lpthread -lgfortran

# Compiler and flags
CC = gcc
CFLAGS_BASE = -no-pie -g -Wall
CFLAGS_OPT = -O3 -march=native -fopenmp -mavx2 -mfma

# Directories
VPATH = ../header:../source:../object
SRC_DIR = ../source
OBJ_DIR = ../object
HDR_DIR = ../header

#------------------------------------------------------------
# Default target
#------------------------------------------------------------
all: header object catcharea

#------------------------------------------------------------
# Main executable
#------------------------------------------------------------
catcharea: $(OBJ_DIR)/catcharea.o $(OBJ_DIR)/trapfloat.o $(OBJ_DIR)/catchment.o \
           $(OBJ_DIR)/file.o $(OBJ_DIR)/matrix_inv.o $(OBJ_DIR)/path_list.o \
           $(OBJ_DIR)/path.o $(OBJ_DIR)/boundary.o $(OBJ_DIR)/geometry.o \
           $(OBJ_DIR)/memory.o $(OBJ_DIR)/matrix.o $(OBJ_DIR)/matrix_multiply_optimized.o \
           $(OBJ_DIR)/co_matrix.o $(OBJ_DIR)/ten_matrix.o $(OBJ_DIR)/scan.o \
           $(OBJ_DIR)/vcalc.o $(OBJ_DIR)/bsolve.o $(OBJ_DIR)/performance_summary.o \
           $(OBJ_DIR)/linear_sys.o $(OBJ_DIR)/terms.o $(OBJ_DIR)/streamline.o \
           $(OBJ_DIR)/area.o
	$(CC) $(CFLAGS_BASE) -o $@ $^ \
	   -lm $(OPENBLAS_INC) $(OPENBLAS_LIBS) \
	   -fopenmp -march=native -O3 -mavx2 -mfma

#------------------------------------------------------------
# Object files
#------------------------------------------------------------
object: $(OBJ_DIR)/file.o $(OBJ_DIR)/matrix_inv.o $(OBJ_DIR)/path.o \
        $(OBJ_DIR)/path_list.o $(OBJ_DIR)/geometry.o $(OBJ_DIR)/boundary.o \
        $(OBJ_DIR)/catchment.o $(OBJ_DIR)/matrix.o $(OBJ_DIR)/matrix_multiply_optimized.o \
        $(OBJ_DIR)/co_matrix.o $(OBJ_DIR)/ten_matrix.o $(OBJ_DIR)/terms.o \
        $(OBJ_DIR)/linear_sys.o $(OBJ_DIR)/bsolve.o $(OBJ_DIR)/scan.o \
        $(OBJ_DIR)/performance_summary.o $(OBJ_DIR)/vcalc.o $(OBJ_DIR)/streamline.o \
        $(OBJ_DIR)/memory.o $(OBJ_DIR)/area.o $(OBJ_DIR)/trapfloat.o \
        $(OBJ_DIR)/catcharea.o

#------------------------------------------------------------
# Core optimized components (with OpenMP, AVX2, OpenBLAS)
#------------------------------------------------------------

# Main program - now with dgemm_type parameter
$(OBJ_DIR)/catcharea.o: $(SRC_DIR)/catcharea.c catcharea.h \
                        boundary_types.h co_matrix_types.h matrix_types.h \
                        ten_matrix_types.h memory_types.h \
                        area.h catchment.h file.h memory.h path.h \
                        scan.h streamline.h trapfloat.h
	$(CC) $(CFLAGS_BASE) -fopenmp -I $(HDR_DIR) -c $(SRC_DIR)/catcharea.c -o $@

# UNIFIED matrix multiply (supports both Hybrid and OpenBLAS)
$(OBJ_DIR)/matrix_multiply_optimized.o: $(SRC_DIR)/matrix_multiply_optimized.c \
                                        matrix_multiply_optimized.h matrix_types.h
	$(CC) $(CFLAGS_BASE) $(CFLAGS_OPT) $(OPENBLAS_INC) \
	   -I $(HDR_DIR) -c $(SRC_DIR)/matrix_multiply_optimized.c -o $@

# Matrix operations wrapper
$(OBJ_DIR)/matrix.o: $(SRC_DIR)/matrix.c matrix.h matrix_types.h file.h \
                     matrix_multiply_optimized.h
	$(CC) $(CFLAGS_BASE) $(CFLAGS_OPT) -I $(HDR_DIR) -c $(SRC_DIR)/matrix.c -o $@

# Matrix inversion (LAPACK)
$(OBJ_DIR)/matrix_inv.o: $(SRC_DIR)/matrix_inv.c matrix_inv.h
	$(CC) $(CFLAGS_BASE) $(CFLAGS_OPT) $(OPENBLAS_INC) \
	   -I $(HDR_DIR) -c $(SRC_DIR)/matrix_inv.c -o $@

# Performance tracking
$(OBJ_DIR)/performance_summary.o: $(SRC_DIR)/performance_summary.c performance_summary.h
	$(CC) $(CFLAGS_BASE) -I $(HDR_DIR) -c $(SRC_DIR)/performance_summary.c -o $@

#------------------------------------------------------------
# Standard components (no special optimization needed)
#------------------------------------------------------------

$(OBJ_DIR)/file.o: $(SRC_DIR)/file.c file.h
	$(CC) $(CFLAGS_BASE) -I $(HDR_DIR) -c $(SRC_DIR)/file.c -o $@

$(OBJ_DIR)/path.o: $(SRC_DIR)/path.c path.h boundary_types.h file.h
	$(CC) $(CFLAGS_BASE) -I $(HDR_DIR) -c $(SRC_DIR)/path.c -o $@

$(OBJ_DIR)/path_list.o: $(SRC_DIR)/path_list.c path_list.h boundary_types.h path.h
	$(CC) $(CFLAGS_BASE) -I $(HDR_DIR) -c $(SRC_DIR)/path_list.c -o $@

$(OBJ_DIR)/geometry.o: $(SRC_DIR)/geometry.c geometry.h boundary_types.h
	$(CC) $(CFLAGS_BASE) -I $(HDR_DIR) -c $(SRC_DIR)/geometry.c -o $@

$(OBJ_DIR)/boundary.o: $(SRC_DIR)/boundary.c boundary.h boundary_types.h file.h path.h
	$(CC) $(CFLAGS_BASE) -I $(HDR_DIR) -c $(SRC_DIR)/boundary.c -o $@

$(OBJ_DIR)/catchment.o: $(SRC_DIR)/catchment.c catchment.h boundary_types.h \
                        boundary.h file.h geometry.h path.h path_list.h
	$(CC) $(CFLAGS_BASE) -I $(HDR_DIR) -c $(SRC_DIR)/catchment.c -o $@

$(OBJ_DIR)/co_matrix.o: $(SRC_DIR)/co_matrix.c co_matrix.h boundary_types.h \
                        matrix_types.h co_matrix_types.h matrix.h
	$(CC) $(CFLAGS_BASE) -I $(HDR_DIR) -c $(SRC_DIR)/co_matrix.c -o $@

$(OBJ_DIR)/ten_matrix.o: $(SRC_DIR)/ten_matrix.c ten_matrix.h boundary_types.h \
                         matrix_types.h ten_matrix_types.h matrix.h
	$(CC) $(CFLAGS_BASE) -I $(HDR_DIR) -c $(SRC_DIR)/ten_matrix.c -o $@

$(OBJ_DIR)/terms.o: $(SRC_DIR)/terms.c terms.h boundary_types.h geometry.h
	$(CC) $(CFLAGS_BASE) -I $(HDR_DIR) -c $(SRC_DIR)/terms.c -o $@

$(OBJ_DIR)/linear_sys.o: $(SRC_DIR)/linear_sys.c linear_sys.h boundary_types.h \
                         co_matrix_types.h matrix_types.h ten_matrix_types.h \
                         co_matrix.h geometry.h matrix.h path.h ten_matrix.h terms.h
	$(CC) $(CFLAGS_BASE) -I $(HDR_DIR) -c $(SRC_DIR)/linear_sys.c -o $@

$(OBJ_DIR)/bsolve.o: $(SRC_DIR)/bsolve.c bsolve.h boundary_types.h \
                     co_matrix_types.h matrix_types.h ten_matrix_types.h \
                     co_matrix.h linear_sys.h matrix.h path.h ten_matrix.h
	$(CC) $(CFLAGS_BASE) -I $(HDR_DIR) -c $(SRC_DIR)/bsolve.c -o $@

$(OBJ_DIR)/scan.o: $(SRC_DIR)/scan.c scan.h boundary_types.h co_matrix_types.h \
                   matrix_types.h ten_matrix_types.h memory_types.h
	$(CC) $(CFLAGS_BASE) -I $(HDR_DIR) -c $(SRC_DIR)/scan.c -o $@

$(OBJ_DIR)/vcalc.o: $(SRC_DIR)/vcalc.c vcalc.h boundary_types.h co_matrix_types.h \
                    matrix_types.h ten_matrix_types.h memory_types.h \
                    bsolve.h catchment.h co_matrix.h matrix.h path.h ten_matrix.h
	$(CC) $(CFLAGS_BASE) -I $(HDR_DIR) -c $(SRC_DIR)/vcalc.c -o $@

$(OBJ_DIR)/streamline.o: $(SRC_DIR)/streamline.c streamline.h boundary_types.h \
                         co_matrix_types.h matrix_types.h ten_matrix_types.h \
                         memory_types.h catchment.h file.h path.h vcalc.h
	$(CC) $(CFLAGS_BASE) -I $(HDR_DIR) -c $(SRC_DIR)/streamline.c -o $@

$(OBJ_DIR)/memory.o: $(SRC_DIR)/memory.c memory.h
	$(CC) $(CFLAGS_BASE) -I $(HDR_DIR) -c $(SRC_DIR)/memory.c -o $@

$(OBJ_DIR)/area.o: $(SRC_DIR)/area.c area.h boundary_types.h matrix_types.h \
                   co_matrix_types.h ten_matrix_types.h memory_types.h \
                   scan.h streamline.h
	$(CC) $(CFLAGS_BASE) -I $(HDR_DIR) -c $(SRC_DIR)/area.c -o $@

$(OBJ_DIR)/trapfloat.o: $(SRC_DIR)/trapfloat.c trapfloat.h
	$(CC) $(CFLAGS_BASE) -I $(HDR_DIR) -c $(SRC_DIR)/trapfloat.c -o $@

#------------------------------------------------------------
# Header file generation (using cproto)
#------------------------------------------------------------
header: file.h path.h path_list.h geometry.h boundary.h catchment.h \
        matrix.h co_matrix.h ten_matrix.h terms.h linear_sys.h bsolve.h \
        scan.h vcalc.h streamline.h memory.h area.h trapfloat.h

# Individual header rules (if cproto is available)
%.h: $(SRC_DIR)/%.c
	@if command -v cproto > /dev/null 2>&1; then \
	    cproto -I $(HDR_DIR) $< > $(HDR_DIR)/$@; \
	else \
	    echo "cproto not found, skipping $@"; \
	fi

catcharea.h: $(SRC_DIR)/catcharea.c
	@if command -v cproto > /dev/null 2>&1; then \
	    cproto -I $(HDR_DIR) $< > $@; \
	else \
	    echo "cproto not found, skipping $@"; \
	fi

#------------------------------------------------------------
# Clean targets
#------------------------------------------------------------
.PHONY: clean distclean rebuild help

clean:
	@rm -f $(OBJ_DIR)/*.o *.o *.d *.out core core.* *~ .*.swp
	@echo "Cleaned object files"

distclean: clean
	@rm -f catcharea performance_results.csv
	@rm -rf benchmark_results_* benchmark_openblas_*
	@echo "Cleaned all build artifacts"

rebuild: distclean all

#------------------------------------------------------------
# Help
#------------------------------------------------------------
help:
	@echo ""
	@echo "BEM Flow Path Analysis - Makefile"
	@echo "=================================="
	@echo ""
	@echo "Targets:"
	@echo "  all       - Build everything (default)"
	@echo "  catcharea - Build main executable"
	@echo "  clean     - Remove object files"
	@echo "  distclean - Remove all build artifacts"
	@echo "  rebuild   - Clean and rebuild"
	@echo "  help      - Show this message"
	@echo ""
	@echo "Usage:"
	@echo "  make                    # Build with default settings"
	@echo "  make OPENBLAS_ROOT=/path/to/openblas  # Custom OpenBLAS location"
	@echo ""
	@echo "After building:"
	@echo "  ./catcharea -h          # Show program usage"
	@echo "  ./catcharea 1.0 100.0 0.001 0 3 64 1   # OpenBLAS optimized"
	@echo "  ./catcharea 1.0 100.0 0.001 0 3 32 0   # Hybrid (custom SIMD)"
	@echo ""
