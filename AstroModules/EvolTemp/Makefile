# =============================================================================
# Makefile for Exact Cooling Integrator
# Based on Zhu, Smith & Hernquist (2017), MNRAS 470, 1017
# =============================================================================

CC = gcc
CFLAGS = -Wall -O2 -std=c99 -fopenmp
LDFLAGS = -lm

# Uncomment to use GetMu library for mean molecular weight calculation
 CFLAGS += -DUSE_GETMU
 GETMU_SRC = getMU.c

# Source files
SRC = exact_cooling_integrator.c
HDR = exact_cooling_integrator.h getMU.h

# Test program
TEST_TARGET = test_exact_cooling
TEST_CFLAGS = $(CFLAGS) -DTEST_PROGRAM

# Library (if building as static library)
LIB_TARGET = libexact_cooling.a

# =============================================================================
# Build targets
# =============================================================================

.PHONY: all test clean help

all: $(TEST_TARGET)

# Build test program
$(TEST_TARGET): $(SRC) $(HDR)
	$(CC) $(TEST_CFLAGS) -o $@ $(SRC) $(GETMU_SRC) $(LDFLAGS)
	@echo ""
	@echo "Build complete: $(TEST_TARGET)"
	@echo "Usage: ./$(TEST_TARGET) <cooling_table.bin>"
	@echo ""

# Build with GetMu support
with_getmu: CFLAGS += -DUSE_GETMU
with_getmu: GETMU_SRC = GetMu.c
with_getmu: $(TEST_TARGET)

# Build static library
lib: $(SRC)
	$(CC) $(CFLAGS) -c $(SRC) -o exact_cooling_integrator.o
	ar rcs $(LIB_TARGET) exact_cooling_integrator.o
	@echo ""
	@echo "Library built: $(LIB_TARGET)"
	@echo ""

# Run test
test: $(TEST_TARGET)
	@echo "Note: Test requires a cooling table binary file."
	@echo "Generate one using cooling_rate_mpi.c first."
	@echo ""

# Clean
clean:
	rm -f $(TEST_TARGET) $(LIB_TARGET) *.o

# Help
help:
	@echo "================================================================"
	@echo "Exact Cooling Integrator - Makefile Help"
	@echo "================================================================"
	@echo ""
	@echo "Targets:"
	@echo "  all         Build test program (default)"
	@echo "  with_getmu  Build with GetMu library support"
	@echo "  lib         Build as static library"
	@echo "  test        Build and show test instructions"
	@echo "  clean       Remove built files"
	@echo "  help        Show this message"
	@echo ""
	@echo "Usage:"
	@echo "  1. Generate cooling table with cooling_rate_mpi.c"
	@echo "  2. make"
	@echo "  3. ./test_exact_cooling cooling_rate_grid.bin"
	@echo ""
	@echo "================================================================"
