# Makefile for Stellar Cluster Supernova Library

CC = gcc
CFLAGS = -Wall -O2 -fPIC
CFLAGS_OPENMP = -Wall -O2 -fPIC -fopenmp
LDFLAGS = -lm
LDFLAGS_OPENMP = -lm -fopenmp
AR = ar
ARFLAGS = rcs

# Library name
LIB_NAME = libcluster_supernova
STATIC_LIB = $(LIB_NAME).a
SHARED_LIB = $(LIB_NAME).so

# Source files
LIB_SRC = cluster_supernova.c
LIB_OBJ = $(LIB_SRC:.c=.o)
LIB_HDR = cluster_supernova.h

# Test programs
TEST_BASIC_SRC = test_cluster.c
TEST_BASIC_BIN = test_cluster

TEST_CONTINUOUS_SRC = test_continuous_sf.c
TEST_CONTINUOUS_BIN = test_continuous_sf

EXAMPLE_SRC = example_usage.c
EXAMPLE_BIN = example_usage

TEST_POPIII_SRC = test_popiii.c
TEST_POPIII_BIN = test_popiii

TEST_POPIII_IMF_SRC = test_popiii_imf.c
TEST_POPIII_IMF_BIN = test_popiii_imf

ALL_TESTS = $(TEST_BASIC_BIN) $(TEST_CONTINUOUS_BIN) $(EXAMPLE_BIN) $(TEST_POPIII_BIN) $(TEST_POPIII_IMF_BIN)

# Installation directories
PREFIX = /usr/local
INSTALL_LIB_DIR = $(PREFIX)/lib
INSTALL_INC_DIR = $(PREFIX)/include

# Default target
all: static shared tests

# Build static library
static: $(STATIC_LIB)

$(STATIC_LIB): $(LIB_OBJ)
	@echo "Building static library..."
	$(AR) $(ARFLAGS) $@ $^
	@echo "Static library built: $(STATIC_LIB)"

# Build shared library
shared: $(SHARED_LIB)

$(SHARED_LIB): $(LIB_OBJ)
	@echo "Building shared library..."
	$(CC) -shared -o $@ $^ $(LDFLAGS)
	@echo "Shared library built: $(SHARED_LIB)"

# Build all test programs
tests: $(ALL_TESTS)

# Build basic test program
$(TEST_BASIC_BIN): $(TEST_BASIC_SRC) $(STATIC_LIB)
	@echo "Building basic test program..."
	$(CC) $(CFLAGS) -o $@ $< -L. -lcluster_supernova $(LDFLAGS)
	@echo "Basic test program built: $(TEST_BASIC_BIN)"

# Build continuous SF test program
$(TEST_CONTINUOUS_BIN): $(TEST_CONTINUOUS_SRC) $(STATIC_LIB)
	@echo "Building continuous SF test program..."
	$(CC) $(CFLAGS) -o $@ $< -L. -lcluster_supernova $(LDFLAGS)
	@echo "Continuous SF test program built: $(TEST_CONTINUOUS_BIN)"

# Build example usage program
$(EXAMPLE_BIN): $(EXAMPLE_SRC) $(STATIC_LIB)
	@echo "Building example usage program..."
	$(CC) $(CFLAGS) -o $@ $< -L. -lcluster_supernova $(LDFLAGS)
	@echo "Example usage program built: $(EXAMPLE_BIN)"

# Build Pop III/PISN test program with OpenMP
$(TEST_POPIII_BIN): $(TEST_POPIII_SRC) $(STATIC_LIB)
	@echo "Building Pop III/PISN test program with OpenMP..."
	$(CC) $(CFLAGS_OPENMP) -o $@ $< -L. -lcluster_supernova $(LDFLAGS_OPENMP)
	@echo "Pop III test program built: $(TEST_POPIII_BIN) (OpenMP enabled)"

$(TEST_POPIII_IMF_BIN): $(TEST_POPIII_IMF_SRC) $(STATIC_LIB)
	@echo "Building Pop III/PISN IMF test program with OpenMP..."
	$(CC) $(CFLAGS_OPENMP) -o $@ $< -L. -lcluster_supernova $(LDFLAGS_OPENMP)
	@echo "Pop III test program built: $(TEST_POPIII_IMF_BIN) (OpenMP enabled)"

# Compile object files
%.o: %.c $(LIB_HDR)
	$(CC) $(CFLAGS) -c $< -o $@

# Run all tests
run-tests: $(ALL_TESTS)
	@echo "========================================"
	@echo " Running Basic Tests"
	@echo "========================================"
	LD_LIBRARY_PATH=. ./$(TEST_BASIC_BIN)
	@echo ""
	@echo "========================================"
	@echo " Running Continuous SF Tests"
	@echo "========================================"
	LD_LIBRARY_PATH=. ./$(TEST_CONTINUOUS_BIN)
	@echo ""
	@echo "========================================"
	@echo " Running Usage Examples"
	@echo "========================================"
	LD_LIBRARY_PATH=. ./$(EXAMPLE_BIN)

# Run individual tests
run-basic: $(TEST_BASIC_BIN)
	@echo "Running basic tests..."
	LD_LIBRARY_PATH=. ./$(TEST_BASIC_BIN)

run-continuous: $(TEST_CONTINUOUS_BIN)
	@echo "Running continuous SF tests..."
	LD_LIBRARY_PATH=. ./$(TEST_CONTINUOUS_BIN)

run-examples: $(EXAMPLE_BIN)
	@echo "Running usage examples..."
	LD_LIBRARY_PATH=. ./$(EXAMPLE_BIN)

run-popiii: $(TEST_POPIII_BIN)
	@echo "Running Pop III/PISN tests with OpenMP..."
	OMP_NUM_THREADS=4 LD_LIBRARY_PATH=. ./$(TEST_POPIII_BIN)

run-popiii: $(TEST_POPIII_IMF_BIN)
	@echo "Running Pop III/PISN tests with OpenMP..."
	OMP_NUM_THREADS=4 LD_LIBRARY_PATH=. ./$(TEST_POPIII_IMF_BIN)

# OpenMP thread safety test
openmp-test: $(TEST_POPIII_BIN)
	@echo "========================================" 
	@echo " Running OpenMP Thread Safety Test"
	@echo "========================================"
	@echo "Setting OMP_NUM_THREADS=4..."
	OMP_NUM_THREADS=4 LD_LIBRARY_PATH=. ./$(TEST_POPIII_BIN)

# OpenMP thread safety test
openmp-test: $(TEST_POPIII_IMF_BIN)
	@echo "========================================" 
	@echo " Running OpenMP Thread Safety Test"
	@echo "========================================"
	@echo "Setting OMP_NUM_THREADS=4..."
	OMP_NUM_THREADS=4 LD_LIBRARY_PATH=. ./$(TEST_POPIII_IMF_BIN)

# Install library
install: $(STATIC_LIB) $(SHARED_LIB)
	@echo "Installing library..."
	install -d $(INSTALL_LIB_DIR)
	install -d $(INSTALL_INC_DIR)
	install -m 644 $(STATIC_LIB) $(INSTALL_LIB_DIR)
	install -m 755 $(SHARED_LIB) $(INSTALL_LIB_DIR)
	install -m 644 $(LIB_HDR) $(INSTALL_INC_DIR)
	@echo "Library installed to $(PREFIX)"
	@echo ""
	@echo "To use the library, add to your compiler flags:"
	@echo "  -I$(INSTALL_INC_DIR) -L$(INSTALL_LIB_DIR) -lcluster_supernova -lm"

# Uninstall library
uninstall:
	@echo "Uninstalling library..."
	rm -f $(INSTALL_LIB_DIR)/$(STATIC_LIB)
	rm -f $(INSTALL_LIB_DIR)/$(SHARED_LIB)
	rm -f $(INSTALL_INC_DIR)/$(LIB_HDR)
	@echo "Library uninstalled"

# Clean build files
clean:
	@echo "Cleaning build files..."
	rm -f $(LIB_OBJ) $(STATIC_LIB) $(SHARED_LIB) $(ALL_TESTS)
	@echo "Clean complete"

# Create documentation
doc:
	@echo "Generating documentation..."
	doxygen Doxyfile 2>/dev/null || echo "Doxygen not installed, skipping documentation"

# Quick test (just continuous SF - most important)
quick-test: $(TEST_CONTINUOUS_BIN)
	@echo "Running quick test (Continuous SF)..."
	@LD_LIBRARY_PATH=. ./$(TEST_CONTINUOUS_BIN) | head -80

# Help
help:
	@echo "Stellar Cluster Supernova Library - Build System"
	@echo ""
	@echo "Available targets:"
	@echo "  all              - Build static library, shared library, and all tests (default)"
	@echo "  static           - Build static library (.a)"
	@echo "  shared           - Build shared library (.so)"
	@echo "  tests            - Build all test programs"
	@echo ""
	@echo "Test targets:"
	@echo "  run-tests        - Run all test programs"
	@echo "  run-basic        - Run basic instantaneous burst tests"
	@echo "  run-continuous   - Run continuous star formation tests"
	@echo "  run-examples     - Run usage examples"
	@echo "  run-popiii       - Run Pop III/PISN tests (OpenMP enabled)"
	@echo "  openmp-test      - Run OpenMP thread safety test"
	@echo "  quick-test       - Quick validation (continuous SF only)"
	@echo ""
	@echo "Installation:"
	@echo "  install          - Install library to $(PREFIX)"
	@echo "  uninstall        - Remove installed library"
	@echo ""
	@echo "Maintenance:"
	@echo "  clean            - Remove all build files"
	@echo "  doc              - Generate documentation (requires Doxygen)"
	@echo "  help             - Show this help message"
	@echo ""
	@echo "Quick start:"
	@echo "  make              # Build everything"
	@echo "  make run-tests    # Run all tests"
	@echo "  make openmp-test  # Test OpenMP thread safety"
	@echo "  make quick-test   # Quick validation"

.PHONY: all static shared tests run-tests run-basic run-continuous run-examples \
        run-popiii openmp-test install uninstall clean doc quick-test help
