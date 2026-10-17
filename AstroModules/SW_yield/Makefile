# ============================================================================
# STELLAR WIND YIELD LIBRARY - Makefile
# ============================================================================
#
# Builds:
#   - libstellar_wind.a  : Static library
#   - libstellar_wind.so : Shared library (optional)
#   - test_stellar_wind  : Test/demonstration program
#
# Usage:
#   make            - Build all targets
#   make lib        - Build static library only
#   make shared     - Build shared library
#   make test       - Build and run test program
#   make clean      - Remove build artifacts
#   make install    - Install to PREFIX (default: /usr/local)
#
# ============================================================================

# Compiler and flags
CC = gcc
CFLAGS = -Wall -Wextra -O3 -fPIC -std=c99
LDFLAGS = -lm
DEBUG_FLAGS = -g -DDEBUG

# Library name
LIB_NAME = stellar_wind
LIB_STATIC = lib$(LIB_NAME).a
LIB_SHARED = lib$(LIB_NAME).so

# Source files
LIB_SRC = stellar_wind_lib.c
LIB_HDR = stellar_wind_lib.h
LIB_OBJ = stellar_wind_lib.o

TEST_SRC = test_stellar_wind.c
TEST_BIN = test_stellar_wind

# Installation directories
PREFIX ?= /usr/local
LIBDIR = $(PREFIX)/lib
INCLUDEDIR = $(PREFIX)/include

# Detect OS for shared library extension
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
    LIB_SHARED = lib$(LIB_NAME).dylib
    SHARED_FLAGS = -dynamiclib
else
    LIB_SHARED = lib$(LIB_NAME).so
    SHARED_FLAGS = -shared
endif

# ============================================================================
# DEFAULT TARGET
# ============================================================================

.PHONY: all
all: $(LIB_STATIC) $(TEST_BIN)

# ============================================================================
# LIBRARY TARGETS
# ============================================================================

# Static library
.PHONY: lib
lib: $(LIB_STATIC)

$(LIB_STATIC): $(LIB_OBJ)
	@echo "Creating static library: $@"
	ar rcs $@ $^

$(LIB_OBJ): $(LIB_SRC) $(LIB_HDR)
	@echo "Compiling: $<"
	$(CC) $(CFLAGS) -c $< -o $@

# Shared library
.PHONY: shared
shared: $(LIB_SHARED)

$(LIB_SHARED): $(LIB_SRC) $(LIB_HDR)
	@echo "Creating shared library: $@"
	$(CC) $(CFLAGS) $(SHARED_FLAGS) $< -o $@ $(LDFLAGS)

# ============================================================================
# TEST PROGRAM
# ============================================================================

$(TEST_BIN): $(TEST_SRC) $(LIB_STATIC)
	@echo "Building test program: $@"
	$(CC) $(CFLAGS) $< -L. -l$(LIB_NAME) -o $@ $(LDFLAGS)

.PHONY: test
test: $(TEST_BIN)
	@echo ""
	@echo "Running test program..."
	@echo ""
	./$(TEST_BIN)

.PHONY: test-quick
test-quick: $(TEST_BIN)
	./$(TEST_BIN) 2 0.02 1e5 1e7

# ============================================================================
# DEBUG BUILD
# ============================================================================

.PHONY: debug
debug: CFLAGS += $(DEBUG_FLAGS)
debug: clean all

# ============================================================================
# INSTALLATION
# ============================================================================

.PHONY: install
install: $(LIB_STATIC) $(LIB_SHARED)
	@echo "Installing to $(PREFIX)..."
	mkdir -p $(LIBDIR)
	mkdir -p $(INCLUDEDIR)
	cp $(LIB_STATIC) $(LIBDIR)/
	cp $(LIB_SHARED) $(LIBDIR)/ 2>/dev/null || true
	cp $(LIB_HDR) $(INCLUDEDIR)/
	@echo "Installation complete."
	@echo "  Headers: $(INCLUDEDIR)/$(LIB_HDR)"
	@echo "  Library: $(LIBDIR)/$(LIB_STATIC)"

.PHONY: uninstall
uninstall:
	rm -f $(LIBDIR)/$(LIB_STATIC)
	rm -f $(LIBDIR)/$(LIB_SHARED)
	rm -f $(INCLUDEDIR)/$(LIB_HDR)

# ============================================================================
# DOCUMENTATION
# ============================================================================

.PHONY: docs
docs:
	@echo "Generating documentation..."
	@if command -v doxygen > /dev/null; then \
		doxygen Doxyfile 2>/dev/null || echo "Create Doxyfile for documentation"; \
	else \
		echo "Doxygen not found. Install with: apt-get install doxygen"; \
	fi

# ============================================================================
# CLEANUP
# ============================================================================

.PHONY: clean
clean:
	@echo "Cleaning build artifacts..."
	rm -f $(LIB_OBJ) $(LIB_STATIC) $(LIB_SHARED) $(TEST_BIN)
	rm -f *.csv
	rm -f *.o

.PHONY: distclean
distclean: clean
	rm -rf docs/

# ============================================================================
# HELP
# ============================================================================

.PHONY: help
help:
	@echo "Stellar Wind Yield Library - Build System"
	@echo ""
	@echo "Targets:"
	@echo "  all        - Build static library and test program (default)"
	@echo "  lib        - Build static library only"
	@echo "  shared     - Build shared library"
	@echo "  test       - Build and run test program"
	@echo "  test-quick - Run quick test with default parameters"
	@echo "  debug      - Build with debug symbols"
	@echo "  install    - Install library and headers to PREFIX"
	@echo "  uninstall  - Remove installed files"
	@echo "  clean      - Remove build artifacts"
	@echo "  help       - Show this help message"
	@echo ""
	@echo "Variables:"
	@echo "  CC         - C compiler (default: gcc)"
	@echo "  CFLAGS     - Compiler flags"
	@echo "  PREFIX     - Installation prefix (default: /usr/local)"
	@echo ""
	@echo "Example:"
	@echo "  make                    # Build everything"
	@echo "  make test               # Build and run tests"
	@echo "  make PREFIX=/opt install  # Install to /opt"
