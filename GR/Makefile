# Z4c DG-FVM Makefile
# Compile: make
# Run: ./bin/z4c_dg -t gauge_wave -n 10 -T 1.0

CC = gcc
CFLAGS = -O3 -Wall -std=c99
LDFLAGS = -lm

all: bin/z4c_dg

bin/z4c_dg: src/z4c_combined.c | bin
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)

bin:
	mkdir -p bin

test: bin/z4c_dg
	@echo "=== Flat spacetime test ==="
	./bin/z4c_dg -t flat -n 5 -T 0.1
	@echo ""
	@echo "=== Gauge wave test ==="
	./bin/z4c_dg -t gauge_wave -n 5 -T 0.2 -A 0.01

clean:
	rm -rf bin output

.PHONY: all test clean
