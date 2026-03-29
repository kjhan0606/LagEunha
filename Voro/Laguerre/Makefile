# Makefile for Laguerre Diagram Library
# ======================================

CC = gcc
CFLAGS = -Wall -Wextra -O3 -ffast-math
LDFLAGS = -lm

NVCC = nvcc
# CUDA Architecture settings:
#   sm_89 : RTX 40xx series, RTX 5000 Ada (Ada Lovelace)
#   sm_86 : RTX 30xx series (Ampere)
#   sm_75 : RTX 20xx series (Turing)
#   sm_70 : V100 (Volta)
#   sm_60 : P100, GTX 10xx (Pascal)
# Override with: make cuda CUDA_ARCH=sm_86
CUDA_ARCH ?= sm_89
NVCCFLAGS = -O3 -arch=$(CUDA_ARCH) --use_fast_math

.PHONY: all cpu cuda clean help

all: cpu

cpu: laguerre.o
	@echo "CPU library built successfully"

cuda: laguerre_cuda.o laguerre.o
	@echo "CUDA library built for $(CUDA_ARCH)"

laguerre.o: laguerre.c laguerre.h
	$(CC) $(CFLAGS) -c laguerre.c -o laguerre.o

laguerre_cuda.o: laguerre_cuda.cu laguerre_cuda.cuh
	$(NVCC) $(NVCCFLAGS) -c laguerre_cuda.cu -o laguerre_cuda.o

example_cpu: example_cpu.c laguerre.o
	$(CC) $(CFLAGS) example_cpu.c laguerre.o -o example_cpu $(LDFLAGS)

example_cuda: example_cuda.cu laguerre_cuda.o
	$(NVCC) $(NVCCFLAGS) example_cuda.cu laguerre_cuda.o -o example_cuda $(LDFLAGS)

run_cpu: example_cpu
	./example_cpu

run_cuda: example_cuda
	./example_cuda

clean:
	rm -f *.o example_cpu example_cuda

help:
	@echo "Targets:"
	@echo "  cpu          - Build CPU library"
	@echo "  cuda         - Build CUDA library"
	@echo "  example_cpu  - Build CPU example"
	@echo "  example_cuda - Build CUDA example"
	@echo "  run_cpu      - Build and run CPU example"
	@echo "  run_cuda     - Build and run CUDA example"
	@echo "  clean        - Remove build files"
	@echo ""
	@echo "CUDA Architecture (default: sm_89 for RTX 5000 Ada):"
	@echo "  make cuda CUDA_ARCH=sm_86   # For RTX 30xx"
	@echo "  make cuda CUDA_ARCH=sm_75   # For RTX 20xx"
