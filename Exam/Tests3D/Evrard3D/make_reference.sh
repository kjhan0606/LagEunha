#!/bin/bash
# Build and run the independent 1D spherical Evrard reference (evrard1d_ref.c).
# Writes ref/evrard1d_N<N>_energy.txt (t Ekin Eth Epot Etot, dt = 0.01) and
# ref/evrard1d_N<N>_prof_t0.800.txt (r rho v P). N=2000 takes ~10 s.
set -e
cd "$(dirname "$0")"
N=${1:-2000}
${CC:-gcc} -O2 -o evrard1d_ref evrard1d_ref.c -lm
./evrard1d_ref "$N" 3.0 "ref/evrard1d_N$N" 0.8
