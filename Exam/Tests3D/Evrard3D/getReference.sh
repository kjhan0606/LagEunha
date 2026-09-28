#!/bin/bash
# Fetch the published t = 0.8 Evrard profile computed with HydroCode1D
# (B. Vandenbroucke; 2000 cells, exact Riemann solver, 2nd order, 3D
# spherical, self-gravity). It is the reference file the SWIFT code uses in
# examples/HydroTests/EvrardCollapse_3D/getReference.sh. Not committed here
# (third-party data); analyze.py uses it if ref/evrardCollapse3D_exact.txt exists.
set -e
cd "$(dirname "$0")/ref"
URL=https://virgodb.cosma.dur.ac.uk/swift-webstorage/ReferenceSolutions/evrardCollapse3D_exact.txt
if command -v curl >/dev/null; then curl -fsSL -o evrardCollapse3D_exact.txt "$URL"
else wget -q -O evrardCollapse3D_exact.txt "$URL"; fi
head -3 evrardCollapse3D_exact.txt
