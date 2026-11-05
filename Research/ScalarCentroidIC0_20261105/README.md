# Common physical version8 and IC(0) evidence

See [the report](../../Docs/UNIFIED_HYDRO_SCALAR_CENTROID_IC0_2026-11-05.md)
for the physical and numerical changes, measurements, and remaining scope.

JSON receipts are plain text. Other retained files are lossless gzip.
The manifest records original paths, byte counts, and SHA-256 hashes.
ELF executables and Python caches are excluded. Build commands, binary
hashes, frozen sources, and accepted/rejected states remain available.

The 52 fixtures and frozen 2D/3D systems are operator checks. Collision
reaches the original t=0.035 target. Blast fails at t=0.03501883433 before
t=0.038. Both continuations preserve physical version8 and resume accepted
checkpoints, with a block-PCG prefix and IC(0)-PCG continuation.
Diagnostic trial checkpoints are not production restart points.

No equal-accuracy speedup is inferred from failed 1D block64 controls.
No GIZMO runtime, production MPI/CUDA, full 2D suite, or manuscript promotion
is established by this archive.
