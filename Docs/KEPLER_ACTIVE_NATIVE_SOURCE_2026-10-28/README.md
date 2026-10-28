# Native Kepler source preservation, 2026-10-28

Reconstruct from LagEunha commit 472e7fcc3419b060d8fa0c74f0961f6c1f9f188b, overlay baseline/, then active_clip/ for the accepted clip-swap source. The baseline is the measured CUDA+MPI+OpenMP build; fixed mass and W=0 are retained. MANIFEST.json records exact source hashes and binary provenance.

The padding_reuse/ wrapper adds only the opt-in GFS_REUSE_TESS_PADDING switch and preserves post-gradient halo order/ID checks. It remains off by default and requires its short exact-operator gate before enablement; this preservation is not a claim of t>6 or GIZMO t=600 stability. Native checkpoints remain on GPFS.
