# Opt-in native fixed-polygon Kepler source, 2026-10-27

This directory preserves the six source files changed for the fixed 1024-gon
W=0 Kepler branch, including the intermediate-corner preclip correction and
the 64-corner host/device polygon ABI. The current candidate preclips 24
additional nearby real bisectors before the fixed wall. Earlier source
receipts are preserved in `MANIFEST_PRECLIP48.json`, `MANIFEST_CORNER64.json`
and `../KEPLER_FIXEDPOLY_NATIVE_PATCH_2026-10-27.patch`.
Reconstruct from LagEunha Git commit
`472e7fcc3419b060d8fa0c74f0961f6c1f9f188b`, apply the previously
[preserved native source overlay](../KEPLER_NATIVE_SOURCE_OVERLAY_2026-10-27/README.md),
then overlay the six files here at their relative `Exam/` paths.
[MANIFEST.json](MANIFEST.json) records source and candidate binary SHA-256
hashes. The candidate was compiled with Intel MPI, Intel OpenMP, CUDA 13.0.2,
`USE_CUDA=1` and `GFS_VARPRESS_W0=1` for sm80/86/89/90.

The code is an opt-in one-rank research branch. It clips CUDA Voronoi cells to
the fixed polygon and adds the static-wall relative-sweep pressure work.
Generator-wall impacts are deliberately fatal because specular reflection
and its wall impulse/torque ledger have not been validated. The current
Slurm test does not establish long-time Kepler stability.
