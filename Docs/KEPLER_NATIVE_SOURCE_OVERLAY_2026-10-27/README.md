# Native Kepler source overlay, 2026-10-27

The active 2D CUDA+MPI+OpenMP source on `syntax` is an isolated Git worktree
under `/tmp/Kepler_sparse_periodic_area_20261027`. This directory preserves
its 42 modified or untracked `Exam/` source files (1,107,880 bytes) so a
node reboot cannot silently remove the code used by the corrected sparse
IC preflights. [MANIFEST.json](MANIFEST.json) records every file's SHA-256,
the base Git commit `472e7fcc3419b060d8fa0c74f0961f6c1f9f188b`, and
the accepted sparse-area binary SHA-256
`eff735174de1b487759f168d75e1e317a201692ef8d6cadba03cd03e49afb4b5`.
An explicit-file rsync with checksum and a second remote dry-run verified
every preserved file byte-for-byte against the syntax worktree.

To reconstruct the source, check out that base commit of `LagEunha`, then
overlay the `Exam/` files at their relative paths. This is a source receipt;
it does not claim a byte-identical binary rebuild without the compiler,
module, make settings, and linked libraries. It also does not promote the
corrected sparse IC: the native preflight and long-time physics gates still
decide that.
