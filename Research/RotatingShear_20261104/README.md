# Geometry reuse and Gresho/KH evidence

See ../../Docs/UNIFIED_HYDRO_ROTATING_SHEAR_2026-11-04.md for scope, results, and unresolved errors.

This is a dated publication capture. A case without result.json is incomplete; its CSV is a partial history. Capturing files does not stop or restart the original calculation.

The native C/OpenMP+MPI+CUDA executable evolves the PDE. Archived Python scripts orchestrate and analyse it. Build commands are site-specific partial relinks, not a portable clean build. Frozen source hashes distinguish the evolution binary from the later read-only acoustic probe.

Large checkpoints and executable binaries remain on GPFS. manifest.json inventories snapshot paths and SHA256 hashes and verifies executed source manifests. The frozen acoustic diagnostic omits the geometric Hessian term. Prior paper results are not reclassified as outcomes of this experimental solver.
