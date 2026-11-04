# Tested native W0 source at the Kepler 20pi endpoint

At commit a1b52477b61570851d1012c3de82274d1a94e455, the 17 files listed in
published_source_manifest.json are byte-identical to the source overlay
used for job 449508. Later commits add an experimental common hydrodynamic
operator to several of the same paths. Recover the historical files from
the stated commit when reproducing job 449508. The overlay is integrated into the
existing optional full-source target Exam/NativeU64, retaining the parent
OpenMP, MPI, CUDA, and U64 particle ABI support.

The qualified executable SHA256 is
b26c3a1140178c973ef29c1af3b3437a87a60af4ac5ff44e9e6283fe6830e667.
The archived build script records the actual site-specific build and reused
object dependencies. The inherited optional full build remains available
through Makefile.native-u64. A clean rebuild of the integrated publication
checkout has not been executed or assigned the historical binary hash.

Job 449508 completed from t=40.3 to exactly 20pi with fixed mass, W=0,
the variational geometric pressure operator, two-level KDK, and exact
adiabatic volume updates. The trajectory before the restart includes RK4.
This is a one-rank fixed-polygon Kepler result. The KDK implementation is a
correctness prototype that still rebuilds global geometry and does not
claim general active-only individual stepping or multi-rank qualification.

The Kepler mode prohibits artificial viscosity and is not a shock solver.
Do not enable the branch unchanged for every hydrodynamic benchmark.
The current expanded research goal is a single consistent operator for
1D, 2D, and 3D, with new shock, contact, shear, and speed verification.
The measured 3.339x GIZMO wall time exceeds the requested 2x target.

Assessment files preserve source and restart lineage, conservation audits,
density fields, and measured throughput. Large snapshots and executables
remain on GPFS. Paths in archived scripts and manifests describe the actual
historical run and are intentionally not presented as portable paths.

The independent 1D contact audit is in Research/Contact1D_CausalAudit_20261104.
The published manuscript snapshot is in Manuscript/Overleaf_20261104.
The contact audit supersedes claims of successful equilibrium preservation
based solely on the earlier manuscript's contact error table.
