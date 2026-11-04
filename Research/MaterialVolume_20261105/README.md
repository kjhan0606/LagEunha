# Material volume consistency and constrained-weight necessary gates

See ../../Docs/UNIFIED_HYDRO_MATERIAL_VOLUME_2026-11-05.md for the mathematical contracts, actual native evidence, limitations, failed attempts, and completed OpenMP stage.

The original C/OpenMP+MPI+CUDA program supplies fluid evolution and GPU geometry. Native C supplies small 1D interval/2D/strict-3D geometry and the bounded coefficient cost screen. NumPy is an independent finite-matrix verifier and analysis tool; it does not integrate a fluid. The exact-circle mapping is a geometry-only control, independently checked against the birth angular velocity.

The candidate's algebra/nonaffine/weighted-geometry necessary gates pass. The current CPU coefficient implementation fails the fixed iteration/tolerance screen even with ACTIVE OpenMP. Production pressure/weights/thermal/time integration remain unchanged, and no candidate PDE or long Kepler trajectory was released. ACTIVE improves the matched original-native short control with bit-identical final fields. GIZMO <=2x wall time remains unverified.

Large native states, per-cell CSVs, executables, and perturbation geometry remain on GPFS with SHA256 inventory. Initial small geometry and compact FD volume samples are archived so the tests can be independently checked. Earlier small-audit source was reconstructed to its recorded hash in verified_source, distinct from a pre-build snapshot. Source snapshots and build manifests distinguish native partial relinks and verifier revisions. Failed technical attempts are preserved rather than counted as physical failures. The manuscript retains the verified earlier baseline.
