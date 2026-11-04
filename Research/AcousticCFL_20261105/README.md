# Complete pressure Hessian and opt-in acoustic CFL experiment

See ../../Docs/UNIFIED_HYDRO_ACOUSTIC_CFL_2026-11-05.md for findings, failures, scope, and remaining work.

All fluid evolution uses the native C/OpenMP+MPI+CUDA program. Python only orchestrates, analyses, or rasterizes exact W=0 Voronoi density regions.

The positive-acoustic block bound is tested algebraically in 1D/2D/3D. Its opt-in operational adapter is currently periodic one-rank 2D. The full geometric Hessian is not bounded, and the default timestep has not changed. Gresho velocity error and kinetic-energy loss remain substantial.

Failed probe/startup versions are preserved separately and excluded from successful evolution claims. Current KH cases all have verified t=2 terminal checkpoints. Large outputs and binaries remain on GPFS with SHA256 inventory. Build manifests and frozen sources identify partial relinks, not clean portable builds. The manuscript retains its previous verified baseline.
