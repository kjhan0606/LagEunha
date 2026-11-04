# Bounded native implicit gate and causal audits

The methods and limitations are documented in
`Docs/UNIFIED_HYDRO_NATIVE_IMPLICIT_GATE_2026-11-04.md`.
This is an experimental candidate, not a completed multidimensional solver.
The existing OpenMP/MPI/CUDA executable completed a bounded 2D Kepler gate.
The persistent numerical heating remains unresolved. Long Kepler and matched
GIZMO performance qualification have not passed. Periodic evolution, full 2D
benchmarks, 3D evolution, and MPI qualification remain outstanding.

The exact-gradient and area-fit alternatives occur only in frozen-state
causal audit programs. The production candidate retains the original common
LSQ fit and constitutive law. Compact results, complete logs, and orchestration
scripts are archived here. These Python scripts do not integrate the PDE.
Large native checkpoints, RHS binary inputs, linked objects, and executables
remain on GPFS at the paths recorded in the reports. Site paths reproduce
actual runs and are not a portable full-source rebuild recipe.

The source manifest in the build directory identifies the frozen v5 native
sources. The publication inventory `manifest.json` records the current source
and compact result files. Historical baseline code remains at a1b5247.
