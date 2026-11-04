# Face-angular fit and actual native 3D evolution

Methods, measured scope, failures, corrections, and remaining work are in
`../../Docs/UNIFIED_HYDRO_FACE_FIT_3D_2026-11-04.md`.
The common face weight is A/r^(dimension-1), an angular proxy rather than an
exact solid angle. Some historical build/run directory names use solid_angle
and are retained to preserve the original provenance.

The original C/OpenMP+MPI+CUDA executable integrates the hydrodynamics.
The 2D branch uses CUDA geometry/pressure and the common OpenMP graph solve.
The 3D branch uses the original CPU/OpenMP linked-cell Voro3D and the same
common operator. Python scripts only orchestrate runs and analyse outputs.
This archive does not qualify GPU 3D or multi-rank evolution.

The compact evidence includes the 2D stage-by-stage causal trace, frozen-fit
comparison, 1D nine-case unchanged regression, 14 current 2D endpoints, actual
3D failure captures, geometry replay, finite differences, endpoint rejection
fixtures, 15 current 3D endpoints, and two timestep-refined 3D shock endpoints.
Historical nonzero failures and loader/build failures are preserved. The first
3D FD range did not meet its unchanged threshold. Two finer perturbations
resolved the derivative limit; the final strict-entry audit uses all six sizes.

The captured 3D geometry inputs are small enough to include directly. Large
native trajectories, raw per-particle traces, objects, archives, and executables
remain on GPFS and are inventoried by path and hash in manifest.json.
Site-specific hybrid build commands do not establish a portable clean rebuild.
Frozen source copies match the executed build manifests. Current source and
historical binaries are distinguished in the inventory.

The stage is incomplete for the full physics suite, new-method long Kepler,
and the matched GIZMO <=2x wall-time target. The manuscript retains the
previous baseline. The next batch is Gresho/KH and geometry/state reuse tests,
followed by stronger shocks and separate force/boundary contracts.
