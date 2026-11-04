# Experimental dimension-independent W0 hydrodynamics

This snapshot is an incomplete research implementation. The qualified native
Kepler 20pi solver remains recoverable at commit
`a1b52477b61570851d1012c3de82274d1a94e455`.
The manuscript snapshot in `Manuscript/Overleaf_20261104` reports the previous
qualified calculation. The manuscript does not validate the new solver.

See `Docs/UNIFIED_HYDRO_PROGRESS_2026-11-04.md` for methods, measured results,
rejected candidates, and outstanding requirements. The shared C/CUDA face and
velocity-fit algebra covers dimensions 1, 2, and 3. Complete time evolution
has been tested only in 1D. The new native CUDA adapter has been evaluated
only as an RHS on saved, fully owned, one-rank, fixed-wall Kepler states.
Periodic images, production KDK coupling, 3D evolution, and multi-rank
qualification remain incomplete. No GIZMO speed requirement has been passed.

`results` contains compact original audit outputs. Explicit shock convergence
outputs precede the implicit candidate and must not be assigned to the
implicit method. `implicit_1d_v1` contains the nine completed N=200 runs,
including original logs and final field files. The original physical cell
masses are fixed during evolution. Inflow boundaries can insert new cells.

The implicit 1D update freezes coefficients per dissipation half-step and
solves the endpoint compression signs. BE dissipation is first order in time.
The separate frozen 2D algebra probe has no endpoint sign solve and can
produce attractive impulses after an edge turns expansive. Neither probe
qualifies coupled multidimensional hydrodynamics. The recorded PCG timing
excludes graph construction, device transfer, and hydrodynamic evolution.
The larger Noh energy error and the incorrect old Shu--Osher reference remain
explicit limitations. Shu norms against the old reference are diagnostic only.

Raw native snapshots and executables remain in
`/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1` and the historical campaign
directories recorded by the existing Kepler audit. Site-specific build JSON
files record actual invocations, not a portable clean build. Archived build
paths can refer to earlier object names superseded by the separate CUDA module.

The source and compact-result SHA256 inventory is `manifest.json`.
Use `Exam/HydroCore/README.md` for the small C checks and 1D build commands.
All calculation directories must use `/gpfs`. Do not use `/scratch`.
