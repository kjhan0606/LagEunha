# Cell-mean KEP native evolution and remaining stress error

See `Docs/UNIFIED_HYDRO_KEP_NATIVE_2026-11-05.md` and the grouped plan
`Docs/UNIFIED_HYDRO_KEP_ALE_PLAN_2026-11-05.md`.

This archive records actual C/OpenMP evolution, not a Python surrogate:
six contact endpoints, Sod200/400, cold64 and its quarter-step control,
and Gresho64. All eleven endpoints complete; both cold physical pressure
gates fail (1.58% and 1.85%, limit 1%). This is not candidate adoption.

- `kep_native_build_v1`: first native cold-birth RHS check.
- `kep_native_build_v2`: full same-core coupling with actual geometry,
  boundary entropy and kinetic ledgers, before any PDE trajectory.
- `kep_native_group_v1`: frozen-settings full endpoints and exact Sod
  actual-cell assessments. Historical launch PIDs are not current status.
- `kep_smooth_reference_v1` and `kep_native_assessment_*`: read-only
  actual-polygon analytic means and physical node reference checks.
- `kep_default_hll_regression_v1`: retained default HLLC contact/Sod
  trajectories and checkpoints remain byte-identical after adding the
  separately selected KEP implementation and observers.

Rebuild the native driver with the ordinary flags in the build receipt,
adding `-DGFS_NATIVE_KEP`. This selects experimental checkpoint version 2;
version-1 late HLLC states cannot be silently resumed under the new core.
Use fresh analytic initial conditions and a new GPFS run directory.
`Exam/HydroCore/run_kep_native_group.py` records the exact test matrix,
fixed law and resource allocation. It returns a failed group status when
the cold physics gate fails, even though individual endpoints complete.

The initial pure optimized C/C++ flux evidence is in
`Research/FVSmooth_20261105/kep_flux_v1`. No actual GPU or MPI FV evolution
qualification is claimed. Original shock cases not yet in this harness,
higher-resolution physical qualification, the full 2D suite and matched
GIZMO timing are still outstanding.

The manifest hashes identify raw and stored bytes. Gzip files are lossless.
Native executables are excluded and identified by build receipts. Their
authoritative location is
`/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1`.
