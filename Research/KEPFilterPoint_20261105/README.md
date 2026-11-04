# KEP stress-filter and point-query development evidence

The scientific scope, failures, version distinctions, and next gates are in
[the report](../../Docs/UNIFIED_HYDRO_KEP_FILTER_POINT_2026-11-05.md).

Original files remain under `/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1`.
JSON receipts are plain files. All other retained files use lossless gzip,
including frozen source, logs, CSV/data, and checkpoints. `manifest.json`
records the original SHA-256, byte count, relative path, and encoding.
Compiled ELF executables and Python cache files are excluded. Executable
hashes and full build commands remain available in the original manifests.

Keep the failed runs and controller exception. Version3 smooth-flow results
are not version4 results. IC-only checks are not evolved trajectories.
`failure_trial.fv.gz` is an unaccepted diagnostic state explicitly rejected
by the production restart reader. `final.fv.gz` contains accepted states.
`push_review_20261105_v1` records the exact failed-trial block-preconditioner
gate, not a successful full strong-shock or Kepler simulation.

The manuscript remains the previously verified baseline. The experimental
candidate is not yet adopted and no matched GIZMO timing is claimed.
