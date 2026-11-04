# Constrained reflecting fields and true PCG stopping

The scope, failures, grouped Astra evaluation, and next unresolved issues
are described in [the report](../../Docs/UNIFIED_HYDRO_WALL_TRUE_RESIDUAL_2026-11-05.md).

JSON receipts are plain text. Every other retained file is lossless gzip,
including frozen source, logs, independent analysis scripts, raw tables,
and accepted/rejected state files. The manifest records original byte counts
and SHA-256. Compiled ELF executables and Python caches are excluded while
build commands and binary hashes are retained.

Both version5 full-target strong tests fail and must not be counted as passed.
The initial 16 and strengthened 22 fixture successes are operator checks, not trajectories. The original
version4 blast replay advances zero accepted steps and preserves its accepted
checkpoint byte-for-byte. failure_trial.fv and retry_trial.fv carry a diagnostic
flag rejected by the production restart loader. They are not restart points.

The stress-length substitution and monotone-limiter observations are frozen
algebraic diagnostics only. Neither candidate has been evolved or adopted.
The final rebuild changes only the native infeasibility message and adds an
audit stress-point CSV. Its wall fixture records match the earlier build.
No successful PDE is inferred from the final build.

The manuscript and accepted Kepler history are unchanged. Production GPU/MPI,
full common 2D tests, and matched GIZMO runtime <=2x remain unfinished.

The default HLLC Sod byte control retains its first roundoff-only comparison
failure and the corrected byte-identical result in separate directories.
