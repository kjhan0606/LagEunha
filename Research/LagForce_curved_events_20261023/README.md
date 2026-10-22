# Kepler curved event research checkpoint (2026-10-23)

This directory preserves the source, runner, and audit programs used for the
accepted-state 256² Kepler curved source-edge experiment. The original
OpenMP+MPI+CUDA application files remain unchanged. The full run history,
physical state identity, proof gates, and outstanding work are in
[`Docs/HANDOVER_CURVED_EVENTS_20261021.md`](../../Docs/HANDOVER_CURVED_EVENTS_20261021.md).

The current accepted physical time is `0.00011263377801805954`. The target
chart uses *trial geometric weights*; the screened fixed-mass law, positive
source/target overlap integration, conservative remap, and later RK4 steps
have not been accepted. These programs must not be treated as a validated
`t=6–10` Kepler evolution.

`material_compact_actual_source_v1/` contains the range producers, native
canonical reconciliation, certified cut ordering, triangle routing,
root-width witnesses, and independent audit programs. The parent directory
contains their bounded geometry and overlap dependencies. `Exam/` retains
four headers required by this source snapshot that were absent from the
tracked application tree. `MANIFEST.sha256` records the exact copied bytes.

The runner scripts retain the cluster's absolute input paths and require the
frozen accepted-state datasets, weight CSV, CUDA/Exam library, and MPI/CUDA
modules on that cluster. Large simulation inputs, CSV event ledgers, compiled
binaries, and Slurm logs stay in the run directory on GPFS; this commit is a
source and handover checkpoint. A new installation should supply equivalent
input manifests and update the runner paths before execution.

The root-width `2^-60` branch is a separate trial. Successful representative
witnesses do not certify every source edge; accept a chunk only after native
exit zero, frozen input hashes, binary128 roundtrips, original-edge presence
and exterior incidence checks, canonical-versus-wire reconciliation, and
parameter-order preflight. The retained `2^-40` chunks remain distinct.
The separate `probe_power_triangle_route_root60_v1.cpp` and matching runner
read the revised canonical ledger for native triangle-boundary routing;
`audit_power_triangle_route_chunks.py` audits both ledger variants. The
full ordered-cut prefix is currently `[0,1450000)`, while native
triangle routing has been independently audited for `[0,750000)` as one
contiguous block. These gates do not construct positive
source/target overlap masses or advance the physical Kepler state.
