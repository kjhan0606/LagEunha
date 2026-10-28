# Halo40 Ada wall measurements and conditional continuation preparation

## Completed native evidence

The four independent 59-step runs restart from the accepted t=2.6549183918680868 checkpoint and reach t=2.66 on the same RTX 5000 Ada, using one MPI rank and two OpenMP threads. Their native physics gates, endpoint fields and printed step traces pass. The candidate keeps halo=40 and vertex certification disabled, and enables padding reuse and local periodic alias refresh.

| Run | Solver wall seconds |
| --- | ---: |
| Canonical A | 256.428 |
| Combined B | 217.385 |
| Combined B repeat | 214.873 |
| Canonical A repeat | 250.816 |

Median wall time is 253.622 versus 216.129 seconds: 1.173475x speed, or 14.79% less time. The fastest canonical / slowest combined ratio is 1.153787, above the predeclared 1.05 threshold. These measurements apply to this GPU and interval.

Three additional endpoint comparisons pass for every stored byte after the native next pointer (record bytes 8 through 575), plus every byte of the 52-byte restart footer. The comparison covers 147456 records per endpoint and retains RK, gradient and auxiliary fields. Pointer bytes 0 through 7 are the only exclusion. This is stored restart equality, not an executed latest-checkpoint restart bridge.

Receipts: `KEPLER_HALO40_NET_WALL_ADA_RESULT_2026-10-29.json` and `KEPLER_HALO40_NET_WALL_ADA_RESTART_BYTES_2026-10-29.json`. The earlier halo20 geometry comparison remains failed; this halo40 evidence does not supersede it.

## Prepared code, not installed

The new staging, orchestration and candidate-chain modules are committed as preparation. At this commit they have not been run to stage or install the new workflow, have not changed the running checker or selected a production successor, and have submitted no new jobs. Existing execution state must be inspected before any later installation; job ID 438957 and t=3 are explicit historical bindings, not reusable defaults.

The proposed finite same-allocation validation requires accepted t=3 native and radial gates, verified cleanup, actual A100 device identity and adequate remaining allocation time. It compares halo40 full geometry/RHS and end states, runs dumps-disabled canonical/candidate ABBA measurements, and executes canonical/off/combined restart bridges using one natural untruncated RK4 step from the accepted checkpoint. It may atomically select the separate continuation only after all proofs pass. The original accepted checkpoint remains the production start; bridge outputs are not promoted. Ordinary validation failure with verified cleanup retains the original successor; uncertain cleanup or transaction state stops advancement.

The separate continuation preserves checkpoints at t=4,5,6,6.5,7,8,9,10 and their native/radial gates. The orchestrator never submits a successor itself.

## Validation and remaining scope

Command from `Docs`: `python3 -m unittest -v test_kepler_halo40_promotion_20261029.py test_kepler_halo40_chain_20261029.py test_kepler_restart_bytes_20261029.py`.

All 17 tests pass. Orchestration tests mock GPU/MPI/Slurm execution and cover cleanup uncertainty, timeout, partial outputs, duplicate/stale receipts, compare-and-swap conflicts and the uncapped natural-step bridge. Chain tests cover original-chain preservation, checkpoint/device binding and mutation detection. Native-byte tests reject physical, auxiliary and footer differences while excluding only the pointer.

These tests do not establish A100 numerical parity or performance, an executed t=3 bridge, t>6 stability, physical spin evolution or GIZMO long-time stability. Those remain open.
