# Native q1 direct-U RK4 arithmetic probe (2026-10-27)

The corrected sparse-q2 native disk has about 52 times less absolute
log-entropy error than full q2 at `t=.825`, while its same-ID density field
is nearly identical. The sparse run took far fewer global steps and has a
different core/outer layout, so that contrast does not identify an error
source. A low-cost arithmetic control uses the **same accepted full-q1
checkpoint** and the **same accepted executable** as the already submitted
energy-mode reference **437231**.

The source's existing guarded `GFS_RK4_DIRECT_U=1` path stores accepted
float `U_n` once in a per-particle double anchor and forms each RK4 stage
from `U_n` with one cast. It retains the same force, EOS, RK4 coefficients,
geometry, fixed mass, W=0, and `GAS entropy_mode=0`; it does **not** make
the accepted internal-energy field double precision. The source guard
requires CUDA W0/av5 Kepler phase 1, no conflicting entropy/weight/stress
mode, and aborts otherwise. This is a narrow test of repeated float update
arithmetic, not a new physical closure.

The native layout check compiled against the accepted 576-byte particle
ABI gives `stress.E_inv_xx` byte offset 472. All 147,456 cells in the
accepted `t=.825` checkpoint have zero there. The
[anchor auditor](audit_kepler_q1_direct_u_anchor_20261027.py) will require
positive finite endpoint anchors, close to final U, plus a passed full
physics gate. This provides a native checkpoint witness that the opt-in
branch actually executed. The already accepted binary contains the source
guard string and has SHA-256
`7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284`.

A source-level trace of the accepted W0 CUDA path found no stage overwrite
of this anchor: `fillParticleSoA` copies it into the device payload, but the
Voronoi tessellation and pressure kernels use that slot as thermal U only
when the mutually exclusive `thermal64` option is on. The ordinary
`writeBackResults_range` writes accelerations, thermal rate and selected
stress diagnostics, not `E_inv_xx`; the function that writes an inverse
moment matrix back is confined to the guarded LagMFM density path. This
supports the opt-in slot contract but does not replace the endpoint anchor
audit, especially across migration and native restart.

The [runner](kepler_q1_direct_u_to0845_20261027.slurm) and
[stager](stage_kepler_q1_direct_u_probe_20261027.py) created a new
`/gpfs/kjhan/LagForce/Kepler_DirectU_H200_20261027` directory from the
byte-verified strict q1 `t=.825` checkpoint SHA-256
`23fa751f5a52f27f8e8d17b51d8c85a8d334ce0e1aa57f56f5fbf2678053f7c7`.
Its staged scripts and inputs passed SHA-256 checks, shell syntax, Python
compilation, and imports. H200:1 job **437239** was submitted for exact
`t=.845` with two OpenMP threads, one MPI rank, 32 GB and a one-hour cap.
The job is pending GPU quota; there is no direct-U physical outcome yet.

The pending bundle was rechecked in place after the unrelated q2c auditor
packaging failure: all `scripts.sha256` and `inputs.sha256` entries pass,
and the strict, radial and anchor auditors all import in the run Python
environment. This checks launch readiness only.

Compare its **incremental** entropy from the shared `.825` state with the
energy reference **437231** and fixed-K arm **437232**, together with E/L,
orbit, density, positivity, injection, wall time, and the stored anchor.
The [same-start pairwise auditor](audit_kepler_q1_direct_u_pairwise_20261027.py)
is ready and checks the baseline/direct checkpoint identity, strict gates,
stored anchor receipt, initial-radius entropy increments, and same-time
disk-field differences. It must run only after both native endpoints exist.
Even a clean `.845` result would only establish a short arithmetic
improvement. The `t>6` and GIZMO-length gates remain separate.

The original H200 direct-U job **437239** is held while its byte-matched
A10 replacement **437268** waits for GPU quota; the energy reference is
A10 **437266**. The fixed-K H200 arm **437232** was canceled before
execution after its missing-CUDA build was detected. See the
[A10 recovery report](KEPLER_Q1_A10_ARITHMETIC_RECOVERY_2026-10-27.md).
