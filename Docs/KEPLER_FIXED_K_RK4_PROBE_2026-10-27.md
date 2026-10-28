# Native q1 fixed-K RK4 attribution probe (2026-10-27)

## Question and status

The active W=0 variational pressure operator has independently passed small-mesh
force/work adjoint checks, but the native q1 `t=.825` disk has a mass-weighted
absolute log-entropy change of `1.53e-4`. A short matched native experiment
tests whether continuously updating float internal energy contributes to that
residual. It is **not** a long-time Kepler result or a general shock solver.

Two H200 jobs were submitted from the **same accepted** q1 `t=.825` checkpoint:

| Arm | Job | Binary | `GAS entropy_mode` | Target |
|---|---:|---|---:|---:|
| Energy reference | 437231 | accepted q1 | 0 | exact `t=.845` |
| Fixed K | 437232 | opt-in K initialization | 1 | exact `t=.845` |

Both use the user's CUDA+MPI+OpenMP RK4 executable, one H200, one MPI rank,
two OpenMP threads, the same 147,456-site q1 IC, fixed mass, W=0, phase-1
variational pressure, no physical viscosity, halo 12, and the same timestep
controls. The run directories and exact staged inputs are in
`/gpfs/kjhan/LagForce/Kepler_EntropyK_H200_20261027/manifest.json`.
The staging code refuses an existing destination and verifies the parent's
complete strict gate and input hashes. The original q1/q2 runs and dependent
long jobs were not altered.

## Source and physical meaning

The accepted energy-mode checkpoint has `stress.K=0` for all 147,456 cells.
Simply setting `entropy_mode=1` would make its EOS use the K floor. The
[opt-in patch](KEPLER_ENTROPY_K_INIT_2026-10-27.patch) initializes each K
once, before the first RK4 geometry/EOS update, from its accepted cell state:

\[
K_i=(\gamma-1)\frac{U_i/V_i}{(m_i/V_i)^\gamma}.
\]

The existing phase-1 force path sets `dK=0`; the entropy EOS subsequently
computes `P=K rho^gamma` and `U=P V/(gamma-1)`. This is an **isentropic smooth
flow closure**. Small measured entropy increments in this arm are expected
from the construction and cannot alone establish a better physical solver.
Its pressure force, orbital/density behavior, total energy, positivity,
floor-injection ledger, and wall cost must be compared with the energy arm.
The closure does not generate shock entropy and is not yet adapted to the KDK
path or individual time stepping.

Provenance:

- Accepted binary SHA-256: `7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284`.
- Experimental binary SHA-256: `6c2ec1ab5b0f1877108010ddb3c838c420712845c4e70c88135da5a7c0a8e26b`.
- Modified `Exam/exam.c` SHA-256: `354938194bae832d3c34903ef2a7bb0ab8e073e114a5d99df8d4d14dec099fdd`.
- Patch SHA-256: `022e716628dbec3fca9a2054958a05cf7939abc8f5b2e45f4a0f5f6d8058ddb0`.
- Shared checkpoint SHA-256: `23fa751f5a52f27f8e8d17b51d8c85a8d334ce0e1aa57f56f5fbf2678053f7c7`.

The source was built in the isolated `/tmp/Kepler_entropy_k_20261027`
copy with the same Intel MPI/CUDA module stack as the accepted branch,
`make -C Exam -j4 libexam.a` and `make -j4 this`. The source patch passed
`patch --dry-run` against the accepted source, the executable linked with no
missing shared libraries, and both staged scripts passed `bash -n` and their
SHA manifests. The C layout check gives a 576-byte particle and byte offset
552 for `stress.K`; direct reading of the parent checkpoint confirms all K
values are zero. This branch remains experimental and is not the active
q1/q2 executable.

After the unrelated q2c post-run packaging failure, both pending `.845`
bundles were rechecked in place: their complete `scripts.sha256` and
`inputs.sha256` manifests pass, and both staged auditor modules import under
the run Python environment. This prevents that known missing-module failure
for these two jobs; it does not predict their physical result.

For that accepted checkpoint, the reconstructed K spans
`9.82806e-7` to `2.15444e-3`; the code's default `K_floor=1e-30` would not
clip any of the 147,456 initial K values. This was checked directly from the
byte-verified native particle payload, not inferred from an idealized IC.

## Decision gate

Each arm must produce a unique complete exact-time checkpoint and pass the
full-state, fixed-mass, W=0, positive-volume/internal-energy, global-area,
`rho=m/V`, zero-injection, orbit, density, and E/L gates. The
[pairwise auditor](audit_kepler_q1_entropy_k_pairwise_20261027.py) computes
**incremental** entropy from the common `t=.825` checkpoint, by initial
radius, rather than comparing accumulated `t=0` errors. It also checks the
stored K payload by cell ID against the accepted state, thermal adiabatic
excess, `U=P V/(gamma-1)` and `P=K rho^gamma` EOS consistency, total-energy
change, angular momentum, and same-time disk field
differences. Only after both native results exist can it issue a comparative
receipt. A clean `t=.845` result would justify a longer controlled test, not
claim stability at `t>6` or at the GIZMO `t=600`/longest-orbit target.

## Superseding CUDA build and A10 jobs

The staged fixed-K H200 executable lacked embedded CUDA device code;
pending **437232** was canceled before running. The same patched source
was rebuilt with `USE_CUDA=1`, and matched A10 jobs **437266** (energy)
and **437267** (fixed K) now use native CUDA. See the
[build diagnosis and replacement manifest](KEPLER_Q1_A10_ARITHMETIC_RECOVERY_2026-10-27.md).
All earlier references to pending 437232 describe the superseded bundle.
