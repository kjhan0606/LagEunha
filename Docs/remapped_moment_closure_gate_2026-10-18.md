# Conservative remap recovery and the remaining closure gate

Read-only analysis of existing t2 conditional positive-density lift and weight
remap outputs; no new trajectory, checkpoint conversion, or production change.
Source: `/gpfs/kjhan/LagForce/Codex_checkpoint_weight_remap_20261018/`.
Evidence: `canonical_remapped_moment_recovery_2026-10-18.json`;
script: `audit_remapped_moment_recovery.py`.

For central spatial, cross and velocity moments S,C,T, recover
G=C^T S^-1 and R=T-G S G^T. The kinetic split is
K_internal=tr(G S G^T)/2+tr(R)/2. Thermal energy is
U=E_hyd-tr(Pi)/2, not E_hyd-|p|^2/(2m).

| Weight displacement fraction | Minimum U | Inner affine share of internal K | Edge affine share |
|---|---:|---:|---:|
| .25 | 3.70952e-8 | 10.26% | 3.59% |
| 1 | 2.64355e-8 | 22.66% | 7.88% |

The minimum eigenvalue of R is -1.235e-17 in both cases: consistent with
positive semidefiniteness at floating precision, not strictly positive.
The kinetic split residual is below 6.62e-24. The unchanged outer patch has
near-zero internal kinetic energy; its affine-share ratio is not meaningful.
Most internal kinetic energy is residual: 77--90% inner and 92--96% edge.
Recovered generator velocities change, but improvement in accuracy is unproved.

## Why retaining R alone does not close transport

The exact Pi flux contains integral rho v_a v_b (v-w).n dA. Its kinetic
energy trace therefore contains third velocity moments. The stored second
moments do not uniquely specify these or the face traces.

A simple normalized scalar counterexample makes this limitation explicit.
Distribution A has u=+1,-1 with probabilities 1/2,1/2. Distribution B has
u=sqrt(2),-1/sqrt(2) with probabilities 1/3,2/3. Both have mean zero and
second moment one, but their third moments are respectively zero and
1/sqrt(2). Thus their kinetic energy advection fluxes differ despite equal
mass, momentum and kinetic energy. These are mathematical distributions,
not a claim that the gas has physical multistream motion. Piecewise spatial
velocity structure also requires additional reconstruction information.

The existing affine flux headers assume an affine trace. They cannot by
themselves transport an arbitrary nonzero residual field. Setting the third
central moment to zero would be a closure assumption requiring validation;
converting residual kinetic energy into thermal energy would additionally
change pressure and is not justified by remap conservation.

## Decision and next gate

Reject affine-only integration of these remapped states. Preserve the existing
original CUDA/OpenMP/MPI runs and initial conditions. Before further component
expansion, compare a donor-subcell reconstruction that retains the remapped
velocity structure against the proposed moment trace on these same saved
states. First establish what face information the saved overlap data retains;
if absent, reconstruct it from the original donor geometry without rerunning
physical evolution. Measure missing momentum/kinetic flux and dependence on
remap size. Do not introduce a residual-pressure or thermalization law merely
to pass positivity. The persistently equal-weight original-model reduction
remains a separate necessary gate, and long-time stability remains unproved.

## Existing long runs, observed 2026-10-17 17:47:30 UTC

The one-shot process checker confirmed executable and cwd for three live runs:
baseline PID2286374 t15.534221; fixed-HG PID2295268 t15.424466;
from-t0 timestep-cap PID2466162 t11.346810. Their latest accepted steps were
finite with zero retries. These observations are not a physical accuracy pass.
No restart, cancellation or new simulation was performed for this audit.
