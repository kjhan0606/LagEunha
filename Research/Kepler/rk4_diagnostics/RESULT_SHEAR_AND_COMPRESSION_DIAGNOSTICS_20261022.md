# Shear and compression diagnostics during the live preprocessing bundle

2026-10-21T23:48:08.962046+00:00. Goal ACTIVE, unchanged.
Previous goal turn made progress: the native paired-rule connection compiled
and the single downstream dependency was installed and verified live. This
turn adds analytical and saved-data evidence without rerunning evolution.

## Actual state and dependency

Original producer launcher2166274, ranks2166305–2166308 and dependency2296727
were revalidated with ps. Producer prefix at report preparation:
222840/1138432triangles (19.57%), all compressed. The one actual
paired RK4 step remains dependent on FULL candidate completion; no partial
rule is adopted. Current accepted state remains four steps,
t=0.000082355880701376864. No old accepted interval or actual k1 force/CG was replayed.
Current native v2 snapshot/runner and original source were not modified.

## New quantified material-shear reference

KEPLER_EXACT_MATERIAL_SHEAR_20261022.md/.csv derive the exact circular-map
singular values from the actual softened IC Omega=(r^2+.047^2)^(-3/4).
Its determinant remains1 while cond(F) grows with differential shear. At t8,
cond(F) is about19.9 at r2, but about127067 near r=.05427 in the central
low-density region. At r=0 the limiting motion is rigid rotation and this
shear/stretch condition number is1. Thus the origin, near-softening shear and
the disk's inner edge are distinct regions. These are analytic flow values,
not the measured conditioning or folding of the actual numerical solution.
The nonpolynomial radial enrichment in st_shape also means a P3-only bound
would omit active functions from a future continuous-map guard.

The actual retained global dt_bound is4.322978134685708e-4. The analytic
origin full-gradient controller estimate is4.3229781401251614e-4, a relative
difference of-1.2583e-9. This supports a central-rotation cost bottleneck but
is not a saved minimum-source/branch witness. It authorizes no dt relaxation.
See KEPLER_ORIGIN_DT_COMPARISON_20261022.json. No physical change is prescribed
from this analytic reference; it is not steering or a background orbit force.

## New existing-state compression decomposition

KEPLER_COMPRESSION_HEATING_IDENTITY_20261022.md derives the positive-rule
Jensen identity for lambda=jhat/jhat0 and the native fixed generating entropy:

    U_min = U0*(V/V0)^(1-gamma)
    U >= U_min

At fixed whole-cell volume, nonuniform subcell relative compression can raise
U while preserving the same entropy law. This gives a mechanism that energy
conservation alone cannot rule out. It does not prove the historical cause.
The diagnostic must never be used to reset physical U to U_min.

The new saved-data-only analyzer read all71080accepted step4cells and the
recorded exact gamma/P0. It completed in2.196s; no point evaluation,
GPU evolution or mass solve occurred. Results in
material_full_compression_heating_audit_v1:

- Max |cell U/U0-1|: 4.92593088333687522282805593649313070386423E-13.
- Max uniform bulk-compression contribution: 4.92591744697980077633865159390978151540968E-13.
- Signed per-cell Jensen remainder range:
  -1.250206032076599846350373640148532045973E-16 to 1.26186322123134215994537015500962026360E-16.
- Global signed remainder/U0: -4.424879982772581339101227557386581859976306938871441548E-17.

At this VERY EARLY time, the observed cell energy changes match tiny bulk
volume changes; the residual beyond uniform compression is about1e-16 and
cannot be distinguished from finite arithmetic/CSV summation noise. This is
not a detection of subcell artificial heating or physical cooling. No signs
were clamped. A future compressed-rule state needs its own generating
positive-measure baseline and the existing quadrature discrepancy ledger;
this original-rule baseline cannot silently be reused across a rule change.

## Resource observation, no submissions/cancellations

Read-only syntax queries: normal QoS MaxTRESPU gres/gpu=8;
h200 partition QoS=part_h200; part_h200 GrpTRES gres/gpu=6.
Node syn104 has8H200 and reports6allocated, CPUAlloc46/64. Thus physical idle
cards alone do not establish immediate allocatability. No other project's
job, QoS or allocation was changed. The Ada producer/dependency remain active.

Continue the existing ASTRA_LONG_EVOLUTION_COST_DECISION_20261022.md bundle.
No routine Astra invocation or new physical method was introduced. Required
next evidence remains the actual same-stage paired step and measured cost;
t6/8/10, continuous injectivity and longest-GIZMO stability remain unproved.

## Latest dispatch update

Ada-only dependency2296727 is superseded; current dependency PID2351718
is material_full_positive_dispatch_v1. See RESULT_H200_PAIRED_DISPATCH_20261022.md
for the H200-first single execution, terminal-pending cancellation/atomic-claim
Ada fallback, and actual verification. Original producer2166274 is unchanged.
