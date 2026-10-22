# Native enriched Euler initialization and fixed-cell operator evidence

Subsequent full initial moving-cell RHS implementation and failed thermal
tangent gate: `ENRICHED_COUPLED_RHS_RESULTS.md`. This document retains the v4
initialization evidence; it is not a time-evolution certificate.

2026-10-21. This bundle advances the active long-time Kepler goal through
representation and operator checks. **No enriched DG timestep has been taken.**
Original native t6-10 failure and GIZMO duration remain unresolved. Current
frozen candidate is `enriched_v4_build`, binary SHA256
`7d6c2824131d70df8cb8f462cc558e8f327c76175c032b6f973e685ce2f0b9f0`.

Repository identity remains the native worktree in `DG_PROJECTION_RESULTS.md`,
HEAD8e25dfe plus accumulated uncommitted driver changes. No original source,
input, checkpoint or output was deleted. No commit/push or shared-context
update occurred. Source manifest was rechecked after runs. Calculations ran
on authorized lageunha CPU/Ada; Python launches/analyzes only. All execution
handles completed; no background monitoring loop or new Slurm job remains.

## Implemented scope and missing runtime

`Exam/gfs_enriched_euler.h` adds the common scalar space P2 plus four disclosed
reference Euler fields, dimension at most10. It includes analytic reference
gradients, physical-coordinate basis functions, complete Euler volume flux and
gravity source, and weak-form application of a supplied canonical ALE face
flux. Every coefficient is an ordinary unknown. Perturbation gives nonzero
polynomial AND enrichment coefficient derivatives. There is no reference
reset, stationary-residual subtraction, independent Pi, thermal refill,
limiter, viscosity, damping, steering, mass or weight change. Spatial momentum
modes carry internal angular information.

`gfs_enriched_gpu.inc` applies projection to actual native CUDA support
polygons. `gfs_enriched_native.inc` exposes it only through
`GFS_DG_PROJECT_ONLY=1 GFS_DG_ENRICHED=1`, order2. Tagged output packets are
representation artifacts. They are not evolved snapshots or restart states.
Default old RK4/SSP physics path passes compatibility replay.

Still absent: full native moving-cell DG RHS with neighboring states and
actual shared-face velocity; stage geometry/mass-matrix/basis update; Euler
mass-flux/weight coupling; DG payload migration/halo/checkpoint/restart;
conservative stage positivity/rollback; shock validation. These fixed-cell
helpers do not implement those features. Projection alone cannot start a
Kepler continuation under the selected Astra plan.

## Conditioning corrections and preserved failed artifacts

- v1 rigid CUDA control stopped on592 negative fan measures. Rebuilding an
  endpoint as a+(b-a) could move an anchor one ulp outside an adjacent edge.
  v2 selects the stored endpoint when projection t=0/1. No finite edge cutoff,
  area tolerance, weight change or mesh change.
- v2 pressure reconstruction is accurate but conditioning fails: maximum
  orthogonality error .99999999988. Independent outer-cell CPU integration
  reproduces it and very large coefficient derivatives. It was never evolved.
- v4 evaluates the density contrast and kinetic reference energy BEFORE
  adding their constant background/thermal parts. Constants are in P2, so
  the mathematical function space is identical. It evaluates the tiny halo
  density component directly. The taper .5(1+tanh z) uses its algebraically
  identical logistic expression with stable exponent sign. Mathematical IC
  and parameters stay the same; pointwise relative rho difference from the
  literal double expression is <=2.45058e-10, exported as metric14.
- Reorthogonalization solves the measured previous-basis Gram instead of
  assuming exact floating-point orthogonality. No ridge or artificial pivot.
  Excluded directions have normalized residual <=1.40808e-14. Their actual
  field and pressure errors are included in the independent assessments.
- v3 derivative audit exposed a factor2 derivative error in the logistic
  rewrite. v4 corrects it and makes the independent derivative gate fatal
  before library compilation. v3 had no native projection run or evolution.
  Failed artifacts are retained.

## Whole native mesh: v4 projection only

Same71080 positions, epsilon=.047, P0=1e-6 and original smooth mathematical
profile. Positive segmented quadrature uses8x8 to build and16x16 independently
to assess, with vertex/face samples. Selected cells use denser rules. All71080
geometry/representation builds complete with no limiter.

| Quantity | Measured value |
|---|---:|
| Nonpositive sampled rho / P cells | 0 / 0 |
| Maximum relative P error | 3.40294e-10 |
| Global volume L1 relative P error | 5.92838e-12 |
| Maximum field error / physical component scale | 9.69277e-14 |
| Maximum independent conserved-integral discrepancy | 3.36037e-13 |
| Maximum quadrature E difference / physical thermal energy | 2.13821e-10 |
| Maximum orthogonality discrepancy, build quadrature | 5.70136e-7 |
| Active scalar dimension9 /10 cells | 15352 /55728 |
| MPI1/4 every exported field and flag | bitwise identical |

Measured build-quadrature Gram is close to identity; future transport must
use its actual stage mass matrix. This is not a moving-stage condition proof.
Cold boosted4096-cell control max relative P error2.43409e-8; rigid4.71269e-10;
varying-rho isobaric1.16361e-10. All are positive.

Native wall times including startup/geometry: actual mesh3.505s MPI1,
5.936s MPI4;4096 controls2.93-3.04s. Four MPI processes shared one Ada for
consistency; this is not a speedup measurement. Authoritative data:
`enriched_v4/{manifest,analysis}.json` and tagged packets.

## Fixed-cell full Euler weak form and perturbation

`enriched_v4_build/contract.json` uses8 actual central polygons:4 touching
the origin and4 historical thermal-depletion cells. Basis rules24x12,48x20,
64x28 are assessed with64x28 and segmented faces. Full HLLC face flux,
softened gravity and Euler volume term are compared with independent analytic
conservative Euler strong derivatives integrated against the same functions.
Cells/faces are static. No stationary RHS is zeroed manually.

- Gradient discrepancy against independent long-double five-point differences
  <=3.785e-11, normalized by1+|derivative|.
- Initial P relative error <=3.645e-12; independent integral error <=2.790e-14
  across rule choices; no initial limiter.
- Equilibrium weak/strong discrepancy <=2.980e-12, normalized by
  area*component_scale*(speed/covariance_length).
- Velocity perturbation (.013,-.007) with fixed potential preserves initial
  P and has weak/strong discrepancy <=4.143e-11 in that normalization. All
  selected cells have nonzero polynomial and enrichment coefficient changes;
  minimum normalized maxima .000132017 and2.60323e-6 respectively.
- Three old ill-conditioned outer cells are rechecked in `outer.json`:
  orthogonality <=5.319e-8, pressure <=1.613e-13 and weak/strong discrepancy
  <=9.379e-9. This finite-precision outer residual remains relevant to future
  long moving-cell calculations.

These are instantaneous spatial identities, not an evolved perturbation,
MPI shared-face RHS check, or shock entropy proof. Separately,
`enriched_legacy_v4/manifest.json` verifies the unchanged old native3-step SSP
replay: all recorded fields and metadata bitwise identical,5.050s.

## Continuous fresh IC versus legacy point-state budgets

`enriched_v4/initial_ledger.json` compares integrated continuous fields with
the old serialized point-lumped state. Their cell momentum/E differ;
projection preserves its continuous integrals. It has not converted an
evolved checkpoint conservatively.

Both total masses are approximately75.9742060105944, with summed cell mass
differences -7.34e-15 and maximum local relative difference5.52e-11. Legacy
gas E8.78999949646204; continuous integral8.78971230175656. Delta is
-2.871947054827e-4, about3.27e-5 of gas E and .3324 of the tiny thermal total.
The separate thermal difference is only -5.99e-12; kinetic difference is
-2.871946994918e-4. Thus resolved continuous versus point-velocity kinetic
integration accounts for the difference. No actual input state was changed.
Fresh initialization must explicitly adopt a consistent ledger before
evolution; this cannot be described as preserving legacy cell p/E or as a
restart from an evolved state. The old failure outputs remain available.

## Immediate next work

Follow `ASTRA_ENRICHED_EULER_PLAN.md` as one bundle. Connect actual moving-cell
face/volume/source RHS and mass flux to the approved weight/mass relation.
Evolve conservative moments; handle actual stage Gram, fixed-within-step
basis, transformations and GCL. Mean gas E minus mean translational kinetic
energy includes resolved velocity variance. Keep the actual thermal/EOS
diagnostic tied to the resolved field instead of feeding that variance into
the old point-cell temperature.

Complete native runtime state, MPI/restart and perturbed evolution, then
same-path Sod/Sedov with enrichment off. Only after those gates: one original
layout t=.2 path followed by t6-10 and the longest-paper duration. Existing
normal restarts remain preserved. No new planning call for each micro-step.
