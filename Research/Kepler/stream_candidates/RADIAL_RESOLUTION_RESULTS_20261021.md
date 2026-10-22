# Smooth inward refinement: native CUDA comparison

2026-10-21. User-approved test, implemented and executed by the driver on
lageunha/RTX5000Ada. The full Kepler goal remains ACTIVE/UNACHIEVED. Both
native trajectories have terminated on strict nonpositive intermediate U;
neither is a t6/8/10 result. The one matched half-step control is complete.

## Frozen implementation and initial state

See RADIAL_RESOLUTION_PLAN_20261021.md for the ONE bundled Astra design.
Native v4 binary SHA256:
`f6f51995a76861b8438ee6649185b860ceb671c3209097005c7fdf6e13cb250f`.
Original RK4, LSF/old limiter, complete raw GPU geometry, W0, fixed individual
masses, Re infinity and no steering during evolution. No intrinsic-spin
evolution or new material/remap operator is claimed for this control.

Both initial meshes use the same non-ring Sobol/CVT algorithm. The graded
target n(r) is smooth with rc.25, central h0/8, and taper2-4. Actual data:

| | Uniform | Graded |
| --- | ---: | ---: |
| particles | 65,536 | 72,518 |
| total mass | 75.97420601059443 | 75.97420601059443 |
| total initial K | 8.789134864276413 | 8.789369513164210 |
| total initial L | 167.50361206293064 | 167.50214019989534 |
| minimum generator r | .055592701 | .005258266 |
| cells with generator r<.047 | 0 | 51 |
| median h in graded r<.047 | — | .011701946 |
| maximum adjacent h ratio | 1.18990 | 1.26150 |
| maximum inertia eigenvalue ratio | 2.24267 | 2.40400 |

N grows10.65%. Mass is integrated from the original continuous field over
actual periodic cells; velocity remains the original pointwise circular
generator value. K/L are not forced to match by changing heat. Initial
cell-mass quadrature changes are8.85e-10/4.17e-9, below the original1e-8
criterion. CVT ran a fixed16 initialization steps, not a proof of an exact
CVT minimizer. The actual mesh properties above are measured.

Authoritative ICs and manifest: `radial_resolution_20261021_v2/`.
Preparation envelope: `radial_launch_iJjAT4En/`; native envelope:
`radial_native_launch_m1jadrdy/`. Source/IC/binary SHA values are recorded.
Existing code, trajectories and restarts were not removed or reset.

## Preparation and audit corrections, preserved

The first preparation attempt (`radial_launch_2t32s3s3/`, exit1) failed its
density integration convergence criterion. No fluid evolution ran from it.
The corrected preparation subtracts the fixed central radial term on every
cell and integrates it by the identity
div[x (r-a log(1+r/a))/r^2]=1/(r+a). Edge quadrature resolves closest-point
structure; the physical density and tolerance are unchanged. Independent
nested integration checks of a central and an adjacent square differ by
3.55e-12/8.45e-12 relative (`radial_ic_integral_check_v2.json`). The first
checker also preserved a JSON serialization failure before that correction.

The first uniform native preflight completed its positive one-step solve,
but the old face-count audit classified all36 records in12 triple groups
as orphans. Inspection proves each group has TWO matched large ends and
ONE extra contact end. The maximum extra length is1.55e-17, extra absolute
force7.95e-23. The new explicit audit retains these extras and rejects any
group with a genuinely unmatched large face or a large split face. It uses
the SAME geometry tolerance2.73e-12. The already computed native preflight
was reused, not rerun. Original audit and exit remain preserved in
`radial_native_launch_q_8hwf3a/` and `.../native/uniform/`.

## Initial native geometry and MPI

- First native volume errors versus independent periodic Voronoi:
  uniform5.95e-8, graded5.93e-8 (float storage precision).
- Full group geometry residual maxima2.55e-15/2.86e-15, and surface closure
  3.93e-17/4.39e-17. Actual absolute force/work sums are in the reports.
- The graded MPI1/4 one-step states are BITWISE equal in every recorded
  field and metadata. This is a one-step layout control, not scaling or
  long-term MPI invariance.
- GPU logs identify the RTX5000Ada. The independent layouts overlapped on
  the same card; timing is not an exclusive-GPU performance comparison.

## First spatial operator

Sum over cells whose initial generators are in r<.5:

| | Uniform | Graded |
| --- | ---: | ---: |
| sum positive Udot | 3.5656144431e-4 | 3.5898616382e-7 |
| max Pstar/P0 in r<.5 | 908.4465 | 3.43159 |

The first positive Udot sum is about993x smaller. It includes reversible
pressure work; it is not a pure irreversible-heating rate or a long-run
convergence order. Actual geometric Vdot/V remains nonzero for the exact
divergence-free input velocity. Work/geometry mismatch is near arithmetic
roundoff, so equating those work velocities does not remove the material
volume issue. `first_face_work_{uniform,graded}.json` retains all bins.

## Actual time evolution and failure witnesses

| | Completed exact checkpoints | Failed attempt |
| --- | --- | --- |
| uniform | .05 | step46, t_start=.05, k2_input, ID35776 |
| graded | .05, .1 | step183, t_start=.110490131230797, k2_input, ID42539 |

Both failures are native exit92. Invalid U values are-5.62054535e-11 and
-1.05409340e-13, respectively. Each retains its NORMAL step-start and
INVALID stage. A failed-stage metadata time is the step-start time, not an
accepted evolved endpoint. No floor or automatic retry is used to continue.
W0 and individual masses remain fixed in the positive checkpoint audits.

At common t.05, graded orbital position/velocity L1 errors in initial
r<.047 are.000621/.000637; entropy-ratio L1 is.24087, with U/Uad range
.7676-2.6158. Thus its small orbital errors do not imply cold accuracy.
At its failure start, the minimum U/Uad in that bin is.01087 and its
entropy L1 is.81842. This path fails the predeclared physical promotion
gate and does not support a t10 launch.

Full native evidence: `radial_resolution_20261021_v2/native_v2/`.
The original abrupt-grid failures near.102/.119 remain separate evidence.
These different layouts do not establish that smooth transitions alone
explain the small difference in failure times.

## Independent density verification at common t=.05

`density_t005.json` integrates the original stationary density over CURRENT
independent periodic Voronoi cells; the reference mass quadrature changes
are below6.18e-9 and native area errors below5.95e-8. The density plot
`density_t005.png` is the actual W0 piecewise constant m/V field, sampled by
periodic nearest-generator queries, with common spatial and color limits.
It does not interpolate site colors. Core color overflow is shown by the
colorbar extensions and numerical extrema are in the JSON.

For fair comparison, `eulerian_density_t005.json` additionally samples BOTH
fields at identical fixed Eulerian pixels, rather than selecting differently
sized cells by their centers. The following are spatial L1 estimates
sum|rho_num-rho_exact|/sum rho_exact, not rigorous integral bounds:

| physical radial bin | Uniform initial | Graded initial | Uniform t.05 | Graded t.05 |
| --- | ---: | ---: | ---: | ---: |
| 0-.047 | .439% | .186% | 5.691% | 11.067% |
| .047-.1 | .110% | .00967% | 10.760% | 14.682% |
| .1-.2 | .0304% | .00301% | 8.604% | 6.220% |
| .2-.5 | .173% | .0390% | 2.148% | 2.223% |

The initial density representation improves, while evolved inner density
noise does not uniformly improve. Therefore the measured initial Udot
reduction cannot be promoted to physical-solution convergence. This is
consistent with a material-volume/strain/transport problem surviving
refinement; this pair alone does not prove a general nonconvergence theorem
or isolate one mechanism uniquely.

## Completed half-step control and next action

ONE effective half-step control reuses the graded positive t1e-7 checkpoint
and targets the SAME t.05. The nominal82 steps have median dt.000617086;
the fixed half ceiling is.000308543. The same binary completed163 accepted
steps versus82 (652 versus328 RK RHS evaluations), native exit0, at the
same exact t.05. Thus this was an effective time-resolution control. Its record is
`radial_resolution_20261021_v2/half_step_t005_v1/manifest.json`, immutable
launcher `radial_half_launch_ams0mjy5/`. Do not repeat it or scan dt/order.

Initial r<.047 entropy-ratio L1 changes only.240869834 -> .240519679;
U/Uad maximum2.615835 -> 2.615329 and minimum.767557 -> .768269.
Position/velocity errors also remain nearly unchanged. This supports
spatial/volume/transport errors dominating this early interval. Two dt
values do not establish a time-convergence order or exclude every temporal
contribution. No half-step trajectory beyond t.05 was run.

All completed positive intervals, including the half control, report zero
floor injection. Half-control individual masses are bitwise fixed and all
weights are0. Preparation/native/half source manifests and the frozen
native binary SHA were checked again after completion and all matched.
All owned calculation and analysis handles from this bundle are terminal.

`entropy_endpoint_ledger.json` and its compressed per-cell CSVs retain the
signed endpoint delta-log-U and (gamma-1)delta-log-V terms. They use saved
positive states only and do not interpolate missing stage rates.

Retain this graded IC and its positive checkpoints as a native regression.
Evaluate the residual material-volume/thermal/strain error using the saved
data and choose the next coherent operator bundle. Do not launch an invalid
trajectory to t10, claim a long-run pass, change eps/P/Re, or return to the
old failed abrupt experiment. One completed end-of-bundle Astra evaluation
selected the next implementation in ASTRA_NATIVE_TRANSFER_PLAN_20261021.md;
no microstep re-review is needed. The full long Kepler goal remains unchanged.
