# Shared Duffy / AD2 / physical-query / target68526 results

2026-10-21. Driver implementation and native execution of the approved
bundle. The long Kepler objective remains **active and unachieved**.

## Results and limits

| Obligation | Result |
|---|---|
| Same canonical function in fixed six vertex Duffy charts, AD2 | PASS, all384 sources |
| MPI1 / MPI4 per-source Duffy ledger | Byte identical |
| Native-cache versus canonical physical all12 point queries | PASS,2304 diagnostic queries |
| Canonical edge / target-plane event completion | Initially11 unproved; exact composition resolves all11 without changing the256-window cap |
| Actual target-domain coverage |344 exclusions +4 whole sources +10 of36 intersection domains proved;26 intersection domains unproved |
| Complete all12 target68526 integral and whole error budget | REJECTED:40 active source contributions remain incomplete |
| Accepted positive target-integration queries / publication |0 / no publication |
| Original coefficient and frozen generating-state SHA | Unchanged |
| Time evolution, RK4 step, historical t6/8/10 stability | Not exercised or established by this bundle |

The current saved finite state is the OLD kinematic control with h=2e-4
and the old acceleration. It is not an accepted evolution step of the new
compatible-pressure operator. No simulation was restarted or reset.

## State and implementation

Authoritative source repository remains
`/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src`,
HEAD `8e25dfe0c2e41a7374c11e2e3102da2a8b028625`,
origin `git@github.com:kjhan0606/LagEunha.git`.
This bundle adds diagnostic implementations under this directory; no source
repository commit/push was made.

Input remains `shared_trace_rhs_v1_mpi4/shared_state.bin`. Original physical
U is read from `geometric_pressure_v2_mpi4/finite/frozen_physical_rank*.bin`.
The full817152-point native-cache/heat proof is retained by protected
source/input/output SHA checks. It was not rerun as a new accuracy proof.
The target call reads those817152 saved points through the native CUDA
operator and uses the CPU native J/grad_v measure to recompute the existing
discrete compatible-EOS context. No new IC, independent heat, or spin
reservoir is introduced. The v2/v3 EOS-context CSV files are byte identical.

The immutable heat recipe applies ONLY to this control:

    U_per_reference_volume = P0/(gamma-1) * J_GENERATING^(1-gamma).

At these queries the generating geometry equals the saved, immutable
current control. This is not a general evolved-thermal restart API.
Physical heat is never assigned from the corrected thermodynamic volume.
The same native operator context supplies fixed per-cell delta and
delta_dot for this call:

    dm = rho0(X) dX
    dU = U_generating(X) dX
    P = (gamma-1) U_generating/(J+delta)
    dUdot = -P [J div(v)+delta_dot] dX.

The push-forward density rho0/J is distinct from EOS rho0/(J+delta).
The physical all12 query includes mass,px,py,angular,K,U,volume,qx,qy,
polar,Udot,potential. The2304 diagnostic queries do not supply a target
quadrature remainder, coverage certificate, or accepted integral.
The new transfer/error calculation is a CPU reference implementation;
CUDA is used for the saved native-field reader. A production CUDA target
transfer certificate has not been established.

## Duffy function and derivatives

`material_shared_duffy.h` implements the same shared enrichment with exact
rational reference nodes, canonical edge functions and stored global
coefficients. Each source has fixed six subtriangles formed by its original
vertices, adjacent edge midpoints and centroid. The original vertex is the
Duffy apex. The opposite-edge factor is differentiated after cancelling
Cartesian1/s algebraically:

    E = s^4 r(t)
    E_i = s^3 [4 s_i r + h_i r']
    E_ij = s^2 [12 s_i s_j r
                 +3(s_i h_j+s_j h_i)r' +h_i h_j r''].

The actual-vertex position/velocity identity is used only at proven exact
radial0, and the canonical endpoint identity only at exact t0/1. Derivatives
retain their bounds. No near-edge flags or epsilon denominators are used.
The same analytic v0-T3 remainder is evaluated with a bounded infinite
binomial tail; the fixed256 cap retains any omitted remainder explicitly.

Native Duffy preflight:

-384 sources,2304 fixed charts;9216 interior,6912 exact apex,
  3456 canonical-edge queries.
-295296 scalar comparisons against the independent original canonical
  first-derivative function and direct, uncancelled AD2 expression.
-0 enclosure mismatches; maximum observed difference/bound .1655359309.
-Minimum real sampled J .999999995876253574; the two other extension
  denominators are bounded below by .5 over each entire real chart.
-`shared_target_v3_mpi4/duffy_preflight` and `shared_duffy_v2_mpi1`
  source CSV files are byte identical. Native times24s MPI4 and93s MPI1
  are for384 diagnostic sources, not production128² Kepler throughput.

Bounds use outward Q128 balls and the inherited <=8ulp elementary-function
engineering model. These are conditional numerical certificates, not
formal proofs of CPU/GPU libm implementations.

## Exact plane composition fixes the event failure

The initial native whole call `shared_target_v2_mpi4` correctly rejected:
11 edge/plane event pairs exhausted the unchanged256-window budget and
30 actual source domains were unproved. The event evaluator was subtracting
ambient coordinates after evaluating the curve.

`material_shared_target_geometry.h` now composes each linear power-plane
with the SAME canonical edge before evaluation. Exact finite expansions
compute the plane values of the stored P3 nodes and contractions of the
stored enrichment coefficients. The resulting 1D function is algebraically
identical to n.chi(t)-b. Position, derivative, curvature and events are
owned by one canonical edge record and reused by the adjacent sources.
Exact vertex/plane zeros have explicit endpoint events.

`shared_target_witness_v1_mpi4` selects ONLY the11 previously unproved pairs:
all11 complete with11-66 windows at the SAME256 cap.165 independent
value/first/second derivative comparisons against n.edge(t)-b agree within
their error intervals, maximum ratio .2324066333. The same run verifies2304
common all12 native/canonical physical queries with positive mass and U.
These files are explicitly named DIAGNOSTIC; they are not target packets.

The final whole call `shared_target_v3_mpi4` uses that composition also for
source-edge fragment signs. It reports0 unproved edge/plane events and26
unproved actual source domains, improving30 to26. Runtime78.05 to37.13s
compares these two same-input reference calls; it is not a full-grid speed
claim.

## Remaining actual failures

The final call processes all384 classifications. Of the36 candidate
intersections,10 have complete closed, positively oriented physical-fan
domains, persisted in `actual_target_nodes.csv`,
`actual_target_boundaries.csv`, and `actual_target_loops.csv`.
The other26 reject on real membership/vertex handling or an insufficient
empty-arrangement witness. An unproved classification is never treated
as an empty contribution.

Old difficult sources:

| Source | Current result |
|---|---|
|118 | Edge events now complete; real membership/Newton still reaches an opposite-vertex Cartesian chart exclusion |
|246 | Actual intersection domain complete;20 attempted sectors reject because the Cartesian complex neighbourhood includes a Duffy apex |
|254 | Edge events complete; strict real inside/outside membership proof remains unclosed |

Final actual-sector remainders:244 attempted,0 certified.148 Cartesian
inverse neighbourhoods include a Duffy apex,86 binomial analytic domains
are not proved,10 J disks are not separated. These are failed bounds/routes,
not evidence of a negative physical Jacobian or absence of a fluid solution.
The clipped reference-Duffy target-boundary route is not yet closed;
physical fans that include a source vertex cannot be promoted to analytic
Cartesian inverse certificates.

The witness run also examines both Cauchy axes for all24 sectors of each
of the four ACTUAL whole-inside sources,192 axis domains. It records86
binomial-domain failures,16 J-disk failures and90 density-square-root
positive-half-plane certificate failures. Scalar sigma/J failed arguments
and their coordinate boxes are retained. Generic sqrt exceptions have
uninstrumented zero argument fields: these zeros are NOT measured sqrt0.

One recorded box has normalized binomial |center|+radius .750981645,
above the unchanged .75 condition. Another J disk has center1.000000016
and radius2.27150218, so that enclosure cannot exclude zero. Real J remains
positive in the real function contract. Independent interval treatment of
affine coordinates and rational factors can overestimate complex ranges;
actual singularity versus enclosure overestimate has not been determined.

## All12 transaction and old physical tolerance

The original whole-target normalizers are retained:
U_lower2.2158233568699930306e-11,
M_lower1.5020370483398437803e-8,
area_lower1.50203704833984374926e-5.
The energy criterion remains U_lower*1e-10, about2.21582336e-21.
No partial allocation is substituted for a whole budget.

All40 active source contributions remain incomplete. Every all12
whole-error field is consequently `unproved`, not0. No positive
target-integration packet is accepted or published. The caller target
packet retains its original sentinel and the saved input coefficients and
immutable heat generator are byte identical to the original state:

`74062ad182ccc5ff1ca0dace272b1b685d78b751af2ea88f8267c9e49def3005`.

Before a future successful publication, the final R80 moment/error
conversion must be included BEFORE the whole tolerance comparison. The
current prototype appends that conversion to the reported error after its
Q128 budget check. This does not change these observed rejections (all40
contributions incomplete), but must be corrected before accepting any
target. This limitation is recorded rather than promoted as a passed gate.

## Files and reproducibility

-`material_shared_duffy.h`, `material_shared_target_geometry.h`,
  `material_shared_target_integral.h`.
-`test_shared_duffy_v1.cpp`, `test_shared_target_v1.cpp`,
  `test_shared_target_witness_v1.cpp`, corresponding frozen-launch runners.
-`shared_target_v1_mpi4` preserves a compile-only missing namespace error;
  corrected in v2. No physical test ran from that failed build.
-`shared_target_v2_mpi4` preserves the first real whole-call rejection;
  `shared_target_v3_mpi4` preserves the composition-improved rejection.
-`SHARED_DUFFY_TARGET_LEDGER_20261021.json`,
  `SHARED_DUFFY_TARGET_V3_LEDGER_20261021.json` and read-only
  `audit_shared_target_v1.py`; integrity PASS is distinct from physics PASS.

All calculations use unique destination directories, frozen source/env/
binary/input/output hashes on authorized lageunha. All native sessions
for this bundle are closed. No Slurm job, background monitor or simulation
is left running by this bundle.

The existing Astra agent completed ONE coherent end review, saved in
`ASTRA_SHARED_DUFFY_REJECTION_REVIEW_20261021.md`. It found no immediate
orientation reversal defect and accepted the rejection as necessary.
The selected next bundle connects certified winding/containment,
coverage independent of star-fan eligibility, ACTUAL clipped reference-Duffy
graphs/event coordinates, cancellation-preserving affine reference
coordinates, and a finite atlas set by actual singularity/event distances.
It also requires an explicit generating/current state identity guard,
reuse of v3 event/domain records by SHA, and final R80 conversion before
the whole gate. The entire ORIGINAL tau vector must be frozen in a manifest;
new sampled velocity/radius scales must not make a criterion more favorable.
One all384 target68526 call then publishes only if every obligation passes.
No additional approval or input is required for this already authorized
work. Long Kepler validation remains outstanding.
