# Kepler periodic-seam control (2026-10-26)

## Why this control is needed

The present native q=1/q=2 Kepler tests evolve a circular velocity field and
central Plummer potential in a **periodic** 6×6 Voronoi box. The disk region
`0.5<=r<2` is well inside the box, but the analytic circular field is not a
smooth, steady field on the periodically identified outer seams. In the
published Hopkins test the effective particle resolution is quoted over
`[-2,2]^2`; §4.2.4 does not establish that the specific 6×6 periodic outer
condition used here matches its boundary treatment.
The paper's resolution footnote describes reducing particle sites in the
inner and outer profile regions and points out that particle methods can
handle vacuum edges. The [isolated-boundary implementation audit](KEPLER_ISOLATED_BOUNDARY_GAP_2026-10-27.md)
shows why the current native W=0 variational CUDA path needs a new boundary
force/work contract before it can replace this periodic control.

With the center at `(3,3)`, at the identified faces `x=0` and `x=6`, write
`s=y-3` and `Q=9+s²+epsilon²`. The initial tangential velocities differ by
`Δv_t=6 Q^(-3/4)` and the normal accelerations by `Δa_n=6 Q^(-3/2)`.
The potential *values* agree at the seam, but their normal derivatives do not.
For `s!=0`, the common normal velocity `-s Q^(-3/4)` also carries gas across
the discontinuity, so the field cannot be a stationary circular-flow IC on
the torus. The same statements hold after swapping x and y.

The [independent t=0 snapshot audit](KEPLER_PERIODIC_SEAM_IC_2026-10-26.json)
checks both 147,456-cell polar layouts. Each seam has 1,150 generator pairs
within distance 0.026 across the periodic identification. Their median
absolute tangential velocity difference is `0.9773`, versus local sound speed
about `0.0127`: median jump/sound speed `76.94` (95th percentile `92.72`).
The analytic midline seam jump is `1.1547`, the sound-speed ratio `92.74`, and
the normal acceleration jump `0.2222`. The finite-distance normal velocity
*jump* is small, while its mean absolute magnitude is `0.245`, so many pairs
advect through the seam.

At native q=1 `t=0.41255`, all 44 short generator pairs were near a box edge,
and cells within 0.2 of an edge had mass-weighted absolute density and
pressure changes `0.0681` and `0.1298` respectively, far larger than the
disk values. See the [snapshot seam audit](KEPLER_Q1_SEAM_T041254_2026-10-26.json).
These observations identify an outer-background defect, **not proof** that
it causes the later disk failure; a native control must decide that.
An [ID-matched radial audit](KEPLER_Q1_RADIAL_T041254_2026-10-26.json) at
the same checkpoint finds mass-weighted absolute radial speed `0.03645`
and absolute density change `0.02543` for initial `r>=2.75`. In the
`1<=r<2` disk they are `9.97e-9` and `0.000517` respectively; the center
`r<0.1` radial speed is `1.97e-6`. Thus the outer disturbance is measurable
while the disk remains close to circular at this early time. The generic
[radial auditor](audit_kepler_native_radial_20261026.py) can apply the same
bins and an energy ledger to the forthcoming baseline/taper sparse outputs.

## Opt-in, mathematically balanced boundary control

An isolated experimental source copy at
`/tmp/Kepler_periodic_taper_native_20261026` adds opt-in variables
`EUNHA_KEPLER_PERIODIC_TAPER_INNER=r0` and
`EUNHA_KEPLER_PERIODIC_TAPER_OUTER=r1`. The proposed diagnostic uses
`r0=2.35, r1=2.75`: the published disk and its near guard retain the
original potential and initial velocity, while the whole box boundary is
outside `r1` and hence static. This is a *different outer-domain problem*;
its role is to isolate boundary contamination, not to pass it off as the
published IC.

Let `z=-(r²+epsilon²)^(-1/2)`, `z0=z(r0)`, `z1=z(r1)`,
`u=(z-z0)/(z1-z0)`, and `T(u)=1-10u³+15u⁴-6u⁵` for `0<u<1`.
Inside `r0`, use the original Kepler potential and speed. In the transition,
set `v_phi=T v_K`, `a=T² a_K`, and
`Phi=z0+(z1-z0) H(u)`, where `H'(u)=T(u)²`, `H(0)=0`.
Outside `r1`, `v=a=0` and `Phi=z0+(z1-z0)181/462` is constant. Therefore
`v_phi²/r=dPhi/dr` in every region with constant initial pressure; the disk
is unmodified, and the periodic boundary becomes smooth. There is no
time-dependent mesh steering or physical viscosity. The W=0,
fixed-mass RK4 scheme is otherwise unchanged in the source variant; the
HLLC label previously used for the main polar runs is superseded by their
[active-operator audit](KEPLER_POLAR_OPERATOR_IDENTITY_2026-10-27.md).

**Stability limit identified in the 2026-10-27 Astra IC review:** radial
force balance does not make the tapered transition a stable equilibrium.
Its specific angular momentum `j=r*v_phi` is positive at `r0` and zero at
`r1`, so `d(j²)/dr<0` somewhere in the transition. Hence the epicyclic
frequency squared `kappa²=r^(-3)*d(j²)/dr` is negative there, violating the
axisymmetric Rayleigh stability condition for this constant-pressure setup.
The taper is a seam-causality diagnostic only; a successful `t=10` endpoint
would not qualify it as a long-time Kepler benchmark boundary. Audit the
transition's `kappa²`, perturbation growth, and radial energy/torque before
interpreting its outer dynamics.

This obstruction is independent of the particular quintic taper. For a
continuous axisymmetric circular field centered in a periodic square,
opposite seam points have opposite tangential components, so matching them
to a static periodic exterior forces `v_phi` to vanish before the seam.
Any smooth profile with positive `j²` in the disk and `j²=0` at a larger
radius has `d(j²)/dr<0` somewhere by the mean-value theorem. Thus no
rotation-to-rest taper with the present constant-pressure, inviscid balance
can be Rayleigh-stable everywhere. Repeating a different taper shape is
not a route to the `t=600` benchmark; a long-time benchmark needs a
different outer-boundary formulation or a separately justified physical
exterior problem.

For the staged `r0=2.35`, `r1=2.75`, `epsilon=0.01175` parameters, evaluating
`j²=T²*r⁴/(r²+epsilon²)^(3/2)` and its analytic derivative on 100,001
uniform radial samples gives `kappa²<0` for approximately
`2.367952<r<2.75`. The minimum is `-0.982773` at `r=2.478104`, corresponding
to a local linear growth timescale `1/sqrt(-kappa²)=1.00873`. This is a
local stability indicator, not a prediction of the finite-domain nonlinear
growth rate. Because `t=10` extends far beyond this scale, pending long
taper job **436966 was canceled before execution** on 2026-10-27. The
preflight **436955** and `t=.825` seam-attribution job **436965** remain.

A standalone C check of the header using `epsilon=0.01175` and these radii
finds maximum absolute `|dPhi/dr-T² r/(r²+epsilon²)^(3/2)|=1.51e-10` across
representative radii, using a centered finite difference of width `1e-5`.
The isolated MPI+OpenMP+CUDA binary built successfully, with SHA-256
`6a735a7c196637b2a33830bc0e73719edc404c3d91f23ea70d51155b28213d6e`;
the [durable patch](KEPLER_PERIODIC_TAPER_NATIVE_PATCH_2026-10-26.patch) and
[build receipt](KEPLER_PERIODIC_TAPER_NATIVE_BUILD_2026-10-26.json) record
the exact change. A native preflight with the same 54,312-site narrow input
as queued job 436940 was submitted as **436955** to H200, requesting one H200
and 32 GiB. It runs to `t=.001` and audits the initial velocity, positivity,
fixed mass, W=0, and orbital accuracy. This checks only the early IC and
code path; same-IC native `t>.4`, `t>6`, and long-time comparison remain.
The running q=1/q=2 jobs and their queued continuations use the original,
untapered binary and inputs.

The audited continuation was originally staged independently: **436965** runs
to `t=.825` only after successful 436955 preflight. The `t=10` job **436966**
was canceled after the Rayleigh stability audit above. Both staged scripts
request one H200 and 32 GiB, use the same restart checkpoint and binary, and
keep the taper variables fixed. The canceled `t=10` audit would have
evaluated total energy with the **tapered** potential as well as the original
Hopkins disk annulus metrics. An additional
[radial audit](audit_kepler_native_radial_20261026.py) records the center,
disk, and outer boundary separately at each endpoint. These jobs are a
controlled attribution experiment; passing them alone does not prove
untapered `t=600` stability.

## Direct seam-geometry evidence at q1 native t=0.691 (2026-10-27)

Untapered full q1 job 436889 stopped at a GPU halo-sufficiency check,
`t=0.691008245`, with four cells whose conservative search radii crossed
the 8-cell (`0.5` length) periodic padding. Its preceding RK4 ledger had
positive internal energies and zero floor/dual-energy injection; the
inner `r<2` disk had no short generator pairs. Independent periodic
Voronoi reconstruction of the saved accepted states finds that the
largest `2R/(0.5+distance_to_physical_edge)` grows from `0.31615` at
t=.2063 to `0.59950` at t=.4125, `0.88180` at t=.6188, and `0.999876`
at the failed step start. Four symmetric outer cells are responsible;
their native saved volumes agree with the independent polygon areas.

The [wrapped exact-orbit control](KEPLER_Q1_BALLISTIC_PERIODIC_HALO_T069101_2026-10-27.json)
predicts **four** 8-cell bound violations at the **same four IDs** at
the same time, without any gas dynamics. This directly supports the
periodic seam geometry mechanism. It is not a valid physical exterior
orbit solution after wrapping. Wider padding can allow tessellation to
continue but does not remove the non-smooth torus initial state or the
large outer-cell deformation. The taper preflight 436955 and short native
continuation 436965 remain the distinct causal test of smooth outer forcing.
The full failure and
recovery contracts are in the [native control](KEPLER_POLAR_NATIVE_CONTROL_2026-10-26.md).

## Superseding sparse-mass audit (2026-10-27)

The original taper preflight **436955** finished and its conditional
continuation **436965** started, but they use the same pre-correction
54,312-site sparse glass and old CPU-initialized fixed masses. Direct
inspection of the taper `khout.000000.dat` confirms the initial Voronoi
areas sum to **30,469.7103055**, rather than the periodic 6×6 box area 36;
the fixed mass sums to **363.2684565**, rather than approximately 13.136.
The glass SHA-256 matches the invalid original narrow candidate. The old
preflight's orbit and positivity checks could not detect this error.

Running job **436965 was canceled** as soon as this was found. Preserve its
output solely for regression diagnosis; neither the taper preflight nor the
partial continuation is admissible evidence about boundary effects. Job
**436966** had already been canceled for the independent Rayleigh-unstable
taper problem. A future boundary control must use the corrected periodic
area sidecar and the new global mass/area gate, and its transition must be
physically stable if a long-time comparison is intended. The
[sparse-mass correction report](KEPLER_SPARSE_IC_MASS_CORRECTION_2026-10-27.md)
contains the valid replacement IC contract.
