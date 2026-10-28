# Kepler–GIZMO comparison contract (2026-10-26)

Source: Hopkins (2015), *A New Class of Accurate, Mesh-Free Hydrodynamic Simulation Methods*, §4.2.4 and Figs. 6–8: https://arxiv.org/pdf/1409.7395

## Published target

- Two-dimensional disk, gamma = 5/3, constant initial pressure P = 1e-6, no gas self-gravity, and an analytic softened external potential Phi = -(r^2 + epsilon^2)^(-1/2).
- Published surface density is Sigma = 0.01 plus a piecewise profile: (r/0.5)^3 for r < 0.5, 1 for 0.5 <= r <= 2, and [1 + (r-2)/0.1]^-3 for r > 2. Epsilon is not numerically specified in the cited section.
- Effective 256^2 initial resolution across -2 < x,y < 2. Figure 6 shows the surface density at t = 120. Figure 7 shows the MFM surface density at t = 600 and says the inner r ~ 0.5 disk has completed >250 orbits without disruption.
- The nearby prose calls t = 600 about 160 inner-disk orbits, while the Fig. 7 caption says >250 at r ~ 0.5. State the reference radius whenever quoting an orbit count. The prose also reports evolution as far as about 1000 orbits, without a corresponding plotted snapshot or an unambiguous reference radius; the final goal must retain this longer claim for explicit interpretation and testing.
- Figure 8 reports annular L1 velocity and density diagnostics at R ~ 1. Its plotted threshold times use L1(v) = 0.01 and L1(rho) = 0.25. These are comparison diagnostics, not an automatic pass by themselves.

### Local continuum equilibrium versus the discrete periodic experiment

The current Hopkins-profile source sets constant `P=1e-6`,
`rho=0.01+profile(r)`, and
`v_phi=r(r^2+epsilon^2)^(-3/4)` in the potential
`Phi=-(r^2+epsilon^2)^(-1/2)`; see
[`Exam/KH/util.c`](KEPLER_NATIVE_SOURCE_OVERLAY_2026-10-27/Exam/KH/util.c)
and [`Exam/exam.c`](KEPLER_NATIVE_SOURCE_OVERLAY_2026-10-27/Exam/exam.c).
For a stationary axisymmetric inviscid flow,
`v_phi^2/r = (1/rho)dP/dr + dPhi/dr`. Here `dP/dr=0` and
`r*dPhi/dr = r^2/(r^2+epsilon^2)^(3/2) = v_phi^2`.
The continuity equation is also stationary for any radial `rho(r)` under
purely azimuthal velocity. Thus the specified smooth fields are an exact
*local continuum* equilibrium for this chosen softening and inviscid Euler
model. The radial density profile by itself does not require an extra
pressure-supported correction to the circular speed.

This local statement does not make the 6×6 periodic-box problem an exact
global equilibrium: the identified opposite faces have mismatched velocity
and potential gradient. Nor does point-sampled `m_i=rho(x_i)V_i`, pressure,
and velocity make every finite Voronoi cell an exact discrete balance. The
polar and corrected sparse layouts are valid discrete samples of the local
equilibrium, but their core resolution, face geometry, timestep and seam
can all introduce evolution. Consequently, native long-time results must be
judged against the same IC and boundary contract before attributing an
improvement to the solver.

## Historical C010 experiment

- Authoritative job 436489 used `EUNHA_KEPLER_RIN=2`, `ROUT=8`, `EPS=0.047`, center (12,12), fixed per-cell mass, W = 0, Re = 0, RK4, one H200, and 71,080 cells. Its accepted original-IC snapshot is `/gpfs/kjhan/LagForce/Kepler_GPUStreamCandidates_H200_20261020/core_refinement/fine_integrated/khout.000001.dat`.
- Direct reading of that t = 0 snapshot finds constant P ~ 1e-6 but a tapered approximately 2/r disk: median density 0.623 for 2 <= r < 4 and 0.317 for 4 <= r < 8. The audit code also identifies the run as a historical tapered 1/r disk. This is **not** the published constant-density annulus. The current density-error metric is relative to this own initial field, not the published GIZMO density norm.
- In its softened point potential, the inner-edge orbital period is T(2) = 2*pi*(2^2 + 0.047^2)^(3/4) = 17.778891990936625. Thus t = 6 is 0.33748 inner-edge orbits, and t = 10 is 0.56246. These are early numerical-stability gates, not long-term disk tests in the paper's sense.
- For orbital-count comparison only, the unsoftened published T(0.5) = 2.221441469079183. Its plotted t = 600 corresponds to about 270.1 inner-edge orbits. The same count in C010's units occurs near t = 4802. Current t = 600 would be only 33.75 inner-edge orbits. If the text's ~1000 orbits refers to the inner edge, C010 would need t ~ 17779. These time mappings do **not** make the two initial-value problems physically identical.
- A genuinely matched GIZMO run needs the published density profile, comparable dimensionless pressure and softening, effective resolution, and density-field outputs. Scaling all radii by 4 at fixed GM changes orbital speed by 1/2 and period by 8; with the same surface-density normalization, the pressure needed to preserve P/(Sigma*v_c^2) scales by 1/4. C010 instead retained P = 1e-6 and has a different Sigma(r).
- The current native binary already contains `EUNHA_KEPLER_PROFILE=hopkins`, which implements the published piecewise profile. To avoid changing the pressure normalization, the new independent preflight uses the paper's original radii (0.5, 2), center (3, 3), box width 6, and 384 cells per side, giving 256 cells across the paper's width-4 effective-resolution domain. It uses epsilon = 0.01175 (same epsilon/Rin as C010); this is an explicit choice because the paper does not specify a numerical epsilon. The later polar q=1/q=2 and profile-sampled experiments use this contract; their full native long-time outcomes remain pending.

The still older 256² HLLC-era jobs **406734** and **407021** are a separate
initial-value and operator comparison, not native evidence for the new W=0
variational polar branch. Their [original job notes](../code_review_grok.md)
identify the radius-2-to-8 disk. Direct re-reading of their saved logs at
`t≈8` gives `floor_cum=5.423e-3` in 406734 and `de_cum=6.017e-3` in 407021;
both have total-energy drift around `7.5e-4`. Job 407021 has
`E_int=1.057826e-2`, versus initial `8.644e-4`, so its roughly 12-fold
thermal rise includes at least `6.017e-3` of explicit dual-energy injection.
The notes locate the overwhelming majority of reset events inside the old
disk's radius-2 hole through `t=7.94`. These observations do not establish
that Kepler shear in the orbiting ring created the injected energy. The new
q1/q2 runs must keep floor/dual/entropy-floor cumulative injection at zero
and compare **signed as well as absolute entropy by initial radius** before
assigning any late thermal change to disk shear; the existing strict and
radial checkpoint auditors record those quantities.

## Operational decision

1. C010 job 436489 ended `FAILED` at t=6.598925 after a valid t=6.500079 checkpoint because a low-density inner-atmosphere cell acquired negative internal energy. Preserve its outputs as historical evidence; an exact t=10 C010 result was never obtained. Its first-annulus density change at t=6.5 is almost entirely reproduced by exact-orbit Voronoi geometry. See [the shell-spacing audit](KEPLER_SHELL_SPACING_DIAG_2026-10-26.md).
2. Continue the native paper-profile polar q=1/q=2/q=2c controls. The q=1 job 436889 is running toward t=.825; q=2 and q=2c jobs are queued. Its t=.206 and t=.4125 checkpoints have positive state and fixed mass, with disk density change dominated by exact-orbit geometry and no demonstrated net heating. The t=.4125 short-range pair activity is confined to the outer background near r=3. See [the native polar control](KEPLER_POLAR_NATIVE_CONTROL_2026-10-26.md). Only later native data can establish whether the polar disk survives the old t>6 failure interval.
3. Use the isolated 54,312-site profile-sampled native preflight 436940 and its validated dependent continuations 436941 (t=.825) and 436942 (t=10) to measure whether removing unnecessary low-density sites makes the same disk tractable. Its disk generators and initial Voronoi areas match the q=2 reference, but the 25-fold orbital-step-bound improvement and 2.715-fold site-count reduction are not measured wall-speedups. See [the step-scaling audit](KEPLER_TIME_STEP_SCALING_2026-10-26.md).
4. After native t>6 accuracy and runtime gates pass, compare fixed-mass W=0 density fields and annular L1 diagnostics with the published t=120 and t=600 figures at equal physical time in the matched-radius problem. Retain the paper's further approximately 1000-orbit prose claim as a separate target requiring an explicit reference radius. If heating or angular-momentum loss remains, diagnose face work and torque before enabling Laguerre weights or cell spin; see [the spin-status audit](KEPLER_CELL_SPIN_STATUS_2026-10-26.md). Re=infinity remains the Kepler physical setting.

## Periodic boundary caveat

The initial centered Kepler velocity has a measured tangential jump across
each periodically identified boundary, with median jump/local sound speed
about 77 in both full polar ICs. The central potential is continuous there,
but its normal derivative jumps. The q=1 `t=.4125` outer background shows
large density and pressure changes and all short pairs lie near an edge.
These facts do not prove later disk contamination, but make the 6×6 periodic
domain a distinct global initial-value problem. The
[seam-control note](KEPLER_PERIODIC_SEAM_CONTROL_2026-10-26.md) derives the
mismatch and records an isolated, balanced outer-potential taper candidate
for attribution. This control leaves the disk IC unchanged but changes
outer-background physics; it must not be presented as the original Hopkins
problem or as a solver improvement.

The 54k/56k profile-sampled ICs also have an unresolved inner-source
discretization: the innermost eight large Voronoi cells have area-averaged
Plummer acceleration about 1.98 times the generator-centered acceleration
used by the native code. This is measured on the initial geometry in the
[sparse-core force audit](KEPLER_Q2_SPARSE_CORE_FORCE_2026-10-26.json), not a
native failure result. Core and seam effects must be distinguished from
disk shear heating in subsequent outputs.

## Native polar status update (2026-10-27)

The earlier operational paragraph above is historical. Both full q1 and q2
now have accepted, complete, fixed-mass W=0 native checkpoints at exact
`t=.825`, with all 15 strict checks passing. The q2 `t=.825` run finished
after a verified halo-12 recovery and shows a 4.49-fold smaller mean absolute
disk density change than q1; independent exact-orbit tessellation predicts
the same geometry-driven difference. See the
[same-time result](KEPLER_Q1_Q2_NATIVE_T0825_RESULT_2026-10-27.md).
The active operator is [variational W0 pressure](KEPLER_POLAR_OPERATOR_IDENTITY_2026-10-27.md)
under RK4; HLLC star pressure is computed but bypassed for pressure force and
work. Neither `.825` result establishes `t>6` stability or GIZMO's
`t=600`/longest-orbit target. The accepted long continuations and 5% site
jitter control remain separate native gates.
