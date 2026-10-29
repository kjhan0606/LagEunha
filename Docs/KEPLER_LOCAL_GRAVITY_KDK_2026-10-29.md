# Local central-gravity KDK prototype — 2026-10-29

## Implemented

The opt-in `GFS_KEPLER_LOCAL_GRAVITY_KDK=1` path in
`Exam/exam.c::exam2d_vph_kdk_int_blend` keeps one globally synchronized
hydrodynamic macro-step. It splits that step as a hydro half-kick, independent
per-generator KDK subcycling in the softened central potential, and a hydro
half-kick. Each generator gets a dyadic substep no larger than
`GFS_KEPLER_LOCAL_GRAVITY_ETA * (r^2 + eps^2)^(3/4)`; the local eta defaults to
0.005 and is separate from `EUNHA_KEPLER_DTETA`, which belongs to the old
globally capped integrator. The global macro-step is reduced only if the
configured maximum bin depth (default 12) cannot cover a local gravity limit.
Both per-generator kick/update loops use OpenMP when this opt-in path is on;
the ordinary integrator path stays serial in those loops.

This does **not** make hydro or face fluxes asynchronous. The GPU hydro force,
Voronoi mesh, pressure work, and all face updates still run at the global
macro-step. The local path rejects XSPH, nonzero constant acceleration, weight
evolution, and centroid steering. Both global-gravity KDK and the opt-in local
gravity split retain the existing phase1 hydro-energy correction. A native
thermal64 path is wired for the one-rank CUDA Kepler W=0 contract: each hydro
half-kick advances authoritative double `U64` with the GPU `Udot64`, refreshes
the legacy float mirror, and re-derives EOS fields from `U64` and the current
double cell volume. Nonpositive or nonfinite `U64` aborts the step. This KDK
mode rejects RK4 transaction and thermal64 initialization-only flags.

The Kepler KDK entry point rejects the generic periodic post-stage callback
and the RK4-only `GFS_KEPLER_WALL_EVENTS` option. The bounded 1024-gon face
clipping and fixed-polygon post-stage callback are now connected to the KH
driver for EVOLMETHOD=1 or 3. A KDK cell that crosses the wall still aborts;
KDK wall-event localization and impulse accounting are not implemented. This
permits an interior bounded-domain comparison and prevents silent periodic
wrapping, but it does not yet qualify trajectories through wall impacts. The
default integrator path is unchanged when the local KDK flag is absent or zero.

The Kepler parameters are parsed and collectively checked once, then copied
from static cached configuration into the per-call locals on every macro-step.
This is required because the integrator function is re-entered each step.
All ranks first verify that the opt-in flag agrees, then check exact agreement
on eta, center, softening, and maximum bin level. Inconsistent environments
abort at startup rather than letting ranks use different forces or diverge at
a feature-specific collective.

## Evidence

`Exam/test_gfs_timebins.c` passes the synthetic 2:1 neighbor-closure and
power-of-two event-schedule checks. `Exam/test_gfs_kepler_kdk.c` integrates a
single circular orbit in the same Plummer potential for 100 periods. Halving
eta from 0.01 to 0.005 reduced endpoint position error from
`7.878e-3` to `1.969e-3` (ratio 4.000); the fine run's relative energy and
angular-momentum errors were `8.148e-14` and `5.663e-15`.

CUDA runtime evidence now exists. Commit 7209164 built on syntax with CUDA 13.0.2 and was exercised on the Ada RTX 5000. The global- and local-gravity KDK paths both stopped at exact t=0.826 from the accepted t=0.825 checkpoint; both passed the 29-check bounded-Q1 audit. Their total-energy changes were about -3.5e-7, with exact masses, W=0, positive volume/internal energy, and no floor or clip injection.

The C helper tests still pass with -std=c11 -O2 -Wall -Wextra -Werror, and
mpiicx -fsyntax-only passes for host and CUDA host-side configurations. The
CUDA objects and executable were actually built with CUDA 13.0.2 on syntax.
One pre-existing format warning remains at Exam/exam.c:8298 (ptrdiff_t passed
to %d). The integrator's MPI startup consistency path has not been exercised
with more than one rank.

### Same-time t=3 comparison

The global-gravity KDK run reached exact t=3 in 895 s but failed one of 29 checks: cumulative KDK floor injection was 4.72e-9 and 17 internal-energy reset events were counted. Its disk gates passed (0.0498% mass-weighted density change; total-energy change -3.07e-7), while its r<0.1 core had 8.63% density change, 0.968 absolute log-entropy change, and 0.00487 mean orbit-position error.

The local-gravity KDK run, with eta 0.005 and maximum level 12, reached the same exact t=3 checkpoint in 616 s and passed all 29 checks. It had zero clip, floor, or pair-pressure events through t=3. Disk density change was 0.0498%; total-energy change was -3.54e-7. Core errors were smaller at this time: 5.03% density, 0.0200 absolute log-entropy, and 0.00185 orbit-position error.

### Local-gravity KDK at t=6.5

The accepted local t=3 checkpoint khout.004447.dat was continued to exact t=6.5 as khout.010970.dat in 1313 s. All 29 bounded-Q1 disk/state checks passed. The run preserved exact masses and W=0 and reported zero negative cells, internal-energy resets, clip injection, and floor injection. Disk density change was 0.0557%, absolute log-entropy change 5.21e-5, orbit-position error 1.05e-6, angular-momentum change -3.31e-12, and total-energy change -3.54e-7.

The inner core remains inaccurate: for original radius r<0.1, density change is 38.4%, absolute log-entropy change 12.5%, and orbit-position error 1.71%. These are comparable to or worse than the existing bounded RK4 t=6.5 result (33.4% core density and 1.70% position error). Local central-gravity subcycling prevents the early KDK floor-reset failure and keeps the disk stable through t=6.5, but it does not resolve the core error.

### Local-gravity KDK at t=10

The accepted local-KDK t=6.5 checkpoint `local/t65/khout.010970.dat` continued to exact t=10 as `local/t10/khout.020441.dat`. The process exited 0 after 2053.356 s. The KDK gate at `local/t10/kdk_gate.json` passes all 29 checks: exact mass and W=0, positive internal energy and volume, no floor/clip injection, bounded polygon geometry, and all configured disk and GIZMO-annulus thresholds. The final ledger is zero for clip, floor, and dual-energy events. Minimum cell volume is positive (`1.636e-5`).

Against the bounded RK4 t=10 checkpoint from the same initial IC, the disk remains effectively unchanged in either method: mass-weighted density change is 0.057431% for local KDK and 0.057430% for RK4. KDK has slightly smaller total-energy drift (`-3.545e-7` versus `-2.229e-6`), while both conserve disk angular momentum to about `1e-12` relative.

The core comparison is unfavorable to KDK. For original radius r<0.1, KDK has 76.8% mass-weighted density change, 0.304 absolute log-entropy change, and 5.36% orbit-position error. The matched RK4 t=10 values are 56.3%, 0.00369, and 5.07%. Thus passing the disk/state gate does not mean that the core is accurate: local-gravity KDK does not solve the late-time core error and is worse than RK4 in these core diagnostics. Keep this implementation experimental; do not promote it to the default integrator or extend it to longer production-like durations on this evidence.

The t=3 global/local full-density and core-zoom comparison uses a shared linear density scale [0,1.05]: /home/kjhan/.codex/visualizations/2026/09/29/01a0eb12-3d2a-7402-8478-d6fa000154a2/kepler_q1_kdk_t3_density_shared.png.

One-rank CUDA+OpenMP runtime was about 10 minutes for t=0.825 to 3 and 22 minutes for t=3 to 6.5 on the Ada RTX 5000. This is not MPI scaling evidence or a t=600 performance qualification.

## Relation to intrinsic cell spin

`GFS_LAGUERRE_ROTATION` in `Exam/exam.c` corrects the **geometric normal
velocity of a moving face**. It is not a transported cell-spin degree of
freedom. The standalone spin reconstruction and face-exchange prototypes in
`Docs/` are not connected to the native integrator, and this local-gravity KDK
prototype adds no spin state.

The saved-state moment audit
(`Docs/KEPLER_Q2_MATERIAL_MOMENT_GAP_2026-10-29.md`) finds that the uniform
polygon geometric covariance rate and affine material covariance rate differ
by about `1.71 Omega` in mass-weighted mean. A freely fitted scalar weight-rate
field leaves `99.6–99.8%` of that mismatch. Separately, the existing native
generator force does not select a unique cell torque. These results rule out
simply adding a spin variable or changing weights as a justified closure.
Before native spin evolution, the model must define a consistent mass
centroid, inertia, face momentum traction/torque, and total-energy transfer,
while retaining the W=0 fixed-mass and `Re=infinity` limits. See also the
mass-centroid conditions in
`Docs/KEPLER_SPIN_CENTROID_CLOSURE_2026-10-27.md` and torque ambiguity in
`Docs/KEPLER_VARIATIONAL_SPIN_PAIR_TORQUE_2026-10-27/README.md`.

## Next validation

The t=10 comparison identifies the inner core—not the disk energy or angular momentum—as the unresolved error. The next useful study should separate inner-resolution effects from timestep effects using matched initial states and explicit core thresholds. Keep the local-gravity KDK path isolated while doing that; the current result gives no basis to add further KDK physics patches or to infer that intrinsic spin is needed. The material-moment and torque-closure objections above still apply. Do not promote any path to t=120/600 until core error, radius/shape sensitivity, and multi-rank behavior are bounded.

KDK still has no wall-event localization or wall impulse/energy ledger, so any wall crossing aborts the run. The current trajectory has remained inside the fixed polygon. Full local hydrodynamic stepping still requires conservative face-flux registers and MPI halo synchronization; the present local KDK only subcycles the central gravity per generator while hydro and face fluxes retain one global macro-step. No multi-rank runtime test has been made.
