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

After correcting the one-time configuration cache, both helper tests pass
with `-std=c11 -O2 -Wall -Wextra -Werror`; the 100-period KDK convergence
numbers above were reproduced. `mpiicx -fsyntax-only` passes for both the
host and `USE_CUDA` host-side configurations after the thermal64 connection.
The compiler reports one pre-existing format warning at `Exam/exam.c:8298`
(`ptrdiff_t` passed to `%d`). No CUDA object was built and no native hydro
simulation has exercised this path; the MPI startup consistency guard has
only been syntax-checked, not run multi-rank. These ODE checks do not establish
disk stability or production speed. This `grammar` shell has neither the Makefile's
`/opt/ohpc/pub/cuda/13.0.2/bin/nvcc` nor a loadable CUDA module, so a device
build must be done in the configured CUDA environment; no build job was
submitted here.

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

The Ada continuation recorded on 2026-10-29 remains an RK4 run and must not be
modified or restarted with this prototype. After the live RK4 continuation
finishes, use a separate run directory and an EVOLMETHOD=3 binary built in a
CUDA-enabled environment. First compare one global-KDK and one local-gravity-
KDK macro-step from the same accepted checkpoint, with RK4-only wall events
disabled; abort-on-crossing keeps this initial comparison bounded. Audit exact
masses and W=0, double `U64`/volume validity, density/entropy, orbit phase,
total energy including the Plummer potential, angular momentum, positive
internal energy, bounded-domain membership, and actual macro/substep counts.
Before extending through any wall impact or toward t>6, implement and validate
KDK event localization and its wall impulse/energy ledger. Full local
hydrodynamic stepping still needs conservative face-flux registers and MPI
halo synchronization; this prototype does not solve that part.
