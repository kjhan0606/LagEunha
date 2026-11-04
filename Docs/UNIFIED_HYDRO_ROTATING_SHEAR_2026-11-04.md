# Geometry reuse and rotating/shearing-flow validation

This experimental continuation of `9182686` uses the original native
C/OpenMP+MPI+CUDA program. Python orchestrates runs and analyses outputs; it
does not evolve the PDE. The run root is
`/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1`.

## Implementation

The 2D global KDK branch retains geometry after the preceding final kick only
when consecutive steps/times, boundary/domain parameters, and exact owned IDs,
positions, masses, volumes, and immutable weights match. Pressure, velocity fits,
and CFL are refreshed from current states. Drift still rebuilds geometry; cold
restarts rebuild initially. `GFS_UNIFIED_REUSE_GEOMETRY=0` is the reference path.

Three native A/B regressions have identical final x, y, mass, vx, vy, W, weight
momentum, thermal energy, and volume:

| Test | Endpoint | Steps | Reused steps | Native phase time, reference/reuse (s) |
|---|---:|---:|---:|---:|
| Moving diagonal contact, 64² | 0.5 | 134 | 133 | 2.203 / 1.387 |
| Axis shear, 64² | 0.25 | 25 | 24 | 0.750 / 0.671 |
| Axis Sod, 64² | 0.2 | 33 | 32 | 1.094 / 0.846 |

These small-run phase timings do not establish matched GIZMO performance.

`unified_initial2d.c` writes native 696-byte ICs for periodic Gresho (unit square)
and smooth KH (1 by 2). Mass and thermal energy are integrated over rectangular
cells; velocities are sampled at generators. KH uses analytic density averages.
Gresho pressure uses 8 by 8 Gaussian quadrature, with a 4 by 4 comparison.
Gamma matches the native float-promoted parameter. Unresolved velocity variance
is not added as thermal energy.

Optional diagnostics record initial pressure acceleration, velocity fit,
instantaneous viscous heat, sparse flow observables, and heat by original
material radius. Initial-only probes update no state and accept no steps.
Exact output clipping requires native KDK, exact final stopping, and a positive
output cadence. These diagnostics do not add a force or thermal correction.

## Results and remaining problems

Gresho 64² and 128² reached t=1. KH 64 by 128 and 128 by 256 reached t=2.
The coarse cases continue accepted t=0.05 preflight checkpoints rather than
restart their ICs. KH 256 by 512 is still running at the publication capture;
captured progress is not a completed t=2 result.

Gresho 64² continuations from the same accepted t=0.05 state give:

| Timestep factor | Relative energy change to t=1 | Volume-weighted vector velocity L1 at t=1 |
|---:|---:|---:|
| 1 | 1.1404e-3 | 0.07656 |
| 1/2 | 1.4954e-5 | 0.07171 |
| 1/4 | 4.0773e-6 | 0.07318 |

Smaller steps control energy error but leave velocity error and kinetic-energy
loss. At 128², default steps give relative energy change 2.7672e-3 and vector
velocity L1 0.03806 at t=1. Endpoint completion does not qualify vortex accuracy.
KH profile differences measure evolution, not error against an exact solution.
The material overlap width measures interleaving rather than molecular diffusion.

A frozen periodic probe evaluates the positive acoustic Hessian part
`M^(-1/2) B^T diag(gamma P/V) B M^(-1/2)`. It omits the pressure-weighted second
derivative of geometry and does not prove full-system stability. Native dt times
the estimated frequency increases from 0.547 at the initial 64² state to 2.001
at t=1; the 128² t=0.5 state gives 2.123. No new timestep restriction or force
correction has been adopted from this probe.

The current face-angular closure retains substantial initial Kepler core heat:
frozen t=0 rates are 3.17197e-3 in the core and 7.57019e-7 in the original disk.
These are instantaneous rates with zero accepted steps, not finite heat
increments or a new long-term Kepler validation.

## Provenance and publication scope

Full Gresho/KH runs use `smooth_native_build_v2`, binary SHA256
`4edfbf8e0e1f532c4fda545a21a72edf92746b38fdf812972c233fa3c20c4f4b`.
The later acoustic probe is
`9f4beb265da2b40e0805e2cd442af57170928d24273ec3388f52a89f95e6ddc0`.
Both are partial relinks, not portable clean builds. Frozen sources and build
manifests distinguish the executed versions; the later addition is opt-in and
read-only. Scripts, measurements, flow histories, source/build provenance, and
a dated partial-run capture are in `Research/RotatingShear_20261104`. Large
snapshots and binaries remain on GPFS with paths and hashes in the inventory.

H200 jobs 449988--449990 remained pending under `QOSGrpGRES`. Only these three
never-started jobs were cancelled before moving the cases to Ada. The unrelated
GIZMO job was preserved. Publication does not stop or restart the active KH run.

The manuscript remains the verified previous baseline (Overleaf `6cb0e1a`).
Experimental results have not replaced its earlier methods or Kepler trajectory.
Full-suite accuracy, new-method long Kepler, vortex error diagnosis, and matched
GIZMO <=2x wall-time qualification remain incomplete.
