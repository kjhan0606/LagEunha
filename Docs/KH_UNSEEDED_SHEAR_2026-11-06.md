# Zero-seed native KH shear diagnostic

Follow-up to KH_FROZEN_RESOLUTION_2026-11-06.md. Same native C/CUDA W=0
fixed-mass variational-pressure KDK, Syntax A10, MPI1/OMP2. No physics solver
change. The native IC generator now accepts an optional final KH seed
amplitude after frozen/dynamic. Default .01 is unchanged: regenerated
32x64 and64x128 default checkpoints are byte-identical to their predecessors.
Invalid/nonfinite/negative amplitudes and use with Gresho are rejected.

The comparison sets only transverse seed amplitude to zero. The corresponding
parameter file also says zero. IC density and vx depend on y alone, vy=0,
P=10, no gravity. This is a stationary solution of the continuous Euler
equations: divergence and advective derivatives vanish and pressure is uniform.
Roundoff perturbations remain possible, so instability of the numerical
trajectory is not by itself an exact linear stability calculation of Euler.

## Actual evolution

Both resolutions reached .1,.6,1,1.5 with individually fixed masses, W=0,
positive cell volumes/internal energies. No clipping or viscosity was added.

| Nx | t | transverse kinetic energy | max abs(vy) | relative pressure RMS |
|---|---:|---:|---:|---:|
|32|1.0|1.42894e-17|2.13592e-8|2.47616e-9|
|32|1.5|1.06850e-9|1.16337e-4|3.16766e-5|
|64|0.6|2.08418e-18|1.29589e-8|4.61473e-10|
|64|1.0|5.20366e-9|5.47355e-4|2.44902e-5|
|64|1.5|1.56337e-2|0.76010|0.04973|

At64,t=1.5 the unseeded transverse energy is42.31% of the seeded result.
This is a comparison of separate nonlinear trajectories, NOT an additive
decomposition assigning42.31% of the seeded energy to error. Final relative
total energy errors are1.68e-13 at32 and2.21e-5 at64. Conserved cell entropy
does not prevent these large pressure/velocity errors.

## Spatial content before large distortion

Direct Fourier projection of vy along x, separately within each y row, then
sum squared coefficients over rows. Confirmed exactly Nx cells per row at
the sampled early times. At32,t=1 and1.5 the largest component is k=14;
at64,t=.6 and1 it is k=32 (two-cell wavelength). More than99.998% of the
nonzero-mode projected power is in k>=Nx/4. The imposed physical seeded run
has initial k=1, so this disturbance is distinct from that chosen seed.
These are direct projections at the actual positions, not exact Parseval
energy fractions on a deformed lattice. No continuum growth-rate comparison
has been made. The immediate useful diagnostic is the grid-scale numerical
response, not another long KH trajectory or a viscosity sensor.

## Time-discretization control

A comparison resumes the IDENTICAL saved64,t=1 checkpoint, reducing Courant
from.4 to.1 and geometric CFL coefficient from.1 to.025, preserving the same
forces and U(V). It saves1.25 and1.5. This tests stricter time limits rather
than assuming each accepted step is exactly one quarter of the old step.
Actual timestep distributions and terminal metrics are recorded in
refined_results.json and refined_metrics.json in the artifact directory.
This cannot remove perturbations already in the shared t=1 state; it tests
their subsequent amplification with finer temporal resolution.

The refinement completed both targets. At t=1.5, Ky is0.01484452 versus
0.01563373 in the original (5.05% smaller), max abs(vy) is0.77753 versus
0.76010, and relative pressure RMS is0.050021 versus0.049730. The energy
change from the shared t=1 state is -5.13e-8, greatly smaller than the old
roughly +2.21e-5, yet the velocity/pressure disturbance remains large.
The old1--1.5 segment median dt is1.91366e-4; refined1.25--1.5 median is
5.58712e-5. Different subintervals mean this ratio is descriptive, not an
exact matched-time refinement factor. First-step maxima on the common
t=1 restart are1.53093e-3 and3.82732e-4, exactly a factor4.

Thus reducing temporal error substantially improves energy accounting but
does not remove this grid-scale shear response. A full time-convergence
study was not done, so do not claim KDK is completely exonerated. The next
useful task is the spatial pressure/geometry response of the short-wavelength
mode on the saved t=1 state; adding a viscosity switch or another long run
would not yet isolate it. This materially changes the previous interpretation:
the seeded density images alone cannot establish physical KH growth.

Artifacts: Research/MaterialVolume_20261105/kh_unseeded_20261106_v1/.
preflight.json, build_command.json/build.log, run.py, results.json, metrics.json,
spectra.json and per-segment lineage/environment/params/log/checkpoint files.
Python runs the native binary and postprocesses; no Python replacement hydro.
No3D run started; no final common algorithm or GIZMO runtime pass claimed.
