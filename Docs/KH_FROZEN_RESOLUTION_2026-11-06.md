# Native fixed-mass KH resolution comparison

## Scope and execution

Executed actual C/CUDA native variational-pressure W=0 KDK trajectories on
Syntax, one NVIDIA A10, MPI rank 1, OpenMP threads 2. This is the smooth
fixed-mass baseline related to the successful Kepler force/volume-work path,
not the new 1D exact-Riemann/PLM material-boundary method and not the rejected
static-condensed weight extension. The parameter file contains historical
av_mode/MUSCL labels; runtime explicitly reports VARPRESS_W0, so those labels
must not be taken as evidence that a shock Riemann path was tested here.

New resolutions: 32x64 and 64x128, box 1x2, saved endpoints .1,.6,1,1.5.
Old 16x32 t=1.5 was reused without rerunning. Each segment resumes the previous
checkpoint. All outputs use /gpfs. No /scratch, new physical parameter,
steering, viscosity, or mass exchange was introduced. Weight is identically
zero, and the existing safety and geometric timestep tests remain active.

ICs are the existing smooth double shear layer: density 1/2, velocities
-1/+1, width .05, pressure 10, transverse seed amplitude .01. Mass and U
are cell integrals; velocities are generator samples. The compiler/native
gamma is the existing float-rounded 5/3. Do not compare against a differently
defined KH reference without matching these details.

## Results at t=1.5

| Nx x Ny | transverse kinetic energy | density range | relative total E error |
|---|---:|---:|---:|
|16 x 32 (existing)|0.06787485|0.87371--2.57165|not recomputed from t=0 here|
|32 x 64|0.05087331|0.84036--2.40565|+1.70577e-5|
|64 x 128|0.03695038|0.83155--2.25626|+3.01434e-5|

Both new runs reach all targets with unchanged individual masses, total
volume 2 to roundoff, positive U and V, W=0. Maximum relative cell entropy
change is 1.16e-14 and 1.74e-14 respectively. Entropy preservation follows
the smooth adiabatic U(V) update; it is not independent shock qualification.

The final .5-time segment needed 947 accepted steps/17 retries at32 and
2597/91 at64. Total native segment wall times are approximately44.90 and
142.28 seconds respectively, including per-segment startup. Neither is a
GIZMO comparison; the requested <=2x GIZMO criterion remains untested here.

Piecewise constant periodic Voronoi density fields are sampled on the SAME
384x768 pixel grid, with no spatial smoothing. Relative L1 differences
between16/32 and32/64 are8.14% and6.27%. These are inter-resolution differences,
not errors against an exact solution. Their decrease alone does not establish
a useful asymptotic convergence rate.

The transverse kinetic energy changes by25.05% from16 to32 and27.37% from32
to64 (relative to the coarser value). Small cell-scale structures are visibly
present in the shear layers. Total E error gets larger rather than smaller
with spatial resolution in this comparison. Therefore reaching t=1.5 and
conserving cell mass are insufficient for a KH accuracy pass. This evidence
also does not isolate a temporal error: timestep histories differ and no
fixed-resolution timestep refinement was performed.

## Decision and next discriminating test

Keep this as a measured baseline, not adoption of W=0 as the final common
hydro algorithm. Before attributing the transverse energy to physical KH
growth, compare a zero-transverse-seed shear trajectory with the same native
operator and geometry. That distinguishes spontaneous mesh/particle response
from the prescribed growing mode without adding a sensor or physical term.
Then a timestep refinement can separate time-integration error from spatial
force error. RT, common shock treatment, general2D weight closure, and later
3D validation remain unfinished. No3D computation was started.

Artifacts: Research/MaterialVolume_20261105/kh_frozen_resolution_20261106_v1/.
Each segment preserves params.dat, environment.json, lineage.json, log and
checkpoint. run.py invokes native executables; analyze.py only postprocesses.
results.json records execution/conservation, metrics.json records state
statistics, density_comparison.json records same-grid differences, and
density_t1p5.png uses common color limits .8--2.6 covering all three states.
