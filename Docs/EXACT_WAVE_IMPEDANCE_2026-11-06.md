# Exact wave impedances and a near-vacuum reference-solver fix

## Common nonlinear face ingredient

Added Exam/gfs_exact_wave_impedance.h, a shared host/device helper. Given a
nonnegative star pressure, it computes the positive secant impedance of the
existing exact ideal-gas wave curve. For a nonvacuum exact fan,

    ZL = (pstar-pL)/fL(pstar),
    ZR = (pstar-pR)/fR(pstar).

The shock branch is evaluated analytically without subtractive division.
The rarefaction branch uses log1p/expm1 near equal pressure and the analytic
rho*c limit at equality. No empirical coefficient or pressure floor is added.
These impedances recover the exact star values through

    w = (ZL uL + ZR uR + pL-pR)/(ZL+ZR),
    pstar = (ZR pL + ZL pR + ZL ZR (uL-uR))/(ZL+ZR).

Thus the previous weak-acoustic face algebra has a nonlinear exact-wave
ingredient. Positivity of the impedances preserves its squared-residual heat
identity. This is not yet a proof of equivalence of the constrained multidimensional
adjoint stencil to the exact1D material solver or a validated shock trajectory.

Vacuum fans have two free edges instead of one shared contact speed. The
helper can evaluate a one-sided wave curve at pstar=0, but callers MUST NOT
infer that it supplies a shared vacuum face. The audit explicitly excludes
vacuum fans from the common-speed claim; no floor/fallback hides them.

## Bug found in the existing reference solver

The initial10000-state face audit failed: maximum scaled speed/pressure
relations differed by0.00216/0.00219 although the heat identity held. The
existing gfs_exact_star stopped when

    hi-lo <= 4e-14 max(1,hi).

That1-pressure-unit floor is dimensionful and accepts an inaccurate pressure
bracket near vacuum. Recorded case3106 stopped at pstar=2.3675e-14 with wave
residual8.858e-4; its corrected star pressure is6.7489e-16. Case4738 had wave
residual0.00663. This is a reference-solver stopping defect, not evidence that
the new secant identity or the native variational KH force caused it.

Changed only the pressure bracket scale in Exam/gfs_exact_riemann_ale.h to
max(DBL_MIN,hi). The velocity-wave residual criterion, Newton/bisection
procedure, iteration cap64 and physical equations are unchanged. The previous
header, failed outputs and example states are retained in the artifact folder.

## Verification

Same deterministic10000 cases, densities1e-2--1e2, pressures1e-6--1e3,
velocities scaled to sound speed, plus explicit uniform/blast/vacuum cases;
all normal faces are2D with normal(.6,.8). Results after the fix:

* 9649 nonvacuum fans checked;351 vacuum fans explicitly excluded;zero exact
  solver failures.
* Maximum scaled recovered speed error1.48e-14, pressure error1.33e-13.
* Maximum heat-identity error3.99e-16; normalized material mass flux6.80e-16;
  normalized material energy-work error2.11e-12.
* Scaling both pressure and density by1e-12 and1e12 in the recorded near-vacuum
  case changes normalized star pressure by at most5.56e-16 and leaves speed
  unchanged. This specifically tests removal of the pressure-unit floor.
* The same19298 one-sided wave inputs evaluated by the C helper and an A10
  CUDA kernel differ by at most4.45e-16 relatively.

Ten initial interfaces from the existing Sod/Lax/collision/Blast/Shu-Osher/
Noh/double-fan/contact definitions are compared before/after in regression.json.
This is an interface regression, not a repeat of complete1D trajectories.
Previously saved figures/results are retained and are not retroactively
reclassified. The native KH evolution binary uses its variational path and
was not rebuilt or rerun as a consequence of this reference-face fix.

## Limits and next integration work

Only the secant helper was evaluated on CUDA; the full iterative exact Riemann
solver is still a CPU reference here. No production C/MPI/CUDA evolution
path enables the new helper. Its coupling to the reconstructed force/work
operator, geometric weight constraints, finite-step positivity and remaining
cold-shear heating still needs verification. Vacuum handling and end-to-end
cost remain open. No3D simulation was started and no <=2x GIZMO claim is made.

Artifact: Research/MaterialVolume_20261105/exact_wave_impedance_20261106_v1/.
It includes original/fixed audits, native C and CUDA tests, unit-scaling and
interface-regression outputs. This is a real source fix and a shared numerical
ingredient, not completion of the common hydrodynamic method.
