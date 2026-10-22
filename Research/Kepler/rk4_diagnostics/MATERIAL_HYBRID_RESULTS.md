# Native material hybrid: face ledgers close; source curvature remains

2026-10-21. Implements `ASTRA_MATERIAL_HYBRID_PLAN.md` on the preserved actual
24-cell failure/neighbor patch. This is an initial smooth operator bundle,
not Kepler time evolution, a new long run, a restart, or a shock result.

## Implementation and trace-space correction

Canonical LagEunha source identity remains HEAD8e25dfe plus existing changes.
Native reference `test_material_hybrid_v1.cpp` reuses the actual polygon/IC
assembly through a guarded hook in `test_material_shared_v1.cpp`.
`Exam/gfs_material_trace.h` supplies C/CUDA positive trace Gram algebra. The
new reference does not replace the user's native OpenMP/MPI/CUDA runtime.

Each cell has five local P1 nodes plus two local velocity-residual coefficients:
168 broken scalar DOFs in total. Actual shared physical faces, not component
membership, impose the constraints. The24-cell patch has28 shared cell faces,
104 raw coefficient trace rows and83 independent rows. Hence the true
conforming kernel has **85**, not the old83, scalar DOFs. On axis-radial
faces one velocity-residual trace is identically zero; tying that coefficient
across the face imposed an extra condition. The native Gram tests verify the
retained rank3/rank4 trace spaces. Non-axis radial faces explicitly reject this
fixture-specific incidence implementation, pending a general trace-rank solver.

The old83-dimensional space is a proper subset of the85-dimensional kernel.
There is no coupling across disconnected components, and cell-center bubbles
are not constrained. The symmetric initial sources happen not to excite the
extra two modes: accelerations agree with the old solve to about1.2e-12.
That agreement is not imposed and is not promised for arbitrary later states.

## Independent solve and actual face traction

Solve block-cell M a=F+B^T lambda, B a=0. F is assembled separately as volume
stress/gravity plus ONLY external pressure traction. Neither old residuals nor
cell energy errors are inputs to this solve. Redundant constraints use a
deterministic canonical-face incidence forest, selected before force assembly.
All nonselected multipliers are zero. The Schur matrix is built as a sum of
positive local Gram forms, not by arbitrary symmetrization.

The83-by83 reference Schur condition number is756.395. Maximum face Gram
condition is24742.3. The small dense reference is not a scalable implementation
for the full71080-cell runtime.

Convert each multiplier to an actual reference traction using the L2-dual
trace basis, G eta=mu. Direct quadrature then checks that its cell load equals
B^T lambda. At finite geometry, use t_physical=t_reference/(dS/dS0), with the
same shared position trace. h=1e-5 is a geometry/Taylor control, not a timestep
of a claimed hydrodynamical integrator.

## Results through v5 (MPI1/4; v1 coarse, v2-5 finer quadrature)

All computations ran natively on authorized lageunha CPU with Intel MPI,
OMP_NUM_THREADS=1, fork bootstrap. Four cases preserve the same rho/v/P and
geometry: gravity off, representable linear gravity, actual Kepler gravity,
and actual gravity with uniform velocity boost(.013,-.007).

| Quantity | Evidence |
|---|---:|
| momentum equation relative residual | <=1.5e-18 |
| direct dual-traction load relative residual, actual gravity | <=8.1e-16 |
| direct dual-traction load relative residual, gravity off | <=1.2e-12 |
| actual-gravity cell numerical energy residual / U | <=2.5e-10/time |
| actual-gravity cell numerical angular residual / 2K | <=4.1e-15 |
| actual-gravity shared acceleration jump | <=5.7e-14 |
| initial pressure derivative / P | <=1.7e-13/time |
| min finite-map volume Jacobian | .9999999793 |

Actual-gravity MPI1/4 broken acceleration differences are below6e-13 and
sampled traction differences below4.5e-16. Increasing quadrature32x20 to64x28
changes acceleration coefficients by at most3.3e-12. Raw sources, binaries,
states, face samples and logs are retained in `material_hybrid_v1` through
`material_hybrid_v5`. `AUDIT_COMPLETE` denotes terminal execution, not a
physics pass.

The original local load residual has kernel projection below5e-16 of the load
scale for actual gravity. Bubbles annihilate at the reported roundoff level.
The independent solve and direct physical-face construction therefore give
the previously implicit intercell transfer an explicit discrete face ledger
on this initial fixture. The thermal equation is unchanged; no residual is
added to Udot.

## Physical pressure and approximation error are NOT removed

The cell ledger using only physical -P n still differs by15.7337 U/time
(15.9180 with the boost), with angular differences about9.4e-5 of2K. Independent
volume integration of rho v dot(a_h-g) reproduces the energy discrepancy.
Thus it remains a source-projection/spatial error; it is not actual heating
or a correction supplied by the hybrid formulation.

The selected multiplier gauge is not a unique pointwise stress. Its traction
L2 difference from -P n is44.93 P even in the zero-gravity control, where the
cell physical ledger and acceleration are essentially exact. There are21
redundant trace rows; the fixed independent-row gauge can retain self-balanced
pointwise components invisible to admissible virtual work. We do not tune
the gauge toward pressure, energy or torque. The direct acceleration/source
errors, rather than that gauge-dependent pointwise norm alone, determine
approximation quality.

For actual gravity, measured cell RMS acceleration error reaches about0.22%
of the exact gravity RMS. This is small compared with gravity but not a cold
thermal accuracy certificate.

## New thermal-curvature gate: first derivative alone missed this

At initial div(v)=0, differentiating the material map and adiabatic Euler
pressure equation gives

    D(div v)/Dt = div(a) - tr[(grad v)^2],
    (D^2 P/Dt^2)/P = -gamma {div(a) - tr[(grad v)^2]}.

For the exact initial softened Kepler field a=g, the two terms cancel. Use
the reconstructed native a_h and its actual basis derivatives instead. The
source projection leaves a nonzero predicted pressure curvature:

- maximum absolute relative curvature error: **1216.89/time^2**;
- maximum per-cell volume L1 curvature error: **255.41/time^2**;
- exact Kepler curvature residual from floating-point evaluation: <1e-9/time^2;
- gravity-off and representable-linear controls: curvature error <3.1e-10/time^2.

These are INITIAL OPERATOR TAYLOR COEFFICIENTS, not measured evolved pressure,
an entropy production rate, or a prediction that can be extrapolated to t10.
They show that Pdot initially zero does not establish cold-flow accuracy of
this representation. v5 independently reconstructs the two-sided maps
x(+/-h)=x +/- h v + .5 h^2 a_h and their volume determinants, then evaluates
P/P_initial=(J/J_initial)^(-gamma). This is an adiabatic geometric derivative
check, not a pressure reset or an evolved solution. Stable determinant
increments and log1p/expm1 avoid cancellation. Its central second-difference
error relative to the analytic numerical curvature is:

| actual Kepler, h | maximum difference / time^2 |
|---|---:|
| 2e-5 | .03044805 |
| 1e-5 | .007612013 |
| 5e-6 | .001903003 |

The error decreases by4 at each halving, while the source curvature error
remains1216.89. The independent geometry calculation therefore confirms the
derivative diagnostic. Gravity-off/linear controls also show the expected
factor4 convergence to their nonzero physical curvature, rather than being
incorrectly required to have stationary Kepler pressure.

Machine evidence and MPI/quadrature differences are in
`MATERIAL_HYBRID_BUNDLE_20261021.json`.

## Promotion status

Discrete face momentum/work/angular representation is established only in the
measured initial smooth scope. Source/thermal-curvature accuracy is unresolved.
Do not start long evolution merely because the numerical traction closes a
ledger. The bundled evidence goes to the existing Astra planner for the next
coherent decision. Curved remap/cold re-expression, density-wave/sliding,
shock stress, actual native scalable CUDA runtime, t6-10 and longest-paper
stability remain required. Original outputs/restarts/runtime sources remain
preserved. No success boundary has been shortened.
