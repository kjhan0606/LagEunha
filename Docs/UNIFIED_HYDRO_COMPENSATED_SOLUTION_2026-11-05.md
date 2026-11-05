# Compensated solution and complete original-test audit

The preceding goal turn committed and pushed authoritative source/evidence.
The present stage diagnoses the late blast failure, implements the grouped
numerical plan, reaches the original blast endpoint, completes N400 strong
refinement, and exposes two failures in the original nine-test suite.
The full common-method, physical-accuracy, 2D, and <=2x GIZMO objective
remains active. No result is promoted into the accepted Kepler manuscript.

The archive is `Research/CompensatedSolution_20261105`. Original runs remain
under `/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1`. Native C/OpenMP evolves
the fluid. Python orchestrates and analyzes. No new production GPU/MPI
trajectory is established by the present CPU prototype.

## The late blast failure is a solution-representation bottleneck

The exact rejected version8 state reproduces the original cap64 failure
1.576629526947e-12. No earlier iteration passes the native 1e-12 criterion.
An independent 60-digit LDL solution rounded to one double per node also
fails at 1.433263661e-12. Evaluating the same rounded vector in extended
arithmetic does not remove the residual.

Rows52/53 have z near 1.273544e-4 with difference 2.217566e-11 and face
Jacobian 2.079473e9. One node ULP contributes approximately 3.0e-12 to the
scaled pair residual. The native matrix has spectral condition about
2.07947e7. Both generators lie in the same physical cell and separate
slightly. The older different-cell reconstruction jump is therefore not
the demonstrated cause at the late failure. The evidence does not prove
the nonexistence of every passing one-double vector.

## A restricted numerical change retains physical version8

The grouped Astra evaluation selects two-word accumulation of the solution
and two-word evaluation of the difference before rounding. Matrix entries,
IC(0), Krylov vectors, and matvec remain double. TwoSum and FMA TwoProd
retain accumulation errors. The standard construction is described by
[Hida, Li, and Bailey](https://www.davidhbailey.com/dhbpapers/qd.pdf).
The build disables floating-point reassociation and contraction.

The native true residual, final conservative face flux, mu*z, and successful
solution copy consume the same two-word solution. Physical equations,
coefficients, birth mu, cap64, and actual residual tolerance1e-12 remain
unchanged. No backward-error allowance is added to the stopping rule.

The 112-byte cell checkpoint remains unchanged. Transient rate/statistics
structures gain rounding diagnostics. The transient rate rounds z to double
and records the rounding error. The harness observes Euler/RK weight-update
rounding against a long-double arithmetic reference. Long double does not
evolve the fluid. Existing stored-W geometry, GCL, and mass-weight checks
still evaluate actual updated states.

## Frozen gates and costs

All 52 physical fixtures pass under the compensated backend. The exact late
blast state converges in two iterations with true residual6.002185e-16.
An independent Decimal reference represented as hi+lo gives6.036879e-16.
Removing low words fails the original criterion. Corrupted vectors and
NaN/Inf are rejected. Conservative pair assembly, KEP, entropy, unchanged
state, and unchanged A/b/scale/geometry/physical-velocity checks pass.
The previous four native 1D/2D/3D systems also pass.

Three interleaved rounds of 100 frozen RHS evaluations give the following
compensated/double median cost ratios. Setup, factorization, compensation,
true checks, final flux, and rate rounding diagnostics are included.
Geometry, I/O, and integrator-level weight-update observations are excluded.

| Frozen system | Solve ratio | Whole RHS ratio |
|---|---:|---:|
| late blast | 0.0530 | 0.2377 |
| earlier blast | 1.2703 | 1.0553 |
| collision | 1.2724 | 1.0352 |
| historical Gresho64^2 | 1.1212 | 1.0915 |
| historical Sod16^3 | 1.1758 | 1.0855 |

The late double control fails cap64 and is not an equal-accuracy speed
comparison. None of the frozen timings establishes GIZMO or GPU/MPI speed.

Two initial audit issues remain archived. The first wall executable was
invoked without the required `fixtures` argument and exited2 without a PDE.
The initial CSR comparison normalized a cancellation-prone constant probe
by the small resulting action. The observed discrepancy3.21e-10 is retained.
The corrected algebra audit checks coefficient-scale backward error and the
native exact constant mode separately. The native root stopping criterion
is unchanged. All final gates use the corrected explicit fixture arguments
and corrected algebra audit.

## Full targets and refinement

Blast resumes the valid version8 accepted state at t=0.03501883433 and reaches
t=0.038 after181 additional accepted steps. Entry checkpoint bytes match the
parent and diagnostic trials remain rejected as production restarts.
The whole physical path combines block PCG, double IC(0), and compensated
IC(0). The mixed path is not a single-backend timing benchmark.
Whole-path maximum conservation error is2.49e-14 and residual-corrected
mass-weight error5.04e-13. The latest segment's sum of accepted weight-update
rounding bounds is1.35e-14 in cell-mass-relative units.

Fresh N400 blast and collision reach their original t=0.038 and0.035 targets
in7135 and4304 accepted steps. Floors and repairs remain zero.
N400 full-run rounding bounds are3.88e-13 and4.34e-14. Existing conservation
and mass-weight checks pass. KEP/entropy checks retain their existing
scale-dependent criteria, not an unscaled absolute threshold.

The independent Athena exact ideal-gas Riemann utility supplies collision
states. Rankine-Hugoniot defects are independently checked. Norms integrate
over actual native intervals restricted to original x in[0,1]. Density,
momentum, and energy comparisons use conservative densities. Primitive
pressure errors are reported separately from conservative energy errors.

| Collision metric | N200 | N400 |
|---|---:|---:|
| density field L1 | 0.1256850 | 0.0627244 |
| relative density field L1 | 1.10156% | 0.54975% |
| relative pressure field L1 | 1.19518% | 0.59143% |
| maximum cell pressure | 2854.986 | 2640.494 |

The exact star pressure is1691.647. Finite-mass shock peaks and neighboring
oscillations remain. Two-level integral-error reduction does not establish
pointwise accuracy or asymptotic convergence. The images and interval tables
are archived with the exact-solution analysis.

## Original nine-test scope is incomplete

The completed physical8 N200 blast and collision paths are reused.
Five additional full original targets pass, namely Sod t=0.2, Lax t=0.13,
Shu--Osher t=1.8 with the original sin(5*pi*x), and stationary/translating
open contacts t=0.25. Contact ROI density/pressure/velocity L1 errors are
at most2.50e-14. The successful endpoint count is therefore7/9, not9/9.
Noh and double rarefaction fail, so the large 2D suite is not released.

Noh stops at t=0.0053020828 of0.6. The last accepted wall cell has
rho252434.9, P1174.877, and volume2.20e-8. The exact finite-pressure mirrored
Riemann star has rho3.999989 and P1.333336. Exhaustive power halfspaces find
positive cells and no hidden cell. At the wall, mean normal jump=-1.9860 and
chi=1, but both constrained reconstructed traces are zero. The additional
reconstructed-compression condition therefore sets artificial pressure to
zero. The next trial also fails the actual residual criterion. Treating
the already unphysical wall pileup as only a numerical residual problem
would miss the physical failure.

Double rarefaction stops at t=0.04991567093 of0.15. The closest generators
have gap2.09753e-7 and gap rate-0.1118284. Their cell means have the same
nonmonotone velocity difference. Centroids remain separated by0.04498163.
The final rejected trial crosses generator order. Exhaustive post-crossing
power cells remain positive at that saved instant. Sorting the arrays alone
would not prove a continuous valid trajectory through coincidence.
The raw gradients are positive while the scalar-limited gradients vanish
at the spurious local extrema. Existing artificial pressure is0.0112834.
An exact rarefaction reference is evaluated at the actual partial time,
not relabelled as the requested endpoint.

## Independent blast reference attempt and next work

Independent Athena PPM/RK3/HLLC reference attempts at N1600/3200/6400 retain
the prior strict abort-on-floor policy. All three abort with
`REFERENCE_FLOOR_WOULD_ACTIVATE` and provide no valid final blast reference.
The new three-region IC integrates cell pressure exactly. The experiment
reuses matching compiled Athena objects, compiles the new IC, and relinks.
No clean-build or successful-reference claim is made.

Next work must evaluate the redundant compression gate and rarefaction
oscillations as common physical issues before promoting the candidate.
Full original tests, padding/reference/refinement qualification, the same
method in complete 2D and3D evolution, production GPU/MPI integration, and
matched GIZMO runtime<=2x remain required.
