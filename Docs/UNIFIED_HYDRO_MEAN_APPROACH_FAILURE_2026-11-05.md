# Mean-approach stress resolves one boundary defect but fails common gates

Physical9 is an experimental opt-in candidate and is not accepted as the
common method. The finite grouped bundle is closed after double-rarefaction
and cold-shear failures. No coefficient, CFL, cap, mu, floor, generator
sorting, or permanent extra state is introduced to pass the failed gates.
Physical8 results remain separate. Full original9, common complete2D/3D,
production GPU/MPI, and <=2x matched GIZMO remain unfinished.

The source flag is `GFS_NATIVE_MEAN_APPROACH_STRESS`. The candidate removes
the reconstructed-jump condition from the artificial pressure. The mean
normal approach, raw chi, positive centroid-normal length, coefficients1/2,
physical reconstruction, geometric face kinematics, and weight law are
otherwise retained. The checkpoint physical version is9. Production loading
of physical8 under the new binary is explicitly rejected.

## Grouped operator gates

All 52 existing wall, frame/contact, and stress fixtures pass. Three extra
1D/2D/3D local affine-compression checks measure linear and quadratic stress
ratios1/2 and1/4 when the physical patch length halves. The local checks are
not PDE convergence. Affine noncompressive checks retain the raw-chi limits.

The saved physical8 Noh state is read only. Mean approach and chi are
unchanged. Removing the reconstructed gate gives positive stress/heating
while reflected wall mass and total-energy flux remain exactly zero.
The actual viscous timestep coefficient and KEP/entropy identities pass.
The saved doublefan closest pair has exactly the same artificial pressure
under both operators. No accepted step advances from either saved state.

## Three independent birth gates

The fixed physical9 binary evolves Noh200, doublefan200, and cold64^2 from
birth. Noh and doublefan use separate single CPU cores. Cold uses eight
OpenMP cores. The cases are independent members of one finite bundle.

| Case | Target | Last accepted t | Status |
|---|---:|---:|---|
| Noh200 | 0.6 | 0.6 | endpoint reached, accuracy still limited |
| doublefan200 | 0.15 | 0.05002352340 | failed before endpoint |
| cold64^2 | 0.25 | 0.25 | pressure gate failed |

Noh completes3182 accepted steps without the previous divergent density
pileup. Finite-pressure exact mirrored Riemann comparison gives density,
pressure, and velocity L1 errors0.0147350,0.00250393,0.00455067 on[0,1].
Maximum density4.06385 and pressure1.36879 are near the exact star4.0 and
1.33334. However, the wall cell density is2.05612. Endpoint completion is
therefore not pointwise Noh accuracy, and no N400 qualification is claimed.

Doublefan again crosses generator order and exhausts the existing geometry
retries. The partial-time density, pressure, and velocity L1 errors are
0.0175871,0.00596996,0.128008. Absolute16/32-point quadrature differences
are8.32e-6,4.39e-6,3.34e-7. The errors are evaluated at the actual partial
time, never at a substituted requested endpoint. Total downward mean-velocity
variation is0.34328 although the exact rarefaction velocity is monotone.
Raw sign-based affected-volume measurements include roundoff signs and
must not be interpreted as the physical extent of oscillations.

Cold reaches0.25 in two accepted SSPRK2 steps but has relative pressure
L1=0.3360922, exceeding the fixed0.01 limit. Integrated stress heat is
5.0415794e-7. The absolute accumulated RK kinetic defect is1.6624026e-11.
With the native gamma1.66666662693, stress heat alone predicts relative
pressure increase0.3361053, close to the measured0.3360922. The dominant
thermal production in the measured trajectory is the stress closure.
Temporal convergence is not claimed, and the failed candidate is not rescued
by reducing the registered timestep.

## Frozen controls isolate smooth-shear dissipation

The same accepted physical9 cold state is evaluated under the physical8
and physical9 operators without advancing a step. Filtered heat rates are
1.38543e-8 and2.35150e-6, a ratio169.73. About99.988% of the physical9 heat
comes from the LSQ residual component rather than the fitted compression
component. The old-operator evaluation is not an old-method trajectory.

An additional manufactured-field control removes accumulated velocity noise.
The saved native geometry is retained, and analytic incompressible
u=(sin(2*pi*y),0), rho=1, P=1e-6 cell means replace Q only in memory.
The analytic cell-mean velocity follows the divergence theorem with exact
cosine integrals over native straight faces. No manufactured checkpoint is
written or accepted, and the manufactured input is not a physical restart.

The resulting heat rates remain1.37351e-8 under physical8 and2.32711e-6
under physical9, a ratio169.43. About99.987% of the new heat is the LSQ
residual contribution. Resolved non-affine incompressible shear therefore
produces artificial heat under the new closure even without evolved velocity
noise. Interpreting every affine-fit residual as unresolved noise is not
sufficient for the cold-shear gate.

The first manufactured-field auditor build had a nested-main naming error
and produced no simulation. The corrected explicit audit entry point,
build commands, frozen sources, and both attempts remain archived.

## Current decision

The blanket reconstructed-gate deletion cannot serve as the final common
method. The branch retains physical9 only as a named reproducible experiment.
The next grouped evaluation must reconcile boundary shock dissipation,
resolved non-affine shear, and rarefaction velocity oscillations. Repeated
constant tuning or importing a successful physical8 endpoint into physical9
would not satisfy the objective. Larger suite release remains withheld.

All source, receipts, failed trials, and analyses are in
`Research/MeanApproach_20261105`. All work uses GPFS and native C/OpenMP.
The accepted Kepler manuscript remains unchanged.
