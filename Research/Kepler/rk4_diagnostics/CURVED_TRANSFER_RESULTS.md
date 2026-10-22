# Actual curved transport estimates, failure of radial allocation convergence

2026-10-21. Native C++/OpenMP/MPI on authorized lageunha. No Kepler time
evolution, new IC, weight Newton or cold compression occurred. The original
finite h=2e-4 actual P3/h2 map/velocity and its mapped generator sites are
retained. Density is transported as rho0/J and physical thermal energy from
the preserved adiabatic finite CONTROL, never assigned from E-K.

## New implemented evidence

`material_curved_trace_bounds.h` bounds the entire shared-edge functional
difference, not only point samples. Each trace is an exact cubic polynomial
plus Omega(r) times a linear polynomial using stored current coefficients.
Interval Bernstein hulls enclose the differences on16 edge subintervals.
All536 shared edges give maximum analytic position gap6.901e-16 and velocity
gap9.342e-15. These are finite mismatch bounds, not exact C0 equality or a
formal native-libm proof. A continuous representative/gluing uncertainty
still has to enter the final geometry certificate.

`material_curved_partition.h` constructs all angular source/target fragments
and retains uncertain bands. `test_curved_transfer_v1.cpp` now integrates
actual m,p,L,K,U,V, position/polar moments, potential energy and thermal rate.
All quadrature weights are positive. A prescribed nonzero weight pattern is
also integrated; it is not a solved dynamical weight state.

`curved_transfer_v3` completes at two quadratures, native MPI1/4 and OpenMP2.
The coarse MPI transfer/uncertainty CSVs bitmatch. Fine W0 saves926330
positive native quadrature records,56 bytes each, with reference source
triangle, target, uncertainty flags, X,Y,positive dX weight and rho0. All
stored records were independently read and checked. Target masses recovered
from stored double weights agree to1.91e-16 relative to original cell mass.
The physical source map remains authoritative; these packets are not a
compressed new thermodynamic state or restart.

Fine source integrals agree with the independent prior native finite control:

- maximum source mass relative difference6.72e-13;
- maximum source K difference divided by cold U2.01e-11;
- maximum source U relative difference6.44e-13.

But target allocation is NOT converged. Raising radial/angular order from
12x10 to24x16 changes target mass by5.038e-5 of original mass, target K by
0.06425 U, and target U by5.028e-5 of original U. This FAILS the cold gate.
Global/source conservation does not repair that per-target allocation error.
Do not increase a fixed quadrature order blindly or promote a weight solve.

Numerical W0 target masses differ from fixed material masses: maximum coarse
relative difference~9.86e-5 and fine~7.24e-5. These establish the need to
examine the finite overlap contract, but are NOT a converged measurement of
the exact drift because target allocation still fails convergence. No mass
was overwritten to impose fixed target masses. A generic W0 finite remap
with simultaneous individual fixed masses/no steering remains unresolved.

`CURVED_TRANSFER_BUNDLE_20261021.json` contains raw results and independent
analysis. Python only reads/analyzes native artifacts; it evolves no fluid.

## Root-isolation recovery and bounded failures

v1 and v2 full coarse executions terminated exit2 after writing complete
diagnostic rows:442 and430 radial nodes hit the angular subdivision budget.
All uncertain intervals were retained; no sliver was dropped. Both outputs
are preserved as failures. Increasing subdivision was not the remedy.

Interval Newton now carries derivative-polynomial coefficients separately.
Differencing uncertain local Bernstein coefficients and dividing by a tiny
interval had created a spurious1/width uncertainty. A monotone root set is
contracted as an interval, and irreducible coefficient uncertainty remains
one explicit band. v3 first reproduces the two troublesome triangles before
the full batch. All full v3 modes have zero angular budget failures.

The initial-velocity remainder in interval evaluation now uses the same
stable binomial-tail identity as the native input helper, with an explicit
infinite value/derivative tail bound. Affine velocity extraction also tightens
the enclosure without changing physical coefficients or the represented map.

`curved_bounds_v4` checks9600 native point/derivative queries, all inside the
new explicit operation-error envelope; all384 triangles remain positive in
the interval check. Maximum native position envelope1.93e-14, gradient
envelope3.65e-12. The position error model separates the single final origin
addition from basis operations, with conservative operation counts and an
explicit <=8ulp native pow assumption. This is an engineering rounding
model, not a formal libm proof. Old raw enclosures are preserved unchanged.

`curved_roots_v2` reruns the five root fixtures plus actual midpoint-radius
partitions with the updated algorithm: no missing known roots, no native
owner/enclosure mismatch, MPI1/4 CSV bitmatch. Fixture root panels1821->297;
actual root panels167494->16594. Root bands remain uncertainty sets, not a
claim of simple-root uniqueness or a radial completeness certificate.

## Actual MPI/OpenMP speed fix

`curved_affinity_v1` directly probes every thread. The inherited
OMP_PLACES=threads implicitly enabled binding even though I_MPI_PIN=0.
Every MPI rank put thread0 onCPU0 and thread1 onCPU64 (each allowed1 CPU).
Four ranks therefore competed on the same hardware threads.

With OMP_PROC_BIND=false and inherited placement variables cleared, threads
can use all128 CPUs. The SAME frozen full coarse4-rank executable/input
cost13.9452s before and3.07795s after, a4.5307x improvement. Transport and
uncertainty CSVs are bitwise identical. This is an observed CPU geometry-stage
speedup, not a whole Kepler/CUDA or four-GPU scaling claim.

The project `run_env.sh` now defaults OMP_PROC_BIND=false while preserving
an explicit user choice. New transfer/topology runners clear conflicting
inherited placements. This does not alter user-wide environment files.

## Current dependency

`material_curved_topology.h` and `test_curved_topology_v1.cpp` implement a
bounded interval-radius root-branch/order audit. Tangencies, endpoint events
and overlapping root-branch ranges are kept with conservative event mass
and absolute K bounds. Priority refinement targets the largest uncertainty;
cap256 radial evaluations per triangle and depth32, no unbounded sweep.
The actual384-triangle W0 execution was terminated after the same method
failed the independent closed tiny-cell fixture; it has no acceptance result.
The frozen diagnostic and signal15 output are preserved in
`curved_topology_v1/TERMINATION_20261021.md`. No Kepler evolution/restart was
stopped or reset. A later exact polynomial reexpression of the unchanged map
passed its algebra/enclosure checks (`curved_compiled_v1`), but the same
tiny-cell fixture still failed: absolute area error1.5154e-11 exceeds1e-11,
and the topology mass/K bounds did not pass (`curved_lens_v3`). The next
geometry method must isolate actual intersection/tangency events and track
resolved branches between them, instead of counting their entire swept bands
as unknown or raising the subdivision cap. Retest the same fixture first.

An adaptive positive integration driver is written in
`test_curved_adaptive_transfer_v1.cpp`. Its launch script refuses execution
unless the radial event guard has completed and passed. It uses8/16 radial
and12/16 angular positive rules on guarded panels, an explicit error estimate
and bounded refinement. Gauss discrepancy is an estimate, not a rigorous
error theorem. It has NOT been executed or promoted merely because source
or global conservation is accurate.

After actual geometry/cold transfer acceptance: coarea weight Jacobian and
one finite solve, retained physical fragments, one admissible cold velocity
compression and reverse remap; then sparse/matrix-free native composition,
sliding/density-wave/shocks, full CUDA/MPI t6/8/10 and longest-GIZMO evolution.
