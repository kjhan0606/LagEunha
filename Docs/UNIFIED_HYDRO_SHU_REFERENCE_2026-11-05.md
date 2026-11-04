# Shu--Osher IC mismatch repaired; independent native reference comparison

The new common 1D initialization silently used `rho=1+0.2 sin(5x)`, while
the original `init_shuosher()` and historical comparison problem use
`rho=1+0.2 sin(5*pi*x)`. This is a separate error from the previously found
upstream-boundary loss in the historical reference. Both invalidate using
that reference to judge the former common Shu result as the same problem.
The 33-control stationary-contact diagnosis remains unchanged.

## Native initialization repair and unchanged-problem regression

`Hydro1DExam/gfs_unified1d.inc` now integrates the original sinusoid with
wave number `5*pi`. This changes only the common Shu initial density
integral, not pressure, viscosity, coefficients, force, or time integration.
Earlier accepted sin(5x) files are preserved as a distinct test variant.

An independent C wrapper includes the actual native initializer and uses
eight-point Gauss integration of `1+0.2 sin(5*pi*x)` on each initial interval.
At N200, before any PDE update:

| Initializer | Density integral L1 error | Maximum density error |
|---|---:|---:|
| Parent sin(5x) | 1.45458047634 | 0.37924986909 |
| Corrected sin(5*pi*x) | 5.29e-14 | 6.88e-14 |

Pressure and velocity retain their intended initial values. Parent and
corrected native executables were built with the same GCC flags and the
same frozen dependencies. Contact, moving contact, Sod, blast, Noh, Lax,
double rarefaction, and strong collision at N200 reach their original
endpoints and have byte-identical parent/current final data. No floor or
step retry activates in these controls.

## Independent Euler reference, same physical problem

The reference is native C++ Athena++, upstream commit
`37869da22cd357c2fbdd056fbb3781cdd4f430da` from
[PrincetonUniversity/athena](https://github.com/PrincetonUniversity/athena).
A pristine committed source archive was copied to a NEW GPFS directory;
the existing Athena++ checkout, including its uncommitted changes, is not
modified or used as the build source.

The reference uses HLLC, PPM, RK3, CFL 0.3, gamma 1.4, domain [-5,5], initial
shock at -4, the original stated left state, and the exact cell-average
sin(5*pi*x) density. Outflow ghosts preserve the constant supersonic left
driver. Every one of the eight density/pressure floor assignment sites in
the adiabatic EOS has an abort precondition: a run cannot conceal a floor
activation. Initial and final primitive tables and history are retained.
Only Python orchestration and postprocessing are used; C++ evolves Euler.

N1600/3200/6400/12800 all reach exactly t=1.8 with positive final density
and pressure, zero floor/reset count, and unchanged upstream pressure
10.33333. The fourth reference level was registered after the first three
revealed appreciable reference uncertainty, before its execution. No native
trajectory was restarted or coefficient changed to obtain this comparison.

The first three-reference pass inadvertently compared to the OLD sin(5x)
native variant. Those raw comparisons are retained with an explicit
`EXCLUDED FROM ACCURACY CLAIMS` explanation. They are evidence of mismatched
inputs, not evidence of native Euler nonconvergence. The valid Athena
trajectories were kept; only newly corrected native IC trajectories enter
the final accuracy table.

## Actual-cell norms and measured reference uncertainty

Integrate absolute differences between piecewise-constant states over the
common interval [-4.9,4.9]. Native W0 faces are neighboring generator
midpoints, not `generator +/- V/2`; the reconstructed widths are checked
against the saved actual V. The complete comparison interval is covered.
These are unnormalised spatial integrals, not mean relative errors.

| Athena reference pair | Density L1 | Velocity L1 | Pressure L1 |
|---|---:|---:|---:|
| 1600 / 3200 | 0.217405 | 0.0388283 | 0.275631 |
| 3200 / 6400 | 0.109839 | 0.0227663 | 0.152695 |
| 6400 / 12800 | 0.0491761 | 0.0107518 | 0.0734945 |

Corrected native W0 fixed-mass variational pressure, existing implicit
viscosity and global KDK give the following comparison to Athena N12800:

| Initial native N | Final cells including inflow | Density L1 | Velocity L1 | Pressure L1 |
|---|---:|---:|---:|---:|
| 200 | 294 | 0.960836 | 0.371140 | 2.016174 |
| 400 | 589 | 0.717536 | 0.267407 | 1.606584 |
| 800 | 1178 | 0.499689 | 0.188728 | 1.184454 |
| 1600 | 2357 | 0.313848 | 0.122776 | 0.793806 |
| 3200 | 4714 | 0.183202 | 0.0748531 | 0.495837 |

All five native runs reach t=1.8, retain positive final states and ordered
generators, and report zero floors/retries. Errors decrease in all three
reported fields. This is now an independent same-IC Euler comparison rather
than only self-convergence. The reference still has finite resolution
error; its inter-resolution differences are not an analytic error bound.
Do not claim the finest native error or formal order is known exactly.

The console's legacy point-interpolated norms use an Athena N6400 table
solely for that diagnostic output. The table above is recomputed from
actual-cell overlaps against N12800. Neither is interchangeable with the
old damaged reference norm or the old sin(5x) self-convergence sequence.

## Scope and reproducibility

Run root:
`/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/shu_independent_athena_v1`.
Compact records, frozen native sources, the upstream problem/EOS patches,
build commands and exact output snapshots are archived under
`Research/ShuReference_20261105`. The Python replay source is
`Exam/HydroCore/shu_independent_reference.py`. Native binaries remain in
the GPFS run directory, with SHA256 provenance.

The new four/five-state FV--ALE candidate is NOT used in these trajectories.
This repair establishes an honest same-problem reference and corrects a
test-definition regression. It does not solve Gresho/Kepler shear,
qualify the remaining 2D/3D suite, or demonstrate GIZMO <=2x runtime.
The full goal remains active. The original manuscript baseline is retained.
