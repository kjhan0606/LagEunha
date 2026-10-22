# Reconstructed HLL volume fluctuation fails its first angular gate

2026-10-21. Implements the single bounded candidate from
`ASTRA_VOLUME_UPWIND_PLAN.md`. Driver runs native C++ on authorized lageunha.
No P4, new Jacobian, runtime integration or evolved snapshot follows this
failure. The full t6-10/longest-GIZMO target remains active.

## Implemented formula and actual state

Native helper `Exam/gfs_sbp_upwind.h` evaluates complete Euler HLL on
reconstructed midpoint traces. It includes endpoint and Roe acoustic speeds,
explicitly rejects invalid traces/intermediate states, and asserts no
positivity theorem for the high-order scheme. There is no tunable factor.

Stored state is the same conserved nodal Euler IC. The same P2 SBP derivative
matrices give each node's primitive gradients; each side extrapolates linearly
to their common physical midpoint. No nodal value or cell energy is replaced
by that trace. All four components use the single correction

    f_new = f_KEP(nodal_L,nodal_R)
              + f_HLL(trace_L,trace_R) - f_KEP(trace_L,trace_R).

The exact same24 actual failure/neighbor polygons and three states are used:
equilibrium, b=(.013,-.007), and b plus prescribed affine mesh motion. All
cases complete with positive traces. Source/executable/outputs are frozen
in `nodal_sbp_upwind_A_v1/`; status is `AUDIT_COMPLETE`, meaning calculation
finished, **not** physics accepted. `analysis.json` contains the measurements.

## A: what holds

15120 affine primitive node-pair checks:

| Contract | Measured maximum |
|---|---:|
| affine trace error (P component divided by1e-6) | 3.50e-11 |
| matching-trace correction / metric magnitude | 1.34e-15 |
| opposite canonical pair correction / metric magnitude | 0 |
| boosted correction / metric magnitude | 1.32e-16 |
| actual patch correction's Euler sum / sum absolute contributions | 2.10e-16 |

Thus pairwise mass, both momentum components and total energy are conserved
to roundoff. This is insufficient for angular weak balance or cold accuracy.

## A: actual angular balance fails

For each cell use the actual quadrature quantity

    Lh = sum_a M_a [(x_a-cx) py_a - (y_a-cy) px_a].

Its derivative includes the node-motion term M_a (wx py - wy px). Compare
with the full numerical angular flux on the physical boundary. Central
gravity has zero torque about this origin, using identical nodal positions.
The volume correction's angular contribution is recorded separately.

Across the24 cells/states, max absolute unexplained angular-rate residual is
**1.11781360e-8**; max residual / sum(M rho v^2) is **.09788694**.
The correction itself reaches .09788723 on the same normalization. These
are rates normalized by an energy/torque scale, not a measured fractional
angular-momentum loss after a time integration. No evolution took place.

Example ID68584, b=(.013,-.007), affine mesh:

- raw angular plus boundary residual -1.08639384e-8;
- correction's angular contribution -1.08639701e-8;
- normalized residual -.09788694.

This violates the explicit first gate. No independent spin or torque repair
is added to absorb the discrepancy. It does not prove that every possible
Euler volume-upwind construction fails; it rejects this concrete SBP-gradient
midpoint reconstruction and correction.

## Simultaneously observed thermal damage

No additional B run was launched. The already evaluated A patch RHS also
provides its pressure derivative, using the same conserved kinetic subtraction:

- max abs(Pdot)/P0 = **33533.4986/time**;
- max volume correction pressure-rate magnitude/P0 = **33536.2902/time**;
- max residual against -gamma P div_h(v) = **33536.2902/time**.

At ID68584 in the affine-mesh case, mean Pdot/P0=2047.04/time and volume
L1 magnitude=2492.25/time. The defect is not an isolated vertex or a small
opposite-sign fluctuation concealed by a mean. The physical initial Pdot is
zero. No heat repair, limiter or physical viscosity was applied.

## Additional baseline angular evidence

`nodal_sbp_angular_v1/` and `nodal_sbp_angular_v2/` evaluate only the previously
unmeasured angular terms; existing density-wave matrices are not rerun.
They use the same selected split flux and full-Euler central diagnostic
control respectively.

| Volume flux | max abs angular residual | max residual / sum(M rho v^2) |
|---|---:|---:|
| original KEP split | 1.16981e-13 | 1.03366e-6 |
| arithmetic full Euler control | 3.16351e-20 | 3.90429e-14 |
| reconstructed HLL correction | 1.11781e-8 | .09788694 |

The original split form is therefore not assumed to inherit the previous
enriched solver's canonical angular ledger. Its exact kinetic/thermal chain
identity does not guarantee the linear-coordinate angular weak identity.
These are cell-local balances; symmetry cancellation in a global Kepler IC
would not establish correct internal angular transport.

## Next planning input

The next bundled planning decision must address actual angular stress/weak
balance together with cold kinetic/thermal compatibility and local linear
stability. A generic pair momentum sum is insufficient. Retain the native
screened mass/weight solver and successful prior canonical face/angular
proof, but connect them only after the volume representation is viable.
Do not reuse this failed correction with a smaller adjustable multiplier.
No current nodal candidate is approved for time evolution.
