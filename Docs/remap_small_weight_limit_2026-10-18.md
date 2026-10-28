# Small weight displacement: a continuous single-remap limit

Reused the original conditional t2 remap, existing compiled CPU/CUDA libraries,
and unchanged initial data in a fresh diagnostic directory. Three new weight
fractions .0625, .015625, .00390625 supplement existing .25 and 1 results.
Only order8 was run for the new cases. No physical trajectory was restarted.

| Fraction | L1 mass change / total mass | Internal K / original thermal U |
|---|---:|---:|
| 1 | .00275031 | .0801869 |
| .25 | .000685254 | .0208548 |
| .0625 | .000171170 | .00533741 |
| .015625 | .0000427836 | .00134683 |
| .00390625 | .0000106953 | .000337638 |

Last-pair displacement scaling exponents are 1.00004 for mass change,
.97086 for mean-velocity change RMS, and .99801 for internal kinetic energy.
All new cases retain positive mass and thermal energy. CPU1/OpenMP4/CUDA
atom outputs are bitwise equal using the existing libraries. Shared-atom
global invariant differences are at most 2.23e-16; comparison to independently
integrated original point invariants is below 9.8e-15. This is a component
check, not full RHS/MPI or timestep convergence evidence.

## Corrected interpretation of the preceding face-trace audit

A large facewise relative difference from a piecewise donor field does NOT
exclude affine reconstruction as a class. Exact preservation of every donor
discontinuity is not a necessary criterion for a convergent numerical method.
The small-displacement data show that the integrated state change and internal
kinetic content tend to zero approximately linearly on this fixed mesh.
The earlier audit supports only the narrower statement that affine traces do
not exactly represent the conditional remapped state. It does not independently
justify rejecting every approximate affine closure.

Nevertheless, discarding or thermalizing an O(delta W) internal kinetic term
at each remap can produce an O(1) accumulated effect when delta W=O(dt) and
the number of remaps is O(1/dt), if no cancellation mechanism is established.
This is an order-counting risk, not measured long-time heating. A single
displacement limit does not prove the approved persistently equal-weight
fixed-mass solver reduction, nor establish consistency of coupled evolution.

Next meaningful test is cumulative remap composition with conserved internal
kinetic bookkeeping, comparing fixed donor-field integration to repeated
projection. Do not add new production moment equations solely to match the
stationary-face diagnostic. Existing long Kepler runs remain the physical
reference and the final stability gate remains open.

Evidence: `canonical_remap_small_weight_scaling_2026-10-18.json`,
`canonical_remap_small_weight_backend_2026-10-18.json`.
Scripts: `run_remap_small_weight_limit.py`, `analyze_remap_small_weight_limit.py`.
Remote: `/gpfs/kjhan/LagForce/Codex_remap_small_weight_limit_20261018`.
