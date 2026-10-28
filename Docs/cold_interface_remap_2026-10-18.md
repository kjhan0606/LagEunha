# Cold density-interface remap: pressure bookkeeping is essential

Repeat the perturbed-grid single-remap diagnostic with rho=.01 at donor
generators r<.5 and rho=1 outside. Uniform thermal energy density e_vol=1e-5
gives uniform pressure for any fixed ideal-gas gamma. The interface follows
donor polygons, not exact circle cut cells. Use the canonical softened circular
velocity samples and unchanged small h^2-scaled weight pattern.

Compare two interpretations of the SAME conservatively remapped state:

    U_true = E_hyd - tr(Pi)/2,
    U_mean_only = E_hyd - |p|^2/(2m) = U_true + K_internal.

Hence the hypothetical pressure excess is exactly K_internal/U_true. This
does not assert which subcell closure is appropriate for a physical gas, nor
that the original point-state solver uses the second interpretation incorrectly.

| Grid | Mixed-density receivers | Max pressure excess if internal K thermalized | Volume RMS excess in mixed cells |
|---|---:|---:|---:|
| 16^2 | 29 | 36.8506 | 17.1451 |
| 32^2 | 67 | 12.5594 | 4.57096 |
| 64^2 | 142 | 4.07256 | 1.40653 |

Entries are relative excesses, not pressure ratios: 4.07256 means a407.3%
increase and resulting pressure5.07256 times the original. With the retained
kinetic moment, constant pressure is reproduced within3.25e-11 relative in
all tested interior cells. The remaining small discrepancy involves subtraction
of thermal energy from much larger kinetic energies and numerical quadrature;
no pressure floor is used. Global moment drift<=1.99e-13 absolute.

These errors decrease with spatial refinement, but can remain large compared
with the deliberately tiny thermal energy at finite resolution. Conserving
total energy alone therefore cannot certify cold-flow pressure accuracy for
a new remapping scheme. Keep internal kinetic covariance separate until a
physical, validated transfer into heat is specified.

## Scope and next priority

This is a synthetic one-remap comparison, not an observed heating event in
the original CUDA solver, not shock/AV evidence, and not a new production fix.
No physical time, generator displacement, pressure force or gravity work is
applied. The original production code has not adopted these experimental
remapping components. Do not conflate this issue with the already observed
geometric volume/shear mismatch in the fixed-mass original scheme.

The next design gate is a consistent physical evolution law for retained
kinetic covariance and interface transport that also satisfies the approved
persistently-equal-weight reduction. Do not launch a new physical Kepler
trajectory just by replacing its energy interpretation. Existing long runs
and their original energy conventions remain the reference.

Evidence: `canonical_cold_interface_remap_2026-10-18.json`.
Script: `audit_cold_interface_remap.py`.
Remote: `/gpfs/kjhan/LagForce/Codex_cold_interface_remap_20261018`.
