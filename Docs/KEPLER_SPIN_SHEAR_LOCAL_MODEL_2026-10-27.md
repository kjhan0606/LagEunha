# Local Kepler shear is not captured by rigid spin alone

This is an analytic local model, not a measurement of the native solver.
It clarifies which velocity components a rigid-cell spin variable could
represent in the user's inviscid Kepler experiment.

At a point `(r, 0)` in an unsoftened circular Kepler flow, write
`Omega = sqrt(GM/r^3)` and use local Cartesian displacements `(dx,dy)`.
The first-order velocity gradient is

```
A = grad(v) = Omega [  0   -1  ]
                    [ -1/2  0  ].
```

Its antisymmetric part is rigid rotation with angular velocity `Omega/4`,
while its symmetric part is a shear with off-diagonal entries `-3 Omega/4`.
Thus even an inviscid Kepler flow contains both local rotation and strain;
neither implies physical dissipation when `Re = infinity`.

For a hypothetical cell with uniform density and isotropic central second
moment `Qxx = Qyy = I/2`, `Qxy = 0`, the local linear velocity is
`v-v_center=A(x-C)`. The cell spin and kinetic energies follow directly:

```
S = integral rho (x-C) cross (v-v_center) dA = I Omega/4
E_rigid = S^2/(2 I) = I Omega^2/32
E_strain = 9 I Omega^2/32
E_linear = 10 I Omega^2/32.
```

In this illustrative isotropic, infinitesimal, unsoftened case, rigid spin
represents one tenth of the linear subcell kinetic energy; strain accounts
for nine tenths. The ratio is not a measurement for the finite, irregular,
softened simulated cells. Their general second-moment tensor and actual
velocity gradient must be used, including any density reconstruction.

It follows that adding only `S^2/(2I)` may represent a physically relevant
component of internal motion, but it cannot by itself represent all local
Kepler shear. Treating the remaining ordered strain energy as heat would
create an artificial thermal contribution. A viable numerical model must
state how it reconstructs and transports that unresolved component, or
demonstrate convergence without explicitly storing it. This observation
does not establish whether the user's native error can be fixed by spin,
nor does it justify adding a physical viscosity to an `Re=infinity` test.

For a general central mass second-moment tensor `Qab` and velocity gradient
`Aab`, the two-dimensional intrinsic spin of the affine field is

```
S = Ayx Qxx - Axy Qyy + (Ayy-Axx) Qxy.
```

The tensor reduces to the scalar inertia `I = Qxx+Qyy` only for the rigid
energy relation. With an anisotropic cell, even the symmetric strain can
contribute to `S`; the actual polygon or density moment therefore matters.

## Measurement on saved native cells

CPU Slurm **437535** completed `0:0`. Its
[diagnostic](audit_kepler_affine_spin_energy_20261027.py) reconstructs the
actual polygon second-moment tensor under a **uniform cell density**
assumption, combines it with the saved native `stress.gUxx/gUxy/gUyx/gUyy`
velocity gradients, and projects the resulting affine field onto the best
rigid spin about the mass centroid. Both snapshots use the same 5%-jitter
initial state at nearly identical `t=1.65137433`:

| Case | Rigid-spin fraction of reconstructed affine kinetic energy |
|---|---:|
| [W0](KEPLER_AFFINE_SPIN_ENERGY_2026-10-27/w0.json) | 11.3299% |
| [Weighted](KEPLER_AFFINE_SPIN_ENERGY_2026-10-27/weighted.json) | 11.3296% |

The saved-gradient RMS error against the analytic softened Kepler gradient
was about `6.1e-4` relative in each case. Each analysis passed the polygon
area, native volume, tensor positivity and rigid-energy bounds. Before the
native analysis, analytic square checks passed rigid rotation, simple shear,
and the `10%` isotropic local-Kepler result.

This is a decomposition of a reconstructed **linear subcell velocity field**
under the uniform-density closure. It is not an evolved spin energy, a
measured heating fraction, or a demonstration that the remainder is lost to
heat. Saved gradients may correspond to a nearby RK stage; nonuniform
subcell density and nonlinear velocity structure are not represented. The
result reinforces the need to account for non-rigid velocity structure when
designing a physically conservative spin model.

To extend this beyond the early checkpoint, CPU Slurm **437539** is queued
with `afterok:437499`. The [runner](kepler_affine_spin_energy_t65_20261027.slurm)
first requires the fresh weighted `t=6.5` 21-check audit, exact binary and
checkpoint hashes, and the same initial state as the accepted jittered W0
`t=6.5` control. It then repeats the analytic square tests and the two
polygon/gradient analyses at the matched historical failure-window time.
This dependency does not delay or change the native `t=10` continuation.
The later superset job **437545** was staged with the same `afterok:437499`
dependency and adds the cellwise thermal-positivity analysis below. The
older **437539** was confirmed `PENDING` with zero elapsed time and then
canceled, so the geometry/gradient analysis will not be duplicated.

## Energy-initialization constraint in the cold disk

CPU job **437543** extended the same two `t≈1.6514` analyses with a
counterfactual positivity calculation. For both [W0](KEPLER_AFFINE_SPIN_ENERGY_2026-10-27/w0_thermal_risk.json)
and [weighted](KEPLER_AFFINE_SPIN_ENERGY_2026-10-27/weighted_thermal_risk.json)
states, the reconstructed rigid-spin energy summed over disk cells is about
61.1% of the current stored native thermal energy; the reconstructed full
affine kinetic energy is about 5.39 times that thermal energy. If one
**hypothetically** held the current discrete total energy fixed while
subtracting just that reconstructed spin energy from the current `U`,
roughly 7,643 of 49,038 disk cells (15.6% of disk mass) would have `U≤0`.

The present native solver does not store this reconstructed subcell kinetic
energy as a separate conserved variable, so this subtraction is not an
initialization prescription and is not evidence that the existing state is
unphysical. It demonstrates why a future spin model must explicitly define
its initial total-energy representation and EOS/pressure treatment rather
than silently transferring energy out of the existing cold gas. The same
uniform-density and affine-field assumptions and t=1.65 limitations apply.

A possible initialization to investigate is to define a **new** conserved
energy representation from the outset, rather than take spin energy out of
the old `U`. For an affine cell field `v(x)=v_C+A(x-C)` about the mass
centroid, its integrated kinetic energy is

```
K_cell = 0.5*m*|v_C|^2 + E_affine,
E_affine = S^2/(2I) + E_nonrigid,    E_nonrigid >= 0.
```

One can keep the intended initial thermal `U` and include these terms in
the newly defined total energy, provided the initial field and mass moments
are physically justified. If the generator and mass centroid differ,
`v_C` must be derived consistently; the generator velocity cannot simply be
assumed to be `v_C`. The non-rigid remainder must be represented or otherwise
closed without silently converting it into heat. This construction changes
the discrete energy reference, so comparisons must use the corresponding
new initial total rather than claiming the old and new totals are identical.
Face fluxes, gravity work, angular momentum, positivity and the user's W0
fixed-mass reduction still need a consistent discrete evolution and native
tests. This is a design option, not an implemented or validated solver.
