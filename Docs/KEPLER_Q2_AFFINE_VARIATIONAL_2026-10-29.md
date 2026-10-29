# Frozen Q2 affine kinetic-energy shape derivative

The earlier [rotor screening](KEPLER_Q2_ROTOR_VARIATIONAL_2026-10-29.md)
represented only the rigid-spin part of an affine velocity field. This
follow-up uses the **full reconstructed affine field**, including its
non-rigid strain, on the same Q2 t30 and t47.5 snapshots. It is a mathematical
counterfactual, not an evolved or proposed-to-be-activated solver.

For each uniform-density polygon, let `S` be the central mass second-moment
matrix, `G` the saved reconstructed velocity gradient, and `A=G S` the
transpose of the centered position–momentum cross moment. Treat `A` as fixed
while differentiating the following kinetic energy with respect to generator
positions:

```
K_affine = 0.5 tr(A S^-1 A^T)
delta K  = -0.5 tr(G delta S G^T),    G=A S^-1.
```

For fixed cell mass and a uniform-density polygon,
`delta S=(m/V) sum_faces integral[((x-C)(x-C)^T-S/m) delta n] ds`.
The implementation evaluates the resulting cubic integrand with two Gauss
points per face. `A` is held fixed only for this **partial derivative**;
the simulation never evolved it. The gradient was compared to independent
finite differences that rebuild the polygons while holding `A` fixed. At
t30, the smaller-step directional derivative differs by 2.34e-10 from
5.44e-6; at t47.5 it differs by 6.33e-10 from 3.30e-5. A second step shows
roughly quadratic error reduction. This is only one direction per state.

| Snapshot | Fixed-affine-moment force projection on main inner-ring L-error mode | Existing reconstructed pressure projection |
| --- | ---: | ---: |
| t=30 | **-6.545 / time** | +0.192 / time |
| t=47.5 | **-6.906 / time** | +0.205 / time |

These are instantaneous projections onto the *existing* error, not mode
growth rates or evidence of damping. This partial derivative points opposite
to the existing pressure feedback, unlike the separate fixed-spin rotor
derivative, and is much larger in this normalization. It motivates studying
the omitted strain coupling, but also raises a risk of excessive forces or
artificial stiffness. No arbitrary scaling or dissipative relaxation is
introduced. A conservative force need not damp a mode and can move it into
oscillations or another instability.

Holding `A` fixed is a derivative convention, not a demonstrated physical
equation of motion. It does not define how `A` evolves or how intrinsic
angular momentum exchanges with orbital momentum; these forces alone are
**not shown to conserve total angular momentum or energy**. Uniform polygon
density and saved velocity gradients are assumptions, and the fitted affine
kinetic energy represents a different energy reference from the current native
cold gas. It cannot be subtracted from the existing thermal energy without a
separate positivity and initialization contract. The `r0<2.35` energy set
has an artificial diagnostic boundary, although all face neighbors needed
by the inner disk are included. There is no native GPU/RK4/MPI integration,
no material-interface closure, no cost measurement and no long-time result.

Reproducible data and scripts: [t30](KEPLER_Q2_AFFINE_VARIATIONAL_T30_2026-10-29.json),
[t47.5](KEPLER_Q2_AFFINE_VARIATIONAL_2026-10-29.json),
[t30 script](audit_kepler_q2_affine_variational_t30_20261029.py),
[t47.5 script](audit_kepler_q2_affine_variational_20261029.py).

## Exact-orbit geometry counterfactual

The negative projection above includes the force of this shape-energy
functional on the unperturbed mesh itself. To distinguish that baseline,
the same particle IDs were placed at their analytic softened circular-orbit
positions at the same time; each site's **actual-state reconstructed `A` was
held fixed** while polygon moments and shape forces were recomputed. Project
both the actual-state torque and this reference torque onto the same actual
angular-momentum error. Their difference is not a growth rate.

| Time | Actual projection | Exact-orbit-position reference, same `A` | Actual minus reference |
| --- | ---: | ---: | ---: |
| 30 | -6.545 | -25.977 | **+19.432** |
| 47.5 | -6.906 | -26.401 | **+19.495** |

Thus the large negative raw projection is not evidence of a restoring
response: in this particular fixed-`A` counterfactual, subtracting the
reference gives a positive projection. The mass-weighted mean magnitude of
the actual affine force in the inner disk is 9.41e-5 and 1.31e-4 of local
gravity at t30 and t47.5, respectively; those small absolute ratios can still
matter for a cold disk. See [receipt](KEPLER_Q2_AFFINE_REFERENCE_2026-10-29.json)
and [reader](audit_kepler_q2_affine_reference_20261029.py).

The reference changes all site positions while preserving the actual-state
`A`, so it is not necessarily a dynamically realizable affine state. The
analytic circular motion is also not the actual outer-boundary model, though
all face neighbors of the inner diagnostic set are included. It is a useful
test of baseline contamination, not a definitive linear-response calculation.
Do not interpret the previous negative numbers as a predicted cure; a proper
coupled perturbation and full angular-momentum/energy closure are still
required before implementing or evolving this model.
