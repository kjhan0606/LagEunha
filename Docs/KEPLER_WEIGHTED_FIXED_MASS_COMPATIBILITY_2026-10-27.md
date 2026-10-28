# Fixed-mass Laguerre volume and energy compatibility

For the current native finite-weight variational branch, let `q` be generator
positions, `v=qdot`, `W` the power weights, and `V(q,W)` the actual cell
volumes. The verified discrete geometry gives

`Vdot = B v + J Wdot`,

where `B` is the generator-motion volume operator and `J` is the symmetric
power-area Jacobian. The native pressure force is `F=B^T P`, and the
fixed-mass adiabatic internal-energy rate is `Udot=-P^T Vdot`. Consequently,
ignoring the separately balanced external gravity,

`d(K+U)/dt = v^T F - P^T Vdot = -P^T J Wdot`.

The present global energy-neutral projection enforces `P^T J Wdot=0`. Exact
volume freezing would instead require `J Wdot=-Bv`. These two requirements
are compatible only when `P^T Bv=0`. They hold together for a spatially
uniform pressure in a closed domain, because total volume is constant. Once
pressure varies, exact cellwise volume freezing and this particular
fixed-mass force/energy closure generally cannot both hold. This is a
compatibility limit of the chosen closure, **not** a proof that a conservative
moving-mesh or weighted hydrodynamic solution does not exist.

The [read-only t=7.426 projection comparison](KEPLER_Q2_WEIGHTED_PROJECTION_METRIC_T7426_2026-10-27.json)
shows the limit numerically. For the independently jittered IC at 64 Jacobi
iterations, the relative volume-rate residual is `0.00832` before the energy
projection, `0.04143` after the current Euclidean projection, and `0.03680`
after a more geometry-aware projection. The latter needs a second elliptic
direction solve; its exact CPU reference took 499 preconditioned CG
iterations. Its predicted one-step disk area-change p99 is `2.74e-6` versus
`3.06e-6` for the current projection. A different projection alone therefore
cannot erase the conflict or justify that cost.

The current finite-weight branch remains useful as a bounded mechanism test:
it maintains fixed masses, positive power cells, and a global pressure-work
identity. The [native `t=10` comparison](KEPLER_Q2_WEIGHTED_NATIVE_V2_T10_RESULT_2026-10-27.md)
found a 27.7% reduction in the jittered disk density error relative to W0,
but still failed the `1e-3` density gate. A general conservative
extension must account for work carried by weight changes or allow
equal-and-opposite mass, momentum, and total-energy flux through the actual
moving power faces. The latter is consistent with the user's conditional
acceptance of weight-linked mass variation. The static `Wdot=0`, `W=0`
Voronoi/fixed-mass mode remains the required reduction case. Any new flux
branch needs a separate local conservation and shock gate before a long
Kepler claim.

## A weight-only conservative remap does not automatically cure density noise

One tempting local closure is to let changing weights sweep mass across the
power face by exactly the volume they exchange. In a constant fluid state,
write the swept cell-volume increment as `delta V_i` and let the remap
transfer `delta m_i = rho * delta V_i` (and the corresponding momentum and
total energy). Then

`(m_i + delta m_i) / (V_i + delta V_i) = rho`

exactly whenever `m_i = rho * V_i`. Thus weight-driven **conservative
relabeling of a constant state leaves density unchanged**. More generally,
to first order the remap's density response is

`delta rho_i = (delta m_i - rho_i * delta V_i) / V_i`.

If the face transfers material at density close to the receiving cell's
existing density, the remap cancels the apparent density benefit from
changing its volume. In particular, applying such a remap *after* the
accepted W0 `t=6.5` jitter history cannot erase the density noise already
present there merely by restoring target power areas. This is an algebraic
screen of a **pure geometric remap**, not a claim that every conservative
weighted hydro flux must fail: pressure, velocity and entropy fluxes may
change the subsequent physical evolution. A new variable-mass branch must
demonstrate a reduction of the measured pressure/orbit defect beyond this
constant-state remap identity. The fresh weighted-from-`t=0` calculation
and same-checkpoint W0/weighted controls are therefore more informative
than implementing a remap on the assumption that mass variation itself
will stabilize the disk.
