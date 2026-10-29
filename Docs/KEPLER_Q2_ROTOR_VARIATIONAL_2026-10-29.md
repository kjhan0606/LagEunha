# Frozen Q2 screening of a conservative rotor counterfactual

An Astra review was requested once for the current W0 fixed-mass, no-steering,
inviscid spin-closure problem. It proposed using the discrete Hamiltonian
`H = H_W0 + sum_i s_i^2/(2 I_i(q))` as a **comparison model** with constant
individual spin `s_i` and geometric inertia `I_i` of a uniform-density cell.
This has a consistent reaction force `F=-grad_q H_spin` and can conserve
energy, linear momentum and total angular momentum at the semidiscrete level
under the appropriate free-boundary symmetries. It is not a derived Euler
closure: spin does not exchange between cells, the true mass centroid/inertia
are not generally the uniform polygon values, and affine strain is omitted.
Conservation alone does not establish physical correctness or usefulness.

For a cell of fixed mass and uniform density, the shape derivative used here is

```
I_i = (m_i/V_i) integral_cell |x-C_i|^2 dA
dot I_i = (m_i/V_i) sum_faces integral_face
          (|x-C_i|^2 - I_i/m_i) w_n(x) ds.
```

The Voronoi normal velocity `w_n` is linear along each straight face and in
the two generator velocities. Two Gauss points per face therefore integrate
this cubic expression exactly for the frozen polygon. Force is obtained as
the negative derivative of `sum s_i^2/(2 I_i)`, keeping each `s_i` fixed.
Cell spins are **reconstructed** from the stored native velocity gradients
and polygon mass moments; they are not evolved spin state. The energy sum is
restricted to cells originally inside r=2.35; all touching neighbors receive
the resulting reaction, so the inner-disk force includes its required local
support. The setup is a diagnostic and is not plugged into native code.

| Existing Q2 checkpoint | Rotor contribution to inner-ring main-mode torque projection | Existing pressure contribution, independently reconstructed |
| --- | ---: | ---: |
| t=30 | +0.6978 / time | +0.1918 / time |
| t=47.5003 | +0.5993 / time | +0.2050 / time |

Positive projection means that this force contribution, at these frozen
states, points toward increasing the squared `0.2<=m/N<0.3` angular-momentum
error amplitude. It is **not** a measured growth rate of the combined or
evolved rotor model. At t47.5 a random-position directional derivative of
rotor energy from independent full polygon recomputation agrees with the
analytic face force to absolute error 3.52e-11 on a derivative 8.92e-6 at
the smaller tested step; t30 has error 5.86e-12 on 8.78e-6. Those are single
directions, not a proof of every force component or topology event.

The result argues against launching a long native run of this **specific**
geometric-inertia rotor candidate without further work: its instantaneous
reaction projects in the same destabilizing direction as the pressure
feedback on both saved Q2 states. Different spin initialization, dynamic
spin exchange, nonuniform mass moments and an evolved trajectory might behave
differently. It does not disprove cell rotation as a remedy or establish
causation. No new steering, viscosity, boundary condition or user simulation
was introduced. The candidate adds a new spin-energy representation; one must
not silently subtract that energy from the original cold gas.

Reproducible receipts and readers: [t30](KEPLER_Q2_ROTOR_VARIATIONAL_T30_2026-10-29.json),
[t47.5](KEPLER_Q2_ROTOR_VARIATIONAL_2026-10-29.json),
[t30 script](audit_kepler_q2_rotor_variational_t30_20261029.py), and
[t47.5 script](audit_kepler_q2_rotor_variational_20261029.py). Pressure
projection is the existing independent W0 face reconstruction, not a new
population-wide native GPU capture. O(faces) arithmetic of the formula is not
an end-to-end speed measurement.
