# Face-midpoint torque versus a linear pressure distribution

Consider a 2D polygon with uniform density and area centroid `C`, and a
linear scalar pressure `p(x)=p0+g·x`. For the normal pressure traction,

```
tau = - integral_boundary (x-C) x (p n) ds
    = - integral_cell (x-C) x grad(p) dA = 0.
```

The result relies on the area centroid and a constant pressure gradient. It
does **not** assert zero torque for nonlinear pressure, nonuniform-density
mass centroids, general stresses or an arbitrary numerical force.

For an oriented straight edge from `a` to `b`, let `e=b-a`, midpoint `c`,
and outward area vector `A=(e_y,-e_x)`. Exact edge integration of the torque
includes a term absent from simply applying the midpoint pressure force at
the midpoint:

```
integral_edge p(x) * ((x-C) x n) ds
 = p(c)*((c-C) x A) + (g·e)*(e x A)/12.
```

The pressure force itself integrates exactly at the midpoint for linear
pressure; its **moment** generally does not. The independent
[triangle test](test_kepler_pressure_linear_torque_20261027.py) agrees with
two-point Gauss integration and shows a nonzero midpoint-only residual in a
skew triangle although the complete linear-pressure torque is zero.

CPU job **437558** applied this diagnostic to 49,038 reconstructed native
disk polygons, using an independent inverse-distance weighted least-squares
fit to saved neighbor pressures. The exact-linear torque residual relative
to the sum of absolute face-torque terms was at most `3.17e-14` for the
[W0 jitter t=6.5](KEPLER_LINEAR_PRESSURE_TORQUE_2026-10-27/w0.json) state and
`2.64e-14` for the [weighted jitter t≈4.1266](KEPLER_LINEAR_PRESSURE_TORQUE_2026-10-27/weighted.json)
state. The RMS midpoint-only residual normalized by `I*Omega^2` was
`1.66e-6` and `8.98e-8`, respectively. The states are at different times,
so these are not a controlled W0-versus-weighted comparison.

This is an independently reconstructed **hypothetical** normal-pressure
traction, not the face force or torque currently implemented by the native
variational solver. The fitted gradient is not asserted to equal its native
Riemann or variational gradient. The result does not measure accumulated
torque, numerical heating, or explain the historical Kepler failure. It
does show that a proposed spin discretization must account for face-moment
accuracy; inserting midpoint forces as torques can introduce a discretization
residual even for a linear pressure field. Any practical correction must
also use equal-and-opposite pair forces and energy-consistent work, and
must be tested with the actual force closure.
