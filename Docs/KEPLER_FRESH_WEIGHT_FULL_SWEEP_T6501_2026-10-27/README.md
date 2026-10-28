# Full weighted-minus-Voronoi face-speed correction at t=6.5

The earlier [weight-sweep estimate](../KEPLER_FRESH_WEIGHT_LOCAL_T6501_2026-10-27/README.md)
used only the changing weight difference. A power face also shifts relative to
the Voronoi face at identical generator positions because the finite weight
difference changes as the generator separation changes. For
`d=|x_j-x_i|`, `ΔW=W_i-W_j`, and `Δv_n=n·(v_j-v_i)`, the normal relative speed
is

```
u_weight - u_Voronoi = ΔWdot/(2d) - ΔW Δv_n/(2 d²).
```

The [full-speed audit](full_speed_sweep.json), CPU Slurm **437597**, includes
both terms. An independent 5,000-case finite-difference check of the radical
plane derivative agrees to maximum relative error `1.35e-6`. On the native
`t=6.5→6.501` endpoint, the velocity-derived and snapshot-secant generator
separation rates agree to RMS `5.84e-9`; midpoint power areas agree with
averaged native volumes to `2.92e-6` relative.

The finite-weight separation correction has RMS `2.23619e-6`, compared with
`8.70862e-5` for the weight-rate term, about `2.57%`. Including both changes
the midpoint arithmetic-density sweep estimate only from `8.68398e-6` to
`8.68095e-6` mean absolute density change; using upwind donor density gives
`8.68096e-6`. All pairwise mass sums cancel to below `2.0e-22`, and all
remapped masses remain positive. Thus the missing kinematic term is real but
does not change the earlier conclusion for this step.

This remains a first-order midpoint geometric estimate. It does not implement
the native RK4 mass, momentum, energy and spin fluxes or establish a
self-consistent variable-mass evolution.
