# Spin closure requires a mass centroid, not just polygon geometry

The isolated face ledger in [gfs_spin_face_exchange_20261027.h](gfs_spin_face_exchange_20261027.h)
does not yet define a full cell evolution. In particular, the polygon moment
prototype assumes a uniform cell density, whereas the ledger evolves the
**mass centroid** independently. Those assumptions need a common closure
before combining them in the native solver.

For a moving cell with boundary velocity `w`, fluid velocity `v`, density
`rho`, mass `m`, momentum `p`, first mass moment `M = integral rho*x dV`,
and centroid `C=M/m`, Reynolds transport gives

```
dm/dt = - integral_face rho (v-w).n dA
dM/dt = p - integral_face x rho (v-w).n dA.
```

The sign uses the outward normal. The `p` volume term is essential: the
face ledger is only a flux substep and intentionally does not contain it.
If the W=0 fixed-mass reduction is implemented with **zero face mass flux**,
then `dC/dt=p/m`. This is a conditional statement about that particular
closure, not a claim that fixed net cell mass alone forces every face flux
to vanish; nonzero face fluxes can cancel in the mass sum while contributing
to the first moment.

If the reconstructed density is uniform in a cell, the mass centroid must
also equal its polygon area centroid `G`. Consequently, zero face mass flux
requires `dG/dt=p/m`. A Voronoi generator moving at `p/m` does not in general
guarantee this relation: evolving faces and neighbors also control `G`.
Equating `C` to `G` after each step without a consistent first-moment update
would silently change the mass moment and, unless compensated, angular
momentum. Thus the fast uniform-polygon inertia calculation alone cannot
close the spin equations under the fixed-mass condition.

Candidate resolutions to evaluate are: an internally nonuniform density or
mass-moment representation compatible with the fixed mass, a face transport
model whose mass fluxes and first moments respect the approved W=0 limit, or
a different explicitly stated rigid-cell closure. Each requires a matched
pressure/work and angular-momentum update and a cost measurement; none is
selected or implemented here. In particular, prescribing viscosity or
friction is not a solution for the user's Kepler `Re=∞` test.

For a specified physical inviscid momentum flux, the total angular momentum
of a cell is `J = C × p + S`, with intrinsic spin `S`. Its boundary flux
must follow from the same momentum traction used for `p`, and its total
energy flux must include the corresponding work. The variational generator
forces alone admit different local traction decompositions, so they do not
uniquely prescribe that spin source. An arbitrary added torque would not
establish a physical solution merely because a global ledger balances.

The already validated two-cell ledger and polygon-moment GPU kernel remain
standalone algebraic and geometry tests. This observation does not alter the
running simulations, which currently do **not** contain this spin model.

CPU-only Slurm **437500** was submitted to measure this condition on two
already audited native snapshots: aligned W0 at `t=20` and 5%-jitter W0 at
`t=6.5`. The [diagnostic](audit_kepler_w0_centroid_kinematics_20261027.py)
reconstructs each Voronoi diagram after small forward and backward
displacements along the saved generator velocities; it compares the inferred
area-centroid velocity with the saved bulk/generator velocity. It also checks
a uniform-translation example and convergence under halved displacement.
The analysis imposes no new forces and does not alter or restart the native
simulations. Any observed mismatch diagnoses the compatibility of the
specified uniform-density/no-face-mass-flux closure; it alone does not prove
that the existing native scheme's long-time error is caused by that mismatch.

Job **437500** completed `0:0` and both tests passed their uniform-translation
and finite-difference convergence checks. The mass-weighted RMS of
`|dG/dt-v_generator|/|v_generator|` was `3.73518e-4` in the
[aligned t=20 state](KEPLER_W0_CENTROID_KINEMATICS_2026-10-27/aligned_t20.json)
and `4.01597e-4` in the
[jitter t=6.5 state](KEPLER_W0_CENTROID_KINEMATICS_2026-10-27/jitter_t65.json).
For initial radii `0.5–0.75` versus `1.5–2`, those RMS values were
`7.4137e-4` versus `2.1144e-4` for aligned and `8.3578e-4` versus
`2.2130e-4` for jitter. These are different snapshots at different times,
so they are not a controlled estimate of the effect of jitter. The result
establishes a nonzero instantaneous kinematic mismatch for these saved
states; it does not establish an accumulated spin error, a heating rate, or
a diagnosis of the native solver's Kepler instability.

CPU job **437508** applied the identical diagnostic to the accepted
[aligned t=6.5 snapshot](KEPLER_W0_CENTROID_KINEMATICS_2026-10-27/aligned_t65.json)
and passed its numerical checks. At the same time the aligned versus jitter
relative RMS mismatch is `3.44700e-4` versus `4.01597e-4` over the disk,
and `7.41931e-4` versus `8.35782e-4` in the inner `0.5≤r<0.75` band.
The jittered configuration is higher by about 16.5% over the disk and 12.6%
in the inner band for these two realizations. This comparison uses the same
time and diagnostic, but two initial conditions do not establish a
statistical trend or causation of heating.
