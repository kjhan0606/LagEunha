# Affine flow: a scalar power weight cannot make every face point material

2026-10-21. Driver's bounded native kinematic check while Astra reviews the
previous rejected spatial bundle. This is a necessary geometry observation,
not a new hydro solver or a claim that general Laguerre ALE is impossible.
Sources, executable, reports and SHA checks are frozen in
`affine_face_kinematics_v1/`; the calculation is terminal. No IC, restart,
generator motion or existing runtime was modified.

## Derivation at W=0

Let fluid velocity be u(x)=B x+b and generator velocities exactly u(q_i),
u(q_j), satisfying the user's no-steering condition. Set q_j-q_i=d n and
m=(q_i+q_j)/2. A point x=m+s t on their initially Voronoi face has n.t=0.
Let z_i=Wdot_i. Differentiating the actual power-bisector support gives

    w_geom,n(x) = v_j.n - (v_j-v_i).(x-q_i)/d + (z_i-z_j)/(2d).

The native screened weight-rate solver already uses this formula. Therefore

    u_n(x)-w_geom,n(x)
       = 2 [(sym B)n].(x-m) - (z_i-z_j)/(2d)
       = 2 [t.(sym B)n] s - (z_i-z_j)/(2d).

The weight-rate difference changes the constant term, not the tangential
slope. For nonzero t.(sym B)n no scalar rate difference can make the entire
planar face material. Pure rigid rotation and isotropic dilation have zero
slope. An individual face aligned with a principal strain direction can also
match; the statement is about general strain/orientation, not every face.

For uniform density on a centered face of length L, optimal constant shift
removes the mean mismatch but leaves

    RMS_min = abs(2 t.(sym B)n) L/sqrt(12).

Even with zero integrated mass flux, spatial transport remains:

    integral_face rho (u_n-w_geom,n) ds = 0,
    integral_face s rho (u_n-w_geom,n) ds
         = rho [2 t.(sym B)n] L^3/12.

Thus a matched cell mass/volume constraint is not proof that all internal
position, momentum, rotation or energy transport vanishes.

## Native evidence

The C++ test independently samples the production fixed lab Euler HLLC/HLL
fan at the geometric normal speed, using identical exact fluid states on
both sides. It checks all four Euler ALE flux components against their
analytic values. A positive four-point face rule exactly integrates these
polynomial examples. Five named controls and2000 random B/orientation/size/
translation cases complete:

| Check | Maximum absolute/normalized error |
|---|---:|
| affine normal-speed formula | 6.11e-16 |
| native complete Euler ALE flux | 3.42e-16 |
| minimum RMS expression | 2.23e-16 |
| integrated mass flux, uniform density | 4.49e-16 |
| tangential first mass-flux moment | 1.39e-16 |

Example: Bxy=1.2, all other B components0, face normal x, L=.4,d=.7,rho=1.3.
Best scalar Wdot difference is0, residual speed RMS=.1385640646, net mass
flux=-3.47e-18, but tangential first mass-flux moment=.00832. Rigid rotation
and isotropic dilation controls match to roundoff.

## Consequence for the next design

A scalar-weight Laguerre ALE method may still transport the required moments
with complete conservative fluxes. The result rejects an assumption of
pointwise material planar faces under generic shear when the generators move
with the prescribed fluid field. It does not reject scalar weights as a way
to control integrated masses or volumes, and does not diagnose the entire
historical Kepler instability by itself.

Any next affine/material-cell formulation must explicitly handle this
nonmaterial-face transport or supply different consistent shape kinematics.
Adding only rotation or pressure work cannot silently drop it. Existing
canonical face moment/angular flux helpers can carry this transport; an
adequate thermal/velocity closure and its actual stability are still needed.
The strict W0 reduction, no steering, Re=infinity, old failure interval and
longest GIZMO objective all remain unchanged.
