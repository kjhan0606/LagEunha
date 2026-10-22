# Exact circular-flow deformation and the current native representation

2026-10-22. This is an analytic diagnostic from the actual IC law, not a new
trajectory, force test, timestep experiment, or prescribed evolution/steering.
No GPU work, accepted state or live build/runner was modified.

## IC and derivation

`material_shared_trace_factory.h::initial_velocity` uses
v=(-y,x)*Omega(r), Omega=(r^2+eps^2)^(-3/4), eps=.047.
For this exact circular solution:

    r(t)=r0, theta(t)=theta0+Omega(r0)*t
    s=r*Omega'(r)*t=-(3/2)*t*r^2/(r^2+eps^2)^(7/4)
    F = R(theta(t)) * [[1,0],[s,1]] * R(-theta0)
    det F = 1
    sigma_max = (sqrt(s^2+4)+abs(s))/2
    sigma_min = 1/sigma_max
    cond_2(F) = sigma_max^2

The two rotations do not change singular values. These formulas follow by
differentiating the circular map in initial/current orthonormal polar frames.
The shear-rate maximum is at r=eps*sqrt(4/3), not at r=0: derivative of
r^2/(r^2+eps^2)^(7/4) vanishes when r^2=4eps^2/3. At the exact origin the
velocity gradient approaches rigid rotation and this shear measure is zero.
At finite small r near the softening length the differential shear is large.

An axisymmetric density transported by this exact map stays unchanged;
detF=1 gives no compressional heating for the inviscid adiabatic solution.
Large physical shear/stretch by itself is not a shock or a reason to add heat.
It does make approximating this map on fixed material triangles difficult.

## Values in the project's current time/radius units

|r|abs(r Omega')|abs(s) at t8|cond(F) at t6|cond(F) at t8|cond(F) at t10|
|---:|---:|---:|---:|---:|---:|
|0|0|0|1|1|1|
|0.054270925|44.55781|356.46248|71476.343|127067.5|198541.84|
|0.1|33.450383|267.60306|40283.412|71613.4|111894.81|
|0.5|4.1778251|33.422601|630.35042|1119.0693|1747.4217|
|2|0.52981794|4.2385436|12.022275|19.915038|30.037414|
|8|0.066287257|0.53029805|1.4846035|1.6892305|1.9180324|

The near-softening region belongs to the low-density central floor in this
setup; density tails/floor remain part of the current operator and timestep
contract. These numbers do not prove which source currently limits the
measured timestep: that requires an actual minimum-bound witness. They are not
GIZMO orbital-normalization claims. Full numerical values are in the adjacent
CSV; Python evaluated only these closed formulas, not fluid evolution.

## Consequence for the pending long-production check

Actual `material_shared_trace_batch.h::st_shape` combines P3 with normalized
radial IC remainders and s^4 edge extensions. Consequently, a polynomial-only
P3 determinant bound would omit active functions. Original positive-node
checks and the new paired quadrature step are necessary numerical comparisons;
they do not prove continuous injectivity between points or indefinite lifetime
of this fixed material representation. DetF=1 also does not control its
condition number. The existing global_map_injectivity_guard_pending stop must
not be removed merely because the positive cubature or finite energy tests pass.

The current adopted cubature/resident-mass cost bundle continues unchanged.
This quantified material-shear reference is evidence for the already pending
mesh-lifetime/representation evaluation after the actual paired step. It does
not authorize another optimization sweep or adopt a new physical method.
Re=inf, no steering, no heat reset and no independent spin reservoir remain.

## Measured global dt is consistent with central rigid rotation

For the same IC, in orthonormal polar coordinates grad(v) has off-diagonal
entries -Omega and Omega+r*Omega'. Its Frobenius norm is therefore
sqrt(Omega^2+(Omega+r*Omega')^2). At r=0 this is sqrt(2)*eps^(-3/2),
although the symmetric strain tensor vanishes there. The native controller
uses eta/max_grad_v with the full Frobenius norm, not just symmetric strain.

With eta=.06, the analytic origin value gives dt=0.00043229781401251614.
The retained actual next high-rule RHS has dt_bound=0.00043229781346857082;
the relative difference is -1.25826527e-09. This very close agreement
is evidence that central rotation is relevant to the cost bottleneck; it is
NOT a recorded limiting-source witness, so it does not independently prove
which triangle and branch won the global minimum. No force or geometry pass
was replayed to obtain this comparison.

This does not justify relaxing dt: explicit native RK4 still integrates the
rapid central material motion and its gravity, even at low density. Nor is
this vanishing central symmetric strain itself evidence of artificial heat.
It separates the cost question from differential shear near the softening
scale and from the artificial-heating question. The adopted dt controller
and current physical run remain unchanged. See
KEPLER_ORIGIN_DT_COMPARISON_20261022.json for the retained numerical comparison.
