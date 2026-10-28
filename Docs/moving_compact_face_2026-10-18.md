# Shared moving weighted-face transport from two density line moments

Use actual t2 generator velocities, finite weights .25*target and1*target,
and explicitly specified instantaneous Wdot=target per unit diagnostic time.
This Wdot is a kinematic probe, not an adopted evolution law. For fixed
instantaneous straight face x(u)=a+u(b-a), derive its normal speed from the
time derivative of the weighted bisector:

    wn(u)=[Xj.vj-Xi.vi-x(u).(vj-vi)+(Wdot_i-Wdot_j)/2]/|Xj-Xi|.

Thus wn=w0+w1*u. The differentiated bisector identity at both endpoints and
midpoint holds within4.35e-17 absolute in this audit.

On a segment with constant donor velocity v, define I0=integral rho ds and
I1=integral rho*u ds. Its mass flux is (v.n-w0)*I0-w1*I1; multiply by v and
|v|^2/2 for momentum and kinetic flux. Thermal energy density is constant in
each original donor and its advective flux is integrated analytically against
the affine relative normal speed. No pressure-work flux is included.

| Weight fraction | Eligible unique internal faces | Retained | Segments |
|---|---:|---:|---:|
| .25 | 3705 | 2760 | 4364 |
| 1 | 3705 | 2750 | 4353 |

Two-line-moment contraction matches direct quadrature within5.20e-19 absolute.
Order8 versus12 differs by at most4.67e-15. Evaluate one flux per retained
face and assign opposite signs to the two receiver cells. Global accumulated
mass/momentum/kinetic sums are at most1.20e-18 absolute; thermal sum is below
2.02e-22. This is shared-owner conservation by construction plus a summation
check, not independently computed MPI-rank agreement.

## Unfinished interface gate

945/955 faces are excluded for ambiguous coincident donor boundaries or
incomplete partition. The diagnostic is not a closed full-mesh RHS. At a face
lying on a donor discontinuity, choosing either donor arbitrarily is not a
physical interface solution. Resolve both one-sided states and their relative
normal propagation, including pressure work and a consistent Riemann/contact
treatment, before claiming full-face coverage or physical conservation.

The compact representation handles known donor advective traces efficiently
in algebraic size, but construction cost, persistent donor history, topology
changes, pressure/thermal closure and approved original-model reduction remain
unproved. No production code, physical trajectory or original energy convention
was changed. No long-time stability or CUDA/OpenMP/MPI parity claim.

Evidence: `canonical_moving_compact_face_2026-10-18.json`.
Script: `audit_moving_compact_face.py`.
Remote: checkpoint remap directory, `moving_compact_face.json`.
