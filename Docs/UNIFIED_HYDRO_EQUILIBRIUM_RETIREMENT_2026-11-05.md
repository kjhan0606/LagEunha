# Uniform-pressure control closes the current strain-weight candidate

The completed mode audit found a large native acoustic frequency for the
candidate `R_i(v)=sum_b psi_i^b G_i(v):S^b`. A final bundled Astra evaluation
selected one minimum equilibrium control before retiring the construction.
The control is complete. The current construction is removed from production
consideration. Source, diagnostics, and accepted original trajectories are
preserved. The result is not a no-go theorem for every `C=B+LR` construction.
The complete user objective remains active and unfulfilled.

## Uniform pressure removes the missing-background-term ambiguity

Use the same 16,384 native cells, original positive fixed masses, and W=0.
Define a diagnostic rest equilibrium with

    v0 = 0.
    P_i = P0 = 1/gamma.
    U_i = P0 V_i/(gamma-1).

The numerical value of P0 is 0.6000000143 because the native gamma is stored
as float32 before conversion to double. Every diagnostic internal energy is
positive. No equilibrium checkpoint or trajectory is written. The original
Gresho state is not changed.

The geometric operator satisfies `C^T 1=0` identically for closed periodic
cells because total volume is fixed and `L 1=0`. Differentiating the identity
gives

    (delta C)^T P0 1 = 0.

At rest the term `delta C v0` also vanishes. The exact inviscid pressure
linearization of the proposed semidiscrete model is therefore

    delta(vdot) = M^(-1) C^T delta(P).
    delta(Pdot) = -H0 C delta(v).
    H0 = diag(gamma P0/V).

The positive acoustic block is the complete pressure linearization at this
equilibrium. The statement excludes numerical viscosity and does not qualify
the nonlinear solver. The measured uniform-pressure force defects are
1.25e-16 for B and 2.74e-16 for C. Relative to the C column norm, the defect is
8.15e-17.

## Conservative inequalities need no further spectrum calculation

Define the rectangular acoustic factor

    T_B0 = H0^(1/2) B M^(-1/2).

Two analytical bounds on the largest B squared frequency are

    lambda_B0,max <= ||T_B0||_1 ||T_B0||_infinity.
    lambda_B0,max <= max_i [|T_B0| (|T_B0|^T 1)]_i.

The second inequality follows from the entrywise absolute-product bound on
`T_B0 T_B0^T`. The earlier B Ritz value is not treated as a certified upper
bound.

The actual-pressure independent C Rayleigh value is q=2.2595997020e10. Since
the original pressure is positive, the equivalent velocity acoustic blocks
obey a positive-semidefinite matrix ordering. Consequently

    lambda_C0,max >= (P0/P_max) q.

This comparison uses the saved Rayleigh value and does not need another
coefficient solve, eigenvector, or eigensolve. A fractional 1e-9 padding lowers
the C estimate and raises the B estimates. The original pressure range is
3.36375--7.20255. The evaluated quantities are:

| Quantity | Squared-frequency estimate |
|---|---:|
| B upper bound, norm product | 1.985780325e6 |
| B upper bound, absolute-product row sums | 1.450670897e6 |
| C lower bound, pressure comparison | 1.882333338e9 |

The tighter B bound gives

    omega_C0,max / omega_B0,max >= approximately 36.0217.

The analytical inequalities are conservative. Their floating evaluation and
padding are reported explicitly and are not an interval-arithmetic certificate.
The difference is already large, so the bounded stage terminates without
additional spectral matrix applications. No native coefficient or PDE solve
was executed.

## Decision and scope

Large candidate stiffness persists in a rest equilibrium with uniform pressure.
The stiffness cannot be attributed solely to the evolved Gresho pressure
deviation or to omitted background geometric force terms. The local strain
construction also worsens measured pressure-force accuracy and fails the
previous bounded CPU coefficient screen. The combined evidence meets the
predeclared retirement condition.

The current `psi_i times fitted strain` construction will not receive a warm
start, a larger iteration cap, a relaxed residual tolerance, an added limiter,
or an implicit-pressure rescue. No long candidate PDE is launched. This closes
one proposed model and preserves the physical constraints and speed objective.

The frequency is not a growth rate or a measured GIZMO runtime factor. Standard
undamped KDK on the frozen harmonic pressure system has a `dt omega <= 2`
condition. The condition is not a full nonlinear/viscous stability theorem.
Future constructions must control the full acoustic norm on deformed cells
with normal volumes, in addition to affine/long-wave consistency and an
affordable coefficient construction. Trimming selected cells does not supply
such a structural argument.

Authoritative evidence is
`/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/material_volume_equilibrium_v1/results/result.json`.
Compact evidence and source are in `Research/EquilibriumBound_20261105`.
The manuscript keeps the verified previous benchmark. A physically qualified
common dimensional scheme, the remaining 2D suite, and matched GIZMO throughput
remain required. None is replaced by this equilibrium diagnostic.
