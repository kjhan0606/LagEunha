# Unlimited-LSQ kinetic recovery fails exact-state pressure admissibility

The preceding net-mass ALE test identified `Kvar_dot` in the first-RHS proxy
`Edot-ubar dot pdot`. A subsequent bundled Astra review selected one inexpensive
pressure-admissibility gate for a specified kinetic-energy reconstruction.
The gate is complete. The specified unrestricted reconstruction is rejected
before any new RHS, gradient derivative, spectrum, or cold-flow replay.
The full physical-accuracy and speed objective remains open.

## The tested definition is fixed

The proposed recovery uses the existing unlimited face-weighted LSQ gradient
on cell-centroid differences and the geometric central second moment Q:

    Krec = |p|^2/(2m) + (rho/2) trace(G Q G^T).
    Prec = (gamma-1)(E-Krec)/V.

No independent gradient/moment state is evolved. This is not the old
13-moment dynamical closure, and the old Jordan-chain failure is not applied
to the new expression. The present issue is simpler. An intercell gradient
can infer variance absent from a valid piecewise-constant cell state.

## An exact positive Riemann state becomes negative before a flux

The native C test uses equal-size 1D Voronoi intervals with rho=P=1 and
Mach number eight. The velocity has one compressive jump from +a to -a with
constant exterior continuation. Each cell is exactly piecewise constant, so
the supplied conservative energy is

    E = P h/(gamma-1) + rho h a^2/2.

The original native `gfs_velocity_fit_add_face` and `gfs_velocity_fit_solve`
produce `|G|=a/h` in the two jump-adjacent cells. Their exact geometric second
moment is `Q=h^3/12`. Consequently the recovery gives

    Prec = P - (gamma-1) rho a^2/24.
    Prec/P = 1-gamma(gamma-1) Mach^2/24.

The stored native gamma is 1.6666666269302368. The independent expression
predicts -1.9629627157140703. Native C evaluation gives
-1.9629627157140677 at both 16 and 32 cells. Exactly two cells are negative
in each fixture. The supplied physical pressure remains one to roundoff.
The prediction and implementation differ by at most 2.7e-15.

The jump-adjacent pressure ratio is independent of h. This is a local
admissibility failure, not a claim that the volume-weighted norm of a fixed
number of bad cells remains constant under global refinement. A single valid
state with negative recovered pressure already fails the declared gate.

Constant-velocity controls recover pressure one. Affine controls supply the
exact affine integrated kinetic energy and also recover pressure one to
roundoff. These controls distinguish a formula/implementation defect from the
admissibility failure. No original Hopkins projection experiment is repeated.

## Scope of the rejection

The supplied m, p, and E are unchanged. No thermal floor, energy reset,
kinetic-energy cap, new limiter, or reservoir is used. The negative recovery
is recorded rather than repaired. The current unlimited-gradient algebraic
deconvolution is removed from production consideration.

This result does not reject every admissible or limited reconstruction.
A different gradient or limiter would define a separate candidate and require
its own conservation, affine accuracy, positivity, physical-flow, and cost
checks. The pressure interpretation cannot be switched after measurement to
count the present candidate as a pass. No Gdot/Qdot compensation or new PDE
is released by this gate.

Authoritative native evidence is
`/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/kinetic_recovery_admissibility_v1/stdout.json`.
Source snapshots, build commands, native-fit hash, and results are archived
under `Research/KineticRecovery_20261105`. The original fluid executable,
accepted states, and manuscript benchmark remain unchanged.
