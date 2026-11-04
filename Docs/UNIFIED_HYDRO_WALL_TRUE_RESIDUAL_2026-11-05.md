# Constrained reflecting velocity and true-residual stopping

## Status

The preceding goal turn made progress by publishing source/evidence and
validating the optional block preconditioner on an exact failed trial.
The present stage completes the conditional block trajectories, identifies
the reflecting-wall failure, implements the grouped Astra plan, and obtains
new full-target strong-shock failure evidence. The common method is still
experimental. The original nine full trajectories, common 2D suite,
production GPU/MPI evolution, and measured GIZMO cost <=2x remain unfinished.

No new result is transferred into the manuscript's accepted Kepler baseline.
No existing Kepler calculation is restarted. All new work is native C/OpenMP
on GPFS. Python orchestrates and analyzes results. The archive is
`Research/WallTrueResidual_20261105` and the original run root is
`/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1`.

## Completed block stage and the actual wall cause

The exact version4 collision-trial block test passed before the two common
block-backend trajectories were launched from birth. Blast stops at
t=0.005558214814 of 0.038 and collision at t=0.01817775922 of 0.035.
Collision reaches the 64-iteration cap with recursive relative residual
1.082392751e-12 but true scaled residual 8.070361495e-14. The exact rejected
matrix has condition 53.8103 under diagonal Jacobi and 34.6526 under block
Jacobi. Both numbers are read-only small-matrix diagnostics, not production
dense solves or trajectory error norms.

The previous diagonal blast accepted checkpoint was used to reproduce only
the next failed step. No accepted step advances and final.fv is byte-identical
to the original accepted checkpoint. All nine trial sizes fail the first
Euler geometry check. An independent exhaustive 1D halfspace audit finds
no hidden cell and positive minimum width 0.000542125. The generator q0
instead crosses the physical reflecting wall. The accepted q0=1.72028e-8
has qdot0=-1.41009, and the last rejected q0=-3.76513e-9 lies outside x=0.

Cell0 has mean velocity 0.00196211097 and centroid 0.00677807538, but the
limited gradient 208.3266203 extrapolates a normal wall trace -1.41009143.
Reflecting ghost means alone do not make the reconstructed wall trace zero.
The reported failure is therefore a physical boundary-field inconsistency,
not demonstrated cell disappearance or owner-transition chatter.

## Grouped plan selected by Astra

The grouped evaluation selects two changes and retains the physical
mass--weight law, birth mu, 64-iteration cap, actual residual tolerance 1e-12,
stress coefficients, CFL, rollback cap, and absence of repairs/steering.

For a reflecting plane x_n=x_wall the constrained affine normal field is

    u_n(x) = mean_u_n * (x_n-x_wall)/(centroid_n-x_wall).

The field preserves the cell mean and vanishes over the entire wall plane.
The normal gradient row is fixed before constructing the stress sensor.
The existing limiter acts only on free velocity rows. Corners constrain
independent rows. A cell touching both opposite walls requires zero normal
mean for an affine field and reports infeasibility otherwise. The point
query evaluates the distance ratio directly. Generators are never clamped
or repositioned. Moving transmissive boundaries do not receive wall constraints.
The physical candidate has a separate checkpoint version5. Version4 accepted
states are read only for causal operator diagnostics and are not production
restarts under the changed reflecting field.

The raw constrained gradient enters the unchanged chi formula. The correct
data identity at constrained cells contains a cross term,

    Edata = Efit + Eresidual + 2*sum omega*(delta_u-G*r) dot (G*r).

The old Pythagorean equality is not asserted at wall cells. The free rows
satisfy their original least-squares normal equations. The filter remains
bounded between zero and one without interpreting the two nonnegative terms
as an orthogonal energy decomposition at a constrained boundary.

For PCG the final criterion is the actual scaled flux residual of the
returned z. The recursive preconditioned relative norm triggers a true
residual check and reliable replacement. A true check also occurs at cap64.
No new iteration, block-size search, RHS mean subtraction, or tolerance
relaxation is introduced. The historical dual stopping rule is explicitly
changed and historical failures retain their original classification.
NaN residuals cannot pass the norm through fmax's NaN suppression.

## Preflight evidence

Sixteen native wall/corner fixtures pass across 1D/2D/3D. The fixtures cover
compatible affine fields, weighted nonaffine fields, free-row normal
equations, mean preservation, wall traces, the cross-term identity, and
kinetic/entropy ledgers. Real strip cells touching both opposite walls test
both feasible zero means and explicitly rejected nonzero means.
Saved-blast operator inspection has wall-trace error 2.05e-19 and unchanged
conserved/birth state. These inspections produce no accepted trajectory.

Two fixed numerical controls pass with unchanged A, b, scales, and cell
diagnostic bytes. The later collision trial returns its actual cap64 z with
true residual 8.070361495e-14 and stop reason true_residual_at_cap. Rebuilding
the earlier blast Euler proposal from the accepted checkpoint and the logged
dt reproduces the rejected candidate residual 2.575026732e-12. The candidate
is not accepted. One reliable replacement produces true residual
7.800464022e-13 after 22 iterations. Independently recomputed residuals of
returned z match the reported values exactly. Nonfinite residual rejection
is checked separately.

The final rebuild adds only a clearer infeasibility diagnostic to the native
operator. Its 16 fixture stdout records are byte-identical to the preflight
records. The final source/build hashes and actual stress-point dump are
retained in wall_velocity_final_build_v1. No successful PDE is manufactured
by repeating the failed full trajectories after the message-only edit.

An additional audit strengthens the affine check beyond the raw gradient.
The check also evaluates the limited gradient and actual generator velocity,
and adds six compatible affine fixtures on weighted geometry. All 22 final
fixtures pass in wall_velocity_affine_audit_v2. These extra tests do not
change the production physical operator or repeat a PDE trajectory.

The default HLLC branch is checked separately at Sod N32, t=0.2. An initial
helper refactor changes final pressure by at most 1.89e-15 and fails the
strict byte-comparison control despite identical step/retry counts. Retaining
the historical trace operation order outside the wall opt-in restores
byte-identical initial/final checkpoints, final table, and node velocities.
Both control attempts remain archived. The default test is not version5
strong-shock validation.

## Full-target version5 trajectories still fail

Both original tests use the same fixed binary from birth and retain their
full original target, guard construction, and ledgers.

| Test | Original target | Last accepted t | Accepted steps | Actual failure |
|---|---:|---:|---:|---|
| blast N200 | 0.038 | 0.006417387664 | 3125 | true residual 1.114485044e-12 at cap64 |
| original collision, physical N200 plus guards | 0.035 | 0.02061232856 | 1159 | true residual 1.168124409e-12 at cap64 |

Blast's failed Euler dt is 3.530697417e-12. The current bound is not relaxed.
Blast's accepted minimum fixed-wall generator margin is 0.00249994734, so
the earlier wall crossing is absent. The exhaustive minimum cell width is
0.000558322 and there are no hidden cells. Collision also has positive
intervals and no hidden cells. Its 273 cap-stop calls legitimately satisfy
the true criterion before the later genuinely failed trial.

The maximum conservation and mass--weight residual-invariant errors remain
small, but the small ledgers do not establish shock accuracy or practical
speed. Both runs exit nonzero. The failed stages are flagged non-restartable
diagnostics distinct from the last accepted final.fv. The conditional N400,
original-nine, cold, and Gresho continuation stages are not launched after
these strong-gate failures.

## New interface and length-scale evidence

Blast's nearest generators 3 and 4 remain in different actual cells and
approach to 2.730828133e-10. Their physical centroid separation is
0.01344604213, about 4.9238e7 times larger. The cell-mean velocities differ
by only -2.325097448e-9, but the two limited reconstructions give shared-face
normal traces 7.393975260 and 6.565504208. The selected generator gap rate
is -0.828471040. Correct containing-cell ownership alone therefore does not
remove the discontinuity between two different containing cells.

The four neighboring mean velocities are 4.596806251, 6.834589103,
6.834589100, and 7.926099242. The central cells form a local maximum and
minimum. Their positive reconstructed slopes survive the current
velocity-vector radius limiter. A standard scalar monotonicity bound would
flatten both slopes in this frozen state. That algebraic observation is not
an evolved control or a selected new limiter. Coordinate dependence,
affine shear, and common multidimensional behavior must be evaluated before
choosing the next method.

There is a second distinct issue to evaluate. The raw LSQ fit and mean
velocity states use physical centroids, but build_points divides the stress
velocity jump by the generator distance f->distance. For the saved pair
chi=0.500000002 and the actual native stress is Pi=0.6429888995. Substituting
only the centroid-normal separation in the frozen stress algebra gives
Pi=1.296172352e-8. The actual native point dump agrees with the independently
computed original Pi. The substituted expression has not been evolved and
no shock improvement is claimed. The power-face motion and screened matrix
still require the true generator distance. A global replacement of
f->distance would be incorrect.

The next grouped evaluation should compare a standard monotone physical
velocity reconstruction and a stress length consistent with the locations
of its velocity states. Any more complex continuous-field alternative needs
a demonstrated necessity and an explicit cost assessment. The next stage
must retain exact affine/boundary checks, frozen causal controls, original
full strong targets, and the unchanged conservation/invariant gates. More
preconditioner, mu, cap, or tolerance sweeps are not a substitute for the
remaining field/length consistency work. The full objective stays active.
