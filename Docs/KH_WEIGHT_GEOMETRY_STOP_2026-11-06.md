# KH dynamic-weight stop: admissibility boundary, not a vanished cell

Diagnostic-only source changes:
- Native CUDA reports weight_clamps separately from other geometry errors.
- Native weight KDK writes the restored step-start state before exhausting
  retries. The rejected trial is never presented as an accepted state.

Replayed only the failed segment from its existing t=1 checkpoint using the
same inputs, physics and timesteps. Last accepted time, energy, Vmin, Umin,
step count and dt match the original failure log. The new result is rc99 at
exactly the same place. All16 geometry hits are weight clamps; other geometry,
volume, velocity-fit and halo errors are zero.

Stored diagnostic state time =1.0731293139821652.
The existing native geometry requires raw face fraction in [.05,.95], where
fraction=.5*(1+(W_i-W_j)/d_ij^2), equivalent to |W_i-W_j|/d_ij^2<=.9.
Saved-state all-pair periodic analysis finds maximum ratio .8999998928338196.
Its outward instantaneous rate is9.33060265318 per unit time. Thus repeated
step shortening approaches a finite geometric admissibility boundary rather
than removing the cause. The rate is evaluated from saved v and pi/mu; it
is not a prediction that the rate remains constant after the stop.

The implicated pair is separated by .03711626. All cell volumes remain
positive (minimum .00375435178243) and thermal energies positive (minimum
.0565111302709). This is NOT evidence of cell disappearance or of failure of
the exact mathematical Laguerre tessellation. The user-stated sufficient
survival condition |Delta W|<d^2 is weaker than the native safety margin .9.
Do not conflate those two bounds, and do not relax .9 just to pass a test:
the native search/geometry implementation uses the margin in its assumptions.

Implication for the physical weight law: the unconstrained canonical
auxiliary dynamics does not preserve the current admissible geometry set.
A physical kinematic closure must incorporate geometric admissibility while
retaining its compatible pressure-work relation. Clipping W after a step or
clipping Wdot while leaving the force unchanged is not a derived conservative
solution. This evidence does not disqualify every evolving-weight method.

Artifact: Research/MaterialVolume_20261105/kh_weight_failure_20261106_v1.
results.json, constraint_metrics.json, replay/log and the native stop snapshot
are retained. analyze_constraint.py examines all minimum-image pairs, not a
certified actual face-neighbor set; the independent device clamp counter
establishes the failing category. The stop snapshot is diagnostic; validate
its reader/restart contract before any continuation. No simulation is live.
