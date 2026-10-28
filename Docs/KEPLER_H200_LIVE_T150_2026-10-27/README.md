# Live 384² periodic W0 snapshot audit

While H200 job **437076** continued, CPU job **437570** audited its complete
`t=1.5000266760004226` snapshot against the original 147,456-site initial
state. All 15 native long-run checks passed. Disk mean absolute density change
was `1.33690e-4`, mean absolute log-entropy change `2.39716e-4`, relative
total-energy difference `-2.13734e-10`, and disk angular-momentum relative
change `-4.77013e-13`. In the `0.8<r<1.2` annulus the literal GIZMO-style L1
density error against 1 is `0.0100000`; against this IC's 1.01 plateau it is
`0.00018124`. The signed state change is not a measurement of irreversible
heating. [Snapshot gate](audit.json), [radial audit](radial.json), and
[source data hashes](data.sha256) record the inputs.

The audit's zero-injection ledger is the most recent cumulative ledger read
from the still-running log, so it can include steps later than this snapshot.
All snapshot quantities refer to the recorded exact time. The 15-check gate
is a preliminary stability check, not evidence of `t>6` or the `t=600`
GIZMO endpoint. The active variational pressure branch and finite halo are
documented [here](../KEPLER_POLAR_OPERATOR_IDENTITY_2026-10-27.md).

The `t=1.50` checkpoint was written at 11:54 KST; the live run reached about
`t=1.666` near 14:01 KST. That observed average is approximately 0.079 time
units per wall hour. A single smaller geometric-CFL step does not predict
the whole run's speed; the instantaneous step size varies. No restart was
launched and the existing computation was not interrupted.

In a later log window spanning 706 steps, `t=1.61821→1.66911`, the median
step was `7.64202e-5` and the 10th percentile `6.76397e-5`, although the
minimum was `1.68989e-7`. This window explains why a single small step is a
poor time-to-completion estimate; it is not a future performance guarantee.
