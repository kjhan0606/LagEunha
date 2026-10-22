# Long Kepler completion requirements

The user goal remains long physical stability in the user's native CUDA+MPI code, beyond the old t6–10 failure and through the longest Kepler duration discussed in the GIZMO paper. Causal short tests are intermediate evidence only.

## Published reference

[Hopkins2015, §4.2.4, Figs6–8](https://arxiv.org/pdf/1409.7395): density images are shown at t120 and t600; the text also reports integrations of approximately1000 orbits. Fig8 measures velocity/density errors in0.8<R<1.2, with comparison thresholds L1(v)=.01 and L1(rho)=.25. These thresholds are paper comparison diagnostics, not universal accuracy guarantees. Figure/text inner-orbit counts differ; retain actual t and explicitly defined orbital periods.

The reference uses gamma5/3, P1e-6, GM1, no self-gravity, and the Eq34 density shape with inner/outer transitions0.5/2. Its effective256² spacing refers to the square of side4, not any arbitrary256² box. The inspected section does not give a numerical epsilon. Record the chosen epsilon and its resolution/inner-radius ratios; do not claim an exact IC match until that ambiguity and boundaries are addressed.

## Required evidence before completion

1. **Native implementation:** frozen source/binary hashes and effective parameters. Evolution uses the existing native solver, not a Python substitute. MPI/halo, volume completeness, face reciprocity and precision checks cover the actual chosen geometry/limiter. Report measured performance with hardware/rank/thread ownership, not inferred acceleration.
2. **Physical model:** no generator steering and Re=infinity for Kepler. W0 reduces to Voronoi and individually fixed masses. If weights/mass transport or intrinsic rotation are adopted, demonstrate their mass, linear/angular momentum, energy and reduction contracts. A reconstruction of a rotating velocity field is not itself intrinsic cell spin.
3. **Old failure range:** retain existing t6/t8/t10 evidence and restarts. An improved trajectory must actually cross the old failing range with local orbital, entropy, density and conservation diagnostics. A small global energy error or a surviving-looking outer ring is insufficient.
4. **Matched paper comparison:** align density profile, gravity/softening, gamma/P, domain/boundary conditions, effective spacing and units. Show actual density fields with common spatial/color scales at the same physical times. Report radial velocity, orbital/phase errors, mass-weighted entropy/heat and global/local angular-momentum transport as well as the paper norms.
5. **Long duration:** cover the published figure times and the stated approximately1000-orbit duration. Fix and publish the reference radius/period before the long run; include inner-edge and R=1 orbit counts to prevent changing the time definition. Do not close the goal at t10, t600 alone, or a short favorable interval. If units are scaled, derive time, velocity and pressure scaling explicitly.
6. **Admissibility and regressions:** audit all RK stages and the final state. No hidden negative-U stage, unreported EOS repair, floor injection, invalid geometry or nonfinite state. Preserve shock accuracy/positivity controls when changing reconstruction or geometry. Reuse valid restarts/checkpoints, with output schedules and ownership recorded.

## Current non-completion evidence

The retained original baseline reached t10.00058 but overheated the core. The vector limiter passes affine tests yet fails actual coarse Kepler/shock gates. Fixed-eps core refinement is a diagnosis, not a new success definition. Intrinsic-spin/weight closure and paper-duration stability are not yet established. The goal remains active.
