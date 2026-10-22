# Long Kepler objective: paper reference and unfinished evidence

Goal remains the user's original CUDA+MPI+OpenMP evolution, including the
historical t6--10 failure and the longest GIZMO comparison. A successful local
transfer or an accepted first RK step does not complete that goal.

## Primary paper checked on 2026-10-22

[Hopkins, arXiv:1409.7395v2, section4.2.4 and Figures6--8](https://arxiv.org/pdf/1409.7395).
The author-hosted PDF timed out in the web reader; the arXiv primary PDF was
read instead. No third-party search summary is used as quantitative evidence.

- Figure6 compares surface density at paper code time120; Figure7 shows600.
- The surrounding text reports integrations extending to approximately1000
  orbits. Thus600 is a required figure milestone, not the longest stated run.
- Figure8 includes velocity and density error measures and their evolution
  with resolution. Density preservation and orbital accuracy both matter.
- The stated setup uses gamma5/3, pressure1e-6, no self-gravity and a softened
  central potential, with an effective256² resolution. The full published
  radial profile must be matched before claiming a direct paper comparison.

The paper's main text associates600 with approximately160inner orbits,
whereas the Figure7caption gives more than250 at radius about0.5. Record
this discrepancy rather than silently using one count as a unique time unit.
An exact endpoint for the stated1000orbits also requires an explicit reference
radius. Keep the1000orbit objective pending while that normalization is
resolved; do not replace it with a convenient shorter time.

## Driver's comparison contract

For a softened potential -GM/sqrt(r²+eps²), derive and report

    Omega(r) = sqrt(GM)/(r²+eps²)^(3/4),
    T(r) = 2*pi*(r²+eps²)^(3/4)/sqrt(GM), N(r)=t/T(r).

This formula is a direct mathematical derivation, not an extra paper result.
Record GM,eps,reference radius and code time together. Existing project runs
may use other length/time normalizations; equal numerical times alone do
not establish a paper-equivalent duration. Preserve all old runs/configs.

First reproduce and pass the original t6/8/10 checkpoints without heat reset,
steering or added physical viscosity (Kepler Re=infinity). Then extend the
same physical model and restart history to the paper figure milestones and
the longest stated orbit duration. If size/IC differs, label that comparison
and retain a matched published-setup run for the final claim.

Required outputs: actual surface-density fields on consistent scales,
radial/annular density and orbit/phase errors, physical U/entropy and their
spatial heating rates, global and annular angular momentum, total energy,
positive-state minima and runtime/cost. Preserve uninterrupted/restart
agreement and the native RK-stage/heat/stress-work ledger. The current
binary80 CPU transfer reference must be connected and compared to the
original CUDA solver before any GPU production stability claim.

2026-10-22 update: RESULT_NATIVE_EVOLUTION_20261022.md proves actual accepted
RK4 and restart continuation on the retained24-cell finite CONTROL patch.
No full-disk t6/8/10 evolution or long solution is proved.

## Saved-input units and metric clarification

See KEPLER_PAPER_TIME_UNITS_AUDIT_20261022.md and its machine-readable JSON.
The current retained disk uses rin=2/rout=8 and a tapered 1/r density profile;
it is not a matched published-setup run. The supplement records actual orbital
periods, the paper's error annulus/thresholds, and why its departure-from-unity
density norm cannot be transferred to the current IC unchanged. It selects no
shorter final endpoint. Full native evolution has since reached accepted step4
at t=8.2355880701376864e-5; this also remains far short of t6/8/10.
