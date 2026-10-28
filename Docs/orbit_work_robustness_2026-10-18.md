# Orbital-work diagnostic: OpenMP and rejected-trial checks

From the preserved checkpoint224, compared two new runs against the existing
same-binary OMP1 instrumented .0025 step:

- OMP4, otherwise identical configuration.
- OMP1, dt cap .005 with forced rejection after final HG, followed by the
  accepted .0025 retry.

Both exited0. At the identical endpoint, all reader fields id,x,y,vx,vy,
mass,volume,den,pressure,ie,w2 match the direct OMP1 run bitwise. All four
accepted stage values wdt,hydro_minus_TCF,TCF,shear,total match exactly.

The rejected attempt produced four stage rows. The analyzer clears pending
stage and HG rows on `[RETRY] step rejected:` and commits only on
`[ACCEPTED_STATE]`. It observes exactly one rejection in the retry case and
none in OMP4, and requires four ordered stages/eight HG-bin rows per accepted
step. Thus the rejected stage rates are excluded rather than summed into
the accepted budget. HG increment numeric parity is not separately asserted
by this analyzer; physical endpoint and accepted RK work parity are.

Scope is one accepted MPI1 step, one forced final rejection and OMP1/4.
This does not validate repeated retries, MPI diagnostic support, CUDA hydro
forces, every snapshot field, or long-time physical stability. The continued
40-step work interval remains a separate run and its completion must be
verified independently.

Evidence: canonical_orbit_work_robustness_2026-10-18.json.
Scripts: run_orbit_work_robustness.py, analyze_orbit_work_robustness.py.
Remote root: /gpfs/kjhan/LagForce/Codex_orbit_work_robustness_20261018.
