# AV compression blend: local encounter interpretation

## Final update

Both branches completed 200 accepted steps to t=6.202110493120931.
Final snapshot/log analysis passed; zero logged floor, dual-energy and entropy
floor injections, fixed mass, zero weights, no retries, all RK4 stage masks 15.
The long RK4E line can be split by runtime output, so injection fields are
checked across the full log with exactly one occurrence per accepted step.
See `canonical_av_branches_final_2026-10-18.json` and corresponding PNG.

Disabling chi decreases inner mean entropy increment from .05165765 to
.03696482, but increases final circular-velocity error RMS from .02058475 to
.02123840 (about 3.18%). Sampled minimum pair gaps are .10688066 and .10725621
initial h; both rebound. Final gap is .13597713 with chi1 and .12647948 with
chi0. This local result does not support adopting chi0 as a stability cure.
Original long runs remain unchanged. The interim observations below are
retained as history; the final analysis is authoritative for this branch.

## Scope

Original RK4/CUDA executable, MPI 1 and OpenMP 1, fixed masses and zero weights.
Both cases restart from the same saved state at t=6.002110493120864.
Only GFS_QBP_CHI differs (1 versus 0). Zero disables the compression blend,
not all artificial viscosity. The target is t=6.202110493120864.
This encounter diagnostic does not replace t=6–10 or long-time validation.

## Evidence through t=6.170110493120920

Authoritative snapshot hashes and 18 common saved frames are recorded in
`canonical_av_branches_partial_003_2026-10-18.json`.

| Quantity | chi1 | chi0 |
|---|---:|---:|
| 803–804 gap / initial material length | 0.1217186312 | 0.1160785705 |
| Separation rate | 0.0118055873 | 0.0082498841 |
| Pressure of 804 / initial pressure | 5.635751718 | 5.033784637 |
| Initial-inner pressure error, mass RMS | 1.428741471 | 1.409204250 |
| Initial-inner thermal energy | 5.83377718e-5 | 5.75617237e-5 |

Both cases have reversed their initial approach and are separating. Removing
the compression blend reduces heating but also leaves this pair closer and
reduces its separation rate. This does not establish a persistent bound pair,
nor improvement of orbital accuracy. Neighbor-face AV pushing toward the
partner is real, but is not sufficient evidence that total AV causes pairing.

At t0, geometric volume compression already produces pressure-rate RMS
0.201314 in the initial inner region, versus 0.00260072 from TCF heating.
These are instantaneous rates, not accumulated errors or additive RMS terms.
AV feedback must be studied together with this upstream volume/pressure error.

## Pending final checks

The final analyzer now includes velocity error relative to local circular
motion, perturbation kinetic energy, radial displacement from original t0,
and logged thermal-floor injection checks, in addition to pair spacing,
entropy, fixed masses, finite states, and RK4 stage masks. The active executor
was not replaced. Analysis-only copy:
`/gpfs/kjhan/LagForce/run_canonical_av_branches_analysis_v2_20261018.py`.

After both exit codes are zero, execute its `analyze` mode on
`/gpfs/kjhan/LagForce/Codex_av_chi_branches_20261018` and inspect `analysis.json`.
No new AV law has been adopted. Short comparison success cannot establish
long-time physical stability or fourth-order temporal convergence.

## Preserved long calculations

Actual process checks at 2026-10-17 16:51:17 UTC found baseline PID 2286374
at t=15.17548, constrained HG PID 2295268 at t=15.23417, and from-t0 timestep-cap
PID 2466162 at t=10.29407, all live. Conservation alone is not an accuracy gate.
