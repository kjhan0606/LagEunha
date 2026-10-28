# Original coupled-response continuation completed

Reference,+,- all exit0 after760 finite accepted steps without retries from
checkpoint224(t2.1029059), ending checkpoint984(t4.002905917763361).
Relative to the original t2.0029059 seed, the duration is2.0. No initial
evolution was rerun. All final snapshots and intermediate results remain.

| Squared response diagnostic | Hydro growth | Independent gravity growth |
|---|---:|---:|
| Raw velocity |250.90734|17.49079|
| Circular-departure velocity |1619.44251|.9043804|
| Omega-weighted displacement |303.19075|18.83011|

Gravity controls use exactly the saved time, with RK4 caps.0005/.00025;
maximum change in growth ratios is2.98e-10 at this endpoint. Gravity-only
orbital phase evolution does not reproduce the hydro departure growth.

The even/odd velocity ratio is.042453 (4.2453%), larger than at earlier
times. Do not interpret this finite-perturbation result as a converged
linear eigenmode or Lyapunov exponent. Half-amplitude validation currently
ends at t3.4; the fine-dt continuation to t3.4 is still running. There is
no hydro refinement claim at this final t4 endpoint.

Reference restart-relative energy drift is-5.44584e-8 and angular-momentum
drift-2.17440e-7. Small conserved-quantity drift does not establish physical
accuracy: the separate unperturbed physical evaluation has inner phase
RMS1.51835deg and relative-pressure RMS.338565 at this time.

Next physical experiment already runs the legacy reference onward from
this exact final checkpoint to t10, alongside the smooth-reconstruction
method starting from the preserved original t2 reference. No additional
copies of the original t2-to-t4 trajectory are needed. Retain both physical
runs, the fine-dt response comparison, and original long-time calculations.
No claim that the new method resolves the physical instability is made.

Evidence: canonical_coupled_t4_response_2026-10-18.json,
canonical_coupled_t4_orbit_control_2026-10-18.json,
canonical_coupled_t4_completion_2026-10-18.json,
canonical_legacy_t4_physical_2026-10-18.json.
