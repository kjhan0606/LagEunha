# Bounded Q1 gravitational timestep diagnostic preparation

This work follows the focused Astra t6.5 assessment. The running bounded Q1 chain and its already submitted post-t10 controls are unchanged. The final objective remains native CUDA+MPI+OpenMP Kepler stability through the GIZMO long-time tests; none of the checks here qualifies t120/t600 or final runtime competitiveness.

## Preregistered core screens and restart selection

`CONTRACT.json` was written before examining the new .825 distributions or the candidate-cap outcomes. It uses original core IDs in r0<.1 and .1–.5, mean plus cell-p99/max density, normalized circular-orbit error and log-entropy criteria. These are numerical diagnostic criteria, not GIZMO-paper tolerances. Regional energy/angular trends are reported rather than assumed separately conserved.

The earliest main-chain checkpoint at exact t=.825 passes the old 28 gates but **fails the added core orbit screen**. For r0<.1 its mass-mean normalized position error is 1.9185%, cell p99 54.9310%, maximum 54.9508%. Normalization is max(original radius, Plummer softening), including the stationary central site. Its density/entropy screens and the .1–.5 screens pass. `RESTART_ELIGIBILITY.json` preserves every outcome; no threshold was changed after observing failure.

A future accurate native branch therefore needs the exact preserved original snapshot, SHA256 `0b7ccf1999a5ec6d22b8b77aeb5935abdda5e442c74aef8daee1fcfdc809f0ca`. It must not inherit the .825 phase error, choose a later favorable checkpoint, regenerate the IC, or repeat the existing coarse branch. The separate Ada capped comparison has since started from a previously accepted t=.0015284035789018555 native checkpoint, with planned checkpoints .4125 and .825; see `ADA_LAUNCH.json`. It is independent of the running t10 chain.

## Independent gravity-only screen

Existing classical RK4 gravity-only ODE code advanced all 3,246 original r<.5 sites from their saved initial state to .825/3/6.5/10. The reference uses their Plummer circular frequencies. It does not evolve gas or update the running native state. The screen reserves 90% of the native orbit-error budget for other errors.

| eta | Global cap | Core normalized position mean at t10 | Core cell p99 at t10 | All-target gravity-only screen |
| --- | ---: | ---: | ---: | --- |
| .2 | 2.54734e-4 | 2.02297% | 59.2667% | FAIL |
| .1 | 1.27367e-4 | 0.064335% | 1.88485% | FAIL |
| .06 | 7.64202e-5 | 0.0050195% | 0.147009% | PASS |

The .1 failure is against the stricter **gravity-only headroom screen**, not a measured native hydrodynamic failure. The largest passing preregistered candidate is .06. This is an eta candidate through t10, not a physically validated long-time choice. The entire ODE diagnostic took 43.5 seconds on Ada. The earlier 95%-attribution test remains a preserved failure; these results do not prove that every core density/entropy error is caused by gravity integration.

## Actual native setup/timing preflight

The original frozen v2 executable read the same preserved t0 IC using native restart argument 1, W0/fixed masses, Re=infinity and the same bounded polygon. Only existing `EUNHA_KEPLER_DTETA=.06` was enabled. No force/integrator/mesh/weights/spin implementation changed. One 52 MB endpoint was requested, with periodic dumps and full-particle TSV disabled.

The actual native CUDA+MPI+OpenMP run on Ada completed 20 full capped intervals plus one 4.34e-19 IEEE endpoint-time tail, at exact t=.0015284035789018555. Maximum accepted h was 7.642017894509277e-5. All 28 original physical/state checks, source/runner hashes, fixed mass, zero injections and exact accepted-clock checks passed. Total native execution was 18.867 seconds; regular RK4 median after two warmup steps was .7446515 seconds. The saved endpoint is **timing/setup evidence only**, never a substitute for long-time testing.

At that frozen early-geometry cost, 130,856 cap-limited steps to t10 would take about 27.1 hours on the original Ada binary. This is an extrapolation, not an observed long run, an upper bound or an H200/cache estimate. It demonstrates why a globally small step cannot simply be adopted as the production solution. The real late-state H200 cache qualification and incremental controls in job 439375 remain an independent data source. A diagnostic capped branch needs an explicit cost bound, early core checks, accepted-restart continuation and comparison against already stored coarse endpoints. Core accuracy must be established alongside acceleration before unqualified t120/t600 expansion.

## Completed same-time capped comparison

The native eta=.06 Ada branch completed exact t=.825 from its preserved t=.4125 predecessor. All 28 original checks and the additional core screen passed; final SHA256 `3a9dc7f01d720c945c3ef21acfe912015adc26bb155bbaf615ef31d4d81f7dc7`. The [same-time comparison](ETA006_VS_UNCAPPED_0825.json) retains identical t0 IDs and masses and independently reconstructed geometry. In the original r0<.1 core, the mass-mean position error normalized by max(r0,eps) is **1.9185%** without a cap versus **0.05487%** with eta=.06; cell p99 is **54.93%** versus **0.2513%**. Relative total-energy change is -3.54e-7 versus -3.17e-11. These are early-time results, not t6.5 or t600 validation.

The cap has a significant countervailing thermal result: mass-mean absolute disk log-entropy change rises from **8.49e-6** without the cap to **1.54e-4** with it, although the underlying smooth trajectory is similar. The two runs used the original single-precision thermal representation; this observation is consistent with the independently measured step-count/thermal-precision sensitivity, but is not proof of a unique cause. The separately staged thermal64 comparison is needed before selecting the cap as a long-time physical solution. The existing saved t=.825 snapshot and all prior Q1/Q2 results are preserved.

For the original r0=.5–.75 inner disk at t=.825, both cases have zero sites whose nearest-neighbor distance is below 80% of the same pair's exact-circular reference distance. The one-percentile ratio is about 0.999999 in both; the RMS difference of the nearest pair's specific angular momenta is 3.50e-9 uncapped and 3.53e-9 capped. Its mass-mean absolute log-entropy change is 3.36e-6 uncapped and 3.92e-5 capped. The [receipt](INNER_PAIRS_0825.json) and [reader](../audit_kepler_q1_0825_inner_pairs_20261029.py) choose each endpoint's nearest neighbor independently; their IDs could differ, and mutual pairs can be counted twice. This is an early-time observation, not a demonstration that Q1 will avoid the close-pair growth seen in the distinct Q2 late run.

## Frozen-state individual-timestep arithmetic

The saved capped legacy-thermal Q1 t=.4125 state contains 90,352 cells. Its
smallest saved final-stage hydro timestep is **6.8301e-4**, while the
eta=.06 softened orbital limit is **7.6420e-5**, about 8.9 times smaller.
One cell occupies the shortest power-of-two timestep bin; seven are within
the first two bins. The initial-radius r0<.1 region has 134 cells; in
r0=.5–.75 the median combined local constraint is **0.03077**. See
[receipt](LOCAL_DT_BUDGET_04125.json) and [reader](../audit_kepler_q1_local_dt_budget_20261029.py).

For a frozen .4125-time interval, counting independently active cell updates
under these stored constraints gives 487,720,096 global updates versus
821,516 idealized power-of-two block updates, a **594-fold arithmetic count
ratio**. This is not a wall-time speedup prediction or a performance upper
bound: the calculation omits same-time face impulses and work, evolving CFLs,
neighbor timestep coupling, mesh reconstruction, GPU/MPI scheduling,
gravity accuracy, impact events and rollback. It does not validate an
asynchronous integrator, and pressure work can itself require short steps.
It does quantify why the user-planned individual stepping is important after
the current (3) development, while the small eta global stepping remains a
diagnostic rather than an AREPO-speed production strategy.
