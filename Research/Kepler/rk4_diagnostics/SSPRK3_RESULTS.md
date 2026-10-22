# Native SSPRK3: verification and limits

2026-10-20. Experimental opt-in only. No long Kepler stability claim.

## Implementation

`Exam/gfs_ssprk3.inc` is included by the user's native C/CUDA solver. `GFS_SSPRK3=1` selects SSPRK(3,3) with stage-wise FE positivity/CFL/orbital bounds and collective rejection. Rejected attempts restore the original full rank payload, ownership/counts, `dtold`, floor/pair counters and clear reusable GPU/halo state. Positions are blended through unwrapped displacements. Individual masses and W=0 are invariant. No steering, physical viscosity, heat top-up or pressure-ratio threshold change is introduced.

The candidate currently supports periodic native phase1/MUSCL/av5, fixed masses and W0. It is not a completed dynamic Laguerre or cell-spin implementation. Default RK4 remains selectable and unchanged in arithmetic.

Current-stage preparation is geometry -> EOS -> halo -> gradient/force. Source inspection found that the old GPU gradient could consume rho/P from the preceding stage before the EOS kernel refreshed them. Thus this comparison includes cache-preparation changes: it does **not** isolate a time-integrator change with identical RHS inputs. The force, limiter, Riemann and pressure/work formulas themselves remain unchanged.

Direct SSP integration of U and velocity does not exactly conserve nonlinear total energy. `ssprk3_steps.tsv` separates semidiscrete hydro work, gravity work and finite-step energy residuals. All attempted RHS evaluations count toward the process budget, including rejected attempts.

## Completed default and rollback checks

`ssprk3_v1/manifest.json` preserves the results. A default RK4 fine preflight matched the frozen native v4 reference in all recorded fields and metadata bitwise.

A uniform 64^2 periodic advection test crossed both periodic and MPI ownership boundaries. The forced discarded FE step was restored before retrying. MPI1 reference versus forced rejection, MPI4 reference versus forced rejection, and MPI1 versus MPI4 all matched bitwise in recorded state and metadata. Advection error <1e-12; energy/mass errors zero; no floor. Reference used six RHS, forced rollback seven.

The first fine interval completed three accepted steps/nine RHS but lacked a final checkpoint because the old output cadence missed the floating endpoint. The output helper then failed. That run is preserved. v2 adds an explicit final output when `GFS_EXACT_STOP=1` reaches the endpoint; only missing interval work was resumed, and the completed default/rollback suite was not repeated.

Frozen v2 binary SHA256: `1d999a23f1dfae0b941133e0db1b6e3a1e1ab292df7a6c2b191643d7f99612c2`.

## Captured intervals

All start from preserved normal states; none reset U or inject heat. Each reaches the same endpoint as its original failed interval.

| Case | Accepted / RHS | min U | Relative total-energy change | Floor |
|---|---:|---:|---:|---:|
| Fine base, t=.10196265 -> .10258819 | 3 / 9 | 7.10093e-12 | 1.6709e-13 | 0 |
| Fine effective half-spacing control | 5 / 15 | 7.20095e-12 | 9.0919e-14 | 0 |
| Sedov base, t=.02458081 -> .02463711 | 2 / 6 | 3.66211e-7 | 1.11083e-6 | 0 |
| Sedov half ceiling | 2 / 6 | see manifest | 8.10058e-7 | 0 |

An initial nominal half-ceiling did **not** bind for fine: positivity already imposed smaller steps, so it reproduced the base bitwise and was not a time-resolution comparison. One additional effective control used ceiling .00014093755556411012 (half the largest actual base step), giving five instead of three accepted steps. No timestep scan was performed.

The effective fine comparison has total-U L1 relative difference 1.00925e-8 but maximum local symmetric relative difference .0138894. ID68526 has U/Uad=.0440095 versus .0446283, a difference .000618789 of its adiabatic reference. Both are still depleted. Whole-domain small differences cannot certify the four affected cells.

For Sedov, sum|delta U|/max(sum|U|) is 1.08136e-4; maximum local symmetric relative U difference .55949. The two ceilings demonstrate admissibility and measurable temporal sensitivity, not a resolved shock convergence study. The finite-step energy residual dominates the semidiscrete net hydro work; positivity alone is not an accuracy pass.

Results: `ssprk3_intervals_v2/manifest.json`, `ssprk3_fine_active_half_v2/manifest.json`, `ssprk3_interval_analysis_v2.json`. Analysis ran with one CPU thread on user-authorized lageunha. Evolution is exclusively native C/CUDA/MPI.

## Fresh fine run: budget stop, severe thermal depletion remains

The same 71,080-cell integrated IC was evolved along one fresh path toward t=.2 on lageunha/Ada, with maximum 3072 RHS/768 accepted steps and eight retries per step. It stopped at **t=.13509103681003889**, native exit99, after515 accepted steps and3072 RHS,584 rollbacks and873.38035 seconds. The final accepted ledger has3070 RHS; the final two discarded evaluations consumed the remaining budget. It did **not** reach t=.2. No further simulation or restart has been launched.

The native solver preserved a normal step-start (~39 MiB): `ssprk3_fresh_v2/fine_fresh/khdiag.step000516.ssp_rejected_start.dat`. The generic native message says “failed stage is not a restart state”, but this particular file is the restored full normal backup, not the rejected trial. Its time equals the last accepted t_start+h and its mass/W0/U checks pass. The launching helper's `missing_checkpoint=true` means no final `khout` was found; it overlooked this differently named preserved normal state. The raw manifest is retained rather than rewritten as a success.

`ssprk3_final_analysis_v2.json` is the final independent analysis, run with one CPU thread on lageunha:

- All individual masses are bitwise fixed; all W=0; no floor injection; minU=3.43356e-12.
- The four focus U/Uad ratios are .0210655/.0194730/.0210763/.0210554, at V/V0~1.26272. The fresh run therefore retains severe thermal depletion even though every committed stage is positive.
- Initial R<.1: mass-weighted |A/A0-1|=.591962, min A/A0=.019473, max4.74045; normalized circular-trajectory position/velocity errors .004067/.003365. Hot and depleted cells coexist. The adiabatic comparison is a diagnostic reference, not a forced general Riemann solution.
- Initial .2<R<.5: entropy error .00649194, normalized velocity error .000180683.
- Independent Voronoi audit certifies69664 interior cells against missing periodic images. R<.5 stored-volume relative error <=5.753e-8, consistent with float storage; stored rho versus m/V <=5.946e-8. Geometry remains accurate at the depleted state.
- Accepted-step semidiscrete net hydro work is1.90856e-16 versus absolute work integral1.26544e-5. Finite-step hydro and potential residual sums are7.25558e-10 and -2.74436e-10. Small global energy residuals do not certify local thermodynamics.

The attempt cost also exposes an implementation inefficiency: selecting a first-stage bound almost exactly can cause a later, marginally smaller orbital bound to reject a trial, while later pressure-switch/positivity restrictions require additional retries. Dense-grid geometry search is much more costly than the force kernel. No speed-up is claimed for this candidate.

**Evaluation:** admissible stages and a full rollback are useful temporal safeguards. This fresh candidate fails the planned practical/thermal gate and does not establish a long-time physical solution. It is not extended toward t6-10 or paper duration on the strength of positivity. Existing outputs remain intact. The complete evidence was sent to Astra in one bundled review before selecting a new spatial implementation.

## User's inner-refinement proposal

The prior fixed-eps comparison already supports the proposal at the initial-RHS level: inner spacing .09375 -> .01171875, total cells 65,536 -> 71,080, and initial R<.5 sum of positive Udot 1.8144767e-4 -> 1.8646217e-7 (about 973x reduction). This includes reversible pressure work and is neither a pure entropy-production rate nor a measured spatial convergence order. Fine old RK4 subsequently encountered negative intermediate U at t~.102.

For the actual Plummer IC, Omega=(r^2+eps^2)^(-3/4), q=-dlnOmega/dlnr=(3/2)r^2/(r^2+eps^2), eps=.047. The strict centre tends to rigid rotation; the coarse h=.09375 did not resolve it. At r>>eps, keeping h/r roughly constant distributes relative geometric resolution more evenly. At r~eps and inward, a finite fine h<<eps and a smooth transition avoid an unbounded central point count. Initial masses can differ by cell, m_i=integral_cell rho_0 dA, while each remains fixed thereafter.

For this potential the actual strain magnitude is s=q*Omega=(3/2)r^2(r^2+eps^2)^(-7/4). Differentiating gives its maximum at r=2eps/sqrt(3), about .05427 for eps=.047; it vanishes as r->0. Thus "inward means ever stronger shear" is incorrect in the softened centre. An initial sampling candidate h(r)=eta*sqrt(r^2+eps^2), smoothly capped by the desired outer spacing, would resolve the softening transition, remain finite at the origin and approach constant h/r outside. This is a proposed initial-resolution profile, not an implemented evolution law or a stability proof. It does not alter weights or steer generators during evolution.

Refinement reduces unresolved curvature and neighbour velocity differences, but q*Omega remains the physical strain rate: cell distortion over a fixed number of local orbits is not removed by shrinking h. The older TCF/ring sliding instability (~cs/h) belongs to a different branch but cautions against treating refinement as a stability theorem for an incompatible operator.

Relevant independent evidence: Zier & Springel, [Simulating cold shear flows on a moving mesh](https://arxiv.org/html/2205.07916), section4, distinguish mesh noise from timestep errors and improve face integration. Their ALE flux mechanism should not be copied as a proof for this fixed-mass pressure-only operator; it motivates checking spatial integration/transport consistency alongside resolution.

## Subsequent completed spatial contract bundle

`ALE_SHADOW_RESULTS.md` records the next native conservative ALE shadow, its
MPI/boost/GCL and observer checks, and the signed work/moment evidence from
the two preserved cold states. That diagnostic does not update mass or evolve
spin and does not make the SSP path a validated long-time model. No further
fresh Kepler evolution was performed on the basis of a diagnostic pass.
