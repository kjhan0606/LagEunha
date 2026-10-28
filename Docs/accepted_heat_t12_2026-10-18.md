# Accepted HG heat through t12

Streamed the existing frozen-consistent HG run log through the last accepted
endpoint not exceeding t12 (11.9997825487). This includes 4113 accepted steps;
222 rejected trials are excluded and the following eight HG records are not
committed. Regional HG sums agree with the accepted hg_dU ledger within
3.83e-17. This is a history count: latest accepted rows reporting retries=0
do not imply no rejection occurred earlier in the trajectory.

For endpoints in (10,12], global delta U=2.3534468560e-4. HG heat is
2.0410413807e-5 in r<.5 and 7.3343417593e-6 in the disk bin, totaling
2.7744811885e-5. HG accounts for approximately 11.8% of net internal-energy
increase in this bookkeeping comparison. This is NOT its fraction of total
dissipation: non-HG delta U=2.0759987372e-4 includes reversible pressure work.
Bins follow instantaneous radius, unlike the initial-ID masks of the t12
pressure/entropy analysis. Do not mix those as equal material budgets.

The from-t0 capped run has no HG_TRANSFER records before the first accepted
step. The regional analyzer rejected that incomplete input; no regional HG
values were inferred for it. No existing simulation was restarted or altered.

Earlier canonical_entropy_step_audit_2026-10-17.json already established
near closure of measured entropy increments by combined TCF+HG heating in
sampled accepted steps (t2/t4/t6), with passive diagnostic on/off equality.
That is combined heating, not a split cumulative t0-to-t12 entropy budget.
Current evidence supports modeled dissipation as the main entropy source at
those sampled states, but does not identify all accumulated t12 contributions.

Next missing measurement: distinct per-cell accepted TCF and HG entropy
integrals along the same trajectory, with rejection rollback and regional
membership handled consistently. Do not equate total-energy conservation,
internal-energy growth, entropy production and orbital stability.

Results: canonical_accepted_heat_t12_2026-10-18.json.
Script: analyze_accepted_hg_heat_t12.py.
