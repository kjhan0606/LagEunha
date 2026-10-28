# Weighted Kepler local pressure work at t≈4.951

The isolated [stage manifest](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_capacity64_v5_fresh_localwork_t4952_ada_v1/stage_manifest.json)
copied the accepted fresh-weighted native `t=4.951043599280444` checkpoint
and the unchanged v5 CUDA executable. The parent endpoint SHA-256 is
`6453ccebd89ea90af807f7ce45f68bde36a7969320e456e0b01adb7290eda63b`.
On an otherwise idle lageunha Ada GPU, one exact RK4 step to
`t=4.952043599280444` took 3.710 seconds. The [native endpoint audit](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_energy_neutral_capacity64_v5_fresh_localwork_t4952_ada_v1/audit_t4952.json)
passes **21/21** checks, including four weighted stages, fixed mass,
positive state, independent power areas (`5.96e-8` maximum relative
error), no injection, and total-energy relative change `4.00e-9` from
the initial state. The endpoint SHA-256 is
`fa45f34188d0b9b1e9c9f67507cdc051014453802a870a9bcae09232237eb716`.
The live long job was not modified or restarted.

The [midpoint face reconstruction](KEPLER_Q2_WEIGHTED_FRESH_LOCALWORK_T4952_2026-10-27.json)
uses the same algorithm and `dt=0.001` as the earlier
[t≈0.825 diagnostic](KEPLER_Q2_WEIGHTED_FRESH_LOCALWORK_T0826_2026-10-27.json).
It reconstructs 148,557 disk-touching power faces; cell and face sums of
pressure-weight power agree within `6.02e-26` absolute power. This is a
secant/midpoint estimate, not a four-stage local RK4 ledger.

| One-step midpoint quantity | t=.825→.826 | t=4.951→4.952 |
|---|---:|---:|
| Cell-local `|P Vdot_weight| dt / U`, p99 | `3.706e-5` | `3.751e-5` |
| Cell-local ratio, maximum | `1.041e-4` | `9.806e-5` |
| Weight-only relative volume change, RMS | `1.403e-5` | `1.426e-5` |
| Actual total relative volume change, RMS | `8.626e-7` | `1.379e-7` |
| Absolute face pressure-difference work / domain internal energy | `7.735e-11` | `5.882e-10` |
| Native maximum global projected stage work | `3.659e-29` | `4.047e-28` |

The weights now cancel geometric volume change by a factor of about 103
in RMS, versus 16 at the earlier state. Nevertheless the p99 **local**
pressure-weight work per cell internal energy remains near `3.7e-5` per
step. The absolute face pressure-difference work is 7.60 times the earlier
estimate as pressure contrasts develop; the signed global work still
cancels to about `6.64e-15` of domain internal energy in this step.
Strong cell-volume stabilization and a small global work sum therefore do
not establish local material transport or a conservative physical
interpretation for moving weights. Conversely, nonzero local work is not
itself a demonstrated violation of the first law: the native update uses
the actual changing cell volume. The independent
[endpoint total-PdV closure audit](KEPLER_WEIGHTED_ONE_STEP_PDV_CLOSURE_2026-10-27.json)
measures `|ΔU + (P_start+P_end)ΔV/2|/U_start` in every active cell. Its p99
is `1.35e-7` in the early step and `8.36e-8` in the late step; the
corresponding one-step entropy-log p99 values agree at this precision.
Thus the observed local weight work is largely counterbalanced by the
generator-motion contribution in the **total** cell-volume and internal-
energy update. This endpoint check supports local adiabatic consistency
for these two smooth steps. It does not verify the four RK4 stages' pairwise
energy and momentum transfer or establish that the moving power faces are
material surfaces. Those remain required before treating this fixed-mass
weighted branch as a general Kepler solution, regardless of its improved
early density and orbit metrics.
