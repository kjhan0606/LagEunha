# Native bounded q2 weighted Kepler comparison at exact t=10

The user's CUDA+MPI+OpenMP solver evolved the accepted weighted `t=6.501`
aligned and 5%-jitter checkpoints to exact `t=10` with one A10 and one MPI
rank per case. The earlier `t=0→6.5` history was the unchanged W=0,
fixed-mass branch; this experiment isolates the effect of enabling evolving
Laguerre weights over `6.501→10`. Physical viscosity and site steering are
off. Cell masses remain exactly fixed in this experimental operator.

Slurm **437384** aligned completed `0:0` in 47m59s. Slurm **437385** jitter
ended `FAILED 1:0` in 47m56s because its endpoint density gate failed; the
complete exact-time snapshot and native audit were retained. The original
job receipts are [aligned](aligned_original_audit.json) and
[jitter](jitter_original_audit.json). The independently run
[corrected equilibrium auditor](../audit_kepler_q2_weighted_equilibrium_20261027.py)
verifies the initial, restart, executable and endpoint hashes, cuts the
stage ledger at each endpoint's own step marker, checks independent power
areas and edge weights, and does not use a W0 Voronoi area as a weighted
physical reference. Its receipts are [aligned](aligned_corrected_audit.json)
and [jitter](jitter_corrected_audit.json).

| Mass-weighted disk measure at t=10 | Original W0 aligned | Weighted aligned | Original W0 jitter | Weighted jitter |
|---|---:|---:|---:|---:|
| Mean absolute `rho/rho(0)-1` | `1.96891e-4` | `1.45954e-4` | `2.37446e-3` | `1.71685e-3` |
| Mean absolute entropy-log change | `1.19396e-5` | `9.06615e-6` | `2.03302e-5` | `1.53116e-5` |
| Mean absolute exact-orbit position error | `3.57456e-7` | `3.35180e-7` | `1.62485e-5` | `1.59889e-5` |

Weighted density drift is **25.87% lower** aligned and **27.70% lower**
jitter than matched W0 endpoints. The aligned weighted case passes all
corrected finite-time gates. The jitter case fails only the strict
`mean |rho/rho(0)-1| < 1e-3` equilibrium-density gate. Its independent
power-area maximum relative error is `5.96e-8`, edge weight difference
over squared site distance is at most `0.0111` against a `0.8` guard, and
total energy drift is `8.08e-9`; these checks locate the remaining density
error outside the final power-area reconstruction. The stage-global
projected weight work is at most `1.11e-26` in the jitter run, but global
cancellation does not establish a local moving-face first law.

The [matched density field](density_jitter_w0_vs_weighted_t10.png) uses
identical pixel positions and color scales for the W0 Voronoi and weighted
power diagrams. Its lower panels show `rho/rho(0)-1` over the initial disk.
The visible high-frequency pattern is reduced in the weighted run but is
still large enough to fail the density gate.

The W0 exact-orbit geometry difference in the corrected weighted receipt is
diagnostic only; [the gate correction](../KEPLER_Q2_WEIGHTED_GATE_CORRECTION_2026-10-27.md)
explains why it cannot define physical success for a power diagram. Neither
this finite-time result nor the one-rank global pressure-work projection
proves pairwise mass, momentum, energy or spin conservation. The IC is the
bounded 55,351-site diagnostic, distinct from the full 384² polar q1/q2
controls and from the paper's exact particle layout. Stability through the
historical `t>6` window improves but remains unproven for jitter and for
the GIZMO `t=600` image or longest-orbit endpoint.

The separately submitted jitter job **437387**, using guard build SHA
`e90a287d25b0d82705f70ad645b5eadade7a18cd3aebff1622606960a333e7a1`,
also reached exact `t=10`. Its [staged audit](jitter_437387_staged_audit.json)
failed the W0 ballistic-area comparison identified above; the
[corrected audit](jitter_437387_corrected_audit.json) instead fails only
the equilibrium-density gate. The 437385 and 437387 endpoints have
bitwise-identical `x,y,vx,vy,mass,volume,ie,den,pressure,w2` arrays even
though their complete file hashes differ. These two submissions duplicate
the same physical trajectory and do not constitute independent evidence of
improved accuracy.

The analogous aligned job **437386** completed `0:0`; its
[corrected audit](aligned_437386_corrected_audit.json) passes, and its same
ten physical endpoint arrays are bitwise identical to **437384**. This is
also a duplicate trajectory, not a new independent physical case.

## Inner-disk pair torque follow-up

The saved aligned W0 and weighted t10 endpoints were re-read using the **same** 3,931 original r0=.5–.75 sites and the W0 endpoint's nearest-neighbor IDs in both cases. Fixed masses and independent power-cell areas were checked. The RMS specific angular-momentum difference is **4.28e-8** for W0 and **3.52e-8** for the weighted case; the mass-mean absolute density change in this region is 3.91e-4 and 2.86e-4. However, the RMS difference of independently reconstructed specific pressure torques is **7.52e-9** for W0 and **1.97e-8** for the weighted case. The instantaneous projection `mean(DeltaL*DeltaTau)/mean(DeltaL^2)` is +0.140 and +0.406 per time unit, respectively. See [receipt](../KEPLER_Q2_WEIGHTED_T10_PAIR_TORQUE_2026-10-29.json) and [reader](../audit_kepler_q2_weight_pair_torque_20261029.py).

This is one time sample with initially small L differences; ratios can be sensitive to the small denominator, and pairs are fixed by the W0 endpoint. It does **not** demonstrate that the weighted trajectory subsequently grows faster or that weights repair or worsen the t47.5 Q2 instability. Weighted evolution was enabled only after t6.501. The weighted reconstruction has not been independently compared to native weighted RHS, and global projected work is not a local first-law proof. The result motivates a matched, time-resolved torque/weight-energy audit before selecting this weight law as a remedy.
