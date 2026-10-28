# Exact W0 Kepler tessellation work on Ada, t=6.5→6.501

The [hash-guarded work receipt](KEPLER_Q2_TESS_WORK_DIAG_V8_ADA_T6501_2026-10-27.json)
uses two diagnostic-only CUDA builds on the same RTX 5000 Ada and the same
accepted 5% site-jitter `t=6.5` checkpoint. Both one-step W0 audits passed
11/11, and each diagnostic endpoint's ten physical fields are bitwise equal
to its corresponding uninstrumented v3 (5×5) or v4 (3×3) control. The builds
add atomic diagnostic counters, so their absolute times should not be used
as production throughput forecasts.

| Per tessellation call, four-call mean | 5×5 seed | 3×3 seed | 3×3 / 5×5 |
|---|---:|---:|---:|
| Directional-preclip bin visits | 21.10 million | 7.676 million | 0.364 |
| Seed-pass bin visits | 21.10 million | 7.676 million | 0.364 |
| Directional-preclip half-plane clips | 896,945 | 885,257 | 0.987 |
| Seed-pass half-plane clips | 470,808 | 463,758 | 0.985 |
| Adaptive bin visits | 784,832 | 1,129,074 | 1.439 |
| Adaptive half-plane clips | 96,444 | 165,493 | **1.716** |
| Input polygon edges processed by all clips | 8.62 million | 9.08 million | 1.053 |
| Median Kernel A time (diagnostic builds) | 87.475 ms | 112.273 ms | **1.283** |

The independent [linked-bin occupancy audit](KEPLER_Q2_GPU_SEED_VISITS_T65_2026-10-27.json)
predicts an average of 400 candidate visits per disk site for 5×5 and 145
for 3×3, consistent with the observed totals after accounting for 682 local
box retries. The older 512-candidate buffer is not exceeded in this state;
its observed maximum seed occupancy is 422.

The 3×3 seed reduces raw bin visits by 63.6% but leaves the complete
polygon less settled before adaptive completion. Adaptive clips rise by
71.6%, and the exact same physics takes 28.3% longer in the instrumented
comparison. The uninstrumented matched comparison showed a similar 29.0%
slowdown. The counters establish more adaptive work; they do not isolate
warp divergence, register spilling or the cost of each clip, so a specific
hardware bottleneck remains unproved. Keep the 5×5 reference for ongoing
physics runs.

An exact next candidate is to retain 5×5 preclipping and skip an entire
seed bin **only when every generator in it is already too far to cut the
current convex polygon**. In the W=0 reduction, a generator at distance
`d` has bisector distance `d/2` from the owner. If a bin's minimum possible
distance `d_min` obeys `d_min² > 4R²`, where `R²` is the largest current
vertex distance squared, all of its candidates would fail the existing
radius rejection. This can reduce the second 5×5 scan while preserving
candidate order and the later adaptive completeness proof. It requires a
separate bitwise endpoint test before adoption. The weighted case needs its
own conservative power-bisector bound; do not apply the W=0 inequality to
arbitrary weights.

The isolated [v9 CUDA patch](patch_kepler_tess_seed_bin_prune_v9_20261027.py)
implemented that strict W=0 bin test with a roundoff margin and excluded
clamped grid-edge bins. Its matched Ada replay passed [W0 11/11](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_tess_binprune_v9_ada_t6501/w0_v9/audit_t6501.json)
and [weighted 16/16](/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/gpu_tess_binprune_v9_ada_t6501/weighted_v9/audit_t6501.json).
All ten physical endpoint arrays are bitwise equal to their respective v3
controls. W0 Kernel A median over four calls was **90.307 ms**, versus
**86.989 ms** for v3 (+3.8%); weighted was 83.818 versus 87.685 ms (-4.4%)
in this small sample. The W0 speed target did not improve, so v9 is not
adopted for long runs. The first `make -j2` attempt hit an unrelated Intel
Fortran module-order error in `bessels.f90`; `make -j1` completed the same
isolated source tree and produced binary SHA-256
`3f18131a31c235703965a094043747d9527fdaa3973d4f71a0c6e7618e465553`.
