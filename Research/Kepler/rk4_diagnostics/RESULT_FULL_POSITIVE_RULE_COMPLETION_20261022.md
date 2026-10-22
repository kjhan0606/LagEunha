# Full positive-rule preprocessing completed; paired step result linked

Latest outcome: the actual single Ada paired step completed and was audited.
See `RESULT_ACTUAL_POSITIVE_TRANSITION_COST_GATE_20261022.md` for accepted
step5, accuracy pass and measured production cost failure. The narrative
below preserves the original preprocessing/fallback observations in time
order; its pending statements are historical, not current execution status.

2026-10-22, 13:03 KST. The t>6 and longest-GIZMO stability goal remains active
and unachieved. This is the completion of the existing bounded preprocessing
bundle, not a new physical trajectory or a stability claim.

## Completed original run

The original producer session 9211, launcher 2166274 and MPI ranks
2166305–2166308 were kept throughout the calculation. At 04:01:50 UTC the
launcher and ranks were gone, and both the native and wrapper terminal files
reported exit 0. Their final records show:

| Quantity | Completed result |
| --- | ---: |
| Original source triangles processed | 1,138,432 / 1,138,432 |
| Compressed triangles | 1,138,432 |
| Original positive quadrature points | 2,040,256,512 |
| Retained points, including any full-rule fallbacks | 249,691,129 |
| Point-count reduction | 8.1711213377 |
| Native preprocessing wall time | 18,966.2571 s, about 5.27 h |
| Four rule files combined | 17,913,969,996 bytes, about 16.68 GiB |
| Signed anchor delta U / U | +2.2939011723e-18 |
| Signed anchor delta Phi / U | -3.8829773771e-14 |

The candidate satisfies the predefined minimum fourfold point reduction and
the anchor integral discrepancy limits of 1e-8 relative to physical U. The
8.17 factor describes point counts; actual RK4 throughput has not yet been
measured with the completed rule.

## Audit scope

`material_full_positive_rule_completion_audit_v1/contract.json` records the
completed producer result, both exit codes, final source/input verification
records, all four 104-byte rule headers and sizes. All of these checks passed.
The native binary was also checked directly against its saved SHA256 manifest
and returned OK. The producer's final records contain 68 source-file and 59
input-file matches.

Each header has the correct rank, four-rank partition, 284608 local triangles,
1138432 global triangles, binary80 ABI, and original accepted anchor time.
The generating-manifest fingerprint agrees across ranks. This audit does not
replace the upcoming native reader's full per-record checksums, positive
weights, moment residuals and state-dependent contractions. No second full
oracle force evaluation or replay was launched for this audit.

The candidate is ready for the actual RK comparison. It is still uncommitted:
all x/v coefficients, fixed original M, generating entropy and the accepted
physical time remain at step 4, t=8.2355880701376864e-5.

## Exactly one downstream execution

The retained dispatcher 2351718 submitted H200 job **433999**. A direct
read-only scheduler query returned:

    433999|PENDING|(QOSGrpGRES)|0:00|gres/gpu:H200:1|16

The existing policy remains H200 first, one card, MPI4/OMP4. If it stays
pending for 600 seconds, the dispatcher requests pending-only cancellation.
Ada fallback is permitted only after terminal cancellation of this same job
and absence of the shared physical-transition claim. There is no parallel
second trajectory, automatic retry after physical execution, or new long-run
release. At the completion audit no physical claim existed.

The next evidence is the actual same-state low/high comparison at all four
RK4 stages, endpoint/cold conservation checks, restart roundtrip and measured
selected-rule cost. Apply the existing accuracy/cost decision and Astra's
conditional map/lifetime guidance to those results. Do not rerun accepted
intervals to repair a wrapper or to obtain a cosmetically clean exit.

The historical t6/8/10 checkpoints, global map injectivity, long material
representation lifetime, W0 dynamic Voronoi reduction and matched-paper long
evolution remain unproved. The future compression-diagnostics build is still
unexecuted and has not been inserted into this frozen paired bundle.

## Same-job cancellation and single Ada fallback observed

At 04:12:54 UTC a direct read-only sacct query of 433999 returned
`433999|CANCELLED by 10396|0:0`. The dispatcher independently retained the
same terminal accounting before fallback. At 04:13:50 UTC it acquired the
shared claim for Ada and launched runner 3026814, start identity 1085719927.
At 04:14:19 UTC the actual runner and its sha256 input-verification child
were live. No H200 run directory existed and no GPU process had started yet.
This is preflight, not an accepted RK step or proof of native CUDA execution.
See the completion audit's separate Ada_fallback_start.json. Preserve this
single runner; do not restart preprocessing or submit a second physical run.

At 04:16:31 UTC native execution was confirmed directly: runner 3026814 had
the retained timeout/mpirun/hydra chain and native ranks 3028401–3028404,
each with process start identity 1085724680. nvidia-smi independently listed
all four native executable paths using approximately 5490–5494 MiB each on
Ada. The wrapper status was RUNNING_ONE_ACTUAL_PAIRED_RK4_TRANSITION. Native
and wrapper terminal results, accepted step5 and selected-rule cost were not
yet present. These observations prove actual native MPI/CUDA execution, not
an accepted next step or t>6 stability.
