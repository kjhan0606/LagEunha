# Q1 short arithmetic comparison: missing-CUDA build caught and corrected

The fixed-K H200 probe **437232** had not started. A device-code audit of
its staged executable SHA-256
`6c2ec1ab5b0f1877108010ddb3c838c420712845c4e70c88135da5a7c0a8e26b`
found no CUDA ELF image, no `gpu_tessellate_and_build_faces` or
`getAccVoro2DBlend_GPU` symbols, and no linked CUDA runtime. It was 7.4 MB,
versus 15 MB for the accepted q1 CUDA executable. The isolated build had
recompiled `Exam/exam.c` without `USE_CUDA`; the existing archive still
contained GPU objects, but the CPU host object did not reference them, so
the linker omitted the native CUDA path. The batch-script syntax and SHA
checks could not detect this build-mode error. **Job 437232 was held and
then canceled before execution**; its staged files remain for provenance.

The fixed-K source content remains SHA-256
`354938194bae832d3c34903ef2a7bb0ab8e073e114a5d99df8d4d14dec099fdd`.
In the same isolated source at `/tmp/Kepler_entropy_k_20261027`, only the
host object was recompiled with explicit `USE_CUDA=1`, then `libexam.a`
was refreshed and the executable relinked with `USE_CUDA=1`. The GPU
`exam_gpu.o` and `exam_gpu_extract.o` hashes match the accepted q1 build
exactly. The corrected executable SHA-256 is
`f836416a58b4743c7b81d348c54a365f03c5bafc0e76c1cc492d8d95182f40c0`;
it contains `sm_86` and `sm_90` cubins, CUDA registration, GPU tessellation,
and variational W0 host calls. `ldd` reports no missing library.

Three new, isolated A10:1, 32 GB, one MPI rank/two OpenMP thread bundles
start from the **same byte-verified q1 `t=.825` checkpoint** and stop at
exact `.845`:

| Arm | A10 job | Binary | Physics change |
|---|---:|---|---|
| Accepted energy reference | **437266** | accepted q1 CUDA | none |
| Fixed-K closure | **437267** | corrected CUDA fixed-K | opt-in entropy mode and K initialization |
| Direct-U arithmetic | **437268** | accepted q1 CUDA | opt-in double RK4 stage anchor |

The [energy/fixed-K manifest](/gpfs/kjhan/LagForce/Kepler_EntropyK_A10_20261027/manifest.json)
and [direct-U manifest](/gpfs/kjhan/LagForce/Kepler_DirectU_A10_20261027/manifest.json)
record all input hashes. Their staged `scripts.sha256` and `inputs.sha256`
checks, shell syntax, and Python auditor imports passed. The new A10
runners differ from the original H200 runners only in the requested
GPU/partition, job name, and device log; physical parameters, RK4 controls,
auditors, and restart step are unchanged.

At the first live check, **437266/437267 were RUNNING on distinct A10 cards
of `syn03`**; **437268 was pending the per-user A10 GPU cap**. Both live
logs show native `[GPU tess]`; fixed-K uniquely reports
`[ENTROPY_K_INIT] rank=0 initialized=147456 retained=0`. The old H200
energy and direct-U jobs **437231/437239** are held to prevent duplicate
evolution; release them only if the A10 replacements fail after diagnosis.
At staging time no `.845` physical or pairwise result was available. The
required full physics, radial, fixed-K payload/EOS, and direct-U anchor
audits are now reported below.

## Exact t=.845 measured result

All three corrected CUDA A10 jobs **437266/437267/437268** finished native
evolution at exact `t=.845` with exit 0. Their full native physics and
radial receipts pass every check, and direct-U's stage-anchor receipt passes
all seven branch-use checks. The independent [fixed-K pairwise]
(KEPLER_Q1_A10_ENTROPY_K_PAIRWISE_T0845_2026-10-27.json) and
[direct-U pairwise](KEPLER_Q1_A10_DIRECT_U_PAIRWISE_T0845_2026-10-27.json)
auditors compare **new error accumulated only over the common `.825→.845`
interval**, not each arm's inherited error before `.825`.

| q1 arm | Disk mass-weighted absolute incremental log-entropy | Signed increment | Native wall |
|---|---:|---:|---:|
| Accepted energy | `7.29314e-6` | `-3.32828e-7` | `623.252 s` |
| Fixed-K closure | `6.71927e-8` | `-2.67925e-9` | `942.163 s` |
| Direct-U stage anchor | `3.63099e-6` | `-2.34256e-7` | `592.672 s` |

Fixed-K suppresses this short disk entropy increment by about **108.5×**,
and direct-U reduces it by about **2.01×** relative to the matched energy
reference. All arms have essentially unchanged disk trajectories and exact
fixed masses; fixed-K versus energy has mass-weighted disk position
difference `3.94e-14` and zero stored density difference at this endpoint.
Fixed-K stores nonzero K for all 147,456 cells; initial K matches the
accepted state to `5.93e-8` maximum relative error and final pressure
matches K-EOS to `1.55e-7` maximum relative error. Direct-U's stored anchor
tracks authoritative final U and is absent at the accepted start as intended.

**Interpretation:** fixed-K is an isentropic closure and is not proof that
the original energy equation's spurious heat source has been repaired; it
cannot serve as the general shock-capable solver. Direct-U changes RK4
floating-point staging while retaining the energy equation, so it is the
more conservative candidate for a later, same-IC long comparison. The
observed wall figures come from separate A10 cards and are not a matched
kernel throughput benchmark. Neither `.02`-time probe tests the historical
`t>6` breakdown. The original held H200 duplicates **437231/437239** were
canceled after the A10 comparisons passed; no accepted output was deleted.
