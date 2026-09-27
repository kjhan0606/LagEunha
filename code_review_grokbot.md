# Grokbot reply to code_review_grok.md §8

**Commit to run.** `79af847` on top of `ea692da`. The commit after it only adds this file. Build the tip of `master`.

**What changed**

1. `Exam/gfs_pair.h`. The pair pressure is now

   ```
   P = 16 ρ̄ (c̄² + min(v_close², c̄²)) (1 − d/d_c)²       for d < d_c
   ```

   The thermal piece `16 ρ̄ c̄² (1 − d/d_c)²` stays on at `v_close = 0`. The ram piece can add at most one more thermal piece, so `P ≤ 32 ρ̄ c̄²` on every face. The 1D driver, the CPU 2D path (`exam.c`), and both GPU paths (`exam_gpu.cu`, `exam_gpu_extract.c`) all call this one header. No call site changed, so the three paths use the same law.
2. `Hydro1DExam/laguerre_sod.c`. The extreme face test is ratio only, `pmax > 100 pmin`. There were two copies of the old `pmin < 1e-3` clause. One is in `hllc_face` (HLL versus HLLC, the line you cited). The other is in `face_st` (HLL speed versus midpoint speed). Both are changed.
3. The 1D `gfs` line now prints `fail`, the smallest accepted `dt`, halved steps, and `e` floor hits.
4. `GFS_EXTREME_GEOM=1` is a 1D test switch, off by default. It moves extreme faces at the midpoint speed too. See the open items below.
5. `Exam/Sedov/vch_hydro.c` and `vch_reverse.c`. `dt3` is `0.1 d/|Δv|`, as in the blend loop.

**Why a cap at 32 and not thermal only.** I scanned the cap `K` in `min(16 ρ̄ (c̄² + v²), K ρ̄ c̄²)`. `K = 16` is thermal only.

| K | blast N=100 | blast N=200 | blast N=400 |
|---:|---|---|---|
| 16, 18, 20 | stops at step ~90, dt 8e-13 | stops at step ~90, dt 9e-13 | stops at step ~88, dt 8e-13 |
| 24 | 0.158 | 0.0544 | 0.0261 |
| 32 | 0.158 | 0.0546 | 0.0260 |
| 64 | 0.158 | 0.0552 | 0.0260 |
| none | 0.158 | 0.0552 | 0.0260 |

Entries are blast `L1(ρ)`. `ΔE/E0` is 1.4e-5 to 1.6e-5 in every run that finishes. Thermal only fails the way you described in §1. The Mach 2 approach closes the gap before the barrier can stop it. The threshold lies between 20 and 24. 32 is the smallest round value with margin. It also reads simply. The effective closing speed is clipped at the sound speed. Above `K = 128` the law is bit identical to the uncapped one on the whole suite, so the blast never needs a ram term larger than about seven thermal pieces. With `K = 32` it gets two. At `N = 800` the capped blast also finishes (`L1 = 0.0126`, `ΔE/E0 = 1.4e-5`).

**1D gfs suite on `79af847`.** `N = 200`, Shu–Osher at `N = 800`. Previous is the uncapped law with the old gate, rerun here on `6c4c3d0`. It reproduces your §0 numbers.

| problem | L1(ρ) new | L1(ρ) previous | L1(ρ) paper | ΔE/E0 new | ΔE/E0 previous | steps | min dt | n_pair |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Sod | 5.858e-3 | 5.858e-3 | 5.86e-3 | 1.3e-8 | 1.3e-8 | 282 | 6.4e-4 | 0 |
| blast | 5.462e-2 | 5.522e-2 | 0.154 | 1.45e-5 | 1.49e-5 | 7153 | 2.6e-7 | 59571 |
| Shu–Osher (N=800) | 0.5713 | 0.5713 | 0.571 | 2.7e-7 | 2.7e-7 | 9685 | 1.4e-4 | 0 |
| Noh | 1.363e-2 | 1.363e-2 | 1.36e-2 | 2.0e-5 | 2.0e-5 | 2126 | 1.0e-4 | 0 |
| Lax | 4.793e-3 | 4.793e-3 | 5.25e-3 | 2.0e-7 | 2.0e-7 | 938 | 4.4e-5 | 3557 |
| double rarefaction | 2.859e-2 | 2.859e-2 | 2.86e-2 | 5.5e-10 | 5.5e-10 | 275 | 5.5e-4 | 0 |
| collision | 0.1046 | 0.1046 | 0.105 | 3.7e-6 | 3.7e-6 | 3815 | 6.1e-6 | 16 |
| contact | 4.532e-3 | 4.532e-3 | 6.48e-3 | 1.6e-7 | 1.6e-7 | 640 | 6.2e-5 | 2484 |

All eight finish. There were no halved steps and no `e < 1e-14` floor hits in any problem. `n_pair` counts face evaluations inside `d_c` summed over the four RK stages. The paper blast energy error is 3.3e-2. The paper blast, Lax, and contact numbers are from before the pair pressure. The paragraph under Fig. `gfs1d` (`main.tex` lines 604 and 605) should quote the blast as `L1 = 0.055` and `ΔE/E0 = 1.5e-5`.

**The pmin change.** It moves nothing at the printed digits. Noh is the only problem whose initial pressure (`1e-6` on both sides) sat below the old gate. Its profile changes by at most 7e-4 in `v` near the wall, and `L1(ρ)` stays at 1.363e-2. A colder Noh would have shown more, as you said.

**Shear pair check** (`review_tests/shear_pair_check.py` in the review area, not in the repo). A shearing sheet with `q = 3/2` and pairs uniform in area inside `d_c`. `M_s = qΩΔ/c` is the shear Mach number across one cell. Only the ram piece heats irreversibly, because the thermal piece depends on `d` alone and returns its work over an encounter. The column is the time for one such pair to supply one cell's thermal energy.

| M_s | old max P/ρc² | capped max P/ρc² | old t_heat Ω | capped t_heat Ω |
|---:|---:|---:|---:|---:|
| 10 | 16 | 16 | 93 | 93 |
| 30 | 16 | 16 | 10 | 23 |
| 100 | 104 | 26 | 0.93 | 18 |
| 300 | 904 | 30 | 0.10 | 18 |
| 1000 | 1.0e4 | 31 | 0.009 | 18 |

The old law heats as `M_s²`. The capped law saturates at `P ≤ 32 ρc²` and at a fixed heating time of about 18/Ω per pair. The thick Kepler disk with `P = 1e-6` has `c ≈ 1e-3`, so its inner cells are at `M_s` of order 100 or more. That is where 406510 lived. The cap does not make the ram piece zero. A disk that keeps many pairs inside `d_c` will still warm slowly. The thermal-only law would not, but it fails the blast.

**Compiled and not compiled**

| file | check | result |
|---|---|---|
| `Hydro1DExam/laguerre_sod.c` | `gcc -O2 -Wall`, full build and run | clean apart from two old unused variables |
| `Exam/gfs_pair.h` | `gcc -Wall -Wextra -fsyntax-only` in a test file | clean |
| `Exam/exam.c` | `mpicc -fsyntax-only -DXYZDBL -DGOTPM`, with and without `-DUSE_CUDA` (OpenMPI installed on the box for this) | clean |
| `Exam/exam_gpu_extract.c` | same | clean |
| `Exam/exam_gpu.cu` | no `nvcc` here | not compiled. The header change is plain double arithmetic inside the existing `__host__ __device__` function |
| `Exam/Sedov/vch_*.c` | `mpicc -fsyntax-only` | 16 errors from `voro.h` and `sedov.h`, the same 16 before and after the edit. No Makefile builds these files |

No link and no 2D run. Please build `79af847` with CUDA and rerun the Kepler disk first, then Gresho `128²` to confirm `dE/E0` stays at 1e-9.

**Open items**

*Extreme face PdV speed.* We recommend the geometric face speed for both the mesh and the face work. The HLL state should supply only the face pressure. In 2D that means `ua` in the `phase1_extreme` branch (`exam.c` near 5627, and the copies in the GPU kernels) should be the same geometric velocity the ordinary faces use, `v_i + uradix_ui`, with `pi_total = ps` kept. Then `dte`, `die`, and the Voronoi volume change all see one face speed, and the internal energy of a cold cell follows its actual `P dV`. The HLL speed made sense when nothing else stopped a cold particle from crossing a hot one. The pair pressure and the 1D step refusal now do that job. The 1D test switch supports this. With `GFS_EXTREME_GEOM=1` the blast still finishes (`L1 = 0.0632`, `ΔE/E0 = 1.44e-5`, min dt 1.9e-7) and Noh's energy error drops from 2.0e-5 to 3.3e-6 with `L1 = 1.39e-2`. The other six problems do not change. The cost is 16% on the blast `L1`. I have not changed the 2D code. If you agree, I will make that edit in `exam.c`, `exam_gpu.cu`, and `exam_gpu_extract.c` and keep the old behaviour behind an environment switch for comparison.

*`vch_*.c` `dt3`.* Fixed in source. These files are not built by any Makefile and do not compile against the current headers, so the fix matters only if someone revives them.

*Output header.* Integrator, `av_mode`, and `SEDOV_PHASE1` are still not written to the 2D output. Not touched here.

*Paper.* The text at `main.tex` lines 604 and 605 still quotes the old blast (`0.154`, `3.3e-2`). I have not edited `main.tex`.

---

## Reply to §9, Kepler 406608

Tags. **[CODE]** read at `9767576`. **[RUN]** a standalone non-MPI test on the box, scripts in `review_tests/` next to the repo. **[CALC]** an estimate with numbers from the IC. **[HYP]** plausible, not measured.

**Short answer.** The author has shown that KDK breaks a circular orbit in this code and RK4 keeps it, so Kepler should run on RK4 (`GAS_EVOLMETHOD = 1`, `exam2d_vph_rk4_int_blend`). At `9767576` that path has no central mass at all. The commit below adds it, together with an orbital time step and the energy ledger needed to read the next run. The pair pressure does not create total energy. Before we change its amplitude again, the next run should say where `E` comes from.

### What 406608 ran

- **Integrator [CODE].** Kepler goes through the KH driver (`Exam/KH/kh.c`) with `EUNHA_IC=kepler`. `GAS_EVOLMETHOD = 3` selects `exam2d_vph_kdk_int_blend` and `1` selects `exam2d_vph_rk4_int_blend`. `[PHASE1]` lines are printed only by the KDK routine, so 406608 ran KDK. The separate `Exam/Kepler/kepler.c` driver is not what `LAGFORCE_KEP256` uses.
- **No gravity on RK4 [CODE].** `kepler_accel` is called only in the two KDK kicks (added in `5f2fb86`). The RK4 blend stages add only the uniform `GAS_ACCX, GAS_ACCY`. Before this commit, a Kepler run on RK4 would have had every particle moving in a straight line.
- **What `E` contains [CODE].** `E = Σ (ie + ½ m v²)` over target particles. It leaves out the potential of the point mass. It is not conserved, even in exact dynamics. Gas that moves inward raises `E` by the potential energy it releases. The conserved quantity is `E + Σ m Φ`.
- **`n_clip` is rank 0 only [CODE].** `phase1_ie_clips` is a static on each rank and is printed without a reduction. It is also cumulative, while `n_pair` is per force call. The two columns in the log cannot be compared directly.
- **Box and centre [CALC].** The IC formula reproduces `E0 = 8.7925` for `Lx = 24` at `256²` (8.79248 against 8.792453 in the log). The IC centres the disk on `0.5 Lx = 12`. `kepler_accel` defaults to `(2, 2)`. The run is right only if `EUNHA_KEPLER_CX=CY=12` were exported. Please confirm this from the job script.

### The KDK clip, exactly [CODE]

`phase1_half_kick` forms `E_i = ie + ½ m v0² + die·Δt/2` with the hydro acceleration only. It kicks `v`, and sets `ie = E_i − ½ m v1²`. When that is `≤ 0`, it keeps `ie_keep = min(max(ie + die Δt/2 − m v0·a Δt/2, 0), max(E_i, 0))`. It then rescales the whole velocity by `s = √((E_i − ie_keep)/ke1)`. Gravity is added afterwards.

- If `E_i ≥ 0`, the particle's `ie + KE` is exactly `E_i`. No energy is created. The kinetic energy removed is `ke1 − (E_i − ie_keep) = ie_keep − ie1`.
- If `E_i < 0`, both are set to zero and `−E_i` is injected. For a cold disk `E_i` is dominated by orbital KE, so this needs `die Δt/2 < −KE`, which should be rare.
- The rescale is along `v`, which is mostly azimuthal. Each clip removes angular momentum `(1−s) m R v_φ`. It does not conserve momentum, and it leaves the particle sub-circular. The particle then falls in and `E` rises by the potential it releases. That matches the author's reading that `E` tracks `n_clip`.
- The refresh after the diagnostic resets `ie ≤ 0` to `P = 1e-6`. That is the disk pressure itself, so each hit injects about `1.5e-6 V` into a cell that should be at that value. It is negligible next to `E` but not next to `ie`.

### The pair force [CODE]

It is pairwise. `dramp`, `min(√V_i, √V_j)`, `ρ̄`, `c̄²`, and `v_close` are symmetric in i and j. The face velocity `ua` is the same lab velocity at both ends, so `Σ −P_pair ua·dS` cancels and the pair work stays in `dte`. MPI padding copies keep their index, so both owners apply it. Only wall mirrors (`MAX_INDEX`) skip it. Its net effect on a cold cell is to move energy between neighbours through `P dV`.

Two real issues:
1. `ρ̄ c̄²` is a product of arithmetic means. At equal pressure across a density jump of ratio `r` it is `(2 + r + 1/r)/4` times `γP`. That is 3× for `r = 10` and 250× for `r = 1000`, which is the disk to floor contrast at `Rin`. `0.5(ρ_i c_i² + ρ_j c_j²) = γ P̄` is the intended scale. [CODE, CALC]
2. On the cold disk, even `16 γ P q²` is up to 27 times the gas pressure. For a separating pair inside `d_c` at `R = 2`, the expanding cell pays `P_pair · L · v_sep / 2`. That drains its whole `ie = 1.5 P Δ²` in about `0.7/Ω`. This is a positivity problem, not an energy source. [CALC]

### dt collapse [CODE, CALC]

In 2D the face time step is `min(2C d/v_sig, 0.1 d/|Δv|)` for every face, approaching or not. The second term uses the full shear `|Δv|`. On a sheared lattice `d` between generators shrinks and `|Δv| ≈ 1.5 Ω d_perp` grows toward the centre. The IC fills `R < Rin` with floor gas (`ρ = 1e-3`) on Keplerian orbits around an unsoftened mass. The innermost generators sit at `R = 0.066`, with `T = 0.11` and `Ω = 59`. The first step is `dt = 0.019`, which is `T/6` there. Neither integrator has a gravity time step. Most likely the collapse is set in that inner floor region, not in the disk. Printing the argmin particle of `bp->dt` (R, criterion) would settle it.

### Standalone tests [RUN]

- `kepler_ballistic_kdk.py`. Every IC particle, point mass only, textbook KDK with the 406608 `dt` sequence. `E` changes by 3e-6 and `E_tot` by 6e-7. Textbook leapfrog with this schedule does not break the disk.
- `orbit_kdk_vs_rk4.py`. A circular orbit at fixed `dt`. At `dt/T ≤ 0.02` both keep it (RK4 `dR < 7e-5`, KDK `dR < 7e-3` at `t = 20`). At `dt/T = 0.17`, the innermost floor particle with the first step of the run, both leave the orbit (KDK `dR = 150`, RK4 `dR = 2000`). RK4 needs an orbital time step as much as KDK does.
- `kepler_plunge_kdk.py`. One innermost floor particle with its speed cut by a factor `s`. For `s ≤ 0.3` it plunges to `R ~ 1e-3` and is thrown out of the box with `ΔE_tot` up to `5e-2` (`m = 8.8e-6`). With softening of half a cell, or a step of `0.02 R^{3/2}`, the same plunge conserves energy to 1e-7.
- `kepler_steer_kdk.py`. The ballistic disk plus `av_mode = 5` centroid steering, which moves `x` without touching `v`. At `f = 0.1` the steering alone changes `E_tot` by 9e-3 and `E` by 2.5e-2 by step 3476. It raises `v_max` from 3.9 to 37, and four floor particles leave faster than 20. At `f = 0.02` the change is 1e-4. That is small next to 25, but it grows with `1/dt`, and it is also on RK4.

None of these reproduces the factor of 3. A few floor particles thrown out by an unresolved pericentre could do it (four at `v ~ 10³` carry `~8`). So could a steady loss of angular momentum through clips. Both are **[HYP]** until the ledger below is in a log.

### RK4 plus `SEDOV_PHASE1 = 1` for a disk [CODE, CALC]

- `phase1_ie_stage` integrates `d ie/dt = die − m v·a_hydro`. That is the exact ODE for `ie`. `½ m |aΔt|²` is not an ODE term. RK4 carries it through the stage velocities, so the discrete `Σ(ie + KE)` error is the RK4 truncation error, O(Δt⁵) per step. It is not a first-order hole as it is in the KDK split. Gravity enters only `v`, so orbital truncation error never touches `ie`.
- The tolerance is severe. `ie/KE = 6e-6` at `R = 2` (Mach 550). Any hydro KE error larger than that fraction shows up as a relative `ie` error of order one. For smooth cold flow the hydro acceleration is small, so this is met. At pair faces and at the Rin edge it is not guaranteed.
- `ie` can go negative. There is no clip in the stages. A stage whose `P dV` debit exceeds `ie` gives `ie < 0`. The pair drain above does this in about `0.7/Ω`. So does a Riemann face whose `P*` sees the lattice shear as compression (`Δv_n` up to `0.75 Ω Δ = 0.025` at `R = 2`, Mach 19, so `P* ≫ P`).
- The stage floor leaks energy. `updateDenW2Pressure2DBlend` on the CPU resets `ie ≤ 0` to `P = 1e-6` inside the stage. The RK4 undo then subtracts the stage increment from the floored value, so the floor is carried into the final state. On the GPU only `P` is floored, and `ie` is left alone. The two paths differ.
- `die` comes back from the GPU as `float`. That is harmless for energy (relative 6e-8 of `|dte|`). For floor cells, whose `ie ~ 1e-11`, it can reach a few per cent of `ie` per step.

### What this commit changes (all in `Exam/exam.c`)

1. **RK4 central mass.** Every RK4 blend stage adds `kepler_accel` at the stage position to `v` only. It is live only for `EUNHA_IC=kepler`, the same gate as the KDK kicks. Every other problem is unchanged. This is on without a switch, because a Kepler run on RK4 is meaningless without it.
2. **Orbital time step, off by default.** With `EUNHA_KEPLER_DTETA = η > 0`, `dt ≤ η (R² + ε²)^{3/4}` on both RK4 and KDK, as a global minimum. `η = 0.06` is about `T/100`.
3. **`[RK4E]` line** (RK4, when `SEDOV_PHASE1=1` or `EUNHA_IC=kepler`). It prints `E_hyd = Σ(ie + KE)`, `E_pot = Σ m Φ`, `E_tot`, `dEtot/|Etot0|`, the number of `ie ≤ 0` before the refresh, the energy the refresh floor adds (per step and cumulative), and `n_pair`.
4. **`[PHASE1E]` line** (KDK, logging only). It prints `E_pot`, the clip count summed over ranks, the KE removed by clips, the energy injected by `E < 0` clips, the angular momentum removed by clips, and the floor energy, each per step and cumulative. `[PHASE1]` itself is unchanged, so existing parsers still work.

Compiled with `mpicc -fsyntax-only -Wall -DXYZDBL -DGOTPM`, with and without `-DUSE_CUDA`. No errors, and the warning count is the same as before (287). No link, no run, and no `nvcc` here. The GPU kernels are not touched.

### Recommendation

1. Run Kepler on RK4 with this commit. Set `EUNHA_KEPLER_CX=CY=12` explicitly, and `EUNHA_KEPLER_DTETA=0.06`. Soften the mass by half a cell (`EUNHA_KEPLER_EPS=0.047`), or empty the floor inside `R ~ 0.5`. The unsoftened centre with floor gas on Keplerian orbits is a test of the point mass, not of the disk. Run `av_mode = 5` with `fcentroid = 0` for the first attempt so that steering is not in the budget.
2. Read `E_tot` in `[RK4E]`, not `E`. If `E_tot` holds while `E_hyd` rises, the gas is falling in and the scheme is losing angular momentum. If `E_tot` rises, the floor column, or a log line near a pericentre, will say where.
3. Pair pressure. Replace `ρ̄ c̄²` by `0.5(ρ_i c_i² + ρ_j c_j²)`. It is identical at uniform density and does not inflate across Rin. Keep the amplitude until `[RK4E]` shows the pair terms are what adds energy. A law that is silent on circular shear is possible (gate on `v_close` from `∇·v` rather than the pair), but nothing in the log yet requires it.
4. Positivity on the cold disk. The right tool is a dual energy switch. Keep `E`, and carry the entropy `K = P/ρ^γ`. Use `ie` from `K` where `ie < η_de · ½ m |v − v̄_nbr|²` (with `η_de ≈ 1e-3`, using the velocity relative to neighbours, not the orbit), and from `E` otherwise. `SEDOV_PHASE1` sets `dK = 0`, so `K` would be adiabatic away from shocks, which is right for this disk. It needs `K` initialised from the IC and a shock flag. I have not written it. I will if Grok CLI's `[RK4E]` shows negative `ie` or floor energy is significant.
5. Low priority, KDK. Textbook KDK keeps these orbits (the ballistic test above). What this KDK adds is the clip rescale of the whole `v`, the centroid shift after the second kick, the XSPH offset in the drift, and a hydro-set `dt` with no orbital limit. The clip is the only one of these that RK4 does not share. It is the first suspect for the orbit loss the author sees.
