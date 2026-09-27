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

---

## Reply to §10–11, work limit and Kepler 406734

Tags. **[CODE]** read at `64ab6c1`. **[RUN]** run on the box (1D, or a standalone script in `review_tests/`). **[CALC]** an estimate from IC numbers. **[HYP]** plausible, not measured.

### Short answer

1. The work limit is correct, and it is silent on the whole 1D suite. All eight problems are bit-identical to cap-only. The blast still finishes (`L1 = 0.0546`, `ΔE/E0 = 1.45e-5`). It still finishes when the budget is squeezed to almost zero. **[RUN]**
2. The 1D call had a units bug. It passed the specific `e` as the budget. Fixed in `716f6bf`. The limit is still silent after the fix. **[CODE, RUN]**
3. In 2D the root, the symmetry, and `nshare = 6` are fine. `GAS_dtold` is the weak point: step 1 uses `dtold = 1e-7`, so step 1 is not limited, and the error scales as `Dtime/dtold` (linear for braking, squared at `q = 0`). **[CODE]**
4. The floor energy is not a slow leak. Each event adds about 8 cells' worth of `ie`, and `floor_cum` is already 3.5× the thermal energy of the whole box. I added a dual-energy switch (`GFS_DUAL_ENERGY`, off by default), `E_int` in `[RK4E]`, and a per-cell floor log (`GFS_FLOOR_LOG`) in `75f00ff`. The switch fixes the value a drained cell gets. On its own it will not remove the energy error. The next run is mainly meant to find where the events happen. **[CODE, CALC]**

### 1D suite with the work limit [RUN]

The section-9 `gfs_pair_*.dat` files are not on this box. I rebuilt cap-only from `79af847` and ran it next to the new code. `N = 200`, Shu–Osher at `N = 800`, run as `GFS_PREFIX=… ./laguerre_sod gfs 200` and `gfs 800 shuosher` from a directory that links the eight reference folders. Outputs are `gfs_pair_*.dat` (cap-only), `gfs_wlimbug_*.dat` (`64ab6c1` as committed), and `gfs_wlim_*.dat` (with the unit fix), all in the scratch area. No `gfs_*.dat` was written.

| problem | L1(ρ) cap-only | L1(ρ) work limit | ΔE/E0 cap-only | ΔE/E0 work limit | steps | dt_min | retries | n_pair | limited faces | e floor hits |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Sod | 5.858e-3 | 5.858e-3 | 1.26e-8 | 1.26e-8 | 282 | 6.4e-4 | 0 | 0 | 0 | 0 |
| blast | 5.462e-2 | 5.462e-2 | 1.45e-5 | 1.45e-5 | 7153 | 2.6e-7 | 0 | 59571 | 0 | 0 |
| Shu–Osher (N=800) | 0.5713 | 0.5713 | 2.75e-7 | 2.75e-7 | 9685 | 1.4e-4 | 0 | 0 | 0 | 0 |
| Noh | 1.363e-2 | 1.363e-2 | 2.00e-5 | 2.00e-5 | 2126 | 1.0e-4 | 0 | 0 | 0 | 0 |
| Lax | 4.793e-3 | 4.793e-3 | 2.05e-7 | 2.05e-7 | 938 | 4.4e-5 | 0 | 3557 | 0 | 0 |
| double rarefaction | 2.859e-2 | 2.859e-2 | 5.5e-10 | 5.5e-10 | 275 | 5.5e-4 | 0 | 0 | 0 | 0 |
| collision | 0.1046 | 0.1046 | 3.72e-6 | 3.72e-6 | 3815 | 6.1e-6 | 0 | 16 | 0 | 0 |
| contact | 4.532e-3 | 4.532e-3 | 1.61e-7 | 1.61e-7 | 640 | 6.2e-5 | 0 | 2484 | 0 | 0 |

All eight finish. The output files are byte-identical (`cmp`) for all three builds. That covers cap-only, `64ab6c1`, and the unit fix, plus the unit fix with `GFS_DUAL_ENERGY=0.5`.

- **The units bug [CODE].** `face_st` passed `P[i-1].e, P[i].e`. In the 1D code `e` is specific (`eos_e = p/((γ-1)ρ)`, `ae = …/m`). The 2D paths pass the cell's total `ie`, which is right. With specific `e` the budget was `1/m` too large (about 200/ρ at `N = 200`), so the 1D limit could never act. `716f6bf` passes `m e`. The `gfs` line now also prints `nwlim` (faces the limit reduced), `nde` and `deinj`.
- **Margin [RUN].** I logged the smallest `Jmax/(P0 A dt)` over every pair evaluation. The blast has 10.7, on a separating face (`q/c = +1.08`). Lax has 17, contact 16, and collision 54. The limit is an order of magnitude away from acting anywhere in 1D.
- **Blast with the budget squeezed [RUN].** I scaled `nshare` so that `B → 0`:

  | nshare | limited faces | L1(ρ) | ΔE/E0 | steps |
  |---:|---:|---:|---:|---:|
  | 2 (default) | 0 | 0.0546 | 1.45e-5 | 7153 |
  | 100 | 158 | 0.0552 | 1.49e-5 | 7154 |
  | 1000 | 4913 | 0.0568 | 2.38e-5 | 7460 |
  | 1e4 | 15955 | 0.0586 | 6.26e-5 | 8465 |
  | 1e6 | 256463 | 0.0633 | 2.62e-5 | 35499 |

  With `B ≈ 0` the blast still finishes, with no retries. Braking is always allowed up to a full reversal, `J = 2μ|q|`. So the limit cannot bring back the thermal-only stall near step 90. That stall came from too small a `P0`, and the limit never touches braking. The same holds at `N = 100, 400, 800` (`L1` 0.158, 0.0260, 0.0126, unchanged).
- `min(ie_i, ie_j)` in place of the sum (see below) is also silent on all eight problems (blast margin 9.6).

### Review of `gfs_pair_work_limit` [CODE]

**Root.** `W(J) = qJ + J²/(2μ)` is the change in relative kinetic energy. `W = B` gives `Jmax = μ(√(q² + 2B/μ) − q)`. That is the positive root, and it is at least `2μ|q|` when `q < 0`. The guards return `p0` when `dt`, `A`, or `m` are not positive. That is right for step 0, but see below.

**Symmetry and MPI.** `p0`, `A`, `q = (v_j − v_i)·e_r`, `μ`, and `B` are all symmetric under `i ↔ j`, and `GAS_dtold` is global. The owner and the ghost copy see the same `ie`, because padding copies the stage state. The CPU stage floor is applied to all `nbp`, ghosts included. So the limited `P` is the same on both ranks, and the force stays antisymmetric to round-off. The GPU kernel and `cpu_reference_force_csr` read `pie_in` / `parts->ie`, which are uploaded from `bp->ie` in the same force call as `v`. So `q` and `ie` are at the same stage level. Two small differences between CPU and GPU:
- The CPU stage floor has already reset `ie ≤ 0` to `1e-6 V/(γ-1)` when the force runs. The GPU floors only `P`, so its `B` term is 0 for that cell. This is tiny.
- `pie_in` and `pmass` are `float`. That is fine for `B`.

**`GAS_dtold` versus the step.** `kh.c` sets `GAS_dtold = 1e-7` before step 1, then `dt` after each step. On the RK4 path the stage-1 force runs before `Dtime` is known. Stages 2–4 use the same `dtold`, even though `Dtime` is known by then. The impulse applied over the step is `Σ w_k P_k A Dtime`, while the cap assumed `dtold`. So:
- `Dtime > dtold` (dt grows): the cap is too loose by `Dtime/dtold`. The work is too large by that factor for braking, and by its square at `q = 0`. After the 4.70e-4 step in §11, a return to 6.4e-4 is 1.37, which gives 1.9× the budget at `q = 0`. `nshare = 6` absorbs this unless one cell has several pair faces at the cap.
- Step 1: `Dtime/dtold ≈ 2e5` on Kepler. The limit is off for that step. This is harmless only because `n_pair = 0` at step 1.
- `Dtime < dtold` (DTETA, a sharp CFL drop): the cap is too tight by the same ratio. The energy is safe, but the barrier is weaker for one step.
- There is no dt growth cap in the blend integrators.

  Fix: add a global `gfs_pair_dt` that the RK4 driver sets to `dtold` before stage 1 and to `Dtime` before stages 2–4. Pass it to the kernel as a new param. Do not overload `params.dtold`, which also feeds `get2dUpqradRk4`. Also cap dt growth at 1.25× per step. Stage 1 carries weight 1/6, so the error that remains is `(Dtime/dtold − 1)/6`.

**`nshare = 6`.** A 2D Voronoi cell has 6 faces on average, and 4 to about 10 individually. Only faces with `d < d_c` draw on `B`. A cell rarely has more than 2 of those, so 6 is safe, even conservative. The real gaps are elsewhere:
- **The budget is summed, but the debit is split.** On the phase-1 path, every face uses `ua = vns` (HLL), so cell `i` pays `J (vns − v_{n,i})`. With a large impedance contrast, `vns` sits near the heavy side, and the light cell pays nearly all of `W`. With `ie_hot > 5 ie_cold`, `B = (ie_i + ie_j)/6` can exceed the cold cell's own `ie`. `B = 2 min(ie_i, ie_j)/nshare` closes this. It is a one-line change in the shared header and is silent in 1D.
- **Only pair faces are limited.** Riemann faces draw on the same `ie` with no cap.

### Where the floor energy comes from [CALC, HYP]

Numbers from §11 and the IC (`P = 1e-6` everywhere, `Lx = 24`, `256²`, `Δ = 0.094`, cell `ie = 1.3e-8`, box `Σ ie = 8.6e-4` for `γ = 5/3`):
- `floor_cum = 3.0e-3` at `t = 5.18` is **3.5× the thermal energy of the whole box**. After 0.3 of an inner orbit, the floor, not the hydro, sets the thermal state of the reset cells.
- Per step, `3.0e-3/7837 = 3.9e-7`, spread over 1–9 cells. That is about `1e-7` per event, or 8 cell `ie`. The cells are driven far below zero within one step. They do not creep down to it.
- Rate: `1.3e-3` per unit time over the last 0.5. Linear to `t = 177.7` gives `dEtot/|Etot0| ≈ 2.6e-2`. More if the pair count keeps rising.

Candidates, with numbers from `review_tests/hll_shear_face.py`. That script uses the same HLL average as `hll_star_state`. With `SEDOV_PHASE1 = 1` and `use_muscl = 0`, it is used on **every** face, with cell-centre states.
1. **Lattice shear seen as compression [CALC].** Keplerian shear gives a normal jump `Δv_n ≤ 0.75 Ω Δ` across a face. For `|Δv_n| ≫ c`, the HLL average has `P* ≈ (γ-1) ρ u c/2` on a diverging face and `P* ~ ρ u²` on a converging face (`u = |Δv_n|/2`). At `u/c = 10`, that is 4.5 P and 109 P. At `R = 2` in the disk, `Δv_n/c = 19`. The diverging face alone drains a cell in `0.9/Ω`. The converging faces heat it about 24 times faster. So a first-order disk should **heat**, turning orbital KE into `ie`. `E_tot` cannot show this, because it is conserved. `E_int` can.
2. **The central floor gas [CALC, HYP].** `ρ = 1e-3` on orbits around the softened mass. At `R = 0.05–0.1`, `Ω dt = 0.02–0.04` and `Δv_n/c = 40–50`. One HLL face changes `ie` by `+1.4…3.3 ie` (converging) or `−0.03…0.06 ie` (diverging) per step. The face set also turns over quickly there. A cell whose faces reorder between RK4 stages can land at `−several ie`. That matches the 8 `ie` per event and the small, fluctuating `n_ie_le0`. Most likely the floor events are at `R < 0.3` and not in the disk.
3. **Pair faces with an uneven split**, as above. This is bounded by `B`, so it is at most 1/6 of the cell `ie` per face and step. It cannot give 8 `ie`.
4. **The RK4 undo carries the stage floor [CODE].** CPU stage floors enter the final `ie` but are not counted in `floor_cum`. The GPU does not floor `ie` in stages. This matters for the ledger, not as a source.

The `[RK4F]` log added below tells 2 apart from 1 and 3 by `R`.

### Fix

**Why dual energy alone will not flatten `E_tot` [CALC].** A cell with `ie < 0` has already gained more kinetic energy than it had internal energy. Refilling it from `K` costs `ie_K − ie`. For this IC, `ie_K` is the same number as the `P = 1e-6` floor unless the density has changed. So the energy added per event is about the same. What dual energy changes is the value: the reset cell gets its adiabatic `P`, not an arbitrary 1e-6 (which is 1000× too hot for a compressed floor cell, and too cold for a compressed disk cell). It also names the energy it adds (`de_cum`). It is the standard fix for thermodynamics. It is not an energy fix.

**What is in `75f00ff` [CODE].** All of it is off by default. Without the two variables, the only change is four extra columns at the end of `[RK4E]`.
- `GFS_DUAL_ENERGY = η > 0`:
  - On the first RK4 blend call, `stress.K = P/ρ^γ` for every particle with `K ≤ 0`. The KH/Kepler IC leaves `K = 0` (the stress block is zeroed). `SEDOV_PHASE1` keeps `dK = 0`, so `K` stays the IC entropy.
  - In the final refresh (host code, so it also acts on GPU runs), `ie < η ie_K` takes `ie_K = K ρ^γ V/(γ-1)`, before the `P = 1e-6` floor.
  - The CPU stage floor in `updateDenW2Pressure2DBlend` uses `K ρ^γ` in place of `1e-6`. The GPU stage path is not touched.
- `[RK4E]` gains `E_int = Σ ie`, `n_de`, `de_inj`, and `de_cum`. The existing keys are unchanged. The ledger is now `dEtot·|Etot0| ≈ floor_cum + de_cum`. With `η > 0`, `n_ie_le0` should be 0.
- `GFS_FLOOR_LOG=1` prints up to 8 `[RK4F]` lines per rank per step, before the change: `x, y, R` (from the Kepler centre), `den`, `ie`, `ie_K`, `ie/ie_K`, `vol`, and `reset`.

**1D smoke test [RUN].** The same rule is in `laguerre_sod.c`, with `K` fixed at `t = 0`:

| η | problems changed | notes |
|---:|---|---|
| 0.5 | none (byte-identical) | `nde = 0` in all eight |
| 0.7 | none | |
| 0.9 | Noh only, ΔE/E0 2.0049e-5 → 2.0052e-5 | |
| 1.0 | six of eight | entropy floor. Blast ΔE/E0 1.45e-5 → 2.6e-3, Sod 2.8e-4. Contacts and rarefactions lose particle entropy numerically, and η = 1 refills it |

Use `η = 0.5`. It never fires on a problem where the scheme works, and it catches a cell that has lost half its adiabatic energy. Nothing adiabatic does that.

**Compile [RUN].** `mpicc -fsyntax-only -Wall -DXYZDBL -DGOTPM`, with and without `-DUSE_CUDA`: `exam.c` 287 warnings both ways (unchanged), 0 errors. `exam_gpu_extract.c` is unchanged (0 and 16 warnings). No link and no 2D run: the full build needs the Intel toolchain and CAMB. `exam_gpu.cu` is not touched. The script is `review_tests/syncheck.sh`.

**The energy fix, design only [CODE].** Limit every face's pressure by the payer's budget. The idea is the same as the pair limit, but it covers all of `pi_total`. It conserves energy and momentum exactly, because the force and `dte` use the same reduced `P`.
```
s_i = (ua − v_i)·n            debit rate of i  = P A s_i   (n outward from i)
s_j = (v_j − ua)·n            debit rate of j  = P A s_j
P  ≤ ie_i /(nshare A s_i dt)  if s_i > 0,   same for j
```
`ua` is the HLL `vns`, the same on both ranks. `ie_j` is the ghost copy. So the factor is symmetric and needs no extra MPI. With `nshare = 6` and forward Euler this keeps `ie ≥ 0`. `dt` should be the `gfs_pair_dt` above. It belongs in `gfs_pair.h` as one `__host__ __device__` function, called from the same three sites as the work limit. I have not written it, because it needs the GPU kernel and `nvcc` to check. If `[RK4F]` shows the events at `R < 0.3`, the cheaper fix is not numerical. Soften harder, or empty or freeze the floor gas inside `R ≈ 0.5`. The disk starts at `Rin = 2`, and gas at `R = 0.05` with `Ω dt = 0.04` and Mach 50 face jumps is not part of the test.

### Recommendation for the next 2D run

1. Keep 406734 running as it is. It is the baseline for the floor rate.
2. Keep 406515–406519 held until a disk passes: bounded `dEtot/|Etot0|`, and `E_int` within a factor of a few of `8.6e-4`, through a few inner orbits.
3. Next Kepler run, when a slot is free: the tip of `master`, the same settings as 406734, plus `GFS_DUAL_ENERGY=0.5` and `GFS_FLOOR_LOG=1`. Please also report `GAS_USEMUSCL` and `SEDOV_PHASE1` from the job, because the HLL-average argument applies only to `use_muscl = 0`.
4. Columns to watch in `[RK4E]`:
   - `E_int`. If it climbs well above `8.6e-4`, the first-order faces are heating the disk from orbital KE (item 1). Then the fix is MUSCL on the disk, not the floor.
   - `n_de` and `de_cum` together with `floor_cum`. `floor_cum` should stop growing, and `de_cum` then carries the injection.
   - `dEtot/|Etot0|`, against `(floor_cum + de_cum)/|Etot0|`.
   - `n_pair`.
5. From `[RK4F]`: a histogram of `R` (split at 0.3, 1.9, and 2.1) and of `ie/ie_K`. That decides between central excision and the face limit.
6. Header follow-ups, not yet committed: `B = 2 min(ie_i, ie_j)/nshare`, and `gfs_pair_dt` together with a 1.25× dt growth cap. Both are small. I will make them together with the face limit once you confirm the GPU side can be built and checked.

### Commits

- `716f6bf` 1D: pair work limit uses `m e`. `GFS_DUAL_ENERGY` smoke switch and counters.
- `75f00ff` RK4 blend: `GFS_DUAL_ENERGY`, `E_int` / `n_de` / `de_inj` / `de_cum` in `[RK4E]`, and `GFS_FLOOR_LOG`.
