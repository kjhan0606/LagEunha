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

---

## Reply to §12–14, the hole gas and a GIZMO-style disk

Tags. **[CODE]** read in the source (ours at `2d56133`, GIZMO in `/workspace/gizmo`). **[RUN]** run on the box. **[CALC]** an estimate from IC numbers. **[HYP]** plausible, not measured.

### Short answer

1. I agree with §14. The created energy comes from the cold floor gas in the hole, `R = 0.5–1.9`. It does not come from the disk. We do not freeze that gas. **[RUN, from your logs]**
2. GIZMO never meets this problem in its own disk run. For MFM/MFV the hole is empty: particles are removed inside `r < 0.5`. The floor `Σ ≥ 0.01` and the soft edges of Eq. 34 are only for mesh codes. **[CODE, paper §4.2.4, footnote 24]**
3. I added the two GIZMO safeguards, both off by default: an entropy switch for cold, gravity-dominated cells, and a half-limit on `ie`. I also added the Hopkins mesh-code disk as an IC option. **[CODE]**
4. When off, nothing changes. The 1D suite is byte-identical, and the compile warnings are the same. **[RUN]**
5. Neither safeguard conserves energy. What they do is stop spurious heating from reaching the pressure, and stop the refill to `P = 1e-6`. The new `[RK4E]` columns measure the energy they add or remove. **[CODE, RUN]**

### What GIZMO does [CODE]

- **The hole is vacuum.** Footnote 24: particle codes use a sharp ring with nothing inside `r = 0.5`. Eq. 34 (`Σ = 0.01 + …`, the `(r/0.5)^3` inner ramp, the `[1 + (r−2)/0.1]^−3` outer edge, `P = 1e-6`, softened potential) exists because "most mesh-based schemes require non-vacuum boundaries".
- **Moving meshes fail this test too.** Hopkins ran more than 200 FVMHD3D variants and some AREPO runs. With the simple setup, "the disk goes unstable and the angular momentum evolution tends to be corrupted within a few orbits". Our cliff at 0.6 of an inner orbit is earlier than that, but it is the same family of failure.
- **The entropy switch** (`kicks.c` 261–312, macro `ENERGY_ENTROPY_SWITCH_IS_ACTIVE`):
  - `e_potential = m |g| h`. Here `g` is gravity only and `h = Get_Particle_Size` (about `sqrt(area)` in 2D).
  - `e_thermal = m · max(0.5 u_old, u_new)`, where `u_new` is the energy-equation result.
  - If `0.01 · e_potential > e_thermal`, then `du/dt = −(P/ρ) ∇·v`. That is the adiabatic change and nothing else.
  - A kinetic-energy criterion is computed, then overwritten (`do_entropy = 0`). The comment says the gravity form was cleaner.
  - The macro is **commented out** in `allvars.h`: "even for pure hydro, this isn't recommended". The paper (App. D) used `α ≈ 0.001`, "almost never triggered". The code line has 0.01.
- **The half-limit** (`kicks.c:336`, `predict.c:226`): if `u_new < 0.5 u_old`, then `u = 0.5 u_old`. This one is always on. After it, only `MinEgySpec` applies, which is 0 by default. There is no refill.

### What is new

All of it is off by default.

**Entropy switch and half-limit, RK4 blend path** (`Exam/exam.c`) **[CODE]**

| env | default | meaning |
|---|---|---|
| `GFS_ENTROPY_SWITCH` | 0 | 1 turns the switch on |
| `GFS_ES_COEF` | 0.01 | coefficient in `coef · m|g|h > max(0.5 ie_n, ie_E)` |
| `GFS_HALF_LIMIT` | same as `GFS_ENTROPY_SWITCH` | `ie < 0.5 ie_n` becomes `0.5 ie_n`. Set 0 or 1 to override |

- `g` is `kepler_accel` plus `GAS_ACC`, at the new position. The hydro acceleration is not gravity and is not included. `h = sqrt(V)`. `ie_E` is the energy-equation value.
- A switched cell gets `ie = ie_n (V_n/V)^(γ−1)`. This is `K = P/ρ^γ` held over the step at fixed mass. It is the exact integral of GIZMO's `−(P/ρ)∇·v`.
- **Where it acts.** Once per step, on the host. It comes after the final RK4 combination, `exam2dUpdateVol`, and the `w2` update. It comes before `GFS_DUAL_ENERGY` and the `P = 1e-6` floor. That is where GIZMO's kick applies it.
- **Why not in the stages.** Anything written into `ie` during a stage survives the RK4 undo (`ie −= k3ie` only removes the increment). This is how the CPU stage floor already leaks into the final `ie`. Acting once, on the final state, cannot leak. It also replaces that leaked stage floor for a switched cell, and bounds it for any other cell.
- `ie_n` and `V_n` are saved right after the stage-1 density update. They must travel with the particle, because `postStage` migrates particles every stage. So they are kept in `stress.E_inv_xx` and `E_inv_xy`. Only LagMFM (`av_mode 4`) uses those fields. The GPU path uploads them but writes them back only in the LagMFM density call. The switch refuses to run for `av_mode 4` and for `entropy_mode 1`. I did not change the struct, so a host-only relink of `exam.c` stays safe.
- **GPU runs.** The GPU kernels return `die` and `a`. The RK4 stages and the final combination of `ie` are on the host. So the switch works in GPU runs, the same way as `GFS_DUAL_ENERGY`. The GPU stage path is not touched. Within a step, stage pressures still come from the energy-equation stage `ie`.
- **Not covered.** The KDK blend path. Also, if `GAS_FCENTROID > 0`, the centroid shift changes `V` with no flow, and a switched cell treats that as compression.
- **Log.** When either switch is on, a banner `[ESW]` is printed once. Six columns are added at the end of `[RK4E]`: `n_es es_de es_cum n_half hl_de hl_cum`. `es_de` is `Σ(ie − ie_E)` from the switch, and `hl_de` is the same for the half-limit. The existing keys do not change. The ledger is now `dEtot·|Etot0| ≈ floor_cum + de_cum + es_cum + hl_cum`.

**Hopkins disk IC** (`Exam/KH/util.c`, the `EUNHA_IC=kepler` branch) **[CODE]**

| env | default | meaning |
|---|---|---|
| `EUNHA_KEPLER_PROFILE` | unset (old IC) | `hopkins` selects Eq. 34 |
| `EUNHA_KEPLER_FLOOR` | 0.01 | floor, relative to the disk density 1 |
| `EUNHA_KEPLER_WIDTH` | `0.05·Rout` = 0.4 | outer edge width |

`ρ = floor + (R/Rin)^3` for `R < Rin`, `+1` for `Rin ≤ R ≤ Rout`, and `+[1 + (R−Rout)/w]^−3` for `R > Rout`. Our radii are exactly 4× Hopkins' (2 and 8 against 0.5 and 2). So `w = 0.1 × 4 = 0.4`, whether it is scaled by `Rin` or by `Rout`. The disk is flat (`Σ = 1`), as in the paper, not `1/R` as in our default. `P = 1e-6` and `v = R (R² + ε²)^−3/4` are unchanged. The old IC has no pressure-gradient term (`dP/dR = 0`), and neither does this one.

**1D** (`Hydro1DExam/laguerre_sod.c`). Same env names. 1D has no gravity, so the switch uses GIZMO's kinetic form (the one GIZMO computes and then discards): `coef · (E_th + ½ m max_nbr |Δv|²) > E_th`. The half-limit is the same rule as in 2D. **[CODE]**

### Tests on the box [RUN]

No 2D run. The full build needs MKL/FFTW-MPI (Intel), CAMB and `nvcc`, and none of them is here. What I did run:

**IC generator.** `review_tests/esw/ic_test.c` includes the `kepler` branch of `util.c` verbatim and evaluates it on the 256² grid of the 24-wide box (`ε = 0.047`). The default IC reproduces the logs: `E_pot = −17.5814`, `E_tot = −8.7917`, `E_int = 8.64e-4`. So the harness is the real code.

| | default | hopkins |
|---|---:|---:|
| gas mass | 75.97 | 209.71 |
| `E_pot` | −17.581 | −42.928 |
| `E_kin` | 8.789 | 21.458 |
| `\|Etot0\|` | 8.79 | 21.47 |
| mass at `R < 1.9` | 0.586 | 3.98 |
| ρ at `R = 0.5–1` | 1.0e-3 – 1.6e-3 | 0.027 – 0.13 |
| ρ at `R = 1.5–1.9` | 0.027 – 0.33 | 0.44 – 0.87 |
| ρ at `R = 2.1–7.9` | 0.18 – 0.80 (`2/R`) | 1.01 |
| ρ at `R = 9–12` | 1.0e-3 | 0.011 – 0.033 |
| max `\|v²/R − g\|/g` | 8e-16 | 8e-16 |

The velocities match the softened `kepler_accel` to round-off in both. Note that the Hopkins hole is not light. Its `(R/Rin)^3` ramp puts 7× more mass inside `R = 1.9`.

**Which cells the switch takes at t = 0.** Same harness, `ie = 1.32e-8` per cell:

| R | default, coef 0.01 | hopkins, 0.01 | default, 0.001 | hopkins, 0.001 |
|---|---:|---:|---:|---:|
| < 0.25 | all | all | all | all |
| 0.25–1.5 | all | all | none | all |
| 1.5–1.9 | all | all | 86% | all |
| 1.9–7.9 | all | all | all | all |
| 8.1–8.5 | none | 87% | — | — |
| > 8.5 | none | none | none | none |

With 0.01 the whole hole and the whole disk are switched. The switch leaves a cell once its `ie` exceeds `0.01 m|g|h`. For the default floor (`ρ = 1e-3`) at `R = 1`, the margin is only a factor of 1.6. So heated floor gas there drops out of the switch quickly. **[RUN, CALC]**

**Unit test of the 2D block.** `review_tests/esw/es_unit.c` includes the new `exam.c` block verbatim, on mock cells. All checks pass.
- 200 steps of compression and expansion (±50%), with the energy equation returning garbage: `K` returns to `K0` within `2e-16`.
- A switched hole cell with `ie_E = −8 ie_n` gets exactly the adiabatic value. `es_de = +9.0 ie_n`.
- The same drain in an unswitched cell (`R = 10`) gives `ie = 0.5 ie_n`, and `hl_de = 8.5 ie_n`.
- A heated edge cell (`ie_E = 10 ie_n`) is pulled back to adiabatic, with `es_de = −9 ie_n`. At `3000 ie_n` it leaves the switch.
- Both off: untouched.
- Half-limit alone, a cell asked for `−2 ie_n` every step: `ie` halves each step, and the total injected over 10 steps is `3.0 ie0`.

What this says **[CALC]**. For one deep event, the switch and the half-limit inject about as much as the `P = 1e-6` refill did (9 `ie_n` against 8–9). That energy went into kinetic energy through the face force before the `ie` update. No `ie` rule can take it back. What changes:
- The refill disappears (`floor_cum`, `n_ie_le0` should stay 0).
- Spurious heating from converging faces never reaches `P` in switched cells (`es_de < 0` there).
- A drained cell keeps its adiabatic `P`. With the half-limit, its `ie` shrinks geometrically, so its own share of injection is bounded while the drain is proportional to `ie`.

**[HYP]** The cliff is a feedback loop: created energy → hotter floor gas → larger face pressures → deeper drains. The switch cuts the first link for the cold gas. If the drain is driven by `ρ Δv²` at shearing faces and not by `P`, the switch will not stop it. Then the energy fix is the face-pressure budget limit (design in the previous reply), not an `ie` rule.

### 1D suite, switch on and off [RUN]

`N = 200`, Shu–Osher also at `N = 800`. Run from `/workspace/lagEunha/scratch_esw1d` with `GFS_PREFIX` set to `base`, `off`, `es`, `half`, `esonly`, or `es1e3`. No `gfs_*.dat` or `gfs_pair_*.dat` was written.

- **Off:** all eight `.dat` files and the summary line are byte-identical to the unmodified source (`cmp`).
- **Half-limit alone** (`GFS_HALF_LIMIT=1`): byte-identical on all eight, `nhalf = 0`. No problem in the suite loses half its `e` in one step.

| problem | L1(ρ) off | L1(ρ) switch 0.01 | ΔE/E0 off | ΔE/E0 switch 0.01 | cell-steps switched | `esde` |
|---|---:|---:|---:|---:|---:|---:|
| Sod | 5.858e-3 | 5.858e-3 | 1.26e-8 | 1.26e-8 | 0 | 0 |
| blast | 5.462e-2 | **0.2159** | 1.45e-5 | **2.75e-2** | 5294 | −7.56 |
| Shu–Osher N=200 | 0.9019 | 0.9019 | 3.11e-7 | 3.11e-7 | 0 | 0 |
| Shu–Osher N=800 | 0.5713 | 0.5713 | 2.75e-7 | 2.75e-7 | 0 | 0 |
| Noh | 1.363e-2 | 1.549e-2 | 2.00e-5 | 5.92e-4 | 2017 | −3.0e-4 |
| Lax | 4.793e-3 | 4.793e-3 | 2.05e-7 | 2.05e-7 | 0 | 0 |
| double rarefaction | 2.859e-2 | 2.859e-2 | 5.5e-10 | 5.5e-10 | 0 | 0 |
| collision | 0.1046 | 0.1046 | 3.72e-6 | 3.72e-6 | 0 | 0 |
| contact | 4.532e-3 | 4.532e-3 | 1.61e-7 | 1.61e-7 | 0 | 0 |

- `nhalf = 0` everywhere, also with the switch on. The result for switch-only is identical to switch plus half-limit.
- The kinetic switch fires where cold gas meets a shock (blast and Noh). There it replaces shock heating by adiabatic compression and destroys energy (`esde < 0`). The blast gets 4× worse. This is the failure GIZMO warns about, and the reason GIZMO dropped the kinetic form.
- With `GFS_ES_COEF=0.001`: the blast is byte-identical to off (`nes = 0`). Noh has `nes = 1757`, and ΔE/E0 goes from 2.0049e-5 to 2.0226e-5. The other six are unchanged.
- The 2D switch uses gravity, not `|Δv|`, so the 1D blast result does not carry over directly. But it is a warning. A disk cell hit by a real shock stays switched until its `ie` passes `0.01 m|g|h`, so its first shock-heating step is lost. `es_cum < 0` in the 2D log would show this.

### Compile [RUN]

`review_tests/syncheck.sh` (`mpicc -fsyntax-only -Wall -DXYZDBL -DGOTPM`):
- `exam.c`: 287 warnings with and without `-DUSE_CUDA`, 0 errors. That is the baseline.
- `exam_gpu_extract.c`: 0 and 16, unchanged.
- `Exam/KH/util.c`: 26 and 27, the same as before the change.
- `laguerre_sod.c` with `gcc -Wall`: 2 warnings, the same as before.

No link. `exam_gpu.cu` is not touched.

### Recommended next Kepler run

Keep 406515–406519 held.

1. **Primary, a GIZMO-style disk.** The binary at the tip of `master`, with 406734's settings (RK4 `way = 1`, `SEDOV_PHASE1 = 1`, `use_muscl = 1`, `av_mode = 5`, `EPS = 0.047`, `DTETA = 0.06`, `P = 1e-6`, `t_stop = 177.72`), plus:
   ```
   EUNHA_IC=kepler  EUNHA_KEPLER_PROFILE=hopkins
   GFS_ENTROPY_SWITCH=1          # half-limit follows (GFS_HALF_LIMIT=1)
   GFS_ES_COEF=0.01              # default; say so in the report
   GFS_FLOOR_LOG=1               # [RK4F] only fires now if ie <= 0
   ```
   Leave `GFS_DUAL_ENERGY` unset so the ledger has one mechanism. Check that `EUNHA_KEPLER_CX/CY = 12` is set. `kepler_accel` defaults to `(2, 2)`, while the IC defaults to the box centre. Check `GAS_FCENTROID` too (see above).
2. **Control, if there is a second slot.** The same, without `EUNHA_KEPLER_PROFILE`. That separates the switch from the IC.

The initial numbers for the Hopkins IC are `E_pot ≈ −42.93`, `|Etot0| ≈ 21.47`, `E_int(0) = 8.64e-4`. Expect `[ESW] entropy_switch=1 coef=0.01 half_limit=1 active=1` once.

**Columns to watch in `[RK4E]`:**
- `dEtot/|Etot0|`, against `(floor_cum + es_cum + hl_cum)/|Etot0|`. They should agree.
- `es_cum`, with its sign. Positive: drained cells are being topped up. Negative: spurious heating is being removed. A large negative value means the disk is being shocked for real.
- `hl_cum` and `n_half`. How much the half-limit adds.
- `floor_cum` and `n_ie_le0`. These should stay 0.
- `n_es`. At t = 0 the switch should take about 25500 cells (every cell inside `R ≈ 8.1` and most of `8.1–8.5`), which is 39% of the 65536 cells. A steady fall means cells are heating out of the switch.
- `E_int`. It should stay within a factor of a few of `8.6e-4`. In 407021 it was 12× by `t = 8`.
- `n_pair` and `dt`. A falling `dt` announced the cliff last time.

**Pass criterion:** `dEtot/|Etot0|` stays on the `10^{-4}` slope (about `1e-4` per unit time, as in 406734 before `t = 9`), with no jump, for several inner orbits: `t > 3 × 17.8 ≈ 53`. The target is 10 inner orbits, `t = 177.7 = t_stop`. Please report the time of the first `|ΔE| = 10^{-3}` and `10^{-2}` crossings, as in §14. Only after a pass should anything go to 406515–406519.

If it still falls off the cliff near `t ≈ 10`, with `es_cum` carrying the jump, then an `ie` rule is not enough. The next step is the face-pressure budget limit (previous reply), which needs `nvcc` to check the GPU kernel.

### Commits

- `aa2d091` RK4 blend: `GFS_ENTROPY_SWITCH`, `GFS_ES_COEF`, `GFS_HALF_LIMIT`, and the `[RK4E]` columns.
- `ee0ea0e` Kepler IC: `EUNHA_KEPLER_PROFILE=hopkins`, `EUNHA_KEPLER_FLOOR`, `EUNHA_KEPLER_WIDTH`.
- `d87f96d` 1D: the same switches, with the kinetic criterion.

The test sources (`ic_test.c`, `es_unit.c`, and their outputs) are in `/workspace/lagEunha/review_tests/esw/` on the box. They are outside the repository.

## Grok CLI께: Laguerre 면 회전 보정 `GFS_LAGUERRE_ROTATION` (기본 꺼짐)

태그는 위와 같습니다. **[CODE]** 소스에서 확인, **[RUN]** 박스에서 실행, **[CALC]** 유도.

### 요약

1. 면 회전 보정(`voronoi_face_rotation`)을 Laguerre 면에도 적용하는 옵트인 플래그를 추가했습니다. 기본값은 0이고, 이때는 아무것도 바뀌지 않습니다. **[CODE]**
2. `w2_i = w2_j = 0`인 면(순수 Voronoi)은 플래그와 상관없이 기존과 비트 단위로 같은 경로를 탑니다. 기존 `0.5` 식을 그대로 두었습니다. **[CODE, RUN]**
3. 따라서 Laguerre 가중치를 쓰지 않는 런(`GAS_Kappa = 0`이면 `w2 = 0`)에는 영향이 없습니다. 가중치가 0이 아닌 런에서만 의미가 있습니다.

### 식 [CALC]

Laguerre 면도 `e = (x_j − x_i)/d`에 수직이고, `de/dt = u_ij,⊥/d`입니다. `get2dUpqradRk4`가 돌려주는 `fact1·u_ij + fact2·e`는 면이 선분 ij를 지나는 앵커 점 `a = x_i + fact1 (x_j − x_i)`의 속도(와 가중치 변화율 항)입니다. `fact1 = (1 + (w2_i − w2_j)/d²)/2`. 그러므로 면 위의 점 `x`에서

```
u_n(x) = u_n(a) − u_ij·(x − a)/d
```

Voronoi 식과 모양이 같고, 기준점만 중점 `m`에서 앵커 `a`로 바뀝니다. `w = 0`이면 `fact1 = 1/2`, `a = m`이 되어 기존 식과 정확히 같습니다. 코드에서는 면 중심 `c_f`(i 기준 국소 좌표, `Voro2D_FindVC`가 코너를 `Vec2DSub`로 저장)에 대해 `offset = fact1·(x_j − x_i)`를 쓰고, `w2`는 `get2dUpqradRk4`와 같은 값(`p->w2`, `q->w2`, GPU에서는 `pw2[]`)을 씁니다.

### 적용 위치 [CODE]

| 파일 | 위치 | 비고 |
|---|---|---|
| `Exam/exam.c` | `voronoi_face_rotation`, 호출 두 곳 (av_mode 5 HLLC 운동량 경로, 에너지/점성 경로) | 플래그는 OpenMP 영역 밖에서 한 번 읽어 인자로 넘김 |
| `Exam/exam_gpu.cu` | `getAccVoro2DBlend_kernel`의 면 속도 보정 | 커널 인자 `laguerre_rot` 추가, `phase1`과 같은 방식 |
| `Exam/exam_gpu.h` | `GPUPhysicsParams.laguerre_rot` | |
| `Exam/exam_gpu_extract.c` | `getAccVoro2DBlend_GPU`, `..._validate`에서 env 읽기; `cpu_reference_force_csr` 두 곳 | CPU 참조 루프도 GPU와 같게 |

운동량과 에너지가 같은 면 속도를 쓰도록 모든 사본에 똑같이 넣었습니다. 켜지면 `exam2d_vph_rk4_int_blend`에서 rank 0이 `[LAGROT] GFS_LAGUERRE_ROTATION=1: ...`를 한 번 찍습니다.

| env | 기본값 | 의미 |
|---|---|---|
| `GFS_LAGUERRE_ROTATION` | 0 | 1이면 Laguerre 면에도 앵커 기준 회전 보정 |

### 검증 [RUN]

- `review_tests/syncheck.sh`: `exam.c` 287/287 경고, 0 오류. `exam_gpu_extract.c` 0/16, 0 오류. 경고 목록은 줄 번호 이동을 빼면 변경 전과 같습니다. `exam_gpu.cu`는 `nvcc`가 없어 컴파일하지 못했습니다. 눈으로 확인만 했습니다.
- 유한차분 검사 `review_tests/lagrot/lagrot_fd.c`: `get2dUpqradRk4`와 새 `voronoi_face_rotation`을 소스에서 그대로 추출해 씁니다. 무작위 p, q, 속도, 가중치 10⁵개에서, 앵커에서 벗어난 면 위 점의 법선 속도를 면 방정식의 중심차분과 비교했습니다.
  - 가중치 고정: 최대 오차 `1.6e-9` (`|u_ij|₁`로 정규화), `|V_fd|` 대비 `5.3e-6` (V_fd ≈ 0 근처의 FD 잡음).
  - `sqrt(w2)`가 선형으로 변할 때(`fact2`의 변화율 항 포함): `1.4e-9`, `|V_fd|` 대비 `1.8e-4`.
  - 비교: 보정 없음은 `|V_fd|` 대비 최대 `10⁵` 이상, 중점을 기준점으로 쓰면 정규화 오차 `0.19`. 앵커가 맞는 기준점입니다.
  - `w = 0`에서 플래그 0과 1의 결과 차이는 정확히 0.
- 1D 스위트(`Hydro1DExam/laguerre_sod.c`)는 `gfs_pair.h`만 포함하고 이 코드를 거치지 않으므로 바뀌지 않습니다.

### 주의할 점 [CODE]

- 메시 생성은 `getwfrac`으로 앵커 비율을 `[0.05, 0.95]`로 자릅니다. `get2dUpqradRk4`의 `fact1`은 자르지 않습니다. `|w2_i − w2_j| > 0.9 d²`인 면에서는 기저 면 속도 자체가 실제 면 위치와 어긋나고, 새 보정은 `get2dUpqradRk4`에 맞췄습니다 (지시대로). 가중치가 이 범위를 넘는 런이라면 알려 주세요.
- GPU는 `w2`를 float로 가지고 있어 `fact1`도 float `w2`에서 계산합니다 (`dev_get2dUpqradRk4`와 같음).

### 가중치를 쓰는 런에서 권하는 시험

`GAS_Kappa > 0`인 같은 설정으로 플래그만 바꿔 두 번 돌려 주세요.

```
GFS_LAGUERRE_ROTATION=0   # 기준
GFS_LAGUERRE_ROTATION=1   # [LAGROT] 줄이 한 번 찍혀야 함
```

짧은 KH 또는 Kepler 몇 궤도면 충분합니다. 보실 것: `dEtot/|Etot0|` 기울기, `[RK4E]`의 `floor_cum`과 `n_ie_le0`, `dt`, 그리고 전단층의 셀 노이즈. 순수 Voronoi 런(`GAS_Kappa = 0`)에서 두 결과가 비트 단위로 같은지도 한 번 확인해 주시면 좋습니다. GPU 런이라면 `getAccVoro2DBlend_GPU_validate`로 CPU 참조와 GPU 커널이 플래그 켠 상태에서 일치하는지 보는 것이 가장 빠른 커널 검사입니다.

### 커밋

- (이 커밋) `GFS_LAGUERRE_ROTATION`: exam.c, exam_gpu.cu, exam_gpu.h, exam_gpu_extract.c, 이 절.

시험 소스(`lagrot_fd.c`, `extracted.h`)는 박스의 `/workspace/lagEunha/review_tests/lagrot/`에 있고, 저장소 밖입니다.

## Grok CLI께: 3D 시험 스위트 `Exam/Tests3D` (3D 드라이버가 생기기 전에는 제출하지 마세요)

요약:
1. **저장소에는 아직 3D GFS 유체 경로가 없습니다** (아래 조사는 저장소 안의 코드만 다룹니다). **[CODE]** 주한님 말씀으로는 예전 3D 드라이버 **Exam3d**가 클러스터에 있고 아직 저장소에 올라오지 않았습니다. **Grok CLI께 부탁드립니다: Exam3d를 저장소에 추가해 push해 주세요.** 그러면 Grokbot이 2D GFS 경로와 비교해 감사하고, Tests3D를 거기에 연결하겠습니다(preflight 표식 확인, IC/에너지 로그 어댑터).
2. 아직 아무것도 이식하지 않았습니다. Exam3d가 올라오면 그것을 출발점으로 삼습니다. 2D 코드는 한 줄도 바뀌지 않았습니다. 이번 커밋은 `Exam/Tests3D/` 아래 새 파일과 이 절뿐입니다.
3. 대신 다섯 시험(음파, Sedov, Noh, Evrard, KH)의 다음을 모두 준비했습니다: IC 생성기, 파라미터 템플릿, `run.slurm`, 해석 스크립트, 합격 기준, 박스 자체 시험(`selftest.sh`).
4. **모든 `run.slurm`은 바이너리에 표식 문자열 `LAGEUNHA_3D_GFS_V1`이 없으면 preflight에서 exit 3으로 멈춥니다.** 지금의 `eunha2.exe`로 제출하면 할당만 받고 곧바로 끝납니다. 3D 드라이버가 들어간 커밋이 생기면 그때 이 절의 명령을 쓰세요.

### 3D 경로 조사 [CODE]

| 항목 | 3D 상태 |
|---|---|
| `ex3d_*`, `nearest3dOpen`, `det3d_dpq(RK4)` (exam.c ~353–1007) | 3D 트리로 탐색 반경(`w2ceil`)만 계산. **어느 드라이버도 호출하지 않음** |
| `Voro/voro.c` `Voro3D_FindVC`, `Voro3D_FaceExtract`, 부피/중심 | 3D Voronoi 기하는 있음. 저장소 안에서는 옛 `Exam/Sedov`만 사용 (Exam3d도 쓸 것으로 추정) |
| `Voro/Laguerre/` 3D (CPU + CUDA `construct_cells_3d_kernel`) | 주 빌드에 없음 (`Voro/Makefile`은 voro.o, voro_eunha.o만 빌드). CUDA 커널은 부피와 꼭짓점 평균 중심만 주고 면 목록이 없음 |
| `Exam/Sedov` (3D Voronoi 유체) | 독립형 OpenMP, MPI 없음. Monaghan AV를 쓰고 HLLC/MUSCL/pair/RK4 없음. 어떤 Makefile에도 없고, 컴파일 오류 12개 |
| HLLC+MUSCL, extreme face (P_max > 100 P_min), pair pressure + `gfs_pair_work_limit`, `voronoi_face_rotation` | **2D 전용** (`getAccVoro2DBlend_impl`, `hllc_face_2d`) |
| 밀도/부피/floor (`updateDenW2Pressure2DBlend`) | **2D 전용** |
| RK4 드라이버, `kepler_accel`, `phase1_ie_stage`, `[RK4E]`, GFS_DUAL_ENERGY, GFS_ENTROPY_SWITCH | **2D 전용** (`exam2d_vph_rk4_int_blend`) |
| rk4 입자의 영역 분할 | `startRkSDD2D`, `MakeDoDeInfo2D`만 있음 |
| `eunha2.c` SIMMODEL | KH/RT/RT_LF/MkGlass2D/Cylinder/Sedov2D만 있음. 3D 모델, IC 읽기, 출력 없음 |
| GPU | `exam_gpu.cu`의 커널은 모두 2D. **3D는 드라이버가 생겨도 당분간 CPU 전용** |
| 자체 중력 | GOTPM(TreePM)은 VORO 입자 질량을 넣지만 주기적 우주론 전용이고 `RunCosmos`는 GFS를 부르지 않음. 2D 드라이버에는 외부 `kepler_accel`과 `GAS_ACC`뿐. **저장소 기준으로 Evrard는 3D 유체와 고립계 자체 중력이 모두 없어 두 겹으로 막혀 있음** (Exam3d의 중력 여부는 감사 후 확인) |

이식 계획, 기대하는 드라이버 인터페이스(LAG3DV1 IC/스냅숏 `snap_%06d.l3d`, `[E3D] step= t= dt= Ekin= Eint= Epot= Etot=` 로그 줄, 제안하는 params 키), 열린 결정은 `Exam/Tests3D/README.md` §1–2, §5에 있습니다.

### 공통 설정

- **빌드**: 3D 드라이버가 들어간 커밋에서 평소처럼 빌드합니다. GPU 플래그는 필요 없고, 3D는 CPU 전용입니다. 제출 전에 표식을 확인하세요.
  ```
  grep -a -c LAGEUNHA_3D_GFS_V1 eunha2.exe   # 1 이상이어야 함
  ```
- **env**: `common/flags_base.env`를 각 `config.env`가 source합니다. 우리가 정한 2D 관례 그대로입니다.
  - `SEDOV_PHASE1=1`, `GFS_FLOOR_LOG=1` (음파와 KH는 0);
  - `GFS_DUAL_ENERGY`, `GFS_ENTROPY_SWITCH`, `GFS_HALF_LIMIT`, `GFS_LAGUERRE_ROTATION`, `EUNHA_KEPLER_*`는 unset;
  - `OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK`.
- **params**: `params3d.dat`. RK4 (`way = 1`), `av_mode 5`, `use_muscl 1`, `kappa 0`, `gpu_enabled 0`, centroid shift 0. pair pressure는 capped law와 work limit를 씁니다. 3D의 `nshare`는 열린 결정입니다.
- **스크립트 공통**: `sbatch`는 시험 디렉터리에서 실행합니다. `LAGEUNHA_SRC`(기본 `$HOME/LagEunha`)의 `eunha2.exe`를 `/gpfs/kjhan/LagForce/<시험>_gfs_…`에 복사하고, `COMMIT`/`ENV`를 기록하고, IC를 만들고, `mpirun` 출력은 `log`에, 마지막 줄 `EXIT:<rc>`를 씁니다. **2D 스크립트와 달리 스크립트가 mpirun의 rc로 끝나므로, 죽은 런이 COMPLETED 0:0으로 보이지 않습니다.**
- **IC 생성**: `LAG3D_PY`(기본 `python3`)에 numpy와 scipy가 있어야 합니다.
- **자원**: 기본 4 rank × 16 스레드이고, 파티션은 비워 두었습니다(`##SBATCH -p`). 클러스터에 맞게 정해 주세요.
- **시간**: 아래 시간은 모두 **추정**입니다. 2D GPU Gresho 처리량(약 2.5×10⁵ 입자·스텝/s)에서 3D CPU를 rank당 약 10⁵으로 가정했고, 3배 정도 틀릴 수 있습니다.
- **공통으로 보고해 주실 것**: 런 디렉터리, `COMMIT`, `ENV`, `EXIT:` 줄, `sacct` 상태, 스텝 수와 wall time, `[RK4F]`/floor 줄 개수, 해석 스크립트의 `*_summary.json` 전문과 PNG.

### 1. SoundWave3D (선형 음파, 수렴과 Galilean 불변성)

| | |
|---|---|
| 설정 | 주기 단위 상자, ρ0 = 1, c0 = 1, A = 10⁻⁶, x 방향, 한 주기 t = 1. cubic 격자, 질량 m = ρ(x) dx³ (2D KH/Gresho IC와 같은 방식) |
| 실행 | `cd Exam/Tests3D/SoundWave3D`<br>`sbatch --export=ALL,N=32 run.slurm`<br>`sbatch --export=ALL,N=64 run.slurm`<br>`sbatch --export=ALL,N=128 -t 12:00:00 run.slurm`<br>`sbatch --export=ALL,N=64,BOOST=1 run.slurm` |
| 해석 | `python3 analyze.py --snaps <N32>/<마지막 snap> <N64>/… <N128>/… --boosted <N64_boost1>/… -o sw3d` |
| 추정 시간 | 32³ 수 분, 64³ 약 3분, 128³ 약 40분 (약 110/215/430 스텝) |
| 합격 | P1: 가장 고운 두 해상도 사이 L1(ρ) 수렴 차수 ≥ 1.7<br>P2: 64³에서 L1/A ≤ 10⁻²<br>P3: L1_boost/L1_rest ≤ 1.25 |
| 보고 | 해상도별 L1, 수렴 차수, boost 비, `sw3d_convergence.png`, `sw3d_profile.png` |

### 2. Sedov3D (γ = 5/3)

| | |
|---|---|
| 설정 | 주기 단위 상자, ρ = 1, P_amb = 10⁻⁵. E = 1을 중심 가까운 64개에 top hat으로 넣음 (cubic에서는 동률 때문에 88개가 되고, 그대로 기록). t_end = 0.05, R_an = 0.3475 (ξ0 = 1.15167) |
| 격자 | 기본 cubic. glass는 `--export=ALL,N=64,LATTICE=glass`. **어느 쪽인지 꼭 보고** |
| 실행 | `cd Exam/Tests3D/Sedov3D`<br>`sbatch --export=ALL,N=64 run.slurm`<br>그다음 `sbatch --export=ALL,N=128 -t 12:00:00 run.slurm` |
| 해석 | `python3 analyze.py --snap <N64>/<마지막 snap> --ic <N64>/ic.l3d [--snap2 <N128>/<마지막 snap>] -o sedov3d` |
| 추정 시간 | 64³ 약 20분 (약 2000 스텝), 128³ 약 6시간 |
| 합격 | 스윕 질량 충격 반경 \|R/R_an − 1\| ≤ 0.03 (128³에서 0.02)<br>\|ΔE/E0\| ≤ 10⁻³<br>ρ_peak ≥ 2.5 (128³에서 3.0)<br>축/대각 비대칭 5% 이내<br>L1(128³) < L1(64³) |
| 보고 | summary.json, `sedov3d_profiles.png`, floor에 닿은 셀 수 |

### 3. Evrard3D (자체 중력 필요, 현재 막힘)

| | |
|---|---|
| 설정 | Evrard (1988), Springel (2010), Hopkins (2015)와 같은 표준 설정. γ = 5/3, G = M = R = 1, ρ = 1/(2πr), u = 0.05, v = 0. E_pot(0) = −2/3, E_tot = −0.6167. 격자를 r < 1로 자르고 r → r^{3/2}로 늘림 (같은 질량). Plummer ε = 0.01 |
| 경계 | 기본 진공. `BACKGROUND=uniform`은 저밀도 배경 가스 (열린 결정) |
| 필요한 것 | `Hydro3D gravity = 1`, 즉 유체 입자와 결합한 고립계 자체 중력. **아직 없음** |
| 실행 | `cd Exam/Tests3D/Evrard3D`<br>`sbatch --export=ALL,N=64 run.slurm` (약 13.7만 입자, t_end = 3) |
| 해석 | `python3 analyze.py --log log --snap <t≈0.8 스냅숏> -o evrard3d`. `[E3D]` 줄이 로그에 있어야 함 |
| 기준 | `ref/evrard1d_N2000_*.txt`: 여기서 만든 1D 구대칭 라그랑주 풀이 (`evrard1d_ref.c`, 에너지 오차 1.5×10⁻⁵). t = 0.8 프로파일은 SWIFT의 HydroCode1D 파일과 log ρ 기준 0.3%로 일치. 논문 곡선을 디지털화한 것이 아님. 기준값: E_kin 최대 0.450 (t = 0.88), E_th 최대 1.758 (t = 1.07), E_pot 최소 −2.543 (t = 1.05) |
| 추정 시간 | n = 64에서 유체만 약 2시간 (약 2×10⁴ 스텝)에 중력 비용 추가 |
| 합격 | E1 (주 기준): max \|ΔE_tot\|/\|E0\| ≤ 1%<br>E2: E_kin 최대 시각 0.88 ± 0.10, E_th 최대 시각 1.07 ± 0.15<br>E3: E_th 곡선의 평균 상대 차이 ≤ 0.15<br>E4: t = 0.8에서 0.05 < r < 0.8 평균 \|log10 ρ/ρ_ref\| ≤ 0.10<br>E3와 E4는 우리가 정한 기준이고 공동체 표준이 아님 |
| 보고 | `evrard3d_energy.png`, summary.json, 중력 해법과 ε |

### 4. Noh3D (구형, 충격 후 밀도 64)

| | |
|---|---|
| 설정 | 주기 상자 L = 6, ρ = 1, P = 10⁻⁶, v = −r̂, t_end = 2. 정답: R = 2/3, ρ_post = 64, P_post = 64/3, 충격 앞 ρ = (1 + t/r)². 주기 경계의 교란이 들어오지 않은 r < 1 안에서만 비교 |
| 실행 | `cd Exam/Tests3D/Noh3D`<br>`sbatch --export=ALL,N=128 run.slurm` |
| 해석 | `python3 analyze.py --snap <마지막 snap> --ic ic.l3d -o noh3d` |
| 추정 시간 | 128³ 약 1.5시간 (1–5시간), 약 1000 스텝 |
| 합격 | \|R/R_an − 1\| ≤ 0.05<br>0.3–0.9 R 평균 밀도가 64의 20% 이내<br>충격 앞 L1 ≤ 5%<br>\|ΔE/E0\| ≤ 10⁻³<br>중심(r < 0.3 R)의 wall heating은 보고만 하고 채점하지 않음 |
| 보고 | summary.json, `noh3d_profiles.png`, extreme face/pair pressure 관련 로그가 있으면 그 개수 |

### 5. KH3D (McNally 2012 층을 z로 확장)

| | |
|---|---|
| 설정 | 1 × 1 × Lz (기본 0.25), P = 2.5, ρ = 1/2, U = ±0.5, L = 0.025. 씨앗 v_y = 0.01 sin(4πx). 격자, 질량 = ρ dV |
| 선형 이론 | 이 층의 압축성 고유값 풀이로 σ = 2.83. 날카로운 경계의 5.92는 비교 기준이 아님 |
| 실행 | `cd Exam/Tests3D/KH3D`<br>`sbatch --export=ALL,N=128 run.slurm` (128×128×32)<br>가능하면 `N=256`도 |
| 해석 | `python3 analyze.py --snaps '<run>/snap_*.l3d' --ic <run>/ic.l3d [--snaps2 '<run256>/snap_*.l3d'] -o kh3d` |
| 추정 시간 | 128×128×32로 t = 2까지 약 50분 (약 2200 스텝) |
| 합격 | K1: n ≥ 128에서 \|σ_fit/σ_lin − 1\| ≤ 0.15 (1.5 M0 ≤ M ≤ 0.06 구간, 4점 이상)<br>K2: 해상도를 올려도 오차가 커지지 않음 |
| 보고 | σ_fit, 맞춤 구간, `kh3d_growth.png` |

### 박스에서 확인한 것 [RUN]

- `Exam/Tests3D/selftest.sh`가 모두 통과했습니다(깨끗한 clone에서도). 검사 내용:
  - C 도구 빌드 (경고 0);
  - 해석해: Sedov ξ0 = 1.15167, Noh 64 / 64/3 / R = 2/3, KH 풀이가 Michalke 0.1897과 일치;
  - 다섯 IC의 질량, 운동량, 에너지 합과 C 읽기 교차 확인 (바이트 단위 왕복);
  - 해석 스크립트를 합성 스냅숏에 적용;
  - `run.slurm`: 표식 없는 바이너리에서 exit 3. 가짜 표식 바이너리와 stub `mpirun`으로 런 디렉터리, IC, 채워진 params, `EXIT:0` 생성. 실제로 제출한 것은 없습니다.
- `review_tests/syncheck.sh`: `exam.c` 287/287 경고, 0 오류. `exam_gpu_extract.c` 0/16, 0 오류. 변경 없음.
- 확인하지 못한 것: LagEunha로 돌리는 모든 부분. 합격 기준은 해석해와 2D 거동에 맞춘 것이고, 3D 런으로 보정한 것이 아닙니다.

### 커밋

- `45c2d67` Tests3D common: LAG3DV1 형식, 격자/glass, Sedov/Noh/KH 해석해, 실행 도우미
- `0ea330b` 다섯 시험의 IC, config, run.slurm, 해석 (Evrard 1D 기준 포함)
- `8573f55` README (3D 조사, 드라이버 인터페이스, 합격 기준)와 selftest.sh
- (이 커밋) 이 절

## Grok CLI께: 3D GFS 경로 `Exam/exam3d_gfs.c` (CPU, 실험 단계)

Juhan 확인: "Exam3d" = `Exam/Sedov`(옛 3D 시제품, 3D에서 실패 후 2D로 내려감). 이번에 주 코드에 새 3D GFS 경로를 넣었습니다. `Exam/Sedov`는 면 루프 틀로만 참고했고 물리는 하나도 가져오지 않았습니다 (ie ≥ 0 clamp 없음, 고정 입자 없음, dt는 합산된 전체 가속도로, RK4). **2D 경로는 비트 단위로 그대로입니다.** GPU 3D는 다음 단계입니다.

### 추가/변경한 파일

- 새 파일: `Exam/exam3d_gfs.c` (약 1500줄), `Exam/exam3d_gfs.h`, `Exam/exam3d_main.c` (독립 실행 `lag3d.exe`), `Exam/gfs_riemann.h` (exam.c의 `hll_star_state`, `hllc_face_2d`, `hllc_face_2d_rest_frame`를 `gfs_*` 이름으로 복사, phase1은 인자).
- `eunha2.c`: `Simulation Model = Hydro3D`이면 `checkarg` 직후 3D 경로로 분기하고 `MPI_Finalize` 후 종료. 다른 모델은 기존 흐름 그대로.
- `Exam/Makefile`: `EXOBJ += exam3d_gfs.o` (libexam.a), `lag3d.exe` 규칙.
- `Exam/Tests3D`: `run3d_common.sh`, 다섯 `run.slurm`, `params3d.dat.template`, `flags_base.env`, `README.md` (§1.3, 새 §1.5, "Exam3d push" 문구 삭제).
- `exam.c`, `gfs_pair.h`, `Voro/`는 건드리지 않았습니다.

### 스위치 (선택 방식)

- 런타임: params의 `Simulation Model = Hydro3D`. 표식 `LAGEUNHA_3D_GFS_V1`을 출력합니다 (`run.slurm`이 이를 검사).
- 물리는 2D의 av_mode 5를 차원만 바꿔 옮겼습니다: 변 길이 → 면 넓이, 변 중점 → 다각형 넓이 중심, sqrt(V) → cbrt(V), 2×2 → 3×3 Green-Gauss 기울기 (Barth-Jespersen). HLLC 정지좌표계 + MUSCL, extreme face (`SEDOV_PHASE1`, 압력비 > 100), Springel 면 회전 (w = 0은 중점, `GFS_LAGUERRE_ROTATION=1`이면 Laguerre anchor), 상한 쌍압력 + `gfs_pair_work_limit`, RK4 (자체중력과 외력 hook 포함), 같은 에너지 변수와 장부.
- `Voro3D_FindVC`에서 중심과 모든 이웃의 w2를 명시적으로 넣습니다 (0 = 순수 Voronoi).
- 환경 변수:
  - 2D와 같은 의미: `SEDOV_PHASE1`, `GFS_FLOOR_LOG` (`[RK4F]`, `[RK4SF]` 줄), `GFS_DUAL_ENERGY`, `GFS_ENTROPY_SWITCH`/`GFS_HALF_LIMIT`/`GFS_ES_COEF`, `GFS_LAGUERRE_ROTATION`, `HYDRO_TSTOP`, `EUNHA_DUMP_DT`;
  - 3D 전용:
    - `GFS3D_PAIR_NSHARE` (기본 16);
    - `GFS3D_STAGE_FLOOR` (0 = 단계마다 P만 바닥값, 2D GPU 경로와 같음 / 1 = 2D CPU의 ie 재설정, `sfl_cum`에 기록);
    - `LAG3D_BC` (periodic|reflect|outflow, 축별 가능);
    - `LAG3D_THETA`, `LAG3D_GRAV_DIRECT`, `LAG3D_MAXSTEPS`;
    - `LAG3D_PM_GM/X/Y/Z/EPS`, `LAG3D_ACC`.
- params 제약: av_mode = 5, entropy_mode = 0, kappa ≤ 0 (kappa > 0은 아직 없음).

### 설계 결정

- **3D work-limit 면 수 (nshare = 16):** 2D의 6은 평면 Voronoi의 평균 면 수(오일러 공식)입니다. 3D 평균 면 수는 Poisson 15.54, 정돈된 glass 약 14.5, bcc 14, fcc 12이고, 입방 격자에서 퇴화하지 않은 면은 6개입니다. 16은 이 값들보다 모두 크거나 같아서 전형적인 셀에서 면 예산의 합이 ie를 넘지 않습니다. 32³ Sedov에서 쌍압력을 탄성으로 바꾸면(nshare = 1e30) dE는 3.2e-3이었습니다 (기본값에서는 1.03e-2, 대부분 바닥값).
- **경계:**
  - periodic은 link-cell 영상;
  - reflect는 거울 영상 (셀이 벽면에서 정확히 잘리고, 벽 면은 반사 HLLC 압력만 받고 에너지 플럭스는 0), 벽을 넘은 입자는 스텝 끝에 반사;
  - outflow는 속도를 복사한 거울 영상, 스텝 끝에 상자 밖 입자를 제거하고 `out_cum`에 기록;
  - 기본값은 IC의 periodic 플래그를 따르고 (아니면 reflect), params의 `Hydro3D boundary` 또는 `LAG3D_BC`로 바꿉니다.
- **중력:** CPU에서 Plummer 연화를 쓴 직접 합(N ≤ 20000이면 기본)과 Barnes-Hut 트리(θ = 0.5, monopole)를 구현했습니다. 고립계만 지원합니다 (periodic 축과 함께 쓰면 오류). Epot = ½ Σ m φ.
- **병렬화:** rank 0만 계산하고 OpenMP를 씁니다. 그래서 `run.slurm`은 `--ntasks=1 --cpus-per-task=32`, `NRANK=1`, 바인딩 끔입니다 (바인딩을 켜면 박스에서 12.2 s, 끄면 1.8 s). MPI 영역 분할은 아직 없습니다.
- 위상 실패 시 1e-9~1e-7 셀 크기의 jitter로 재시도합니다 (`njit`로 기록, 지금까지 시험에서 0).

### 빌드와 실행

- 클러스터 (기존 Makefile 흐름): 평소처럼 `make` 하면 `exam3d_gfs.o`가 libexam.a에 들어가고 eunha2가 `Hydro3D`를 처리합니다. 독립 실행 파일은 `cd Exam && make lag3d.exe` (MPI/FFTW 불필요; gcc 박스에서는 `make lag3d.exe CC=gcc OPT="-O2 -fopenmp"`).
- 실행 (Slurm 제출은 Grok CLI/Juhan이 합니다. 저는 제출하거나 취소하지 않았습니다):
  - `cd Exam/Tests3D/SoundWave3D && sbatch --export=ALL,N=32 run.slurm` (N = 32, 64, 128, boost 런 포함);
  - `cd Exam/Tests3D/Sedov3D && sbatch --export=ALL,N=64 run.slurm` (그다음 128);
  - `Evrard3D`, `Noh3D`, `KH3D`도 같은 방식입니다 (각 README 절 참고);
  - 바이너리 선택: `LAG3D_BIN=<경로>/eunha2` (기본) 또는 `LAG3D_BIN=<경로>/Exam/lag3d.exe`, 실행기 `LAG3D_LAUNCH=mpirun|direct`, 추가 옵션 `LAG3D_MPIRUN_OPTS`.
- 로그: stdout에 `[E3D] step= t= dt= Ekin= Eint= Epot= Etot= dE_rel= floor_cum= sfl_cum= de_cum= es_cum= hl_cum= out_cum= npair= N= Rmax= njit=`, stderr에 `[RK4E]`. 스냅숏은 `snap_%06d.l3d` (LAG3DV1, VOL 포함, 중력이 있으면 POT 포함)로 t0, dump 간격마다, t_end에 씁니다. 각 `analyze.py`가 그대로 읽습니다.

### 박스에서 확인한 것 (측정값만; 8 스레드, gcc 14, `SEDOV_PHASE1=1`, C = 0.3, 입방 격자)

- **2D 불변:**
  - `syncheck.sh`: exam.c 287/287 경고 0 오류, exam_gpu_extract.c 0/16 0 오류. 전과 같습니다.
  - exam.o md5 `a7d63a46…` 전후 동일.
  - es_unit, ic_test, lagrot, `hll_shear_face.py`, `shear_pair_check.py`, `voronoi_face_velocity_check.py` 출력이 바이트 단위로 같습니다. `laguerre_sod_diag`는 실행 시간(ms) 표기만 다릅니다 (저장소 코드를 쓰지 않는 독립 프로그램).
  - `exam3d_gfs.c`는 `gcc -Wall -Wextra`에서 경고 0.
- `Tests3D/selftest.sh`: ALL SELFTESTS OK.
- **SoundWave3D (A = 1e-6, t = 1):**
  - 16³/32³/64³에서 L1(ρ)/A = 3.99e-2 / 8.64e-3 / 2.12e-3, 차수 2.21 / 2.02 (맞춤 2.12);
  - 진폭비 0.9665 / 0.9962 / 0.9996, |dE/E| ≤ 1.2e-13;
  - analyze PASS;
  - 32³ boost (v = 1): L1_boost/L1_rest = 0.9999998, |dE/E| = 1.5e-14.
- **Sedov3D (t = 0.05, ξ0 = 1.15167 → R_an = 0.3475):**
  - 16³: dE/E0 = 9.8e-5 (바닥값 없음), R_peak 0.357, R_meas 0.418, ρ_peak 1.40;
  - 32³: dE/E0 = 1.03e-2 (그중 floor_cum 1.01e-2, 장부는 닫힘), R_peak 0.331, R_meas 0.404, ρ_peak 1.77, 비등방성 0.985, 368 스텝, 150 s;
  - 32³ C = 0.15: dE 2.6e-3, ρ_peak 1.64;
  - 32³ glass: dE 1.2e-3, R_meas 0.394, ρ_peak 1.79;
  - 16³ glass: dE 4.1e-5;
  - 32³ `GFS3D_STAGE_FLOOR=1`: dE 2.2e-2 (sfl 2.19e-2).
  - 원인: 충격파 바로 뒤에서 ie ≤ 0이 되어 스텝 끝 바닥값이 들어갑니다. 입방 격자에서는 뜨거운 입자 약 90개가 격자 방향으로 약 3셀 앞서 나가서 R_meas가 커집니다 (glass에서는 3개).
- **Evrard3D:**
  - n = 16, 균일 배경 (5968 입자): 직접 합 Epot(0) = −0.6594462 (IC 값과 같음), 트리 θ = 0.5는 −0.659604;
  - t = 0.8에서 dE 2.4e-3 (바닥값 1.4e-3, R ≈ 1.14의 배경 셀), 596 스텝, 292 s. analyze E1, E3 참, E2, E4 거짓 (t = 0.8까지만 돌림);
  - 진공 + 반사 상자: 40 스텝 (t = 0.41), dE 6.5e-6.
- **Noh3D, KH3D 16³:** 30 스텝 스모크 이상 없음 (Noh dE 2.9e-6).
- OpenMPI mpirun으로 `run.slurm`을 끝까지 돌렸습니다 (SoundWave3D, Sedov3D N = 16, EXIT:0, Slurm 없이 직접 실행).

### 클러스터에서 해야 할 것

1. icx/MKL/CUDA 환경에서 전체 eunha2 링크 확인 (박스에는 nvcc/MKL이 없어서 3D 경로는 `lag3d.exe`로만 빌드했습니다. eunha2.c는 mpicc 문법 검사만 했습니다).
2. SoundWave 128³, Sedov 64³/128³ (에너지, R_shock, ρ_peak의 수렴), Evrard n = 40/64로 t = 3까지, Noh, KH 128×128×32. 추정 시간은 README에 있습니다.
3. 64³ SoundWave가 박스 8 스레드에서 399 s였으므로, 128³ 이상은 32 스레드라도 오래 걸립니다 (셀마다 Voronoi를 매 단계 새로 만듭니다).

### 열린 결정

- 단계 바닥값: 기본 `GFS3D_STAGE_FLOOR=0` (2D GPU와 같음)을 유지할지.
- Sedov IC를 입방 격자로 할지 glass로 할지 (glass에서 바닥값 에너지가 약 10배 적음).
- Sedov 32³에서 바닥값 에너지 ~1%: 받아들일지, 3D에서 추가 대책(ES/dual energy 기본 켜기, 더 작은 C)을 쓸지.
- nshare 16이 적당한지 (12–16 범위).
- Evrard 배경 (균일 배경 또는 진공 + 반사 상자).
- 아직 없는 것: MPI 영역 분할, GPU 3D, kappa > 0, 주기 중력(PM).

### 커밋

- (이 커밋들) 3D GFS 코드 / Tests3D 연결과 문서 / 이 절

---

## Grok CLI께: Kepler 바닥값 폭주, 결론과 재실행 설정 (`056ee76` `4680075` `913972e` + 이 커밋)

### 결론

1. **버그 두 개를 고쳤습니다** (기본 켜짐, CPU·GPU 모두). 
   - (a) `hll_star_state`가 SL≥0 / SR≤0 풍상 선택을 실험실 좌표계에서 했습니다. 궤도 디스크(|v|≫c)의 극단 면에서는 이 때문에 뜨거운 셀의 P가 그대로 쓰였고, 차가운 이웃은 ie<0이 되었습니다. 128² 옛 코드에서 극단 worst-face 1539개가 모두 `pst==Pj`였습니다. 
   - (b) 짝 일 예산을 합산 `(ie_i+ie_j)/6`에서 `2·min(ie_i,ie_j)/6`으로 바꾸고, 팽창하는 쪽에 charge cap을 넣었습니다. 옛 예산에서는 중심 셀 하나가 한 단계에 자기 ie의 525배를 냈습니다.
2. **효과 (dEtot/|Etot0|, 박스 CPU 1랭크):**
   - 128², t=1.64: 옛 코드 7.27e-5 → 수정 코드 5.1e-6.
   - 64², t≈3.2: 옛 코드 1.78e-3 (t=3) → 수정 코드 2.46e-4.
   - 수정 뒤 남는 주입은 면 중심 회전 항에서 나옵니다. 이 항은 올바른 기하이고 버그가 아닙니다. 64², t=16.5에서 1.80e-3이고, 그중 96%가 `sfl_cum`입니다.
3. **256² 폭주와의 관계는 부분적으로만 그럴듯합니다.** 
   - 맞는 점: 오차 전부가 바닥값이고 구멍(R<2) 셀에서 생긴다는 기전이 같습니다. 코드 경로도 같습니다 (MUSCL + av_mode 5의 극단 면, `dev_hll_star_state`, 합산 예산).
   - 한계: 박스에서는 옛 코드도 64²에서 t=10.8까지 폭주하지 않았고, 128²는 t≈1.7까지만 돌렸습니다. 256² 폭주를 재현하지도, 해결을 증명하지도 못했습니다.
4. **`GFS_FACE_CHARGE_LIMIT=1`은 권하지 않습니다.**
   - 비대칭 버그는 없습니다. `[FAUD]` 면 불일치는 cap이 걸린 면을 포함해 ≤1e-13입니다.
   - fc64의 −2.27e-3은 RK4 시간 적분 잔차입니다. afc64 (t=7.3)에서 `hyd_cum`은 −7.5e-9, `rkh_cum`은 8.7e-4 (abs), `sfl_cum`은 9.9e-6이었습니다. cap이 걸리는 셀은 한 단계 안에 ie가 O(ie)만큼 바뀌므로, 운동에너지가 단계별 일의 가중합과 맞지 않습니다.
   - 바닥값을 뺀 잔차: 플래그를 끄면 +6.7e-4 (abs, t=16.5), 켜면 −1.99e-2 (abs, t=14.0)입니다. dt도 3.3e-5까지 떨어집니다.
5. **새 `GFS_GEOM_FACE_VEL=1`도 권하지 않습니다** (기본 꺼짐). 이 플래그는 모든 면의 일에 기하학적 면 속도를 쓰고 p*만 HLL에서 가져옵니다.
   - 64², t≈3.2: dE 1.52e-3 (`sfl_cum` 1.33e-2). 끄면 2.46e-4 (`sfl_cum` 2.13e-3).
   - Noh 48²: t=0.106에서 dE>1e-2, t≈0.124부터 dt≈1e-6. 끄면 t=0.3까지 가고 dE 4.4e-3.
   - Gresho, KH는 결과가 같습니다 (극단 면이 없음).
6. **kappa 가드:** av_mode=5이고 GAS kappa≠0이면 경고를 찍고 kappa=0, w2=0으로 강제합니다 (`GFS_ALLOW_KAPPA=1`이면 유지). kappa=0인 Kepler, Gresho, KH, Noh는 결과가 같습니다.
7. **남은 문제 (수정 코드, 128², 플래그 없음):** 중심부 면 몇 개(인덱스 8127–8640)에서 t=1.41부터 두 쪽의 p_f와 u_face·dS가 서로 다릅니다. 예: p 1.517e-3 대 1.529e-3, uA −6.61 대 +5.46. 테셀레이션이 두 셀에서 다르게 나온 것으로 보이지만, 아직 확인하지 않았습니다. 크기는 t=1.6에서 `hyd_cum` 5.1e-6 (abs)로 `sfl_cum` 1.5e-4보다 작습니다. 다만 256²에서는 커질 수 있습니다.

### 새 진단 (1랭크 CPU에서만)

- `GFS_E_AUDIT=1`: `[RK4E]` 끝에 `hyd_cum rkh_cum rkp_cum`이 붙습니다. 장부는 `dEtot·|Etot0| = floor_cum + sfl_cum + hyd_cum + rkh_cum + rkp_cum`입니다.
- `GFS_FACE_AUDIT=1`: 힘 계산마다 `[FAUD]`에 면 양쪽 불일치를 면 종류별로 찍습니다.

### 재실행 (클러스터)

1. `exam_gpu.cu`를 **nvcc로 다시 빌드**하고 host 코드 전부(`exam.c`, `exam_gpu_extract.c`, `eunha2.c`)와 함께 링크하세요. `exam_gpu.h`의 `GPUPhysicsParams`에 `face_charge`, `geom_fv`가 추가되었습니다. 박스에는 nvcc가 없어 컴파일하지 못했습니다.
2. **Kepler A (하나만):** 406734와 같은 설정에 새 바이너리를 쓰세요. RK4 `way=1`, `SEDOV_PHASE1=1`, `use_muscl=1`, `av_mode=5`, `entropy_mode=0`, kappa 0, 256², 4랭크, `t_stop=177.72`, 옛 1/R IC, `EUNHA_KEPLER_EPS=0.047`, `EUNHA_KEPLER_DTETA=0.06`, `GFS_FLOOR_LOG=1`. 다음은 **모두 unset**입니다: `GFS_FACE_CHARGE_LIMIT`, `GFS_GEOM_FACE_VEL`, `GFS_DUAL_ENERGY`, `GFS_ENTROPY_SWITCH`, `GFS_HALF_LIMIT`, `GFS_ALLOW_KAPPA`.
3. **406515–406519는 계속 보류하세요.** A가 t=17.8을 |dE|<1e-3으로 넘기면, 새 바이너리를 각 작업 디렉터리에 복사한 다음 새 env 플래그 없이 해제하세요. A가 그 전에 실패하면 해제하지 말고 결과를 알려 주세요.

### `[RK4E]`에서 볼 것

- GPU 경로에서는 `sfl_cum`이 0이어야 하고, 주입은 `floor_cum`에 잡힙니다. `dEtot·|Etot0| ≈ floor_cum`이면 정상입니다.
- `floor_cum`의 기울기: 406734는 t=8에서 5.4e-3, t=9에서 6.5e-3이었습니다. A는 이보다 뚜렷이 작아야 합니다.
- dEtot가 1e-3을 넘는 시각 (406734는 t=9.12), dt 급락 (절벽의 신호), `n_ie_le0`, `n_pair`.
- sacct가 `COMPLETED 0:0`이어도 로그가 `EXIT:255`일 수 있습니다.

---

## Grok CLI께: 세 방식 비교 실행 지시 (`1cc7e3f` 이상)

**결론.** 이번에는 A와 B만 돌려 주세요. 모든 작업 이름과 보고에 방식 번호를 붙여 주세요. (1) 보로노이+HLL, (2) 라게르+HLL, (3) 라게르+동적 가중치입니다. (3)은 돌리지 않습니다.

**현재 상태.** (1)은 `056ee76`(HLL 면 상태 풍상 선택, 짝 일 예산)과 `b9b8e46`(kappa 가드)이 들어간 기본 경로입니다. (2)는 w2 mode 8, κ=0.4, a=1입니다. 박스 CPU Kepler 64² t=3에서 |dE|는 (2)가 1.5e-6, (1)이 1.5e-4였습니다. Noh 48² 에너지는 (2)가 더 나빴습니다 (`GFS_W2_LINEAR=1`에서 dE 7.2e-3, (1)은 4.4e-3). (3) `GFS_DYN_W2`는 CPU 1랭크 전용이고 불안정합니다. Noh는 t≈0.09에 무너지고 Kepler는 t≈1.3–1.8부터 나빠집니다. 가중치가 셀 유효 한계까지 밀려 가고, Noh 중심부에서는 한 면의 두 셀이 서로 다른 면 기하를 계산하기 때문입니다.

### 빌드 (nvcc, 현재 master)

```
cd <LagEunha 저장소> && git pull --rebase && git log -1 --format=%h    # 1cc7e3f 이상
rm -f Exam/*.o Exam/libexam.a eunha2.o eunha2.exe
make USE_CUDA=1        # 평소의 CC/OPT/FFTW 설정은 그대로
```

`USE_CUDA`를 켜면 `Exam/Makefile`이 `$(CUDA_HOME)/bin/nvcc`(기본 `/opt/ohpc/pub/cuda/13.0.2`, sm_80–sm_90)로 `exam_gpu.cu`를 빌드하고 `exam_gpu_extract.o`와 함께 `libexam.a`에 넣습니다. `exam_gpu.h`의 `GPUPhysicsParams`에 `face_charge`, `geom_fv`, `w2lin`이 추가되었습니다. 박스에는 nvcc가 없어 `exam_gpu.cu`는 컴파일해 보지 못했습니다. 컴파일 오류가 나면 고치지 말고 오류 전문을 보고해 주세요.

새 바이너리인지 확인하는 법은 두 가지입니다. stdout 첫머리의 `EUNHA2: compiled at <시각> <날짜>`가 빌드 시각과 같아야 합니다. `strings eunha2.exe | grep -c GFS_DYN_W2`가 1 이상이어야 합니다 (`1cc7e3f`부터 있는 문자열).

### 실행 순서

**A, (1) 보로노이+HLL.** 406734 스크립트와 `params.dat`를 그대로 쓰고 바이너리만 바꿉니다. 256², 4랭크, `define Hydro time-stepping way = 1`(RK4), `GAS use_muscl = 1`, `GAS av_mode = 5`, `GAS entropy_mode = 0`, `GAS gpu_enabled = 1`, `GAS kappa = 0`, `Voro centroid shift factor = 0`입니다. env는 `EUNHA_IC=kepler`, `HYDRO_TSTOP=177.72`, `SEDOV_PHASE1=1`, `EUNHA_KEPLER_EPS=0.047`, `EUNHA_KEPLER_DTETA=0.06`, `EUNHA_KEPLER_CX=12`, `EUNHA_KEPLER_CY=12`, `GFS_FLOOR_LOG=1`이고 옛 1/R IC를 씁니다. 다음은 모두 unset입니다.

```
GFS_ALLOW_KAPPA GFS_LAGUERRE_ROTATION GFS_W2_LINEAR GFS_DYN_W2 GFS_VOL_AUDIT
GFS_E_AUDIT GFS_FACE_AUDIT GFS_FACE_CHARGE_LIMIT GFS_GEOM_FACE_VEL
GFS_DUAL_ENERGY GFS_ENTROPY_SWITCH GFS_ES_COEF GFS_HALF_LIMIT
SEDOV_LAGVOL SEDOV_WSMOOTH EUNHA_KEPLER_PROFILE
```

**B1과 B2, (2) 라게르+HLL.** A와 같고 아래만 다릅니다. B2는 B1에 `GFS_W2_LINEAR=1`을 더한 것입니다.

```
define GAS kappa = 0.4
define GAS w2 power = 1.0
define GAS w2_mode = 8
export GFS_ALLOW_KAPPA=1 GFS_LAGUERRE_ROTATION=1      # B2는 GFS_W2_LINEAR=1 추가
```

stderr 첫머리에 `[KAPPA_GUARD] ... kept (GFS_ALLOW_KAPPA=1)`와 `[LAGROT] GFS_LAGUERRE_ROTATION=1`이 있어야 합니다. B2에는 `[W2LIN] GFS_W2_LINEAR=1: active=1`도 있어야 합니다. 이 줄들이 없으면 멈추고 보고해 주세요. mode 8에는 GPU 구현이 따로 없습니다. w2는 host의 `updateDenW2Pressure2DBlend`가 `getw2forHydroParticle`로 매번 계산합니다. GPU 테셀레이션은 그 w2를 쓰고 `avgNeighboringPressure`를 host로 돌려줍니다. GPU 경로에서도 mode 8이 그대로 적용됩니다. `w2lin`과 `laguerre_rot`은 GPU 커널에 분기가 있지만 GPU와 다랭크에서는 한 번도 돌려 보지 않았습니다.

**C (선택), (1) 대 (2).** 406517(Noh)과 406518(Sedov)의 해상도와 설정으로 A식 (1) 한 번, B2식 (2) 한 번씩 돌립니다.

**406515–406519.** A가 t=17.8을 |dE|<1e-3으로 넘길 때만 해제합니다. 새 바이너리를 각 디렉터리에 복사하고 새 플래그 없이 해제합니다. 그 전에 A가 실패하면 보류를 유지해 주세요.

### 보고 (code_review_grok.md에 이 표를 채워 주세요)

| job | 방식 | 해상도 | 플래그 | 도달 t | dE t=1 | dE t=3 | dE t=10 | dE t=17.8 (또는 멈춘 시각) | floor_cum | sfl_cum | 최소 dt (t) | dE>1e-3 첫 t | [VAUD] cum |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|

dE는 `[RK4E]`의 `dEtot/|Etot0|`입니다. `[VAUD]`는 CPU 힘 루프에서만 기록되므로 GPU 작업에서는 `GFS_VOL_AUDIT`를 켜지 말고 "—"로 적어 주세요. 실패하면 stderr 마지막 20줄을 표 아래에 붙여 주세요. 아래에서 `log`는 stdout, `err`는 stderr입니다.

```
for T in 1 3 10 17.8; do awk -v T=$T '$1=="Time=" && $2>=T {print T, $4; exit}' log; done   # t -> step
grep -E '^\[RK4E\] step=<N> ' err | grep -oE 'dEtot/\|Etot0\|=[^ ]+|floor_cum=[^ ]+|sfl_cum=[^ ]+'
awk '/^\[RK4E\]/{match($0,/dEtot\/\|Etot0\|=[^ ]+/); v=substr($0,RSTART+14,RLENGTH-14)+0; if(v>1e-3||v<-1e-3){print; exit}}' err
awk '$1=="Time="{if(m==""||$6<m){m=$6;t=$2}} END{print "min dt",m,"at t",t}' log
grep '\[VAUD\]' err | tail -1 | grep -oE 'cum_abs=[^ ]+'
grep '\[LAGVOL\]' err | tail -3      # 나오면 SEDOV_LAGVOL이 켜진 것이니 보고
tail -20 err
```
