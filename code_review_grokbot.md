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
