# LagEunha GFS: review of the tree at 5f2fb86

**Scope.** The geometric-face hydro path after commit `5f2fb86` (`master`, 2026-09-27). Files: `Exam/exam.c`, `Exam/exam_gpu.cu`, `Exam/exam_gpu_extract.c`, `Exam/gfs_pair.h`, `Exam/KH/util.c`, `Hydro1DExam/laguerre_sod.c`, `Colorize/color.c`.
**Relation to `code_review.md`.** That note reviewed `75aaeed` (2026-09-26). This note says what that review got right, what `5f2fb86` already changed, and what is still wrong. It does not repeat the Laguerre-weight solvability argument. That argument stands: one weight per particle cannot place every 2D face on a chosen Riemann speed.

**Evidence.**
- **[CODE]** read in `5f2fb86`.
- **[RUN]** a finished run of this scheme, or of the binary copied into a job before `5f2fb86` but after the pair-pressure change. The distinction is marked on each run.
- **[OPEN]** plausible, not measured on this commit.

---

## 0. Summary

1. The 1D GFS suite with the bounded pair pressure finishes all eight problems. Sod, Noh, the double rarefaction, and Shu–Osher (`N=800`) are unchanged to the published digits. The blast `L1(ρ)` fell from 0.154 to 0.055 and its energy error is `1.5e-5`. Lax and the contact improved. **[RUN, pair-pressure binary, before the centroid term]**
2. 2D Gresho at `128²`, `t=1`, finished with `dE/E0 = 3.9e-9`, `n_neg = 0`, `dt = 1.9e-4`. **[RUN, same binary]**
3. The thick Kepler disk on 4 H200s (job 406510, that same binary, point-mass `v_φ = R^{-1/2}`, `P = 1e-6`) is not a disk test anymore. By `t = 2.58` (step 12955) `dE/E0 = 6.1e4`, `ρ_max = 61`, `dt = 4.2e-7`, `n_pair ≈ 210`. The thermal-plus-ram pair pressure is creating energy inside the cutoff. **[RUN]**
4. `5f2fb86` fixes three items from the previous review: the Voronoi centroid term on zero-weight faces (`voronoi_face_rotation`, `exam.c` 4168), the CPU `dt3` dimensions (`exam.c` 5671), and the RK4 stage update when `die` means `dE/dt` (`phase1_ie_stage`, `exam.c` 4185). None of these are in the binary that job 406510 is running.
5. What is still open, in the order that matters for the runs now queued: the pair-pressure normalisation, the extreme-face HLL speed versus the geometric face, the 1D cold-pressure gate, and the old `dt3` left in `Exam/Sedov/vch_hydro.c`.

---

## 1. Pair pressure

**What is in the code.** `gfs_pair.h`:

```
d_c = 0.2 * min(length_i, length_j)
P   = 16 * ρ̄ * (c̄² + v_close²) * (1 - d/d_c)²    for d < d_c
P   = 0                                         otherwise
```

In 2D the length is `sqrt(cell area)`. In 1D it is the cell length. `v_close` is the approaching speed along the pair and is zero when the particles are separating or already at rest. The pressure is added to the face pressure, so the force is pairwise and the work is booked in `dte`.

**Why the factor is 16.** A form `4 ρ c² (1-d/d_c)²` stays finite as `d → 0` and matches `ρ c²` at `d = d_c/2`. On the 1D Woodward–Colella blast that barrier did not stop a Mach-2 approach: the gap fell to `3e-11` while the closing speed was still 9, and the face CFL then drove `dt` to `1e-12`. Replacing `c²` by `c² + v_close²` and raising the coefficient from 4 to 16 let the blast finish. The coefficient is an empirical fix for one shock, not a derived stopping condition. **[RUN]**

**What that does to a disk.** On a circular shear flow, neighbours inside the cutoff have a nonzero `v_close` even when no particle is falling in. `16 ρ v_close²` is then a large face pressure on an ordinary shear pair. Job 406510 shows the consequence. Energy is already wrong before the timestep collapses:

| step | dt | dE/E0 | ρ_max | n_pair |
|---:|---:|---:|---:|---:|
| 2000 | 1.4e-4 | 0.14 | 1.4 | 84 |
| 8000 | 1.0e-4 | 1.3 | 14 | 114 |
| 10000 | 1.9e-5 | 14 | 18 | 110 |
| 12000 | 1.2e-6 | 3.8e3 | 73 | 194 |
| 12955 | 4.2e-7 | 6.1e4 | 61 | 210 |

`n_neg` at the last step is 12 out of 65536. The smallest cell area is still of order the mean (`Vmin*Nx² = 1.6`), so this is not a collapsed cell. It is the ram-pressure term on ~100 faces. **[RUN, binary without the centroid term]**

A stationary pair (`v_close = 0`) is still pushed by `16 ρ c² (1-d/d_c)²`. That part is what KH and the central Kepler cavity need. The `v_close²` piece is what the blast needed and what the disk cannot afford.

**Fix.** Keep the thermal piece. Drop `v_close²` from the pressure, or gate it so it cannot exceed the local thermal pressure by a set factor, and retune the blast on that law. Do not turn the term off for subsonic approach: the faces that killed the first blast were Mach 2.1–2.6, and the pairs that stall KH have already reached `Δv = 0`.

---

## 2. Face velocity

**Ordinary faces.** For `w² = 0` the code now adds the centroid correction (`exam.c` 4168, applied at 5089 and 5503; the same block is in `exam_gpu.cu` and `exam_gpu_extract.c`):

```
Δu_n = −(u_q − u_p)·(c_f − m) / d
```

`c_f` is the midpoint of the two Voronoi corners. The correction is along the generator direction `er`. Both ends of the face compute the same lab-frame velocity, so the pairwise energy cancelation is unchanged. The hard-wall branch still removes the normal velocity after this correction. Extreme faces still replace the face velocity by the HLL normal speed (`exam.c` 5583 in the previous layout; the `phase1_extreme` branch after 5503).

**What this does not fix.** The extreme-face branch is still the inconsistency described in `code_review.md` §1.5. Geometry moves at the Voronoi speed. `PdV` on that face uses the HLL speed. A cold cell next to a strong shock can lose internal energy that does not match its volume change. The 2D trigger is now only `P_max > 100 P_min` (`exam.c` 5169). The old `P_min < 1e-3` clause is gone, which is required for a uniformly cold Noh inflow. The 1D driver still has the old clause (`laguerre_sod.c` 1503: `pmin < 1e-3 || pmax > 100 pmin`). The 1D Noh at `N=200` did not show it, because that run's pressures after the first step are not all below `1e-3`. A colder rescaling would.

**Float equality.** The rotation term is skipped unless both weights are exactly zero. `Kappa = 0` writes that zero. A roundoff weight silently drops the correction.

---

## 3. Time step

**2D KDK.** `exam.c` 9841–9845 takes the minimum of the uniform-grid CFL and the face CFL returned by the force routine. The previous review's statement that `dt_dummy` is discarded is no longer true. The face CFL is why one close pair freezes every particle. That is intentional. It is not a cure: if nothing reverses the approach, `dt` still goes to zero. Job 406510 is that case.

**CPU `dt3` in the blend loop** is now `0.1 d / |Δv|` (`exam.c` 5676). The old expression had units time²/length. The same old expression is still compiled in `Exam/Sedov/vch_hydro.c` 169 and `vch_reverse.c` 169. Those are not the 2D GFS path. They will bite any run that still calls them.

**1D.** `laguerre_sod.c` refuses a step whose positions cross, halves `dt`, and also limits the step by the closing speed and by the pair acceleration. That is why the blast no longer dies in one stage. It does not make a too-weak barrier strong. It only stops the run from accepting a crossed state.

---

## 4. Energy variable

`SEDOV_PHASE1=1` stores `dE/dt` in `die`. The KDK path already knows this (`phase1_half_kick`, `exam.c` 4193). Gravity is applied after that fix, so the potential, which is not part of `E`, does not wipe `ie`.

The blend RK4 path used to add `die*dt` to `ie` and also kick `v`, which counts `ΔKE` twice. `phase1_ie_stage` (`exam.c` 4185) now subtracts `m v·a` before the stage increment. This is the first-order split only. The `½ m |a Δt|²` piece is still handled only in `phase1_half_kick`, and only by the clip that rescales velocity. RK4 stages do not do that clip. Production runs use KDK (`GAS_EVOLMETHOD=3` with `av_mode=5`), so this RK4 path is not what jobs 406510–406519 execute. It is still the path the paper text calls GFS. The environment variable remains the switch. There is no `GAS_SCHEME=GFS` parameter, and output headers still do not record the integrator.

The KDK clip (`exam.c` 4191) keeps total energy and drops momentum when `ie` would go negative. `n_clip` on 406510 grew from 178 at step 2000 to 6142 at step 12955. That counter is the run saying the closure failed.

---

## 5. Boundaries and images

A missing `kh.sao` used to segfault inside `readsao` because `fopen` was not checked. `color.c` 290 now fills a grey image and returns. Flat fields no longer divide by `ccdmax-ccdmin = 0`. The hydro step is no longer killed by the plotter. Gresho's first attempt died this way at `t=1` with an otherwise exact energy; the rerun exited 0.

The cylinder wall still reflects a particle that steps inside and drops the inward normal velocity without putting `½ m v_n²` into `ie` (`Exam/Cylinder/cylinder.c`). That is unchanged. The open cylinder's reported energy difference is still not a closed budget unless inflow and outflow fluxes are counted. The 1D code already does this with `bc_eng_net`. The 2D runs do not.

---

## 6. What the previous review measured, and what changed

| Previous claim | On `5f2fb86` |
|---|---|
| 1D blast energy error is the `e < 1e-14` floor, `L1(ρ) = 0.154` | True for the `wtrack` driver they ran. The current `gfs` blast is `L1 = 0.055`, `ΔE/E0 = 1.5e-5`. The floor is not an immutable GFS error. |
| Paper 1D numbers do not match the `vr` driver | True, and expected. `vr` is CD10 viscosity without the extreme-face HLL. The paper numbers match `./laguerre_sod gfs`. |
| KDK ignores the face CFL | False now. See §3. |
| Extreme if `P_min < 1e-3` | False in 2D. Still true in `laguerre_sod.c` 1503. |
| Midpoint velocity omits the centroid term | Fixed for `w² = 0` in the blend and GPU kernels. Not measured again on the perturbed lattice. |
| RK4 double-counts `ΔKE` under `SEDOV_PHASE1` | The linear term is removed in the blend RK4. The quadratic term is not. Production uses KDK. |
| `dt3` has the wrong dimensions | Fixed in the blend loop. Not fixed in `Exam/Sedov/vch_*.c`. |

The RT error `7.6e-4` was guessed to be the RK4 double count. The production RT log is KDK. That guess should not be kept.

---

## 7. Order of work

1. Rerun the Kepler disk with `v_close²` removed from `P_pair`, or with that term capped by the local `ρ c²`. The thermal piece stays, including at `Δv = 0`. The 1D blast has to be rerun before accepting the coefficient. Job 406510 is not a usable disk.
2. Remove `pmin < 1e-3` from the 1D extreme-face test so it matches 2D.
3. Decide whether an extreme face's `PdV` speed is the HLL speed or the geometric speed, and use one of them for both the flux and the mesh motion. Using both is the remaining energy source at shocks.
4. Replace `dt3` in `Exam/Sedov/vch_hydro.c` and `vch_reverse.c`, or stop calling those files.
5. Record integrator, `av_mode`, and `SEDOV_PHASE1` in the output header. Do not describe RK4 as the GFS integrator while the runs use KDK.

The least-squares weight solve in `code_review.md` §2c is a different scheme. It is not the next patch.

---

## 8. Handoff, 2026-09-27

Two sessions share this repository. The user carries messages between them. Code, comments, and instructions live only on `kjhan0606/LagEunha`.

| Who | Writes | Does |
|---|---|---|
| This session | `code_review_grok.md` | MPI+CUDA runs on Slurm, and the comment on those runs |
| Grokbot | `code_review_grokbot.md` | Coding, 1D and other non-MPI tests, and replies to the instructions below |

Do not edit the other session's review file. Read it, then answer in your own file.

### Instruction for Grokbot

The next code change is the pair pressure in `Exam/gfs_pair.h`. Keep the thermal piece `16 ρ̄ c̄² (1-d/d_c)²`, including when `v_close = 0`. Remove `v_close²`, or cap the whole pressure by a small multiple of `ρ̄ c̄²`. The uncapped ram term is what job 406510 used.

After that change:

1. Rerun the 1D `gfs` suite, including the Woodward–Colella blast at `N=200`. The blast must finish. Report `L1(ρ)` and `ΔE/E0` for all eight problems against the numbers in §0.
2. Remove `pmin < 1e-3` from the 1D extreme-face test (`Hydro1DExam/laguerre_sod.c` around line 1503) so it matches the 2D ratio-only test.
3. Push to `master` and say what commit to run.

Do not start Slurm jobs. Do not cancel 406510 or 406515–406519. This session will rerun MPI+CUDA after the commit is on `master`.

### State of the queue at handoff

406510 (`LAGFORCE_KEP256`) was stopped by this session after 57 minutes. Final state: `t = 2.58`, `dt = 8.4e-12`, `dE/E0 = 2.5e14`, `n_neg = 63/65536`, `n_pair = 437`. It is not a result. The binary is the pair-pressure build from before `5f2fb86`.

406515 (KH), 406516 (RT), 406517 (Noh), 406518 (Sedov), and 406519 (cylinder) are held (`JobHeldUser`) so they do not start on that same binary. Release them only after a new binary is copied into each run directory. Do not cancel them.

1D `gfs` and 2D Gresho `128²` `t=1` (`dE/E0 = 3.9e-9`) already passed on that older binary. Do not repeat them unless the pair-pressure formula changes. If it changes, the blast is the required rerun.
