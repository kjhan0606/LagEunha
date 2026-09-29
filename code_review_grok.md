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


### Kepler 406608 stopped

Stopped at about 9.5 minutes. Last sample in `debug/2026-09-27_kepler_gresho.txt` is already decisive: `t=0.65`, `dE/E0=2.81`, `n_clip=2333`, `n_pair=136`. The cap `min(v_close^2, c^2)` did not keep the disk energy flat. Do not release 406515–406519 onto this binary. Gresho 406611 was still running on syn09 when Kepler was stopped.

---

## 9. Run report, capped pair pressure (`79af847`)

Binary built from `master` at `79af847` (`P ≤ 32 ρ̄ c̄²`). H200 was not available, so both runs used H100. Sampled lines are in `debug/2026-09-27_kepler_gresho.txt`. That file was taken while both jobs were still running. The Gresho numbers below are the finished state.

### Gresho — pass

Job 406611, `LAGFORCE_GRES128`, syn09, 2× H100, 4 ranks (two ranks per GPU). Directory `/gpfs/kjhan/LagForce/Gresho_gfs_Nx128_cap`. `128²`, `t = 1`, wall time 4 min 21 s, `EXIT:0`.

| | step 1 | step 993 | step 1985 | step 2977 | step 3969 |
|---|---:|---:|---:|---:|---:|
| t | 0.0005 | 0.370 | 0.577 | 0.792 | 1.000 |
| dt | 5.0e-4 | 2.5e-4 | 2.1e-4 | 2.0e-4 | 2.4e-4 |
| dE/E0 | 0 | −1.0e-6 | −1.2e-6 | −1.2e-6 | −1.3e-6 |
| ρ_max | 1.001 | 1.241 | 1.176 | 1.086 | 1.085 |
| n_pair | 0 | 0 | 0 | 0 | 0 |

`n_neg = 0/16384` and `n_clip = 0` for the whole run. `n_pair` was 0 except for a few steps where it reached 2 or 4. The density image was written. No segfault. This is the same quality as the earlier Gresho on the uncapped binary (`dE/E0 = 3.9e-9`).

### Kepler — fail, stopped

Job 406608, `LAGFORCE_KEP256h4`, syn08, 3× H100, 4 ranks (ranks 0 and 3 shared device 0). Directory `/gpfs/kjhan/LagForce/Kepler_gfs_disk_h200_Nx256_cap`. Point-mass disk, `256²`, target `t = 177.72`. Stopped by request after 9.5 min. It was not going to recover.

| step | t | dt | dE/E0 | ρ_max | n_clip | n_pair |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.019 | 1.9e-2 | 0 | 0.80 | 0 | 0 |
| 236 | — | 5.9e-4 | 0.010 | 0.80 | 171 | 14 |
| 732 | — | 3.9e-5 | 1.01 | 0.80 | 602 | 40 |
| 2001 | — | 9.5e-5 | 2.57 | 1.53 | 1951 | 130 |
| 3476 | 0.647 | 8.8e-5 | 2.81 | 5.45 | 2333 | 136 |

`n_neg` stayed at 0 or 1. The gas energy grew from 8.79 to 33.5 while the density was still near the initial disk value, so this is not a collapsed cell. Capping `v_close²` by `c²` did not stop the heating. The thermal piece `16 ρ c² (1-d/d_c)²` is still on for every face inside the cutoff, and ~100 such faces are enough to triple the energy before one inner orbit (`T(R=2) = 17.8`).

406515–406519 stay held. Do not copy this binary into them. The next Kepler run needs a pair pressure that is silent on a circular shear flow. Gresho shows the gate can stay at `n_pair ≈ 0` when the mesh does not pair. The disk does pair, and the present amplitude is too large.

---

## 10. Pair impulse limited by internal energy

`P0` is the same capped law as `79af847`. The pressure added to the face is smaller when the discrete work would overdraw the internal energy:

```
μ = m_i m_j / (m_i + m_j)
q = (v_j − v_i)·n          (q < 0 approaching)
B = max(ie_i, 0)/nshare + max(ie_j, 0)/nshare
Jmax = μ (sqrt(q² + 2B/μ) − q)
P = min(P0, Jmax / (A Δt))
```

`nshare` is 6 in 2D and 2 in 1D. The work being limited is `W(J) = q J + J²/(2μ)`, so a stationary pair is limited by `sqrt(2 μ B)` and an approaching pair may still brake when `B = 0`. The function is `gfs_pair_work_limit` in `Exam/gfs_pair.h`. It is called from `exam.c`, `exam_gpu.cu`, `exam_gpu_extract.c`, and `Hydro1DExam/laguerre_sod.c`. On the 2D paths `Δt` is `GAS_dtold`, the previous accepted step. On the 1D path it is the `dt` of the current `rk4_step`. The 1D suite has not been rerun with this limit. **[CODE]**

### Kepler 406734 — running

406654 (`LAGFORCE_KEP_RK4`, RK4 with the cap but without this limit) was stopped after it had already gone to `E_hyd = nan` and `floor_cum = inf`. `mpirun` returned 255. The batch script still exited 0, so `sacct` says `COMPLETED`.

406734, `LAGFORCE_KEP_W`, syn08, 3× H100, 4 ranks. Directory `/gpfs/kjhan/LagForce/Kepler_gfs_rk4_budget_Nx256`. Same disk as 406654: RK4, `EPS = 0.047`, `DTETA = 0.06`, `P = 1e-6`, `t_stop = 177.72`. At `t = 1.08` (step 1244):

| | value |
|---|---:|
| dt | 6.3e-4 |
| dEtot/\|Etot0\| | 1.8e-6 |
| floor_cum | 6.5e-6 |
| n_ie_le0 | 0 |
| n_pair | 8 |

406654 was still mild at `t ≈ 4.3` and crossed `|ΔE| = 0.01` at `t = 4.59`. This section is the launch, not a pass. 406515–406519 stay held. Do not copy this binary into them until this disk stays bounded through that time.

---

## 11. Report to Grokbot: Kepler 406734 through the old blow-up

406734 is still running (syn08, 3× H100, 4 ranks, `/gpfs/kjhan/LagForce/Kepler_gfs_rk4_budget_Nx256`). This is the same RK4 disk as 406654 (`EPS = 0.047`, `DTETA = 0.06`, `P = 1e-6`, `t_stop = 177.72`) with `gfs_pair_work_limit` in the binary. Wall time about 41 min when this was written. No NaN. **[RUN]**

406654, without the work limit, crossed `|ΔE| = 0.01` at `t = 4.59`, reached `0.1` at `t = 5.01`, and was at `~10^3` by `t = 5.84` (`floor_cum = inf` when it was stopped). 406734 went through that interval as follows.

| t | step | dt | dEtot/\|Etot0\| | floor_cum | n_ie_le0 | n_pair | E_hyd | E_pot | E_tot |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 4.00 | 5945 | 6.44e-4 | 1.69e-4 | 1.46e-3 | 1 | 32 | 8.79003 | −17.5802 | −8.79017 |
| 4.59 | 6889 | 6.13e-4 | 2.52e-4 | 2.17e-3 | 4 | 50 | 8.79086 | −17.5803 | −8.78944 |
| 5.00 | 7546 | 6.49e-4 | 3.34e-4 | 2.85e-3 | 4 | 52 | 8.79072 | −17.5794 | −8.78872 |
| 5.18 | 7837 | 4.70e-4 | 3.57e-4 | 3.03e-3 | 2 | 64 | 8.79086 | −17.5794 | −8.78852 |

Over the run so far the largest `|dEtot/|Etot0||` is `3.6e-4`, the largest `floor_cum` is `3.0e-3`, `n_pair` peaked at 74, and `n_ie_le0` peaked at 9. In the last `Δt = 0.5` (`t = 4.68` to `5.18`) `floor_cum` rose from `2.38e-3` to `3.03e-3`. `E_pot` has moved from `−17.581` to `−17.579`. The gas-energy error tracks `floor_cum / |Etot0|` (`|Etot0| ≈ 8.79`). The cells that hit `ie ≤ 0` are still being replaced by `P = 1e-6`, and that replacement is the whole energy error.

So the work limit stopped the `t ≈ 5` runaway. It did not stop the slow floor. `t = 5.18` is about 0.3 of an inner orbit (`T(R = 2) = 17.8`). Do not call this a pass, and do not copy this binary into 406515–406519. Leave 406734 running.

Please rerun the eight 1D GFS problems with this limit. `Hydro1DExam/laguerre_sod.c` already calls `gfs_pair_work_limit` (`nshare = 2`, `Δt = gfs_1d_dt` set at the start of `rk4_step`). Do not overwrite `gfs_*.dat`. Compare with the cap-only `gfs_pair_*.dat` from section 9, especially the blast, where the old closing faces were Mach 2 and the work of braking is negative, so `B = 0` must still allow `P0`. On the 2D paths the `Δt` passed to the limiter is `GAS_dtold`, the previous accepted step, not the RK4 stage fraction.

---

## 12. Kepler 407021 started; MUSCL is on

406734 is still the baseline and was left running. The diagnostic run is 407021, `LAGFORCE_KEP_DE`, syn09, 3× H100, 4 ranks, `/gpfs/kjhan/LagForce/Kepler_gfs_rk4_dual_Nx256`. The binary is `fe28651` with the host `exam.c` relinked. `exam_gpu.cu` is unchanged. Same disk as 406734. **[RUN]**

Flags from the job, not from a guess:

| | |
|---|---|
| `SEDOV_PHASE1` | 1 |
| `GAS use_muscl` | 1 |
| `GAS av_mode` | 5 (MUSCL + HLLC) |
| `GAS gpu_enabled` | 1 |
| time-stepping `way` | 1 (this is the RK4 path; the params comment still says KDK) |
| `GFS_DUAL_ENERGY` | 0.5 |
| `GFS_FLOOR_LOG` | 1 |

The cell-centred HLL overwrite in `exam.c` runs only when `sedov_phase1_on() && !use_muscl`. This job does not take that branch. The face pressure is the MUSCL/HLLC state from `av_mode` 5.

Step 1 of 407021: `E_int = 8.644e-4`, `n_de = 0`, `de_cum = 0`, `floor_cum = 0`, `dEtot = 0`, `n_pair = 0`. That `E_int` matches the thermal energy of the initial box. No `[RK4F]` line in the first steps. 406515–406519 stay held.

---

## 13. Baseline runaway at t = 9.5, and the first [RK4F] radii

Both jobs were still running when this was written. 406734 was not stopped. **[RUN]**

### 406734, the work-limit baseline

The slow floor lasted past the old blow-up and then failed the same way. 406654 crossed `|ΔE| = 0.01` at `t = 4.59`. 406734 did it at `t = 9.54`.

| t | dEtot/\|Etot0\| | floor_cum | n_ie_le0 | n_pair | dt |
|---:|---:|---:|---:|---:|---:|
| 6.00 | 4.4e-4 | 3.7e-3 | 6 | 62 | 6.2e-4 |
| 8.00 | 7.5e-4 | 5.4e-3 | 4 | 92 | 6.4e-4 |
| 9.00 | 9.0e-4 | 6.5e-3 | 8 | 138 | 5.3e-4 |
| 9.50 | 6.6e-3 | 2.6e-2 | 8 | 132 | 6.1e-4 |
| 9.54 | 1.0e-2 | 3.4e-2 | — | — | — |
| 9.60 | 1.39 | 11.8 | 11 | 172 | 7.7e-5 |

`|ΔE| = 10^{-3}` was crossed at `t = 9.12`. From `t = 9.50` to `t = 9.60` the floor jumped by about 400× and `dt` dropped. `E_pot` was still `−17.577`. This is the pressure-floor runaway, delayed, not an orbit decay. The useful floor-rate baseline is the interval before `t ≈ 9.1`.

### 407021, dual energy, t = 0.83

527 `[RK4F]` lines, every one with `reset = 1`. None are in the disk.

| region | lines | R |
|---|---:|---|
| R < 0.3 | 105 | 0.243–0.264 |
| 0.3–1.9 | 422 | 0.336–0.399 |
| 1.9–2.1 | 0 | |
| R ≥ 2.1 | 0 | |

`ie_K > 0` on all of them. 24 lines have `ie ≤ 0`; the rest are positive and sit at `ie/ie_K = 0.13–0.50`, which is the `η = 0.5` threshold and not an 8-cell-energy drain. `den` is `6.6e-4–2.5e-3`, the central floor gas (`ρ_floor = 1e-3`). `floor_cum = 0` and `n_ie_le0 = 0` at `t = 0.83`: the dual-energy reset is catching these cells before the `P = 1e-6` floor. `de_cum = 4.8e-6`, about `1e-8` per logged event, and `dEtot/|Etot0| = 1.2e-6`. The ledger matches that small injection.

`E_int` went from `8.64e-4` at `t = 0` to `1.13e-3` at `t = 0.83`. `E_hyd` and `E_tot` did not move with it, so that rise is a conservative exchange, kinetic energy into internal energy. `de_cum` is 50× too small to be the source. The log does not split `E_int` by radius, so this does not yet say the disk is heating. It does say the cells that lose half their adiabatic energy are at `R = 0.24–0.40` and not at `Rin = 2`. Too early to call the `t ≈ 9` runaway fixed. `n_pair = 8`.

---

## 14. Kepler closed: both runs hit the same cliff

Both jobs are finished. Neither is a pass. 406515–406519 stay held. Do not copy either binary into them. **[RUN]**

Shared setup: RK4 (`way = 1`), `SEDOV_PHASE1 = 1`, `GAS use_muscl = 1`, `av_mode = 5`, `EPS = 0.047`, `DTETA = 0.06`, `P = 1e-6`, annulus `R = 2–8` in a box of side 24 centred at `(12, 12)`, `t_stop = 177.72`. H100×3, 4 ranks. Initial `E_pot ≈ −17.581`, `|Etot0| ≈ 8.79`, initial thermal energy `E_int(0) = 8.644e-4`.

| job | what | node | wall | end |
|---|---|---|---:|---|
| 406734 `LAGFORCE_KEP_W` | work limit, no dual energy | syn08 | 1 h 37 min | log: signal 11, `EXIT:255`. `sacct` says `COMPLETED 0:0` because the batch script does not fail on `mpirun` |
| 407021 `LAGFORCE_KEP_DE` | `fe28651` plus `GFS_DUAL_ENERGY=0.5`, `GFS_FLOOR_LOG=1` | syn09 | 2 h 00 min, ended 2026-09-28 02:43 | `exam_gpu.cu:3428` CUDA `invalid argument`, then signal 9, `EXIT:255`. Same `sacct` `0:0` |

406654, the RK4 disk without the work limit, crossed `|ΔE| = 0.01` at `t = 4.59`. The work limit moved that crossing to `t = 9.537`. Dual energy moved it to `t = 10.311`. The slope before the crossing is the same, and the jump after it is the same.

### 406734

`floor_cum` is the whole energy error. `E_pot` is still `−17.577` at the last finite step.

| t | step | dt | dEtot/\|Etot0\| | floor_cum | n_ie_le0 | n_pair | E_pot |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 6.00 | 9141 | 6.21e-4 | 4.41e-4 | 3.67e-3 | 6 | 62 | −17.579 |
| 8.00 | 12381 | 6.38e-4 | 7.53e-4 | 5.42e-3 | 4 | 92 | −17.578 |
| 9.00 | 14024 | 5.31e-4 | 8.96e-4 | 6.54e-3 | 8 | 138 | −17.578 |
| 9.123 | 14224 | 4.40e-4 | 1.02e-3 | 7.62e-3 | 11 | 152 | −17.577 |
| 9.50 | 14858 | 6.14e-4 | 6.64e-3 | 2.63e-2 | 8 | 132 | −17.577 |
| 9.537 | 14918 | 6.13e-4 | 1.03e-2 | 3.37e-2 | 12 | 138 | −17.577 |
| 9.586 | 15016 | 3.51e-4 | 0.102 | 0.545 | 16 | 159 | −17.577 |
| 9.597 | 15075 | 1.04e-4 | 1.32 | 9.55 | 13 | 172 | −17.577 |

NaN starts after `t = 9.60`. The log then advances to `t = 17.68` with `dt ≈ 0.8` and `floor_cum = inf`. The record that can be compared is `t ≤ 9.1`.

### 407021

`floor_cum` stays 0 through the finite log. `n_ie_le0` stays 0 because that count is after the reset. The reset sets `ie` to `ie_K`, and `de_cum` adds `ie_K − ie`. A cell that has gone largely negative therefore injects its whole deficit. That is the same created energy as replacing `ie` with the `P = 1e-6` floor. Early on, the logged cells were only shallow (`ie/ie_K = 0.13–0.50`, about `1e-8` each). The cliff is the deep deficit.

| t | step | dt | dEtot/\|Etot0\| | de_cum | E_int | n_de | n_pair | E_pot |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 6.00 | 9294 | 6.20e-4 | 4.48e-4 | 3.72e-3 | 6.39e-3 | 3 | 62 | −17.579 |
| 8.00 | 12636 | 6.24e-4 | 7.34e-4 | 6.02e-3 | 1.06e-2 | 3 | 126 | −17.577 |
| 9.00 | 14270 | 6.14e-4 | 9.04e-4 | 7.32e-3 | 1.32e-2 | 7 | 168 | −17.576 |
| 9.304 | 14778 | 6.13e-4 | 1.00e-3 | 8.27e-3 | 1.42e-2 | 4 | 168 | −17.575 |
| 10.00 | 15895 | 6.58e-4 | 1.60e-3 | 1.40e-2 | 1.71e-2 | 3 | 140 | −17.570 |
| 10.311 | 16349 | 6.93e-4 | 1.24e-2 | 8.54e-2 | 7.72e-2 | 6 | 110 | −17.566 |
| 10.413 | 16500 | 5.20e-4 | 0.103 | 1.18 | 0.679 | 10 | 134 | −17.563 |
| 10.432 | 16578 | 2.34e-4 | 1.07 | 5.07 | 6.30 | 11 | 132 | −17.562 |

At `t = 8` the two jobs match: `dE ≈ 7.5e-4` and about `6e-3` of created energy, booked as `floor_cum` in 406734 and as `de_cum` in 407021. `E_int(8) / E_int(0) ≈ 12`. Of that rise, `6.0e-3` is the injection and the rest is a conservative exchange of kinetic energy into internal energy. By the cliff, `E_int` and `de_cum` move together. `E_pot` at the last finite step is `−17.556`. NaN starts at `t = 10.441`. The last line in the log is `t = 10.93`, `dt = 0.099`, `de_cum = inf`.

### Where the resets are

Through `t = 7.94` (`dE = 7.3e-4`, still on the slow slope) there were 35200 `[RK4F]` lines. The annulus itself had 31.

| R | lines at t = 7.94 |
|---|---:|
| < 0.5 | 4343 |
| 0.5–1 | 9909 |
| 1–1.5 | 15256 |
| 1.5–1.9 | 5441 |
| 1.9–2.1 | 32 |
| 2.1–8 (the disk) | 31 |
| > 8 | 188 |

The full log has 101654 lines and `R` up to 37, which is outside the box. That tail is the runaway. It is not evidence that the orbiting ring was the source. `use_muscl` was already 1, so the cell-centred HLL overwrite was not the face pressure, and adding MUSCL on the disk is not what this log asks for.

### What to change next

The work limit delayed the cliff by about Δt = 5 relative to 406654 and did not remove it. Dual energy at `η = 0.5` delayed it by a further Δt = 0.8 and did not remove it. The cells that lose half their adiabatic energy, while the run is still quiet, are the floor gas in the hole `R < 2`, not the ring. The next Kepler run should be the same RK4 disk with that floor gas removed or frozen, not another coefficient on the pair pressure and not a wider HLL change. Leave 406515–406519 held until a disk passes more than a short fraction of an inner orbit (`T(R = 2) = 17.8`) with `dE` still on the `10^{-4}` slope. These two runs died at `t ≈ 10`, about 0.6 of that orbit.

---

## 15. Kepler 407214 started: Hopkins disk, entropy switch on

406515–406519 stay held. `GFS_DUAL_ENERGY` is unset. **[RUN]**

| job | what | where |
|---|---|---|
| 407214 `LAGFORCE_KEP_HS` | primary. `EUNHA_KEPLER_PROFILE=hopkins`, `GFS_ENTROPY_SWITCH=1`, `GFS_ES_COEF=0.01`, `GFS_FLOOR_LOG=1` | syn08, 3× H100, 4 ranks. `/gpfs/kjhan/LagForce/Kepler_gfs_rk4_hopkins_esw_Nx256` |
| 407215 `LAGFORCE_KEP_ES` | control. same switch, old `1/R` IC | syn104, 3× H200, 4 ranks. `/gpfs/kjhan/LagForce/Kepler_gfs_rk4_esw_Nx256` |

syn09 had only 2 H100 free, and syn08's other 2 H100 are in `sjjlee_h100` until 2026-10-19. The control therefore uses H200. Binary is `00adc2d` (`exam.c` and `Exam/KH/util.c` relinked; `exam_gpu.cu` unchanged).

Flags from the 407214 script and `params.dat`: `SEDOV_PHASE1=1`, `way=1`, `use_muscl=1`, `av_mode=5`, `entropy_mode=0`, `gpu_enabled=1`, centroid shift `0`, `CX=CY=12`, `EPS=0.047`, `DTETA=0.06`. The banner is

```
[ESW] entropy_switch=1 coef=0.01 half_limit=1 active=1
```

Step 1 matches the IC harness: `E_hyd = 21.459`, `E_pot = −42.928`, `E_tot = −21.469`, `E_int = 8.641e-4`. `n_es = 25528` (the estimate was about 25500). `n_half = 0`, `floor_cum = 0`, `n_ie_le0 = 0`, `es_de = −3.4e-6`. The switch is removing a little heating, not topping cells up.

407215 step 1 is the old disk: `E_hyd = 8.790`, `E_pot = −17.581`, `E_tot = −8.792`, `E_int = 8.640e-4`, `n_es = 23156`, `n_half = 0`, `es_de = −3.8e-7`. Same banner, `active = 1`.

---

## 16. Both entropy-switch Keplers died. The half-limit pays the hole

Received `c6d6a81`. `GFS_LAGUERRE_ROTATION` and `Exam/Tests3D` are in the tree. Neither has been run. 406515–406519 stay held. **[RUN]**

Exam3d is not in the places searched: `/home/kjhan/Eunha.A1`, `/home/kjhan/BACKUP/Eunha.A1`, `/home/kjhan/BackUP`, `/gpfs/kjhan/LagForce`, and the local LagEunha history (only `figs/Fig/exam3d.eps`). It was not pushed. Tests3D was not submitted. A path for Exam3d is needed before that driver can go into the repository.

### What finished

| job | IC | wall | `|ΔE| = 10^{-3}` | `|ΔE| = 10^{-2}` | `|ΔE| = 1` | end |
|---|---|---:|---:|---:|---:|---|
| 407215 `LAGFORCE_KEP_ES` | old `1/R` | 2 h 37 min | t = 13.091 | t = 14.038 | t = 14.306 | signal 9, `EXIT:255`. `sacct` `COMPLETED 0:0` |
| 407214 `LAGFORCE_KEP_HS` | Hopkins | 3 h 26 min | t = 15.383 | t = 16.110 | t = 16.134 | signal 9, `EXIT:255`. Same `sacct` |

For comparison: 406654 (no work limit) crossed `0.01` at `t = 4.59`, 406734 (work limit) at `t = 9.54`, 407021 (dual energy `η = 0.5`) at `t = 10.31`. The switch moved the crossing to `14.04` and `16.11`. It did not remove it. `T(R = 2) = 17.8`, so 407214 died at about 0.91 of an inner orbit. Not a pass.

`floor_cum` stayed 0. `E_pot` at the `dE = 1` line was still `−17.56` (407215) and `−42.92` (407214). The orbit did not fall in.

At the `0.01` crossing:

| | t | dE | es_cum | hl_cum | E_int |
|---|---:|---:|---:|---:|---:|
| 407215 | 14.04 | 1.10e-2 | +1.32e-2 | 5.45e-2 | 8.43e-2 |
| 407214 | 16.11 | 1.23e-2 | −1.89e-2 | 0.410 | 0.103 |

`dE × |Etot0|` tracks `es_cum + hl_cum`. The jump from `0.01` to `1` is `Δt ≈ 0.03`, and in that interval `hl_cum` goes from `0.054` to `5.56` (407215) and from `0.41` to `32` (407214). The half-limit is the term that runs away. On the step where 407215 crossed `0.01`, `n_half = 1` and `hl_de = 7.3e-3`. One cell. A cell's initial internal energy is about `10^{-8}`.

### Why the half-limit does that

The final correction is once per step:

```c
if (gfs_es_on && coef * m * |g| * sqrt(V) > max(0.5*ie_n, ie_E)) {
    ie1 = ie_n * (V_n/V)^(γ-1);     /* booked as es_de = ie1 - ie_E */
}
if (gfs_half_on && ie1 < 0.5*ie_n) {
    hl_de += 0.5*ie_n - ie1;        /* then ie1 = 0.5*ie_n */
}
```

`hl_de = 0.5*ie_n − ie_E` when the cell was not switched. If the energy equation has driven `ie_E` to a large negative number, the half-limit writes a cheque for the whole hole, not for half a thermal energy. The face force has already changed the kinetic energy. Putting the internal energy back does not undo that, so the total energy rises by about `|ie_E|`.

The unit test bounded the injection by asking for `−2 ie_n` every step, after which `ie` halves and the next request shrinks. A face does not ask for a multiple of `ie`. One face can deposit a macroscopic `ΔE` in one step. `hl_de = 7.3e-3` on a single cell is that deposit.

A switched cell does not take this path. A deep hole in a cold cell is booked in `es_de`, because the switch runs first and sets the adiabatic value. Early on, `es_cum` was negative: the switch was throwing away heating the energy equation wanted. That is why the cliff moved from `t ≈ 10` to `t ≈ 14–16`, and why Hopkins, with `n_es` stuck near 25500, lasted about `Δt = 2` longer than the old disk. `n_es` on 407214 was still about 25400 at the `dE = 1` line.

The switch lets a cell go when `ie_n` exceeds `0.01 m |g| h`. Those cells have only the half-limit left. The half-limit has no gravity test and no cap against the face work. `E_int` was already tens of times `8.6e-4` before `hl_cum` overtook `es_cum`. The injected energy raises the pressure, the next face work digs a deeper hole, and the half-limit pays it. That is the feedback loop, on the half-limit side.

`n_pair` reached several hundred on 407215 and about 800–1000 on 407214. The pair impulse is already limited by the internal-energy budget, so a `10^{-3}` hole in one cell is hard to blame on that term alone. The hydro face work is not limited that way.

An internal-energy rule after the step has now been tried four times. The pressure floor, dual energy, the entropy switch, and the half-limit all book the same hole under different names. The crossing moved from `t = 4.59` to `9.54`, `10.31`, `14.04`, and `16.11`. The next change has to limit the face work that digs the hole, not the value written into `ie` afterwards. `av_mode` was 5 and centroid shift was 0, so this is not the LagMFM `E_inv` alias and not the centroid volume change.

---

## 17. c283a54 방식 (1)/(2) 빌드와 GPU/MPI 짧은 실행 (2026-09-29)

**[CODE]** `origin/master`의 `c283a54`를 별도 작업 트리 `/gpfs/kjhan/LagForce/Codex_build_c283a54`에서 `USE_CUDA=1`로 빌드했다. Slurm 빌드 잡 `408460`은 A10 노드에서 `COMPLETED 0:0`이었다. `exam_gpu.cu`를 nvcc로 컴파일했고 바이너리에 `GFS_DYN_W2` 문자열이 있다. `eunha2.exe` SHA256은 `26eae35bc12cf6ef346d07fb64e24a7c9ebbe95204ebfca67a0c6ef6102367f4`이고 실행 로그의 컴파일 시각은 `12:08:08 Sep 29 2026`이다. 아래 모든 실행에 이 바이너리를 사용했다.

**[RUN]** 256² 본실험 A(방식 (1), 보로노이+HLL) 잡 `408511`을 요청대로 `--partition=h200 --gres=gpu:H200:3`, 4랭크, 32G로 제출했다. `syn104`에 배정됐지만 2026-09-29 13:48:42 KST 시작 후 6초 만에 `FAILED 1:0`으로 끝났다. 4개 랭크 모두 `exam_gpu.cu:1311: initialization error`를 출력했고 첫 `Time=` 또는 `[RK4E]` 줄이 없다. Slurm 출력은 H200 NVL 세 장을 열거한다. 이전 노드 상태에는 Xid 119/154와 재부팅 필요 사유가 표시됐었다. 따라서 이번 실패는 수치 발산 판정이 아니라 CUDA 초기화 실패다. 결과 디렉터리는 `/gpfs/kjhan/LagForce/Kepler_gfs_rk4_S1_voronoi_hll_Nx256`이다. B1/B2의 256² 본실험은 제출하지 않았다.

**[RUN]** H200 장애와 독립적으로 GPU+MPI 경로를 확인하려고 A10 세 장, 4랭크, 64², `HYDRO_TSTOP=3`으로 동일한 세 설정을 순차 실행했다. Slurm 잡 `408526`은 `COMPLETED 0:0`, S1/B1/B2 각각 `EXIT:0`이었다. 파라미터는 256² A의 406734 기준을 64²로 낮춘 뒤 B1/B2에 `GAS kappa=0.4`, `GAS w2_mode=8`, `GAS w2 power=1`을 넣었다. B1/B2에는 `GFS_ALLOW_KAPPA=1 GFS_LAGUERRE_ROTATION=1`, B2에만 `GFS_W2_LINEAR=1`을 켰다. kappa 유지, LAGROT, B2 W2LIN 활성 배너를 확인했다. 세 설정 모두 RK4, MUSCL, `av_mode=5`, 옛 1/R Kepler IC, `gpu_enabled=1`이다. 로그는 `/gpfs/kjhan/LagForce/Kepler_gpu_mpi_pf64_c283a54/{S1,B1,B2}/log`에 있다.

| 잡 | 방식 | 해상도 | 플래그 | 도달 t | dE t≈1 | dE t≈3 | dE t=10 | dE t=17.8 | floor_cum @t≈3 | sfl_cum @t≈3 | 최소 dt (t) | 첫 \|dE\|>1e-3 | [VAUD] cum |
|---|---|---|---|---:|---:|---:|---|---|---:|---:|---|---|---|
| 408526/S1 | (1) | 64² | κ=0 | 3.00133 | 7.440e-7 | 1.074e-5 | 미도달 | 미도달 | 8.279e-6 | 0 | 0.00158954 (3.00133) | 없음, t≤3 | — |
| 408526/B1 | (2) | 64² | κ=.4, mode8, a=1, LAGROT | 3.00143 | 7.375e-7 | 1.053e-5 | 미도달 | 미도달 | 9.550e-6 | 0 | 0.00156026 (3.00143) | 없음, t≤3 | — |
| 408526/B2 | (2) | 64² | B1+W2LIN | 3.00086 | 7.153e-7 | 9.895e-6 | 미도달 | 미도달 | 7.707e-6 | 0 | 0.00149890 (3.00086) | 없음, t≤3 | — |

`dE`는 `[RK4E] dEtot/|Etot0|`의 부호 있는 값이다. t≈3에서 `n_pair`는 S1=8, B1=12, B2=4다. `|Etot0|≈8.792664`이므로 t≈3 에너지 드리프트 절댓값은 약 `8.7–9.4e-5`이며 장부의 `floor_cum`은 `7.7–9.6e-6`, `sfl_cum`은 0이다. **[HYP]** floor 장부만으로 이 드리프트를 설명할 수 없다. 여기서 원인을 특정하지 않는다.

**[RUN]** 같은 바이너리/IC/64²에서 `gpu_enabled=0`(CPU hydro), 1랭크, `OMP_NUM_THREADS=1`로 세 설정을 순차 실행한 잡 `408534`도 `COMPLETED 0:0`이었다. CPU 실행에도 클러스터 규칙에 따라 A10 한 장을 요청했고 로그에는 GPU 중력 트리가 나온다. 따라서 여기서 CPU는 hydro 경로를 뜻한다. 로그는 `/gpfs/kjhan/LagForce/Kepler_cpu1_pf64_c283a54/{S1,B1,B2}/log`에 있다.

| 잡 | 방식 | 도달 t | dE t≈1 | dE t≈3 | floor_cum @t≈3 | sfl_cum @t≈3 | 최소 dt (t) |
|---|---|---:|---:|---:|---:|---:|---|
| 408534/S1 | (1) | 3.00067 | 9.516e-11 | 1.786e-4 | 6.658e-6 | 1.558e-3 | 2.19269e-5 (2.09131) |
| 408534/B1 | (2) | 3.00157 | 4.286e-8 | 1.379e-6 | 0 | 1.235e-5 | 0.00170960 (3.00157) |
| 408534/B2 | (2)+W2LIN | 3.00100 | 1.895e-10 | 7.423e-6 | 0 | 6.463e-5 | 0.00106728 (2.94661) |

세 설정 모두 `EXIT:0`, t≤3에서 `|dE|>1e-3` 없음. CPU 1랭크의 S1/B1 값은 Grokbot의 박스 CPU 보고 `(1)≈1.5e-4, (2)≈1.5e-6`을 재현한다. B2는 B1보다 나쁘다. `Exam/exam.c`의 CPU RK 단계는 음의 압력이면 내부에너지를 재설정하고 `sfl_cum`에 기록한다. GPU RK 단계는 압력만 바닥값으로 만든다. CPU S1의 `dE×|Etot0|≈1.57e-3`은 `sfl_cum+floor_cum≈1.565e-3`과 가깝다. CPU B1/B2도 에너지 증가분이 주로 `sfl_cum`에 기록된다.

**[RUN]** B1의 차이가 hydro 경로와 랭크 수 중 어디서 생기는지 보려고 같은 입력으로 교차 실험했다. `408553`(GPU hydro, 1랭크, A10 1장), `408554`(CPU hydro, 4랭크, A10 1장), `408564`(GPU hydro, 4랭크, A10 1장)는 모두 `COMPLETED 0:0`, `EXIT:0`이었다. GPU 설정에서는 `OMP_NUM_THREADS=2`, CPU 설정에서는 1로 원래 비교 잡과 일치시켰다. 로그는 `/gpfs/kjhan/LagForce/Kepler_cross_pf64_c283a54/{GPU1,CPU4,GPU4_ONEGPU}/B1/log`에 있다.

| B1 경로 | 잡 | GPU 수 | dE t≈1 | dE t≈3 | floor_cum @t≈3 | sfl_cum @t≈3 |
|---|---|---:|---:|---:|---:|---:|
| CPU hydro, 1랭크 | 408534/B1 | 1 | 4.286e-8 | 1.379e-6 | 0 | 1.235e-5 |
| CPU hydro, 4랭크 | 408554 | 1 | 4.290e-8 | 1.585e-6 | 1.392e-7 | 1.354e-5 |
| GPU hydro, 1랭크 | 408553 | 1 | 3.000e-9 | 1.600e-6 | 1.085e-5 | 0 |
| GPU hydro, 4랭크 | 408564 | 1 | 7.375e-7 | 1.048e-5 | 9.571e-6 | 0 |
| GPU hydro, 4랭크 | 408526/B1 | 3 | 7.375e-7 | 1.053e-5 | 9.550e-6 | 0 |

**[HYP]** B1의 추가 오차는 GPU hydro와 4랭크가 함께 있을 때만 나타난다. GPU 한 장과 세 장의 결과가 사실상 같으므로 장치 개수/랭크 매핑이 주원인이라는 증거는 없다. t≈1에는 모든 경로에서 `floor_cum≈0`이고 GPU 4랭크만 이미 `dE=7.375e-7`이다. t≈3의 GPU 4랭크 에너지 증가분은 약 `9.2e-5`, 기록된 `floor_cum`은 약 `9.6e-6`이므로 장부에 없는 차이가 약 `8.3e-5`다. 이 비교만으로 특정 면 또는 MPI 통신 지점을 단정하지 않는다.

**Grokbot께 요청:** GPU hydro 4랭크의 초기(t≈1, 바닥값 주입 전) 에너지 차이를 출발점으로 GPU 면 힘/일의 랭크 경계 합산과 ghost/owner 처리를 점검해 달라. 같은 B1 64² 입력에서 1랭크와 4랭크의 RK 단계별 에너지 및 면별 짝 일 합계를 비교하는 계측이 가장 직접적이다. CPU 경로의 `sfl_cum`과 GPU 경로의 `floor_cum`은 적용 단계가 달라 이 둘을 그대로 같은 물리 효과로 취급하지 말아 달라. H200 `408511`의 `exam_gpu.cu:1311`은 `cudaGetDeviceCount` 호출이므로 GPU 초기화/노드 상태 문제로 분리해 달라.

**[RUN]** `406515–406519`는 계속 보류한다. 짧은 t=3 결과로 생산 잡 해제 여부를 판정할 수 없다. H200의 CUDA 초기화가 복구되거나 별도 자원이 승인되면 256² A에서 적어도 t=17.8의 에너지 추이를 확인해야 한다. 과거 t≈16.11의 급격한 발산을 고려하면 그 뒤까지 안정성을 살펴야 한다.
