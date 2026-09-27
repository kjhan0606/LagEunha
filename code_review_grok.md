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
