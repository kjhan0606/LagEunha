# LagEunha GFS hydro: code review of error sources and physical limitations

**Scope.** The geometric-face-scheme (GFS) hydro path used for the paper's 2D tests (`Exam/exam.c`, `Voro/voro_eunha.c`, the test drivers in `Exam/{RT,KH,Cylinder}`) and the 1D suite (`Hydro1DExam/laguerre_sod.c`).
**Sources.**
- Code: `kjhan0606/LagEunha`, branch `master`, commit `75aaeedc52a2183077ad45f73f4a1c1471b11a6d` (2026-09-26 16:29 KST). Cloned read-only to `/workspace/lagEunha/code`. Nothing was written to GitHub.
- Paper: `/workspace/lagEunha/overleaf/main.tex` (Section 3.4 and Section 5).
**Evidence labels.**
- **[CONFIRMED-CODE]** read directly in the source.
- **[CONFIRMED-RUN]** reproduced by running code.
- **[HYPOTHESIS]** plausible but not verified.

**What was run** (about 15 minutes total; the 2D MPI/FFTW build was not attempted):
- `Hydro1DExam/laguerre_sod.c`, built with gcc from a scratch copy in `/tmp`. I added an instrumented copy that counts internal-energy floor events: `review_tests/laguerre_sod_diag.c`.
- A Python/scipy check of Voronoi face kinematics: `review_tests/voronoi_face_velocity_check.py`.

---

## 0. Executive summary (most important first)

1. **GFS is not one well-defined code path.** It is spread across `av_mode=5` + `use_muscl=1` + the environment variable `SEDOV_PHASE1=1` + the choice of integrator. With `SEDOV_PHASE1=1`, the per-particle `die` means **dE_total/dt**. Only the KDK integrator knows this. The RK4 integrator, which the paper says GFS uses, would add the total-energy rate to `ie` and also kick `v`, so every ΔKE is counted twice. Without the environment variable, the extreme-face HLL fallback is silently replaced by a non-physical "Fix C" (P* = min(P_L,P_R)). **[CONFIRMED-CODE]**
2. **In the 1D blast, all of the energy error is internal-energy floor clipping.** It is not truncation error. 0.900 of 275.02 energy units (3.27e-3) is injected by `e<1e-14 → 1e-14`, which matches the measured error of 3.29e-3 exactly. The "hllface" variant (Laguerre face moved at the HLL speed, weights reconstructed, which is your original design in 1D) has no floor events and an energy error of 1.1e-6. The "wtrack" variant reproduces the paper's blast numbers (L1 = 0.154, ΔE = 3.8e-2), and 100% of its error is also floor energy. **[CONFIRMED-RUN]**
3. **On ordinary faces, the P dV face velocity is not the true Voronoi face velocity.** The code uses (u_p+u_q)/2 for the whole face. The true normal velocity of a Voronoi face has a rotation term −(u_q−u_p)·(c_f−m_pq)/d. Measured on a perturbed lattice, the mismatch in dV/dt is 50% of the true rate in shear flow and 31% for grid-scale noise. With the rotation term included, the match is exact to 1e-10. This does not break total-energy conservation, but it makes internal energy inconsistent with the geometric volume. That produces spurious heating or cooling in shear (KH, Gresho, Kepler) and feeds positivity fixes. **[CONFIRMED-RUN]**
4. **The KDK time step ignores all local stability limits.** It uses Δt = C·(L_x/N_x)/max(v_sig), with the initial uniform spacing. The per-face CFL, acceleration CFL and dK CFL computed in the force routine are thrown away (`dt_dummy`). In compressed gas the local Courant number is too large by about √(ρ/ρ̄): about 3× in the cylinder bow shock and 4–5× at the Noh centre. The paper's cylinder time step (≃8e-3) equals this global formula. **[CONFIRMED-CODE, strong circumstantial evidence from paper numbers]**
5. **Extreme faces break geometric consistency by design.** On extreme faces the P dV uses the HLL normal speed, while the Voronoi face moves at the mid-plane speed. Cells then do P dV work that does not correspond to any volume change. Near cold gas this drives `ie` negative, and the result is clipped or floored, which creates energy and breaks momentum conservation. In 1D the same situation produces particle overtaking (a generator outside its own cell) and a 20,000-step Δt collapse. **[CONFIRMED-RUN (1D), CONFIRMED-CODE (2D)]**
6. **Boundary and postStage operations are non-conservative and badly placed inside RK4.**
   - Wall clamps zero the velocity.
   - The cylinder wall removes the inward normal velocity and reflects positions without any P dV.
   - These run between RK4 stages, and the "undo K4 shift" arithmetic then subtracts from an altered state.
   - The cylinder has an open inflow/outflow boundary with pulsed injection, so its "−1.0% energy" is not a conservation diagnostic.
   **[CONFIRMED-CODE]**
7. **Scale and unit dependence.**
   - Absolute thresholds: P<1e-3 for "extreme", P floor 1e-6, 1e-10 floors, the Fix C guard 0.01.
   - A dimensionally inconsistent time-step criterion: Δt₃ = 0.1·d/|Δv/2|², which has units of time²/length.
   Results therefore change if the problem is rescaled. **[CONFIRMED-CODE]**
8. **Independent per-cell tessellation can make face terms asymmetric.** Each particle rebuilds its own Voronoi cell from a limited neighbour stencil (3×3 cells, with an adaptive radius up to R=4 and a "band-aid fallback"). If the cells of p and q disagree about the shared face, then F_pq ≠ −F_qp and the pairwise energy exchange does not cancel. **[HYPOTHESIS; cheap to test]**
9. **The paper's numbers and the current code disagree.** The 1D `vr` driver, which writes `vr_*.dat` and is presumably the "one-dimensional problems" figure, uses CD10 artificial viscosity and no extreme-face HLL, and it crashes on the blast. Its L1 values differ from the paper (Sod 5.63e-3 vs 5.86e-3; Shu–Osher 0.780 vs 0.571; Noh 1.90e-2 vs 1.36e-2). The paper numbers were likely produced by uncommitted settings or drivers. **[CONFIRMED-RUN]**

---

## 1. Findings in detail

### 1.1 GFS depends on an environment variable, and the meaning of `die` changes with it (Critical)

**What.**
- `sedov_phase1_on()` (`Exam/exam.c` ~l.4072) reads `getenv("SEDOV_PHASE1")`. It gates three things:
  - the extreme-face HLL fallback in `hllc_face_2d` (l.4208),
  - the Riemann-speed P dV on extreme faces (l.5080),
  - redefining the particle rate at l.5605–5608 as `ibp_rk4->die = dte` (total energy), with `stress.dK = 0`.
- The KDK integrator `exam2d_vph_kdk_int_blend` handles this through `phase1_half_kick` (l.4130–4160): E = ie + KE + die·Δt/2, then ie = E − KE.
- The RK4 integrator `exam2d_vph_rk4_int_blend` (l.9223–9595) always does `k*ie = die*Dtime` and `k*v = a*Dtime` (l.9322, 9373, 9424, 9475). It has no phase-1 branch.

**Mechanism.** With RK4 and `SEDOV_PHASE1=1`, each step gives ΔE_total ≈ 2·ΔKE + (work). All kinetic–thermal exchange is counted twice.
- Where KE barely changes (KH, Gresho) this is invisible.
- In RT, where potential energy turns into a small KE, it would give ΔE/E of order KE/E. **[HYPOTHESIS]** The paper's 7.6e-4 could have this origin.
- In Noh, where all KE is converted, it would be O(1). The Noh error is therefore probably not from RK4 + phase 1, which suggests the runs used KDK.

The paper says GFS uses fourth-order RK. Commit `51cc2f1` says "The cylinder uses the same KDK integrator". Without the environment variable, `hllc_face_2d` instead applies "Fix C v2" (l.4214–4225): if P_min > 0.01 and P_max > 100 P_min, then P* = P_min and v* = mean. That is not a Riemann solution, and it underestimates the face pressure at strong interior jumps.

**Fix.**
- Make GFS an explicit, named parameter (for example `GAS_SCHEME=GFS`) and remove the environment gating.
- Carry one energy variable consistently in every integrator. Either integrate E_total per particle (recommended), or integrate `ie` with die_int = dte − F·v. Never switch meanings.
- Log the run configuration (integrator, av_mode, MUSCL, environment variables) into every output header.

**Test.** Run RT with `EVOLMETHOD=1` and `SEDOV_PHASE1=1`, and log Σ(ie+KE) per step. It should drift at about the rate of ΣΔKE.

### 1.2 Energy error comes from positivity fixes, not from time integration (High)

**What.**
- 1D (`laguerre_sod.c` l.916, 946): `if(P[i].e<1e-14) P[i].e=1e-14`. This runs in every RK stage and at the end of the step.
- 2D final update (`exam.c` l.9564–9566, and the same at 9891, 10138): if pressure ≤ 0, set P = 1e-6 and rebuild `ie`. This creates energy.
- KDK phase 1 (`phase1_half_kick`, l.4138–4152): if E < KE, velocities are rescaled so KE = E and ie = 0. Energy is kept, but **momentum is not**, and the next step's P = 0 then hits the 1e-6 floor.

**Evidence** (1D blast, N=200, GFS-like settings: Voronoi, MUSCL + HLLC, `extreme_hll=1`, no AV). This is the `voronoi` row of `./laguerre_sod blast 200`. E₀ = 275.02.

| variant (all w=0 initially) | L1_ρ | ΔE/E₀ | floor events | floor energy / E₀ |
|---|---|---|---|---|
| voronoi (GFS-like) | 0.0717 | 3.29e-3 | 732 | 3.27e-3 (≈100%) |
| hllface (face moved at HLL v*, weights rebuilt) | 0.0874 | 1.1e-6 | 0 | 0 |
| wtrack (RK-tracked HLL faces) | 0.1545 | 3.77e-2 | 4207 | 3.77e-2 (≈100%) |
| wsmooth | 0.0689 | 3.29e-3 | 594 | 3.27e-3 |

- Reducing the CFL from 0.3 to 0.1 and 0.03 does **not** reduce the error (3.37e-3 and 3.36e-3). Those runs also failed before t_end (step cap reached, ρ_max up to 1.7e4). The error and the failures are therefore not temporal truncation.
- Only 183 of the 732 floor events are next to a face that is extreme at the time of the event (instantaneous criterion). The rest follow particle overtaking and tiny cells (Section 1.5).

**Fix.**
- Do not floor. Detect negative internal energy and treat it as a failure of the face-velocity or pressure closure.
- A dual-energy switch (entropy for cold, supersonic cells, with a total-energy remainder) is the standard remedy.
- If a floor is unavoidable, take the deficit from the neighbours' total energy and report it.
- Replace the KDK velocity rescaling with an energy-consistent pressure limiter at the face.

### 1.3 The ordinary-face P dV velocity omits the Voronoi face rotation term (High)

**What.** With w=0, `get2dUpqradRk4` (`Voro/voro_eunha.c` l.69–104) returns u_face − u_p = ½(u_q−u_p). The face velocity is then the midpoint velocity ū, used for the whole face in `die_rev_face = −p·(uradix·dS)` (`exam.c` l.5463) and in `dte` (l.5510). The exact normal velocity of the Voronoi bisector at a point x on the face is

  w_n(x) = n·ū − (u_q−u_p)·(x−m_pq)/d_pq.

Integrated over the face this gives L[n·ū − (u_q−u_p)·(c_f−m_pq)/d]. The second term vanishes only if the face centroid c_f sits at the pair midpoint m_pq or the relative velocity is purely radial. This is the same correction Arepo uses (Springel 2010, face-velocity correction for Voronoi faces).

**Evidence** (`review_tests/voronoi_face_velocity_check.py`: 1024 points, perturbed lattice at 30% of Δx, periodic, finite-difference dA/dt):

| velocity field | rms(dA/dt) | rms(mean-face − dA/dt) | rms(exact − dA/dt) |
|---|---|---|---|
| uniform translation | ~0 | ~0 | ~0 |
| shear u = sin 2πy | 8.0e-4 | 4.0e-4 (50%) | 6.5e-11 |
| compression u = −sin 2πx | 4.4e-3 | 1.9e-4 (4%) | 6.7e-11 |
| grid-scale noise | 2.9e-3 | 9.1e-4 (31%) | 7.4e-11 |

**Mechanism.**
- Total energy is still pairwise conserved, because the same face velocity enters both cells.
- Internal energy, however, is updated with a P dV that differs from the actual geometric dV, and P is then recomputed from ie and V_geo. The result is spurious heating or cooling of either sign, of order P·|∇u|·Δx relative error per step, largest in shear.
- This is a direct source of pressure noise in KH, Gresho and Kepler. In the Keplerian ring it fits the reported secular spreading and the 12 slow particles at R=1.89 **[HYPOTHESIS]**.

**Fix.** Use the face centroid in `get2dUpqradRk4`, which is cheap and keeps pairwise antisymmetry:

  u_face = ū − [(u_q−u_p)·(c_f−m)/d] n.

The Laguerre case (w≠0) needs the same term.

### 1.4 The KDK time step ignores local CFL, acceleration and entropy limits (High)

**What.** `exam2d_vph_kdk_int_blend` (l.9697–9715) sets Dtime = Courant·dx_uniform/max_global(vsig_max), with dx_uniform = L_x/N_x fixed. The Dtime returned by `getAccVoro2DBlend` includes:
- the per-face limit 2C·max(d, ¼√V)/v_sig,
- the viscous limit,
- the acceleration limit (`ACC_CFL_FRAC`),
- the dK limit.

That returned value is stored in `dt_dummy` and discarded (l.9695, 9769).

**Mechanism.** The cell size scales as (m/ρ)^{1/2} in 2D, so the local Courant number exceeds the nominal value by about √(ρ/ρ̄):
- cylinder, ρ_max = 9.4: ≈3×
- Noh, ρ ≈ 16–30: ≈4–5.5×

This is a direct route to overshoot, particle interpenetration and positivity failures.

**Evidence.**
- The paper's cylinder Δt ≃ 8e-3 with C = 0.3 and Δx = 0.05 implies v_sig ≈ 1.9, which is exactly the global formula.
- The per-face limit in the bow shock would be about 3× smaller.

**Fix.** Use Δt = min(global, local) = min over faces and cells of the limits already computed, and apply it within the step.

### 1.5 Extreme-face HLL speed versus geometry: negative ie and particle overtaking (High)

**What.**
- In `av_mode=5` with phase 1, an extreme face (P_min<1e-3 or P_max>100·P_min, from reconstructed states, l.5076–5086) uses the HLL normal speed v* in `dte` (l.5500–5505).
- The internal-energy part is therefore −P*·(v*−u_p)·dS.
- The geometric face, however, moves at ū.

**Mechanism.**
- For the cold cell next to a strong shock, (v*−ū) is O(Δv), and P* is huge compared with that cell's own P.
- The cell loses or gains P*·(v*−ū)·A·Δt of internal energy without a matching volume change. With a tiny initial `ie`, this goes negative within one step and is then clipped or floored (Section 1.2).
- In 1D the same inconsistency lets a particle overtake its neighbour. In the diagnostic dump at step 40, x₂₀ = 0.11310 > x₂₁ = 0.11216. Particle 20 then lies outside its own cell [0.10996, 0.11263], and cell 21 shrinks to V = 9.8e-5 with P = 3.8e3.
- From t ≈ 0.0304 to 0.0307 the run takes about 20,000 steps at Δt ≈ 4e-9, with cells of V ~ 1e-7 and P ~ 5e4.

**Fix.** Make geometry and face velocity the same thing (Section 2). As a stop-gap, bound the extreme-face work so that no cell's ie drops below a fraction of its value in a step, and redistribute the remainder symmetrically.

### 1.6 Boundary operations inside RK stages, and non-conservative wall and cylinder treatments (Medium–High)

**What.**
- The RK4 integrator calls `postStage` after every stage and then "undoes" the K4 shift with `x -= k3x; v -= k3v` (l.9480–9488).
- For walls (`wallx_postStage_blend` l.9136–9152, `walls_xy_postStage_blend` l.9095–9129), `postStage` reflects positions, flips velocities, and **clamps** particles to half a grid spacing from the wall with v_n = 0. The undo arithmetic then subtracts from a modified state, so the RK4 combination is not RK4 and energy is lost.
- `cyl_postStage` (`Exam/Cylinder/cylinder.c` l.154–165):
  - reflects the position (r → 2R−r), which changes the volume without P dV;
  - removes the inward normal velocity, destroying ½m v_n² without putting it into `ie`.
- The cylinder domain is open in x: injection when `inject_accum ≥ dmean` (pulsed, l.~335–395) and outflow deletion (`cyl_outflowCleanup`). The paper's "energy differs by −1.0%" does not subtract inflow and outflow fluxes, so it is not a conservation measure.
- Non-target particles get v = 0 at the end of each step (l.9572).

**Fix.**
- Apply boundary operations only at full steps.
- Convert any removed kinetic energy into internal energy of the same particle.
- Report an energy budget with boundary fluxes, as the 1D code already does through `bc_eng_net`, including injected energy.
- Prefer mirror-ghost Riemann walls over position clamps everywhere. The flat-wall ghost path still uses M(n,m) plus Monaghan AV (l.~4847–4875), not GFS.

### 1.7 Pairwise antisymmetry is not guaranteed with independent cell construction (Medium, HYPOTHESIS)

**What.** Each particle builds its own cell (`Voro2D_FindVC`, l.~4630) from the 3×3-cell neighbour list, retrying up to R=4, then falling back with a "band-aid" message (l.~4650–4680). Forces and energy are gathered only onto `ibp`.

**Mechanism.**
- If the two cells disagree about the existence, length or normal of the shared face, the face terms do not cancel. Causes include a missing neighbour in one list, near-cocircular points, or roundoff.
- MUSCL states are also evaluated at the face midpoint expressed in each generator's local frame.

**Test.** Log |ΣF_i| and Σ dE_i/dt over all particles each step. Both should be at roundoff level. Also count faces seen by only one side.

**Fix.** Build a global Delaunay or Voronoi structure per stage, or compute each face once (ordered pair p<q) and scatter ±.

### 1.8 Scale and unit dependence and dimensional errors (Medium)

- **Δt₃** (l.5545–5551): 0.1·d/|Δv/2|² has units time²/length. It is dimensionally wrong and changes with the units. It should be something like 0.1·d/|Δv|.
- **Absolute thresholds:**
  - "extreme" if P_min < 1e-3 (l.5078, 4208);
  - Fix C guard 0.01 (l.4176);
  - pressure floor 1e-6;
  - floors of 1e-10 on reconstructed ρ and P (l.~5062).
  Use ratios to a local or reference scale, for example P_min < ε·P_max or ε·ρc².
- The extreme criterion uses reconstructed P_L and P_R. Different limiter outcomes can flip the classification between stages, which gives time-discontinuous forcing inside RK4.

### 1.9 Other issues (Low–Medium)

- **Centroid steering** (`av_mode=5 && GAS_FCENTROID>0`, l.9522 and KDK path) moves generators after the step with no remap or P dV. V, ρ and P jump while `ie` stays fixed. It is not Lagrangian and it changes entropy.
- **KDK second kick** evaluates velocity-dependent HLLC pressures at v_{n+1/2}, not at a predicted v_{n+1}. This is first-order in the dissipative part.
- **Sound speed in HLLC** (`av_mode=5`) uses cell values, not the reconstructed states, so wave-speed estimates are slightly inconsistent.
- **`get2dUpqradRk4`** divides by `dtold`. With nonzero weights and `dtold=0` on the first step this gives NaN.
- **Weight time integration** (nonzero-weight case, see Section 2b):
  - w² is reset to its start-of-step value at each RK stage (l.9340, 9391, 9442, 9494).
  - The face speed uses a *lagged* finite difference (w−w_old)/dt_old clamped to c_s (`voro_eunha.c` l.87–96).
  - w² is then changed at the end of the step (`applyW2Controls`, l.9543).

---

## 2. Your design intent, "Voronoi + HLL = Laguerre with designed weights"

### (a) How much of the energy error comes from the geometry/face-velocity mismatch?

- **Directly, none of the total-energy error.** In both the 1D code and the 2D phase-1 KDK path, the same face velocity enters both neighbours. Σ dE/dt then cancels pairwise for *any* face velocity, whether mid-plane, HLL or wrong. What the mismatch breaks is thermodynamic consistency: internal energy is updated as if the cell volume changed at the rate implied by the face velocities, while density and pressure are recomputed from the re-tessellated volume.
- **Indirectly, most of it.** The inconsistency drives ie < 0 in cold cells next to extreme faces, and in 1D it drives particle overtaking. The energy error then enters through the positivity fixes (floors and clips). In the 1D blast, 100% of the measured error is floor energy (Section 1.2).
- **When face velocity equals the actual face motion, the floors disappear.** "hllface" moves the Laguerre faces at v* and rebuilds the weights. It is exact in 1D and has zero floor events and ΔE = 1.1e-6, compared with 3.3e-3 for GFS-like Voronoi. The limitation is that its faces had to be clamped 17,185 times to stay between their particles, and that clamping is a geometric change with no P dV.
- **The 2D phase-1 path should behave the same way** **[HYPOTHESIS]**. The paper's RT comparison supports this: the contact speed on every face gave a 10⁷ energy increase, which a pairwise-conservative scheme can only produce through such fixes.
- **Ordinary faces add a separate, always-present mismatch.** Even on ordinary faces the mid-plane velocity is not the true Voronoi face velocity (Section 1.3: 50% of dV/dt in shear). This is a second, smaller but pervasive source of spurious entropy.

### (b) Why the nonzero-weight attempts were unstable in 2D

1. **Topological over-determination** (the paper's Section 3.4 argument, confirmed). For a power diagram, the offset of face pq from the midpoint is s_pq = (W_p−W_q)/(2d_pq), with W = w². The face kinematics give

   Ẇ_p − Ẇ_q = 2d_pq (v*_pq − n·ū) + (W_p−W_q)(n·u_qp)/d_pq ≡ b_pq.

   In 2D there are about 3N faces but only N−1 independent Ẇ. The system B·Ẇ = b is solvable only if b sums to zero around every Delaunay triangle, which is about 2N cycle constraints. Riemann speeds do not satisfy this. Geometrically, three faces must keep meeting at one Voronoi vertex, so independent face offsets are impossible. In 1D the face graph is a chain with no cycles, which is why "hllface" works there.
2. **The weights in the code are not solutions of that system.** `getw2forHydroParticle` (`exam.c` l.68–190) sets w² from local state:
   - P^{(γ−1)/γ} (mode 0),
   - entropy (1, 2),
   - pressure × volume equalization (3),
   - c_s (4–6),
   - dK (7).
   At a shock these give large ΔW over one face, so the face jumps toward a particle and cells degenerate (hence `w2ceil`, rate limiters and floors in `applyW2Controls`). Nothing aims the face at v*.
3. **The weight time integration is inconsistent with the geometry.** Within an RK4 step W is frozen (reset each stage), but the P dV uses a lagged, clamped Ẇ from the previous step. W is then changed discontinuously at the end of the step. The faces move without P dV at the jump, and the P dV includes motion that did not happen during the step. This is exactly the mismatch of (a), now on every face, and it grows with the weight variation that shocks produce. It fits the Sedov blow-up "within a few tens of steps".
4. **The rotation term is missing**, as in the zero-weight case (Section 1.3).

### (c) Can a weight construction give "Laguerre face velocity = HLL speed" consistently?

Exactly, no, because of (b)1. Approximately and consistently, **yes, if the realised face motion is what enters P dV**:

1. **Target.** At each stage compute b_pq from the desired speed v*_pq, which is the HLL speed on extreme faces. On ordinary faces use the mid-plane speed, so b = 0 up to the ΔW term.
2. **Weighted least squares for Ẇ.** Minimise Σ_f a_f (Ẇ_p−Ẇ_q−b_f)². This is a graph-Laplacian solve, L_a Ẇ = Bᵀ A b: sparse, symmetric positive definite, a few conjugate-gradient iterations per stage, with Ẇ fixed up to a constant. Take a_f large on extreme faces and small on ordinary faces, so the unavoidable residual r_f is pushed into ordinary faces where it is smallest. The residual is a measurable, tunable modelling error, not a conservation error.
3. **Use the realised face velocity in the energy equation**, including the rotation term: w_n = n·ū + Ṡ_pq − (u_q−u_p)·(c_f−m)/d with the solved Ẇ. Do not use v*. Geometry and P dV then agree identically, total energy stays pairwise conserved, and ie can no longer go negative through bookkeeping.
4. **Integrate W as a state variable** in every RK or KDK stage (W ← W + Ẇ·Δt at each stage), not by finite differences between steps and not with end-of-step jumps.
5. **Keep cells valid.** Require |W_p−W_q| < (1−δ)d²_pq for all neighbours. This keeps each generator inside its own cell and each cell non-empty. Enforce it with a Δt limit, Δt ≤ ((1−δ)d² − |ΔW|)/|ΔẆ|, or with bound-constrained least squares (active-set). The global shift in W is free, so W ≥ 0 is trivial.
6. **Cost and benefit.** One extra sparse solve per stage. Around a shock front the extreme faces form a chain, and the cycle constraints close through ordinary faces, so the error is spread over the neighbouring faces rather than concentrated in the cold cell.

**Alternative worth knowing.** Cell-centred Lagrangian schemes with nodal Riemann solvers (GLACE: Després & Mazeran 2005; EUCCLHYD: Maire et al. 2007) solve the same consistency problem at the vertices. The vertex velocities come from a nodal solver, the geometric conservation law holds exactly, and total energy is conserved. They are the mature form of "geometry moves at the Riemann velocity". They give up the Voronoi or Laguerre structure and need rezoning (ALE) in shear.

---

## 3. Fundamental physical limitations of the approach

- **Pure Lagrangian mesh in shear.** A Voronoi mesh is always valid, but cells elongate and neighbour lists change quickly. Without explicit mass exchange there is no sub-cell mixing, so contact interfaces stay sharp forever. That is good for contacts, but KH and RT mixing layers cannot mix at the molecular or sub-grid level: entropy of mixing is zero, and metal mixing in galaxies would need a separate sub-grid diffusion model.
- **Resolution follows mass.** Low-density regions (voids, hot haloes, the bow-shock stand-off and the cylinder wake) are poorly resolved. The two-dimensional Noh corner particles with V = 432 are an example. Refinement or derefinement (splitting and merging with conservative redistribution) is needed for cosmology.
- **The Galilean invariance of the face velocity choice is the key trade-off.** The mid-plane velocity is Galilean invariant and geometric. The HLL speed is invariant only if computed in the face frame (`hllc_face_2d_rest_frame` does this). Mixing them breaks consistency (Section 2).
- **Dissipation.** A Riemann solver on a Lagrangian contact gives no dissipation at contacts. The code comment at l.~4878 notes that pure HLLC failed on KH. Shocks still need entropy generation. In GFS it comes only from the HLLC/HLL star pressure. At strong shocks with P*·(v*−ū) mismatches, that is where negative `ie` appears. Some explicit, energy-consistent artificial viscosity or a proper entropy fix at shocks remains necessary.
- **Particle overtaking and mesh tangling.** In 1D the generator order is not enforced, and in the blast it is violated. In 2D the Voronoi mesh cannot tangle, but generators can come arbitrarily close, which gives sliver faces, tiny d_pq in the time step (hence the ¼√V floor), and noisy gradients. Mesh regularisation (steering) is required for long runs and must be accounted for (Section 1.9).
- **Time integration.** Explicit RK4 and KDK do not conserve energy for the nonlinear v² term. That error is small compared with the fixes found here, but it is nonzero. Only a total-energy formulation with pairwise fluxes makes the integrator conservative by construction.

---

## 4. Suggested order of work

1. Make GFS explicit, remove the environment gating, and use one energy variable (E_total) in all integrators. Log the configuration in outputs.
2. Use local and global minimum time steps in KDK. Move boundary operations to full steps.
3. Add the Voronoi face rotation term to the face velocity.
4. Remove floors and clips in favour of a dual-energy or positivity-preserving face limiter, and count and report every fix.
5. Build the least-squares "designed weights" of Section 2c, integrated as a state variable, with realised face motion in P dV.
6. Add diagnostics each step: ΣF, Σ dE, floor count and energy, one-sided faces, the geometric-versus-flux volume residual Σ|ΔV_geo − Σ face flux·Δt|, and a boundary energy budget.
7. Regenerate the paper's tables with committed parameter files.

---

## Appendix: reproduction commands

```
# 1D (scratch copy; the repo was not modified)
cp -r /workspace/lagEunha/code/Hydro1DExam /tmp/h1d && cd /tmp/h1d
gcc -O2 -o lsod laguerre_sod.c -lm
./lsod vr 200          # current 'vr' driver (CD10 AV, no extreme HLL); blast stops at dt=1e-12
./lsod blast 200       # voronoi / laguerre / hllface / material / wtrack / wsmooth / pgate
gcc -O2 -o lsod_d /workspace/lagEunha/review_tests/laguerre_sod_diag.c -lm
./lsod_d blast 200 voronoi | grep floor   # floor_n=732 (183 next to an extreme face), floor_eng=0.900 of E0=275.02
# 2D face-velocity check
python3 /workspace/lagEunha/review_tests/voronoi_face_velocity_check.py   # needs numpy and scipy
```
