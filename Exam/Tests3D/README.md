# Exam/Tests3D: 3D test suite for LagEunha GFS

This directory has five 3D verification problems: IC generators, parameter
templates, Slurm scripts, analysis scripts against analytic or reference
solutions, and pass criteria.

| # | test | directory | reference |
|---|------|-----------|-----------|
| 1 | linear sound wave, convergence + Galilean invariance | `SoundWave3D/` | analytic |
| 2 | Sedov–Taylor blast, γ = 5/3 | `Sedov3D/` | exact similarity solution (ξ0 = 1.15167) |
| 3 | Evrard adiabatic collapse | `Evrard3D/` | 1D spherical reference computed here + HydroCode1D profile |
| 4 | spherical Noh | `Noh3D/` | exact (ρ_post = 64) |
| 5 | Kelvin–Helmholtz, McNally 2012 profile in 3D | `KH3D/` | compressible linear theory solved here |

> **Status: the repository has no 3D GFS hydro driver yet** (see "3D path
> audit" below, which covers only what is in the repo). Juhan Kim reports that
> an older 3D driver, **Exam3d**, exists on the cluster but has not been
> pushed. **Grok CLI: please add Exam3d to the repo and push it.** Grokbot will
> then audit it against the 2D GFS path and hook Tests3D up to it (preflight
> marker check and an IC/energy-log adapter). Until then every `run.slurm`
> stops at a preflight check (exit 3) unless the binary contains the marker
> string `LAGEUNHA_3D_GFS_V1`. Nothing here was run with LagEunha. The IC
> generators, analytic solvers, and analysis scripts were verified on
> synthetic data (`./selftest.sh`).

---------------------------------------------------------------------------

## 1. 3D path audit (Sep 2026, master after 187550c)

### 1.1 What exists in 3D

| piece | where | state |
|-------|-------|-------|
| 3D tree / neighbour helpers | `Exam/exam.c` `ex3d_findCentroid`, `ex3d_FindCellSize`, `ex3d_idivision`, `nearest3dOpen`, `ex3d_dist`, `det3d_dpq`, `det3d_dpqRK4` (~l.353–1007) | Only compute the search radius (`w2ceil`) from a 3D tree walk. **No driver calls them.** |
| 3D Voronoi cell primitives | `Voro/voro.c` `Voro3D_FindVC`, `Voro3D_FaceExtract`, `Voro3D_Volume_Polyhedron`, `findPolyhedronCentroid`, `Voro3D_norm_polygon`, `getAvgPressureOnSurface3D`; type `Voro3D_GasParticle` in `voro_eunha.h` | Cell, face, and volume geometry. In the repo, used only by the legacy `Exam/Sedov` code (presumably also by Exam3d). |
| 3D Laguerre library | `Voro/Laguerre/` (`constructLaguerreCell3D` CPU; CUDA `construct_cells_3d_kernel`) | CPU gives 2D+3D cells. The CUDA 3D kernel returns volume and a **vertex-average** centroid (not the true centroid) and **no face list**. **The main build does not use it:** `Voro/Makefile` builds only `voro.o` and `voro_eunha.o`. |
| legacy 3D Voronoi hydro | `Exam/Sedov/` (`vch_hydro.c` …) | Standalone, OpenMP, non-MPI. Monaghan-type artificial viscosity (`alphavis`/`betavis`). No HLLC, MUSCL, pair pressure, or RK4. Not in any Makefile, and **does not compile** (12 errors in our check). |
| particle types | `treevorork4particletype` has x,y,z / vx,vy,vz / ax,ay,az | Storage is 3D-capable. |
| 3D domain decomposition | `Cosmos/cosmology.c` `MakeDoDeInfo3D` (for `voroparticletype`/`treevoroparticletype`) | Cosmology only. |
| gravity with gas mass | `pmmain.c` `Get_TSC_Den` and the tree correction include VORO particle mass | Periodic cosmological TreePM (GOTPM) only. `RunCosmos` calls no GFS/Voronoi hydro. |

### 1.2 What the 2D GFS path has and 3D does not

All of these are 2D-only in `Exam/exam.c` / `Exam/exam_gpu.cu`:

1. **Force/flux routine**: `getAccVoro2DBlend_impl` (~l.4776–6093, ~1300
   lines), which contains:
   - HLLC face solver `hllc_face_2d` (~l.4440) + MUSCL reconstruction;
   - extreme faces `phase1_extreme` (P_max > 100 P_min, ~l.5045/5704);
   - pair pressure with the capped law and `gfs_pair_work_limit` (`gfs_pair.h`, ~l.5726);
   - face-centroid rotation correction `voronoi_face_rotation` (~l.4312; now also the `GFS_LAGUERRE_ROTATION` flag, 187550c).
2. **Density/volume/pressure/floor update**: `updateDenW2Pressure2DBlend` (~l.2776, ~800 lines).
3. **RK4 driver** `exam2d_vph_rk4_int_blend` (~l.9527, ~600 lines), with:
   - central forces `kepler_accel` (~l.4220) and uniform `GAS_ACC`;
   - `phase1_ie_stage` (SEDOV_PHASE1) and the `[RK4E]` energy ledger;
   - GFS_DUAL_ENERGY and GFS_ENTROPY_SWITCH.
4. 2D linked-list cell search, padding, `postStage`, periodic/reflective boundaries.
5. Domain decomposition for rk4 particles: `startRkSDD2D` / `MakeDoDeInfo2D`.
6. SIMMODEL dispatch in `eunha2.c` (KH, RT, RT_LF, MkGlass2D, Cylinder, Sedov2D only), params model, IC reader, output.
7. **GPU**: every kernel in `exam_gpu.cu` is 2D (`getAccVoro2DBlend_kernel`, `voronoi_tessellate_kernel`, …). There is no 3D GPU kernel, so a 3D port would be CPU-only until one is written. The 2D GPU path also returns `die`/`pie_in` as float, which limits a 10⁻⁶ sound wave in precision if a 3D GPU path copies that.
8. **Self-gravity coupled to non-cosmological hydro particles**: the 2D drivers have only the external `kepler_accel` and uniform `GAS_ACC`. The GOTPM gravity is periodic/cosmological and is not called from the GFS drivers. **Evrard is blocked twice: no 3D hydro and no isolated self-gravity.**

The paper draft (`main.tex`, section "HLLC + Cullen–Dehnen on 3D Laguerre
Cells") describes 3D equations that the code does not implement.

### 1.3 Decision

Nothing was ported. Exam3d (on the cluster, not yet in the repo; see the
status note at the top) is the starting point. Once it is pushed, Grokbot
will audit it against items 1–8 above and hook Tests3D up to it. Every 2D
code path is untouched (this commit adds files only under `Exam/Tests3D/`).

### 1.4 Port plan (for whoever writes the 3D driver)

- **Geometry**: take faces from `Voro3D_FindVC`/`Voro3D_FaceExtract` (CPU) or `constructLaguerreCell3D`, which needs a face list (area, normal, face centroid). Use the true polyhedron centroid (`findPolyhedronCentroid`), not the vertex average.
- **Fluxes**: `hllc_face_3d` with two tangential components. MUSCL gradients by least squares over the face neighbours (2D uses the face-weighted form; generalise it). Keep `phase1_extreme` and the P_max > 100 P_min trigger unchanged.
- **Pair pressure / work limit**: `gfs_pair_work_limit` assumes `nshare = 6` (hexagonal 2D cells). A 3D Voronoi cell has ≈15.5 faces on average. The value is an **open decision** (below); keep the capped law P ≤ 32 ρc².
- **Rotation correction**: the 2D scalar torque becomes a vector r_fc × F. The face centroid is in the face plane.
- **RK4 driver**: copy the 2D driver, generalise the state to (x,y,z,vx,vy,vz,ie), keep `phase1_ie_stage`, the floor, and the `[RK4E]` ledger. Add a hook for the gravitational acceleration and potential (Evrard) next to `kepler_accel`.
- **Gravity for Evrard**: isolated (non-periodic) Barnes–Hut on the hydro particles with Plummer softening ε, or direct summation for N ≲ 3×10⁴ (enough for a first test).
- **DD**: 3D rk4 domain decomposition, or 1 rank + OpenMP at first.
- **Interface the scripts here expect** (§2).

---------------------------------------------------------------------------

## 2. Interface between this suite and a 3D driver

- **Marker**: the binary embeds `LAGEUNHA_3D_GFS_V1`, e.g. `__attribute__((used)) static const char lag3d_marker[] = "LAGEUNHA_3D_GFS_V1";`.
- **Parameters**: `common/params3d.dat.template`. Keys marked `[2D key]` already exist in `Params/params.h`. The `[proposed 3D key]` keys are:
  - `Simulation Model = Hydro3D`;
  - `Hydro3D IC file`, `Hydro3D t_end`, `Hydro3D dump dt`;
  - `Hydro3D gravity` (0/1), `Hydro3D softening`.
  - Fixed: RK4 (`Hydro time-stepping way = 1`), `av_mode 5` (MUSCL + HLLC), `kappa 0` (Voronoi faces), `gpu_enabled 0`, centroid shift 0.
- **Environment**: `common/flags_base.env`, the established 2D names:
  - `SEDOV_PHASE1=1`, `GFS_FLOOR_LOG=1`;
  - `GFS_DUAL_ENERGY`, `GFS_ENTROPY_SWITCH`, `GFS_HALF_LIMIT`, `GFS_LAGUERRE_ROTATION` and `EUNHA_KEPLER_*` unset.
  - Per test, `config.env` overrides these; the run script also exports `HYDRO_TSTOP` and `EUNHA_DUMP_DT`.
- **Snapshots**: LAG3DV1 files `snap_%06d.l3d` every `dump dt` plus one at `t_end`. The analysis scripts accept any LAG3DV1 file names.
- **Energy log** (one line per step or per dump, to stdout):
  `[E3D] step=<n> t=<t> dt=<dt> Ekin=<> Eint=<> Epot=<> Etot=<>`
  (Epot = 0 without gravity). `common/lag3d_io.py:parse_energy_log` reads it.

### LAG3DV1 format (`common/lag3d_io.py`, `common/lag3d_ic.h`)

Little-endian. 128-byte header:
- `char magic[8]="LAG3DV1\0"`, `int32 version=1`, `int32 flags` (bit0: pot present);
- `int64 np`, `double time`, `double gamma`, `double box[6]` (xmin,xmax,ymin,ymax,zmin,zmax);
- `int32 periodic[3]`, pad, `double G`, `double softening`.

Then structure-of-arrays of `np` doubles:
- x y z vx vy vz mass u rho vol;
- optionally pot;
- then `int64 id[np]`.

`u` is the specific internal energy. In ICs, `rho`/`vol` hold the intended
values (the code recomputes them from the tessellation). Every IC also gets a
JSON sidecar (`<file>.json`) with the test parameters and the totals.
`common/l3d_check <file>` prints the totals from C, and `-w <out>` rewrites
the file through the C writer (byte-identical round trip verified).

---------------------------------------------------------------------------

## 3. Tests

Common to every test: γ = 5/3, RK4, MUSCL+HLLC, Voronoi (κ = 0), CPU-only,
4 MPI ranks × 16 OpenMP threads (placeholder; set the CPU partition in
`run.slurm`).

To launch (Grok CLI, **only once a 3D binary exists**):

```
cd $LAGEUNHA_SRC/Exam/Tests3D/<Test>
sbatch --export=ALL,N=<n>[,…] run.slurm
```

Each script:
1. creates `/gpfs/kjhan/LagForce/<Test>_gfs_…`;
2. copies the binary and records `COMMIT` and `ENV`;
3. generates the IC with `$LAG3D_PY` (numpy+scipy);
4. fills `params3d.dat`;
5. runs `mpirun`, writing `log` and `EXIT:<rc>`;
6. exits with mpirun's rc. Unlike the 2D scripts, a crash is therefore not reported as COMPLETED 0:0.

It ends by printing the analysis command.

Runtime estimates below are **estimates only**. They extrapolate 2D GPU
Gresho throughput (≈2.5×10⁵ particle-steps/s) to an assumed ≈10⁵
particle-steps/s per CPU rank for 3D cells (≈15 faces instead of 6), and
could be off by 3× either way.

### 3.1 SoundWave3D (test 1)

- **Setup**:
  - periodic unit cube, ρ0 = 1, c0 = 1, A = 10⁻⁶;
  - right-moving adiabatic wave along x (`--dir diag` for k ∥ (1,1,1));
  - one period, t = 1.
  - Cubic lattice with mass m_i = ρ(x_i) dx³, the same construction as the 2D KH/Gresho ICs (`--mode displace` gives equal masses). A glass is refused below A = 10⁻³ because its volume noise swamps the wave.
- **Boosted run**: `BOOST=1` adds v_bulk = 1 along x.
- **Runs**:
  - `N=32`, `N=64`, `N=128` (`-t 12:00:00`);
  - `N=64,BOOST=1`.
- **Analysis**: `analyze.py --snaps r32/<last> r64/<last> r128/<last> --boosted r64b/<last> -o sw3d`
  - L1(ρ) against the IC profile advected by (c0 k̂ + v_boost) t;
  - writes the convergence order, `sw3d_convergence.png`, `sw3d_profile.png`, `sw3d_summary.json`.
- **Pass**:
  - P1: order ≥ 1.7 between the two finest runs;
  - P2: L1/A ≤ 10⁻² at 64³;
  - P3: L1_boost/L1_rest ≤ 1.25.
- **Estimate**: 32³ a few minutes, 64³ ≈ 3 min, 128³ ≈ 40 min (≈110/215/430 steps).

### 3.2 Sedov3D (test 2)

- **Setup**:
  - periodic unit cube, ρ = 1, P_amb = 10⁻⁵;
  - E = 1 as a top hat on the 64 particles nearest the centre (as in the 2D Sedov IC; on the cubic lattice ties make it 88, which is reported and used);
  - t_end = 0.05, R_an = 1.15167 (E t²/ρ)^{1/5} = 0.3475.
- **Lattice**: cubic by default; `LATTICE=glass` for the glass (`common/lattice.py`). **Record which one.**
- **Runs**: `N=64`, then `N=128` (`-t 12:00:00`).
- **Analysis**: `analyze.py --snap r64/<last> --ic r64/ic.l3d [--snap2 r128/<last>] -o sedov3d`. It reports:
  - the swept-mass shock radius, the density-peak radius, ρ_peak, L1(ρ);
  - energy error, axis/diagonal anisotropy;
  - radial profiles vs the exact solution.
- **Pass**:
  - |R/R_an − 1| ≤ 0.03 (0.02 at 128³);
  - |ΔE/E0| ≤ 10⁻³;
  - ρ_peak ≥ 2.5 (3.0 at 128³);
  - anisotropy within 5 %;
  - L1(128³) < L1(64³).
- **Estimate**: 64³ ≈ 20 min (≈2000 steps), 128³ ≈ 6 h.

### 3.3 Evrard3D (test 3)

- **Setup** (Evrard 1988, MNRAS 235, 911; as used by Springel 2010, MNRAS 401, 791 and Hopkins 2015, MNRAS 450, 53):
  - γ = 5/3, G = M = R = 1, ρ = 1/(2πr) for r < 1, u = 0.05, v = 0;
  - E_pot(0) = −2/3, E_tot = −0.6167.
  - Particles: cubic lattice cut at r < 1 and stretched r → r^{3/2} (equal masses).
  - Plummer ε = 0.01.
  - Outer boundary: vacuum (default) or `BACKGROUND=uniform` low-density gas (open decision).
- **Needs** `Hydro3D gravity = 1`, i.e. isolated self-gravity coupled to the hydro particles, which does not exist yet (§1.2 item 8).
- **Run**: `N=64` (≈137k particles), t_end = 3.
- **Reference**, not digitised from any paper:
  - `ref/evrard1d_N2000_*.txt` come from `evrard1d_ref.c` (1D spherical Lagrangian, compatible energy, exact enclosed-mass gravity, |ΔE/E| = 1.5×10⁻⁵). Regenerate with `make_reference.sh`.
  - Its t = 0.8 profile matches the HydroCode1D reference (SWIFT's `evrardCollapse3D_exact.txt`, fetched by `getReference.sh`, not committed) to 0.3 % in log ρ.
  - Reference numbers: E_kin max 0.450 at t = 0.88, E_th max 1.758 at t = 1.07, E_pot min −2.543 at t = 1.05.
- **Analysis**: `analyze.py --log log --snap <snapshot at t≈0.8> -o evrard3d`. It reads the `[E3D]` lines and writes the energy-curve plot and a profile plot vs the 1D reference.
- **Pass**:
  - E1 (primary): max |ΔE_tot|/|E0| ≤ 1 %;
  - E2: t(E_kin max) within ±0.10 of 0.88 and t(E_th max) within ±0.15 of 1.07;
  - E3: mean |E_th − E_th,ref|/max E_th,ref ≤ 0.15;
  - E4: mean |log10 ρ/ρ_ref| ≤ 0.10 over 0.05 < r < 0.8 at t = 0.8.
  - E3 and E4 are thresholds we chose, not community standards.
- **Estimate**: ≈2 h hydro at n = 64 (≈2×10⁴ steps) plus the gravity cost. **Blocked**.

### 3.4 Noh3D (test 4)

- **Setup**:
  - periodic cube L = 6, ρ = 1, P = 10⁻⁶, v = −r̂, t_end = 2;
  - exact: R = t/3 = 2/3, ρ_post = 64, P_post = 64/3, pre-shock ρ = (1 + t/r)².
  - Only r < L/2 − t = 1 is uncontaminated by the periodic edge, so the analysis uses that sphere only.
- **Run**: `N=128` (R = 14 dx).
- **Analysis**: `analyze.py --snap <last> --ic ic.l3d -o noh3d`.
- **Pass**:
  - |R/R_an − 1| ≤ 0.05;
  - plateau (0.3–0.9 R) density within 20 % of 64;
  - pre-shock L1 ≤ 5 %;
  - |ΔE/E0| ≤ 10⁻³.
  - Wall heating at r < 0.3 R is reported but not graded.
- **Estimate**: 128³ ≈ 1.5 h (1–5 h), ≈1000 steps.

### 3.5 KH3D (test 5)

- **Setup** (McNally, Lyra & Passy 2012, ApJS 201, 18, extended uniformly in z):
  - box 1 × 1 × Lz (Lz = 0.25 default), P = 2.5;
  - ρ = 1/2, U = ±0.5, smoothing L = 0.025;
  - seed v_y = 0.01 sin(4πx); `--znoise` optional.
  - Lattice with mass = ρ dV.
- **Linear theory**: `common/kh_linear.py` (compressible eigen-solve of the exact profile) gives σ = 2.83 for the seeded k = 4π mode. The sharp incompressible value, 5.92, is **not** the right reference.
- **Runs**: `N=128` (128×128×32); optionally `N=256` for K2.
- **Analysis**: `analyze.py --snaps 'snap_*.l3d' --ic ic.l3d [--snaps2 …] -o kh3d`. It computes the McNally mode amplitude M(t) and fits ln M over 1.5 M0 ≤ M ≤ 0.06.
- **Pass**:
  - K1: |σ_fit/σ_lin − 1| ≤ 0.15 at n ≥ 128;
  - K2: the error does not grow with resolution.
- **Estimate**: 128×128×32 to t = 2 ≈ 50 min (≈2200 steps).

---------------------------------------------------------------------------

## 4. Verification done on the authoring machine

`PY=/workspace/venv/bin/python ./selftest.sh` checks:

- **C tools**: `common/l3d_check.c` and `Evrard3D/evrard1d_ref.c` build with `-Wall -Wextra`, no warnings.
- **Analytic solvers**:
  - Sedov ξ0 = 1.15167 (γ = 5/3), 1.0328 (γ = 1.4); cylindrical γ = 1.4 matches α = 0.984; energy integral 1.000001; post-shock ρ = 4.
  - Noh ρ_post = 64, P = 64/3, R(2) = 2/3, ρ(1, 2) = 9.
  - KH solver σL/U0 = 0.18968 vs Michalke (1964) 0.1897 at kL = 0.4446; McNally σ converges (2.8295 at n = 512, 2.8227 at n = 256).
- **IC generators** at small N (all five, including boost, diag, displace, glass, and background variants):
  - finite values, m > 0, u ≥ 0, inside the box, unique ids;
  - total mass, momentum, E_kin, E_int (and E_pot for Evrard) printed and cross-checked by the C reader;
  - C writer round trip byte-identical.
- **Analysis scripts**: each `--selftest` builds synthetic snapshots from the analytic solution (with controlled errors) and checks the script recovers them:
  - sound-wave order 2.00;
  - Sedov R ratio 0.996;
  - Noh plateau and radius;
  - Evrard 1D reference vs HydroCode1D, and energy drift 1.5×10⁻⁵;
  - KH fitted σ.
- **Run scripts**:
  - with the real-style binary (no marker) each script exits 3;
  - with a fake marked binary and a stub `mpirun`, each builds the run directory, IC, fully filled `params3d.dat` and `log` with `EXIT:0`.
  - Nothing was submitted.

**Not verified**: anything involving LagEunha itself. There is no 3D
driver, and this machine has no nvcc, MKL, or CAMB. The pass thresholds
are calibrated on analytic expectations and on the 2D behaviour of the code,
not on a 3D run.

## 5. Open decisions

1. **Evrard outer boundary**: vacuum (particle-code convention; needs an unbounded-cell treatment in the tessellation) or a uniform low-density background (`--background uniform`, which perturbs the energy budget slightly and is reported separately).
2. **Softening** for Evrard: ε = 0.01 (default), or scale it with the particle spacing.
3. **Work-limit `nshare`** in 3D: 6 (the 2D value) or the 3D mean face count ≈ 15.
4. **CPU partition, ranks, and threads** in `run.slurm` (placeholders).
5. **GPU precision**: if a 3D GPU path is written, `die`/`pie_in` must be double for the 10⁻⁶ sound wave.
6. **KH Lz**: 0.25 (default, cheap; the linear mode is z-independent) or 1 (full cube).
7. **Glass vs lattice**: lattice by default for every test (exact initial density). A glass is available for Sedov/Noh/Evrard to expose lattice imprint.
