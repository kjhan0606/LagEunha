# Kepler integrator and original-reference audit, 2026-11-04

**Latest user criterion:** compare with GIZMO's error at the same time and
resolution; no need to suppress valid errors far below that level. Read
[the comparison criteria](/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs/KEPLER_GIZMO_COMPARISON_CRITERIA_2026-11-04.md).
The old1% original-radius cohort screen below is now an internal diagnostic,
not an independent accuracy rejection. Historical numbers/flags remain.

**Latest result:**449508 reached exact20pi and COMPLETED0:0 at16:36:32KST.
Root and Astra verified practical disk stability under the updated criterion.
Read [the final result](/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs/KEPLER_20PI_RESULT_2026-11-04.md)
and its completion audit. The notes below preserve earlier interim
observations; their running/pending status statements are historical.

Work is on LagEunha with shared
/gpfs access; Slurm commands go through syntax.kias.re.kr. No /scratch,
no new simulation starts/restarts in this audit, no production binary/source
changes, and no Git commit/push. This Docs directory is not a Git worktree.

## Live handles at 16:01:50 KST

- Native baseline 449508 RUNNING on syn104, one H200/one rank/eight CPUs.
  Last exact audit t=60.4304666507566, continuation-relative energy drift
  +6.44593e-9, invalid=clips=floors=0. Target remains exact 20pi.
- GIZMO reference 449138 RUNNING; do not cancel/restart it.
- Short 4A coupling pilot 449658 COMPLETED 0:0 at 15:33:05. No pending
  candidate job exists and no duplicate is needed.

Baseline package:
/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/native_u64_kdk_individual_20261104_v1/kdk_block2_wallkdk_v1

Main run: run_wall_lts_kdk_t4030_to20pi_v2. Latest saved checkpoint at this
poll is khout.219476.dat (t~60; footer not yet audited here).

## User's KDK/RK4 experience

Do not select KDK from bounded endpoint energy alone. The user observes
plain KDK losing circular orbits, RK4 preserving circular motion better,
and RK4 deteriorating after eccentricity develops. Existing read-only
orbit tests separate radius/phase/energy/angular momentum and pericentre
resolution. See analysis_tools/orbit_audit_20261104_v1 and
analysis_tools/gravity_compose4_pilot_20261104_v1; their source is bound to
the actual native gravity operator. They do not replace full hydro tests.

## New completed audit and lineage correction

Authoritative results and full explanation:
[full-reference README](/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/native_u64_kdk_individual_20261104_v1/kdk_block2_wallkdk_v1/analysis_tools/full_reference_audit_20261104_v1/README.md).

Current trajectory is RK4 t0→40.000091, global KDK→40.2, then block2 KDK.
The exact adiabatic-volume map starts at t=40.02. Five producer/consumer
checkpoint pairs are SHA-identical; hashes/contracts are preserved. It is
not a pure KDK run from t0. The EOS Gamma is float32(1.666667), or
1.6666669845581055. Entropy checks must use that configured exponent.

Use the already verified native t0 file:
/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/native_u64_w0_initial_20261101_v1/initial696.bin

The new read-only audit.py reuses the established ordered(), total_energy(),
weighted_mean(), BINS and Gamma. Identity checks give zero changes and the
prior RK4/circle endpoint's energy/density metrics reproduce within 1e-14.
No legacy reader or physical threshold was replaced. Wall and full-lineage
injection gates remain unqualified rather than fabricated from KDK logs.

At saved t=59.000378, total energy relative to original t0 is -2.44648e-6,
disk mass-weighted density change .95028%, entropy-log change 7.00464e-7,
and position error 1.03558e-4. The .5≤R0<.75 and 1≤R0≤2 density changes
exceed the original 1% cohort thresholds. All finite/mass/positive-state
checks pass. Do not call this full physical stability from small energy.

The old 17 large-amplitude cells inside the *restart-radius* bin1.2..2
all originate at t0 radii2.47756..2.62605: outer buffer, not original disk.
Stable IDs are joined explicitly. The t0-based orbit diagnostic has one
undefined epicycle/phase value at the exact r=0 site; this is not an invalid
simulation state. Details and return codes are in orbit_receipts.json.

## New geometric attribution, reusing the prior t45 method

Independent actual Qhull volumes calibrate against native t59 geometry;
selected disk cells are all inside the physical polygon. Actual max rho/rho0
is3.02299 with112 cells above1.5. Replacing only original outer positions
by exact circular reference positions, solely in a non-evolving geometric
diagnostic, eliminates all112 peaks and lowers density change from.95028%
to.52637%. Replacing inner positions alone gives.94765%; replacing both
backgrounds gives.52373%. All-exact geometry gives.06307%.

No simulation positions were steered. This extends the established t45
outer-support intrusion evidence to the current branch, without proving
the complete force history or a conservative cure. Residual error remains
with the actual disk positions. Metadata inherited by two extra control
outputs has an explicit separate erratum; numerical arrays/results stand.

## Short forward-4A coupling completed, not promoted

Candidate sibling package gravity_forward4_candidate_20261104_v1, job449658:
same GPU, same native t55 start, exact t55.0502026 endpoint for both cases.
Baseline47.15384s/151steps/451meshes; candidate48.44532s/152steps/454meshes.
Maximum measured relative energy excursions8.16370e-9 vs7.90226e-10.
One pair shows+2.74% runtime and10.33x smaller energy excursion, but both
gravity map AND cap/subdivision differ. No repeated timing or attribution
to 4A alone is established. Codec/mass/positivity/no clipping checks pass.

Original-disk endpoint velocity difference is only2.44e-8 max, position
1.68e-8, density relative L2 4.03e-8. Large all-cell velocity differences
are in the low-density core. Disk-density improvement is not demonstrated.
Two ordinary global-KDK fallback steps remain per case; hydro remains
second order. Pure-gravity 108-case convergence and eccentric-orbit results
are in the candidate README. No extra mesh per gravity microstep is added.

## Next actions

1. Preserve 449508 until its exact20pi endpoint. No new long candidate run
   is justified by the short energy result alone. Monitor sparsely.
2. At completion verify codec/footer/hash, compute original-t0 metrics
   with the new adapter, and render actual mass/volume density using the
   requested linear .5..2 GIZMO scale. Assess failed density cohorts and
   inherited/missing wall qualification explicitly.
3. Do not claim an integrator-only comparison with RK4/circle endpoint
   448969: it changes boundary type after t40. The paper sharp-annulus
   vacuum-boundary setup remains distinct and unresolved.
4. Keep no steering, Re=infinity, fixed mass/W=0 baseline scope and actual
   single-rank/block2/all-cell-geometry performance limitations visible.
