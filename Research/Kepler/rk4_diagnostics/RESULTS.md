# Native RK4 stage admissibility audit

Current next implementation is fixed in `ASTRA_CURVED_BRANCH_PLAN.md`:
common smooth branch representatives and separate root/representation,
event and target-specific gluing budgets, first on triangle225/137 and the
same tiny-cell fixture. Recovered complete source-complement artifacts in
`curved_complement_v2_first_completed_recovery` give344/384 quadrature
estimates passing and40 failing, worst ratio42.0019. Positive packet mass
roundtrip is9.322e-20; this is NOT full cold-remap acceptance. First native
exit status/JSON were lost by a driver Bash-runner editing error that caused
a duplicate; the duplicate was terminated, the first complete artifacts
independently recovered, and immutable exec runner/env snapshots now prevent
that failure. No full rerun was launched to replace the JSON. No new Kepler
time evolution, weight Newton or compression is promoted.

Latest event isolation: `CURVED_EVENTS_RESULTS.md` and
`ASTRA_CURVED_EVENT_PLAN.md`. Explicit source-edge crossings, power vertices,
radial tangencies and chart crossings recover the formerly missed closed
tiny cell: area error~2e-21 at unchanged1e-11 tolerance, event mass/K bounds
pass. Corrected exact-vertex domain predicates pass all384 actual triangles
for W0 and the preserved nonzero weight pattern. These are event certificates,
not full cold-transfer or evolution acceptance. First positive outgoing plus
source-complement run failed on a zero quadrature weight: an exactly zero
reference vertex had been expanded by interval arithmetic into an artificial
radius~1e-2475. Outputs remain frozen; exact zero arithmetic is now preserved,
with a bounded origin regression preceding any replacement full run.

Latest actual transfer work: `CURVED_TRANSFER_RESULTS.md` and
`CURVED_TRANSFER_BUNDLE_20261021.json`. Entire shared-edge functional gaps
are bounded, and native positive curved m,p,L,K,U,V quadrature now runs on
all384 triangles with MPI bitmatch and a preserved926330-point packet.
Independent source integrals agree, but target radial convergence FAILS:
12x10->24x16 target K changes0.06425 U. No weight solve/compression/evolution
is promoted. The old bounded interval-radius guard was terminated after its
small-cell fixture failed; its data are preserved and it will not be resumed.
The replacement event/complement work is described above. CPU binding was proved incorrect and
fixed: the same4-rank geometry calculation is4.53x faster with bitwise equal
results. The original historical/longest-paper Kepler target is unchanged.

Latest actual curved geometry: `CURVED_GEOMETRY_RESULTS.md` under
`ASTRA_CURVED_REMAP_PLAN.md`. Native interval bounds cover all384 triangles,
with J lower bound0.891667 and separated exterior components. Global
injectivity remains conditional on exact shared-trace gluing. Interval
degree6 angular root/partition code runs OpenMP2/MPI1/4, including actual
finite intercell transfer and small-lens/endpoint/tangency fixtures; CSVs
bitmatch. Radial topology, full overlap moments, cold error and finite weight
solve remain incomplete; no new Kepler evolution follows these controls.

Latest actual finite-state bundle: `MATERIAL_FINITE_RESULTS.md` and
`MATERIAL_FINITE_BUNDLE_20261021.json`. The12-entry material map, physical
stress/work, Piola traction and actual centered spin/rotation moments now
pass CPU/OpenMP/MPI1/4 and real CUDA at two quadratures. Same-backend MPI
states/forces bitmatch; rigid rotation generates negligible physical heat.
Actual deformed x/v are retained. The native curved packet's degree6
power-plane ray algebra is independently verified against the map. These
are finite RHS/geometry controls, not a time step or long Kepler result.
Next: actual curved overlap/weight remap and admissible cold re-expression.

Latest spatial accuracy gate: `MATERIAL_P3_RESULTS.md`. Fixed complete P3 plus
one internal h/2 refinement preserves the24 physical cells and continuous IC.
Generic physical-trace hybrid solves pass native MPI1/4 and the predeclared
curvature gate: max chi6.2912e-4, max cell mean3.9603e-5; curvature maximum/L1
improve8.50x/7.77x. Stable initial-residual arithmetic is independently checked
at80 digits and on actual CUDA queries. This is an initial spatial gate only.
The12-entry finite material map is now tested in the bundle above; curved
Kepler remap and cold re-expression remain. No time evolution or long
stability result follows this gate.

Latest independent hybrid bundle: `MATERIAL_HYBRID_RESULTS.md`. Actual physical
trace constraints give85 conforming scalar DOFs instead of the old83. The
independent multiplier solve and direct dual-face traction close numerical
cell momentum/work/angular ledgers in native MPI1/4. Physical pressure/source
approximation errors remain. A new initial material-pressure curvature check
finds1216.89/time^2 maximum relative error despite initially tiny Pdot. This
requires bundled evaluation before promotion; no new time evolution occurred.

Latest shared-field gate: `MATERIAL_SHARED_RESULTS.md`. Native globally
continuous kinematics/mass assembly now run MPI1/4 on the24 actual cells.
Global energy/angular ledgers and shared traces hold, but the original
physical per-cell Euler work/torque ledger differs by up to15.918 U/time.
Gravity-off/linear-source controls and the exact M_i*a-f_i factorization
identify consistent-mass projection/intercell transfer. This is not Pdot.
The local transfer structure needs bundled evaluation before runtime use.

Latest structural work: `MATERIAL_ALE_RESULTS.md`, machine summary
`MATERIAL_ALE_BUNDLE_20261021.json`. Compatible material-map local RHS preserves
initial Kepler thermal/torque identities and agrees with prior native IC
integrals. Finite affine common-refinement/weight remap controls and the CUDA
core are checked. This is a partial initial bundle: global shared kinematics,
actual curved Kepler remap/compression, composed density-wave/sliding response,
shock stress and actual native runtime remain unimplemented/unverified. No new
Kepler evolution or long stability certificate follows these controls.

Additional kinematic evidence: `AFFINE_FACE_KINEMATICS_RESULTS.md`. Native
affine controls verify that scalar Wdot differences shift a power face's
normal speed but cannot remove its tangential strain slope. Zero net mass
does not eliminate position/momentum moment transport. A material-cell
redesign must handle that actual ALE transport. No new evolution was run.

Latest follow-on: `VOLUME_UPWIND_A_RESULTS.md`. The single reconstructed-HLL
volume candidate fails the first angular gate and produces large artificial
thermal rates on the same24 actual cells. Pair Euler conservation alone is
not enough. No P4, follow-on Jacobian, nodal runtime or evolution was started
from that candidate. Failed native sources/packets are frozen and retained.

Latest spatial bundle: `NODAL_SBP_RESULTS.md`. Native positive-diagonal nodal
P2 SBP on24 actual failure/neighbor cells satisfies the thermal chain identity,
but has large cold-flow divergence error. An advected density-wave Jacobian
has a verified growing mode1.810393/time, despite the stationary translating
subcase passing. This candidate is excluded from runtime/long evolution.
Details and all source/matrix records remain for the next bundled review.

Latest bundle: `ENRICHED_COUPLED_RHS_RESULTS.md`. Native initial moving-cell
Euler RHS and screened weight-rate solve now run on71080 cells, with canonical
mass/angular ledgers and MPI consistency. A velocity-perturbation thermal
tangent gate fails: independent projection reproduces the large spurious
pressure derivative, including volume L1 errors. No enriched time evolution
is started. The operator representation requires a revised bundled review.
Historical failures and the actual t6-10/longest-GIZMO goal remain intact.

Latest implementation/evidence: `ENRICHED_EULER_RESULTS.md`. Enriched Euler
initialization runs in native CUDA:71080 cells, max relative cold-P error
3.403e-10, no nonpositive samples, MPI1/4 all exported fields bitmatch.
Conditioning is measured/corrected and fixed-cell Euler perturbation helpers
are checked independently. Full moving-cell DG transport, weight/mass
coupling and runtime state remain absent. The fresh continuous IC integrals
also differ from legacy lumped cell momentum/E; the explicit initial ledger
must govern future initialization. No enriched timestep has occurred.

Previous spatial evidence: `DG_PROJECTION_RESULTS.md`. Native P2/P4/P6 Euler
projection is implemented; same-IC cold-pressure reconstruction fails in1552,
4 and4 cells respectively. Original IC contains an independently verified
central density feature on scale1e-6, distinct from gravity softening. No IC
was changed. Next selected plan: `ASTRA_ENRICHED_EULER_PLAN.md`, with common
Euler reference-field enrichment and freely evolving coefficients. No DG time
evolution exists yet. Prior affine initialization/MPI/restart results remain
in `MOMENT_INITIAL_RESULTS.md`; its thermal closure and independent
collisionless kinetic-tensor evolution are rejected. t6-10 remains unverified.

2026-10-20. Physics remains native C/CUDA. Experimental checks default OFF.

## Confirmed

- `GFS_EXACT_STOP=1`: caps only the last selected RK4 step at the requested endpoint. `GFS_DIAG_DT_MAX` supplies an optional additional diagnostic ceiling. The native physical CFL remains an upper bound.
- `GFS_FAIL_NONPOSITIVE_U=1`: checks every stage input/RHS and final combination before a negative U can be concealed by the EOS pressure floor. Defaults add no MPI collectives.
- `GFS_RK4_AUDIT_STEP_MIN/MAX` and optional `GFS_RK4_AUDIT_STAGE` select read-only native stage rows.
- `GFS_RK4_FAILURE_SNAPSHOT=1`, one rank only, preserves the normal step-start and invalid stage using the native full-stride format. Only `.start.dat` is a restart state. No rolling checkpoint or old run is overwritten.

430381 and the v3 repeated one-step preflight verified that default/strict+audit produced identical recorded fields and metadata to the frozen prior native binary. Exact stop reached t=dt=1e-6. Initial incomplete build430372 omitted rebuilding the KH library; it was never used for physics. Build430375 explicitly rebuilt both libraries. v3 SHA and logs are retained.

## Hidden intermediate Sedov failure

The old-limiter raw-geometry Sedov run's final ledger first reported negative U around step458. The strict replay instead caught **step256, k4_input**, at start time.024580813969805604:

- ID1699: U=-4.240326688886853e-6.
- ID1829: U=-9.677064554125536e-6.
- Normal start of that step is retained in `validation/sedov_first_failure/khdiag.step000256.start.dat`.

430394 replays only this captured interval, h=5.629704901366495e-5, endpoint.02463711101881927:

| Replay | Result |
|---|---|
| original full step, strict | Same k4 failure; captured x/y/v/U/V/m match original failed state exactly |
| original full step, observe | Final U positive, but k4_input and k4_rhs contain2 negative-U cells |
| two capped half steps, strict | First half passes; second k4_input has1 negative-U cell, minU=-1.246182108616267e-7 |

Thus final-state positivity misses an inadmissible RK4 intermediate state. Halving once is insufficient. This evidence does **not** prove that an arbitrary CFL reduction repairs the spatial work/flux operator.

## Variable contract

In `SEDOV_PHASE1`, `die` stores a total hydrodynamic energy rate, but `phase1_ie_stage` evolves the internal-energy rate `die-m*v·a`. The RK4 stages do not directly update conservative total E and then subtract kinetic energy. The logged conservative-total-energy Forward-Euler positive-step root is explicitly a **counterfactual diagnostic**, not the exact stage update. The actual RK4 Butcher combinations also differ from FE; a bound based only on each stage's instantaneous U cannot by itself guarantee all stage inputs are positive.

All original failure outputs remain. No new positivity method or integrator has been adopted from this audit.

## Fine cold interval and signed work ledger

430482 replayed the normal fine-old step-start at t=.10196264987779348. All four original failure IDs67935,68001,68526,68592 recur. The full replay has bit-identical U/V/m and only roundoff differences in x/v. One full step fails at k2_input; two half steps fail at the first k4_input. Four quarter steps reach the same endpoint .10258818651626195 without a floor. This is timestep sensitivity, not a physics pass: ID68526 ends with U/U_adiabatic_from_IC=.04059 (start .04462).

The read-only ledger ran on user-authorized lageunha (1 CPU thread), after checking all four native source hashes against frozenv5. Pending Slurm analysis430502 was owner/workdir checked and cancelled because QOSGrpGRES prevented it starting. Source git metadata resides in a syntax-local /tmp worktree, unavailable remotely; exact frozen source hashes established the shared source identity instead. All evolution was the native CUDA solver; Python here only analyzes preserved outputs.

`fine_cold_interval/ledger/results.json` and `selected_faces.tsv` reconstruct the first RHS. Maximum force difference1.06e-22, Udot/Edot difference7.01e-15 (stored rate precision), surface closure4.44e-17. Define s_geom=(v_geom-v_i)·S and s_used=(v_work-v_i)·S. The signed identity is

R_i=Udot_i+P_i Vdot_i_geom
   =sum_f(P_i-Pstar_f)s_geom_f
    -sum_f Pstar_f(W_used_f-W_geom_f)
    +P_i v_i·sum_f S_f.

For ID68526 at this state:

| Rate | Value |
|---|---:|
| Native Udot | -2.3690151e-8 |
| Reversible -Pi Vdot_geom | -2.9313445e-11 |
| Signed pressure contrast residual | -2.3660843e-8 |
| Work velocity residual | 3.32e-21 |
| Pair pressure contribution | 0 |

The other three failed cells agree by rotation to the precision of the actual states. Thus this interval's excessive cooling is assigned to pressure-contrast work; the current pair protection and geometric/work velocity mismatch are not its source. This does not prove that every Riemann local thermal exchange must equal its adiabatic reference, but a fixed-mass inviscid closure still needs a compatible entropy/work contract.

The actual fixed-weight Voronoi Jacobian reference F_i=sum_j P_j dV_j/dx_i, Udot_i=-Pi Vdot_i has total pressure power residual -1.84e-20 against sum|Pi Vdot_i|2.89e-6. Independent Qhull tessellations at x±delta v for certified interior cells verify Vdot: halving delta1e-5→5e-6→2.5e-6 gives L2 relative error3.94e-7→9.83e-8→2.45e-8, as expected for centered differences. At1.25e-6 it is8.73e-9. No evolution using this reference has been performed.

The quarter stages repeatedly change Udot sign (~-2.3e-8 to+4.7e-9). Current CUDA source has a cell-pressure ratio100 switch that changes both pressure solver and work speed. It is a specific suspect requiring direct branch evidence, not yet a proven cause of the entire cooling history. Do not silently alter that threshold.

The older user handover `/home/kjhan/BACKUP/Eunha.A1/Claude/LagForce/Docs/HANDOVER_2026-10-09.md` separately reports a conservative sliding-cell instability with growth ~cs/h in different TCF/ring tests. Increasing resolution reduces initial curved-flow reconstruction errors but can increase the rate of that unresolved mode at fixedcs. Therefore a consistent volume/force identity or a short positivity pass cannot alone certify long Kepler stability. Old weight-DAE geometric feasibility failures and staggered-remap limitations remain relevant; no old flag has been reactivated.

## SSPRK3 implementation and fresh test

See `SSPRK3_RESULTS.md` for the subsequent opt-in native SSPRK3 implementation, verified MPI1/4 full rollback, successful stage-positive captured intervals, and the fresh fine run's thermal/practical failure. Default RK4 matches the frozen reference bitwise. This candidate also refreshes current-stage geometry/EOS before halo gradients and therefore does not isolate a time-method change with identical RHS cache inputs.

The fresh 71,080-cell path stopped at t=.13509103681003889 on its3072-RHS budget,515 accepted steps,584 rollbacks,873.38seconds. Every committed state stayed positive without heat top-up, but the four focus U/Uad ratios fell to .01947-.02108. The preserved normal restart state is `ssprk3_fresh_v2/fine_fresh/khdiag.step000516.ssp_rejected_start.dat`. It did not reach t=.2 and was not extended to t6-10. Independent current-volume error is at float precision; individual masses and W0 remain exact. Final analysis: `ssprk3_final_analysis_v2.json`. SSP safeguards alone have not repaired local thermodynamics.
