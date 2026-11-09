# Face-local common scalar reconstruction: periodic gates,2026-11-09

This is the continuing candidate from `COMMON_RECONSTRUCTION_CONTINUITY_2026-11-09.md`, not an adopted default. Original successful Kepler, current cylinder and RT simulations are untouched. No3D evolution has been started. Full-method accuracy, common boundary reconstruction and ≤2xGIZMO runtime remain unproven.

## Full-period native2D wave results

All8predeclared trajectories completed: axis/diagonal16²,32²,64², plus diagonal32² half-amplitude and half-timestep controls. Each uses the existing native C/CUDA execution path and immutable candidate SHA256 `aae8f1f6c426ddd025fafb7c6a575c7cb1476ae94eed33e4a5da6020bddb881a`. Archived IC checkpoints and original endpoint/CFL/dt/volume-inverse settings were retained; runtime CPU/CUDA raw-gradient comparisons passed. Fixed masses, positive states, volume closure and total energy checks passed.

|Direction,N|Control densityL1/epsilon|Candidate densityL1/epsilon|
|---|---:|---:|
|Axis16|.09728173|.04599282|
|Axis32|.04003256|.01951584|
|Axis64|.01162539|.00613015|
|Diagonal16|.15619249|.04847881|
|Diagonal32|.05702602|.01991671|
|Diagonal64|.02049637|.00616720|

These errors are normalized by initial wave amplitude epsilon=1e-5, not by background density. Actual-cell integration and sine/cosine fitting are saved in `common_face_local_native_20261109_v1/wave_interpretation/report.json`; `convergence.png/pdf` were visually inspected. Axis and diagonal at equal N do not have equal cells per wavelength.

Density convergence orders for16→32,32→64 are1.237,1.671(axis) and1.283,1.691(diagonal). The candidate reduces error but **second-order evolution has not been demonstrated**. At32² diagonal, halving dt reduces density/momentum errors by2.674%/2.759%; halving amplitude changes them by−0.001261%/−0.001364%. Thus neither finite-amplitude contamination nor the tested time increment explains most of the remaining error. This does not isolate a particular remaining limiter or geometric term. At64² amplitude ratios are.998919(axis) and.998856(diagonal), phase errors−.00385678/−.00394658rad. Both damping and phase/shape error remain.

## Actual1D-to2D operator reduction

`face_local_candidate_v3/dimensional_reduction` embeds the actual initial and final states of all9candidate1Dproblems into exact weighted rectangular strips. For an interval of width h_i with generator at its midpoint, W_i=h_i²/4 places each shared Laguerre face at the prescribed interval boundary. Four identical periodic y rows form the2D state. The unchanged candidate1D driver RHS and candidate2D free-RHS are evaluated, and force, volume rate and total energy rate are compared after dividing by strip height.

All18states pass: maximum scaled discrepancy below4.6e-14, transverse force exactly zero. Both means and gradients evolve substantially in the shock snapshots, so this is stronger than constant-state algebra. The outer two x cells on each side are excluded because the1D driver zeros outer slopes while the2D fixture wraps periodically. It verifies the planar **semi-discrete interior** reduction, not the full moving mesh/weight inverse, time integration, boundaries or3D adapter. No other hydro modifications were made for this audit.

## Ongoing periodic shear and vortex gate

`common_face_local_native_20261109_v1/flows` runs sequentially in the existing third Ada slot, CPU12/13:

1. Current frozen common-method64² low-pressure shear control through t=.1,.4,.7,1.
2. Same IC and settings with the face-local candidate through the same times.
3. Candidate64² Gresho through t=.002,.3,.6,.9,1.2, compared with the already completed current common Gresho control.

The new shear control is necessary because archived shear runs used older volume coupling; they are not silently relabelled as the current control. Its IC is the existing periodic rho1,P1e-6,vx=.1sin(2pi*y),vy0 checkpoint, verified field-by-field. Initial Gresho13fields match the prior common control exactly. CFL.4,dtcap.00025,endpoint-volume tolerance1e-11 and original physics are unchanged. Each segment restarts only from its own accepted preceding endpoint; no running calculation is restarted for diagnostics. Runtime scalar parity checking is disabled after the wave gate, identically for control/candidate. Maximum39checkpoints, about106MiB, free space recorded in plan. Native and supervisor handles/status/results are saved. Hard invariant failure stops the batch without automatic retuning.

Current control shear completed t=1 with relative thermal increase/volume-mean pressure error4.84235e-4 and exact fixed masses; candidate and vortex results must be compared when available. Conservation and positivity alone will not decide adoption. Wall/exterior reconstruction still uses the prior route and is not ported in this experimental binary; no RT/cylinder test may use it until that work is done.

## Other tests and speed qualification

Existing high-resolution cylinder was observed live near t=3.0515 with stationary-wall power zero; target3.50025 is not yet reached. Its current frozen model is separate from the face-local candidate.

Syntax Slurm was rechecked:456829 pendingQOSGrpGRES,456838 and456949 pendingDependency,454560 pendingPriority, each H200:1. No jobs were cancelled or resubmitted. The≤2xGIZMO requirement is still unproven; an isolated CPU free-RHS speed reduction is not a complete-solver benchmark.


## Completed shear comparison and partial vortex evidence

Both matched64² shear trajectories completed t=1. All masses are bitwise fixed, states remain positive, and energy residuals are below4.6e-15. Relative thermal increases are0.00048423507112(control) and0.00048423504601(candidate), effectively unchanged. Accepted-step ledger integration attributes the increase to acoustic heat:7.26352641e-10 versus7.26352597e-10. Pressure work is only6.47e-18 versus8.71e-18; the SSPRK2 compatible correction is about−1.30e-27. Thermal-budget residuals are below6.5e-18. This is attribution in the discrete ledger, not proof of the initial cause of the residual velocity error. Scalar reconstruction continuity alone does not remove this smooth-shear heating.

At Gresho t=.3, the candidate accumulates9.996% less acoustic heat and9.027% less net thermal increase than the frozen control, but its generator velocity error is9.237% larger. Independent exact-field integration on the actual saved convex Laguerre polygons confirms deterioration: momentum L1 is0.017880454 versus0.016703370; pressure L1 is0.00045145185 versus0.00038345407; density L1 is0.00029650928 versus0.00027176570. Quadrature orders16/32 agree far more closely than these differences, and the recovered volumes agree with the saved ones within3.2e-14 for the candidate. Thus the result is not explained solely by generator sampling or lower dissipation. It is a partial-time result, not the final t=1.2 gate.

Evidence: `common_face_local_native_20261109_v1/flows/budget_comparison_v2/report.json`, `gresho_candidate/segment_01/cell_assessment`, and `interpretation_t03`. The v1 budget-report velocity ratio used an incorrect fallback for a nested control field; it is explicitly superseded. Its energy/thermal calculations were unaffected. The v2 code requires the correct field and contains no missing-field substitute.

No candidate is adopted. Continue the already-running vortex trajectory to its predeclared endpoints and compare shape and cell-integrated error, rather than stopping at lower heat or adding a Gresho-specific correction. Wall/exterior consistency, full second-order convergence,3D and matched ≤2xGIZMO runtime remain unproven. Previous successful Kepler and current cylinder/RT runs are preserved.


## Reproduced affine velocity-trace discontinuity (no production change)

`common_face_local_native_20261109_v1/velocity_affine_audit_v1` compiles the frozen native velocity reconstruction functions unchanged. An8² periodic jittered Voronoi fixture supplies actual polygon volume centroids. Select an interior face by its tangential offset, then impose a local affine divergence-free shear whose normal means at the two owners coincide. Both complete reconstruction stencils stay away from the periodic discontinuity. The physically exact common face-normal velocity is−0.023448735562234754.

Raw LS reconstruction yields−0.023448735562234778 and−0.023448735562234705 (maximum error4.86e-17). Their tiny difference is nonzero in floating point, while the two owner means are exactly equal. `gfs_material_ordered_jump_alpha` therefore sets alpha=0; the current mapping scales BOTH complete increments, returning zero on both sides and discarding the finite common face motion. Perturbing one owner by±1e-14 switches alpha to1 and restores the correct face velocity to within9.4e-15. This is an actual frozen-function reproduction on valid irregular cell geometry, not merely an invented pair of trace numbers.

The faulty condition is suppression of the common reconstruction component when only the jump needs limiting. A tolerance around zero would hide this example rather than give a continuous map. A possible algebraic route is to preserve the raw common component while clipping the reconstructed jump to the interval between0 and the owner jump. HOWEVER, native volume constraints, force/work maps and their frozen linear actions must all use the same resulting coupled map. Changing only a Riemann input or only the alpha formula would not satisfy that contract. Such an operator has not been implemented or qualified; shock bounds and1D reduction would also need to be retained. No new physical parameter or experiment-specific sensor is justified by this result.

This defect does not yet prove the initiating cause of Gresho error, residual acoustic heating, wave convergence order, or earlier Kepler failures. Current scalar-candidate trajectories continue unchanged to their declared endpoints. Keep this read-only fixture as a reproduction test for the next common velocity-map investigation.


## Completed periodic flow batch: candidate not adopted

All13 planned segments completed, with both shear cases reaching1 and Gresho reaching1.2. Supervisor637235 and final child662451 are absent. No trajectory was interrupted or restarted from its IC. The final budget report is `flows/budget_comparison_final/report.json`; final actual-cell shape/error comparison is `flows/interpretation_t12` and its plots.

|Method|Integrated acoustic heat|Kinetic loss fraction|Cell momentum L1|Cell pressure L1|Density L1|
|---|---:|---:|---:|---:|---:|
|Control|0.003813917241|0.04333629011|0.0408223137|0.00101346514|0.001034765641|
|Candidate|0.0038060648|0.04323593941|0.04315202791|0.001039033937|0.001050320556|

At t=1.2 the candidate's acoustic heat reduction has shrunk to0.206%, while the cell-integrated momentum error is larger. The early t=.3 heat reduction did not persist as a material overall gain. The scalar candidate repairs a proven vanishing-face discontinuity and improves acoustic-wave errors, but does not qualify as a generally better common solver. Retain it as a research branch; do not adopt it for RT, cylinder, Kepler or3D. The original baseline remains intact.

**Historical correction:** The unbounded common-speed route described below was already rejected in `D/continuous_jump` and repeated in `S/continuous_normal_v2`. It is not a new or approved next method; see `COMMON_BOUNDED_VELOCITY_2026-11-09.md`. The following paragraph is retained as the superseded proposal.

The next core implementation target is the independently reproduced common-mode loss in the existing velocity trace. For raw trace rows T_L,T_R and owner selectors E_L,E_R, investigate preserving C=(T_L+T_R)/2 while limiting the jump J=(E_R-E_L)+alpha[(T_R-T_L)-(E_R-E_L)]. The paired traces are C−J/2,C+J/2. This is a candidate algebraic definition, not a verified solver. Row sums, exact affine reproduction, face reversal, boost invariance, force/work adjoint identity, nonnegative dissipation where claimed,1D reduction and shock behavior must be checked with the same map. Do not add a tolerance-based switch or problem-specific sensor to suppress the reproduced failure. CPU/CUDA and geometric constraint applications must stay consistent.

No complete-solver ≤2xGIZMO performance result has been obtained.3D remains deferred until2D qualification. Existing cylinder and RT/GIZMO comparison processes are independent controls, not runs of this candidate. No manuscript or remote repository was changed in this stage.
