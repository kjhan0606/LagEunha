# Hydro operator continuity audit, 2026-11-06

Scope: existing native 2D C material RHS/PCG with CUDA Laguerre geometry, fixed cell masses, no steering. 2D must finish before3D. Goal remains common1D/2D/3D physics and runtime<=2xGIZMO; neither completion nor final method selection is established.

## Demonstrated defects and corrections

1. Velocity-limited reconstruction used in reversible pressure-force transpose: at fixed initial Sod geometry, velocity noise1e-18 changes acceleration~89.2. Energy still conserves. Separating geometric reversible C/C^T from limited dissipative K removes this finite jump (response1.47e-16).
2. Equal angular LS weights for every retained face: position perturbation1e-12 introduces sliver faces and acceleration change~53. Face-measure LS reduces it to3.70e-8 and scales with perturbation. Match limiter gradient and trace construction metric.
3. Reusing pressure limiter in reversible C instead (experimental split mode2) is rejected: relative pressure noise1e-12 on finite cold shear gives Vdot change1.55e-4; geometric reversible control gives1.13e-16. Moving a homogeneous nonlinear limiter to the conjugate field is not a cure.
4. Dissipative velocity limiter on nonzero cold sine shear: transverse noise test did NOT show a finite jump; force differences tend tozero. Do not generalize defect1 beyond demonstrated scope.

Evidence directories under Research/MaterialVolume_20261105: material_sod_map_sensitivity_20261106_v1, material_geometry_continuity_20261106_v1, material_pressure_reversible_20261106_v1, material_transverse_noise_20261106_v1. All fixtures are distinguished from actual evolution.

## Time-evolution qualification still fails

Periodic planar Sod64x16, gamma float1.4, t=.05: max transverse velocities original limited angular.09799, face measure only.10531, split/angular.03498, split/face measure.05230. Exact fixed masses and energy residual~1e-15 do not establish accuracy. Central point-reference plots compare cell averages at generators; they are not cell-integrated exact-error estimates. Early symmetry loss precedes weight constraints.

Next simplification test: keep density/P limiting, leave velocity reconstruction unlimited (GFS_MATERIAL_THERMO_LIMIT_ONLY), retain acoustic damping/heat, face-measure geometry, existing shared operator. This omits both velocity limiter computation and a distinct reversible trace policy. It is a comparison, not an accepted solution; shock stability/positivity and shear accuracy must be demonstrated.

## Implementation and speed evidence

First CFL retry uses computed stage2 bound instead of automatically halving; all stage checks and full rollback retained. Short8192-cell restart1->1.005:55->28steps,33.6->17.5seconds under concurrent load. Forced rollback matches reference within roundoff; constrained continuation to1.28 matches old state exactly when timestep unchanged. Not a GIZMO benchmark.

Seven C/CUDA assembly fixtures for split reversible operator match to4.20e-17; this does not qualify full GPU-resident hydro/constraint evolution or multi-rankMPI. Running8192-cell KH continuation uses immutable earlier binary, not later experimental changes. Query live process/log for its time; do not restart from beginning.

Kepler baseline retained separately: W0 fixed mass, variational pressure/adiabatic geometric volume evolution, mixed integrators through20pi. It motivates geometric consistency but does not qualify new shock closure. RT gravity/boundary integration, broader2D suite,3D qualification and runtime target remain open. No new empirical damping/threshold/symmetry filter has been accepted to hide the planar failure.

Update: thermodynamic-only limit comparison completed Sod t=.05 with max|vy|1.14e-10 and column rho/P spreads3.18e-11/4.79e-11. This is a concrete planar-symmetry improvement by removing velocity limiting, not all-shock accuracy qualification. Pressure undershoot and resolution errors remain; see material_thermo_only_20261106_v1/metrics.json.

Resolution update: identical thermodynamic-only binary at128x32 completes Sod t=.05. Point-reference density/pressure errors1.638%/1.969%, versus2.983%/3.603% at64x16; transverse velocity5.10e-10. Pressure undershoot remains. Controlled8192-cell KH restart1->1.005 with the same face-measure binary: enabling only density/P maps changes heat by+0.0953%, mass identical and total energy equal to printed precision. This is short restart evidence, not full KH/RT qualification. See MATERIAL_RT_PORT_AUDIT_2026-11-06.md for verified RT gaps.

Longer controlled KH continuation launched: material_thermo_only_kh_longcompare_20261106_v1 resumes each branch's own verified t=1.02 checkpoint to1.1, immutable same thermodynamic-only candidate binary, sequential default/limited comparison on syntax GPU1. provenance.json stores checkpoint hashes and controller handle. At launch follow-up, controller2690286 and child2690625 live; first branch reached1.0369346. Original production2348638 independently remains live near1.2558 toward1.3. These are distinct binaries/experiments; do not combine as one trajectory. Final comparison pending.

KH controlled continuation t1.02->1.1 completed for both immutable thermo-only branches:516accepted/508retries each, cell masses unchanged, energy relativechange-1.13e-16, minimumP9.75469/9.75494. Limited/default interval heat ratio1.0007551, pressure relativeL1 difference2.4146e-6. This is still a late-time restart comparison, not a full growth convergence result. Evidence material_thermo_only_kh_longcompare_20261106_v1.

Speed candidate: optional accepted-CFL midpoint-slope extrapolation reduces t1.02->1.03 KH retries57->12 with58accepted steps. All original CFL/geometry/positive-state/rollback checks retained; mass exact, energy<=2.27e-16relative, pressure difference5.99e-11. No hydro/weight equation change. Disabled for external kicks; exact time-continuity guard on history. Same-binary timing rerun pending, not a GIZMO<=2x claim. Evidence material_cfl_extrapolate_20261106_v1.

Same-binary native timing follow-up completed: KH8192,t1.02->1.03,2CPU threads/GPU2 sequential off/on,37.169s->22.976s (1.618x locally,38.18% time reduction),58accepted both,57->12retries. Off endpoint bitwise matches pre-extrapolation off reference. Original old KH stayed onGPU0; host variability remains, no GIZMO benchmark claim. Prediction remains opt-in pending difficult constraints and shock qualification. Evidence material_cfl_extrapolate_benchmark_20261106_v1.

CFL shock regression: native2D Sod64x16,.005->.05,231accepted each, retries48->12, fixed masses/energy maintained, pressureL1 difference4.773e-8. Constrained KH32x64,1.271825->1.28,active<=8:82steps,no retries, on/off bitwise equal, but dtcap1e-4 dominates, so predictor is not materially exercised. A variable-CFL cap1e-3 follow-up was launched rather than treating this as full active-case qualification. Evidence material_cfl_shock_constraint_20261106_v1 and material_cfl_active_variable_20261106_v1. Existing physical closure concerns remain.

Active variable-CFL follow-up completed: same constrained32x64 restart withcap1e-3,52accepted steps both,active<=8,retries51->24. Pressure relativeL1 difference7.20e-10, heat ratio1.0000001145, exact masses, energy<=2.27e-16relative. Unlike fixedcap fixture, timestep proposal genuinely changes intervals here. Supports opt-in use in later controlled2D runs, not physical closure/default/MPI/3D/GIZMO qualification. Evidence material_cfl_active_variable_20261106_v1.

Latest Gresho refinement launched with saved identical old ICs, face-measure/thermo-only maps/constraint projection and CFL predictor, n32/n64 target.3; preserves old results. n32 completed: radialRMS.0062742 vs old.00419094 (worse), tangentialRMS.0150290 vs.0157279 (better), heat.000507126 vs.000572048 (lower), pressureL1.00423747 vs.00383515 (worse). New dt differs, so a same-cap1e-4 n32 temporal control was launched in material_thermo_gresho_time_20261106_v1. n64 remains running; no unmatched-time convergence or overall accuracy improvement claim. Evidence material_thermo_gresho_20261106_v1.

Latest Gresho64 completed t=.3: radialRMS.00247942,tangentialRMS.00701428,densityL1.000156902,pressureL1.00150149,heat.000107337. Compared with latest32, radial/tangential errors fall by2.53x/2.14x and heat by4.72x; dt also changes so no pure spatial-order claim. Compared with old64, radial/pressure errors are larger (old.00159003/.00105197), heat nearly unchanged. Temporal-control32 has now completed; its quantitative comparison follows. Do not claim uniformly better physics or adopt on conservation alone.

Gresho32 temporal control cap1e-4 completed t=.3. RadialRMS0.00627426003 versus adaptive0.00627420166; relative difference-9.30363e-06. Thus the larger radial error versus old method is not removed by matching its timestep cap. No correction fitted to Gresho was added. Evidence material_thermo_gresho_time_20261106_v1.

Same-state Gresho isolation: final-state force change due geometry metric RMS.4678 versus thermo limiter.0003566; initial regular-grid metrics coincide. Reverting metric does not uniformly improve frozen total force, and pressure/dissipative errors partially cancel. A dual-face-volume metric A*d gave mixed force/shear results, not adopted; removed its temporary production macro, retained reproducible diagnostic snapshot. All unrestricted metrics fail the frozen P1e-6 high-frequency cold shear mode4 (face RHS status2 vacuum). This is not a Gresho trajectory failure but a remaining robustness limit. No floor/extra sensor added. Evidence material_gresho_operator_isolation_20261106_v1.
