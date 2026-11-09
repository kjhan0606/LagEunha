# Common face-tangent reconstruction and scalar QR evidence — 2026-11-09

## Purpose and preserved controls

Develop one reconstruction across dimensional reductions without new physical fields, shock-specific parameters, or a test-specific force. Existing Kepler, stationary-cylinder, RT and completed control data remain unchanged. These experiments are isolated native periodic2D paths; no3D evolution, MPI qualification or <=2xGIZMO speed certificate is claimed.

The previous actual-neighbor velocity extrema produced a finite force jump as an edge vanished. The log-scalar parent fixed positive scalar representation but inherited that velocity defect. The new face-tangent map uses only the two owner values and their two tangentially translated values for its common-shift interval. Raw jump is clamped between zero and the owner jump, then common shift is bounded. The same signed-coordinate helper handles log(rho), log(P), and normal velocity. It is a continuous sequential map, not a joint Euclidean projection.

Tangential velocities use the same face-measure LS inverse already used for velocity reconstruction, evaluated at r_t=r-(r dot n)n. One geometric coefficient per CSR entry is added; no physical state is added. Continuity depends on nonzero particle separation and a bounded fit inverse. This does not certify colliding sites or rank-deficient geometry. A linear field is preserved locally when the face separates the centroids and the fit is nonsingular.

Scalar LS is performed in log coordinates. The pressure-offset reconstruction uses rho/P as its coefficient, following grad(P)/rho=grad(log(P))/(rho/P). This remains a reconstruction hypothesis; it is not an entropy inequality, exact balance across a shock, or a new constitutive law. Log values are exponentiated and nonpositive/nonfinite numerical outputs are rejected; no physical pressure floor is added.

## Face-tangent parent gates and first native runs

Artifacts: `/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1/Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/common_face_tangent_pair_20261109_v1`.

All9 1D endpoints complete with fixed cell masses; Shu–Osher800, other physical cells200. The18 initial/final states extruded into exact weighted strips agree in force, volume and work to5.20e-14 relative. This is a free-RHS reduction excluding the two boundary cells on either end, not full2D motion. Local affine/boost/scale and4096 signed-coordinate symmetry checks pass. Native C/CUDA RHS disagreement is1.61e-15; cached pressure-response discrepancy6.59e-14. Rebuilding the native executable required exam_gpu_extract.o, material_cuda_rhs.o and material_cuda_project.o; the scalar object must not be omitted.

The previously finite velocity-topology jump (~8.1%) now decreases with geometry perturbation: maximum relative force difference1.515e-4,1.515e-6,1.515e-8 for perturbations1e-4,1e-6,1e-8. This conditional continuity test is not a complete positivity or time-evolution theorem.

Both32-square one-period native waves completed. Density error divided by input wave amplitude is0.0198795(axis) and0.0207928(diagonal). Compared with the bounded-velocity parent0.0197526/0.0200255, these are slightly larger; no across-the-board accuracy improvement is claimed. Fixed masses, endpoint geometry and energy checks pass. Shared-machine elapsed times are not performance evidence.

The same candidate failed glass Sedov at t=0 before any accepted step: CUDA_PROJECT_FAIL cell913, phase gradient_project. The declared regular case was not launched. Halving dt repeats the same initial-state failure. Parent failure and all diagnostics are retained. The overall parent method is not adopted.

## Identified numerical failure and same-objective correction

`sedov/glass/initial_matrix_audit` reconstructs the actual initial power cells independently; area error4.42e-16 relative. The CPU replay also fails, excluding a CUDA-only explanation. Cells913 and2589 have valid geometric fit condition numbers1.11 and1.08, but pressure-offset design matrices have condition numbers1.3428e7 and2.8064e7. Squaring the design into normal equations gives an ill-conditioned Gram system and its rank check rejects these cells.

The QR candidate keeps every original row, weight and right-hand side. The LS rows are sqrt(w)*x_effective and RHS sqrt(w)*Delta(log q). Streaming Givens rotations solve that objective directly. The small shared kernel supports dimensions1..3; only1D/2D hydrodynamics is exercised. No row normalization, new diffusion, pressure floor or physical retuning is applied.

The numerical rank criterion DOES change:128epsilon is now relative to QR factors rather than the squared Gram matrix. Retaining the constant does not mean retaining the old acceptance threshold. Truly rank-deficient synthetic rows are rejected. The two previously rejected actual systems agree with70-digit QR to2.96e-15 and zero reported relative difference, respectively. The full saved initial CPU RHS succeeds with zero summed energy-rate residual. The older saved Sedov failure also succeeds, with summed energy residual1.60e-14.

QR artifacts: `/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1/Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/common_face_tangent_qr_20261109_v1`.

## QR validation and remaining limits

Same QR kernel used for scalar reconstruction in both1D and2D. All9 1D runs were recomputed from the original byte-identical ICs, with the same exact Riemann solver, SSPRK2, CFL and boundaries. Density L1 values remain equal to the parent at numerical precision; they are not newly improved by changing the LS solver.

| Case | Time | Density relative L1 | Maximum density | Energy balance error |
|---|---:|---:|---:|---:|
|sod|0.2|0.0037420156|1|8.07e-16|
|blast|0.038|0.034394656|8.1493168|1.03e-15|
|noh|0.6|0.0044067245|4.9730888|-9.78e-16|
|shu|1.8|0.0031983262|4.7338356|2.66e-15|
|lax|0.13|0.0031122348|1.3537716|-3.64e-15|
|doublefan|0.15|0.023003504|1|-2.32e-15|
|collision|0.035|0.0036417891|33.312427|2.25e-15|
|contact_open|0.25|2.5895335e-13|1|0|
|moving_contact_open|0.25|1.0991208e-13|1|4.18e-16|

Blast and Noh density overshoots remain:8.1493 versus reference6.4417, and4.9731 versus3.99999. Endpoint completion and low integral error do not certify peak accuracy.

The18 strip comparisons pass to5.20e-14; local affine and edge-disappearance tests are retained. Nineteen standalone C/CUDA fixtures include the newly identified initial Sedov failure. Exact objects linked into the native executable also pass19-fixture CPU/CUDA scalar+RHS and cached-response tests. The scalar QR change requires rebuilding only exam_gpu_extract.o and material_cuda_project.o relative to the face-tangent parent; material_cuda_rhs.o is byte-identical.

The native QR binary SHA is `0be5ac55098a166887ce86e4e35bcb28ba85476ad2fdc16da50969292c0746c1`.

Two native batches were frozen and launched:32-square axis/diagonal waves for a complete period, and64-square glass then regular finite-disk Sedov to0.1. Same archived ICs and physics settings as the controls; only CPU allocation, failure-dump paths and diagnostic CPU/CUDA slope auditing differ. All plans, effective params, commands, native PIDs and sources are recorded. At launch the total checkpoint bound was16 including archived inputs,28.51MB, with149TB free. Waves use CPU12/13 and Sedov4/5; existing RT stays on14/15. No live process was stopped or restarted.

At the original launch checkpoint, QR evolution assessment was pending; the completed results and decision appear below. The batch stopped on a failed gate without automatic tuning. Active constraint behavior, strong-shock morphology, low-pressure shear, KH/RT accuracy, stationary-wall integration, full MPI, Kepler compatibility and matched-accuracy GIZMO runtime remain unproven. Any future boundary/3D port must preserve the same mathematical operator;3D evolution waits for2D completion.

## Completed native comparison and adoption decision

The QR waves, shear64 through t=1, and Gresho64 through t=1.2 completed with fixed masses, positive states and invariant gates. No live run was interrupted. Glass Sedov first reached t=.1 but was correctly marked failed because one reported start-volume-normalized defect was1.015340978e-11, just above the declared1e-11 gate. The old record was not relabelled as a pass.

The inverse had divided its residual by target volume; the final audit divided by initial volume and used a different subtraction ordering. At the same early step the inverse reported9.867523526e-12. This was a convergence-contract mismatch, not proof of physical blow-up. `common_face_tangent_qr_volume_20261109_v1` changes ONLY the inverse's convergence metric to the exact audit arithmetic and start-volume denominator, retaining1e-11. The Newton RHS, geometry, pressure force/work, timestep and all scalar/velocity reconstructions are unchanged. A10000-case contract audit includes an expanding-cell counterexample; native rebuild changes only exam_gpu_extract.o.

With that contract aligned, both original glass and regular64-square Sedov ICs reach t=.1 and pass the unchanged invariant thresholds. This is a subsequent binary; the Q flow results and V Sedov results must not be represented as a complete suite on one executable. Images use actual power polygons, equal0..5 colors, periodic clipping and no smoothing. Visual inspection finds stronger local shock-edge structure in glass and markedly stronger axis-aligned lobes/low-density pockets on the regular grid.

| Sedov IC | Density L1 control → candidate | Pressure L1 control → candidate | rho=2 radius RMS control → candidate |
|---|---:|---:|---:|
|glass|7.2095% → 10.1722%|4.4504% → 6.7517%|0.7635% → 1.1051%|
|regular|10.9362% → 10.7284%|9.5776% → 9.8774%|0.9040% → 2.1223%|

These norms still use the ideal point-source solution for a finite-disk IC; angular contours are morphology diagnostics, not exact shock-position errors. The glass maximum density rises4.328→5.531, with four cells above the plotted upper color5. Regular density L1 changes little despite visibly worse morphology, demonstrating why integrated norms alone cannot decide adoption.

At shear t=1 the pressure relative L1 is0.000484235 (0.0484%); cumulative acoustic heat agrees with the unchanged common control within6.1e-8 relative. Thus the new map does not remove the remaining shear heating. At Gresho t=1.2 acoustic heat is18.69% smaller, cell momentum L1 changes4.0822%→4.1694%, pressure L1 changes0.10135%→0.08696%, and the minimum material entropy proxy change becomes more negative(-0.000712→-0.003074). Reduced total heat is not proof of more accurate dynamics or a cellwise entropy property.

Decision: do not adopt this reconstruction as the common production method. Preserve the reusable QR and inverse-metric corrections and all positive/negative results. Before another physical modification, separate the scalar-coordinate and normal-velocity reconstruction effects on the saved shock states. No extra viscosity sensor or shock-specific parameter is justified by this batch. The complete comparison and hashes are in `common_face_tangent_qr_volume_20261109_v1/cross_test_review.json`.

The matched-accuracy GIZMO<=2x requirement remains open. The older shared/different-hardware ratio1.61291 is not a certificate for these binaries. Reuse the existing H200 ABBA harness only after a common candidate has matched accuracy; include initialization, geometry, weight inversion, rejected steps and I/O.


## Final image review and interpretation limits

The authoritative comparison images are `common_face_tangent_qr_volume_20261109_v1/figures_v2/glass_density.png` and `regular_density.png` (PDFs alongside). Both were visually inspected after reducing each panel title to one line, exposing the regular candidate's upper axial protrusions. Image hashes and review findings are recorded in `figures_v2/manifest.json`; earlier images remain preserved.

The inverse metric correction aligns both the denominator and subtraction order with the existing final audit. It is not uniformly stricter for every cell. The Gresho 18.69% reduction is specifically the acoustic-dissipation ledger; the material entropy proxy is not a cell-integrated thermodynamic entropy certificate. Its initially generator-sampled velocity also limits the exact continuum comparison. Neither the scalar log coordinate nor the normal-velocity map is separately established as the cause of the Sedov regression by this combined comparison.

Current decision remains not to adopt the combined reconstruction. The completed glass-layout control in `COMMON_SEDOV_GLASS_2026-11-09.md` instead changes only the IC using the frozen control operator. Its 128-square comparison must not be attributed to this newer candidate. Saved-state ablations and the separately frozen RT halo diagnostic are next; no production replacement, 3D evolution or speed certification follows from these results.


## Subsequent mechanism separation

See `COMMON_RECONSTRUCTION_ABLATION_2026-11-09.md`. Six frozen operator combinations isolate the current log pressure-offset fit as a large instantaneous force change while the velocity maps preserve star pressure and change common boundary speed. The actual initial hot-side fit becomes nearly constant before limiting. This strengthens the non-adoption decision without declaring the prior primitive interval admissible; it still generates vacuum faces.
