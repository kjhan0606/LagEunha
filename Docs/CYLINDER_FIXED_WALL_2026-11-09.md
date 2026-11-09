# Cylinder: fixed wall geometry, 2026-11-09

## Authorized scope

The user requested correction of the moving Cylinder boundary in the existing boundary operator and a new SU2 comparison, without adding empirical patches. The user subsequently explicitly approved representing a boundary cell by several convex integration pieces while retaining one material ID and one mass/momentum/energy state. This is an exception for boundary geometry, not permission for mass flux or a change to the interior hydrodynamic model.

## Established cause

The production radial-image power boundary rotates as its generating fluid site moves tangentially. At the saved t=3.50025 endpoint, the 252 actual wall faces have maximum normal speed 0.860357693826. The 122 self-image faces already move (maximum 0.146787434701), although their weight-rate contributions cancel exactly. Thus changing a weight limiter would not fix this boundary motion. The existing C/CUDA geometry functions reproduce recorded speeds within 3.33e-16. Their C/CUDA difference is at most 5.72e-13.

The existing stationary-wall flux gives exactly zero volume transfer and total-energy work. Setting this speed to zero while leaving the moving radial-image geometry would be inconsistent; the geometry itself must change.

## Current implementation

Private source and build: `Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_stationary_wall_20261109_v1/native`.

- The fluid integration region is the ordinary real power cell minus the fixed 128-facet SU2 obstacle, using its actual saved vertices. This finite polygon is not an exact circle.
- Exposed power-face segments and stationary solid facets are emitted. Internal interfaces of a convex-piece decomposition are not physical flux faces. Signed boundary integrals compute the aggregate cell volume and centroid without evolving multiple states.
- Existing one-sided stationary Riemann traces and C/CUDA pressure-force/pressure-work operators are reused. Stationary walls are excluded from the weight Jacobian; actual clipped fluid-face lengths enter it.
- Fixed-wall normal velocity, volume flux, and total-energy work are zero. The internal-energy contribution `-v_bulk dot F_wall` must remain: zero wall work does not mean zero change in a cell's bulk kinetic energy.
- Radial images are removed only in the new compile-time branch. Periodic aliases, inlet/outlet samples, fixed mass, SSPRK2, second-order reconstruction, and existing contact kinetic-energy accounting remain.
- Point contact uses the same fixed polygon, not the circumscribing circle. The existing removal of a remaining inward normal velocity is retained and reported, not hidden as conservative wall work.
- `GFSFIXS1` checkpoint manifests include a fingerprint of the polygon vertices; ordinary `GFSSTRM3` restarts are rejected. New uniform initial states are initialized in native C/CUDA from actual cut volumes at t=0. Evolved checkpoints are never reinitialized.

## Checks completed so far

Independent polygon-difference tests: 2,018 configurations, including touching, holes, fully covered and random cells; relative area error <=1.39e-16 and C/CUDA measure difference zero. This tests geometry, not fluid accuracy.

The first native initialization exposed a missing boundary in one nearly degenerate power cell. Inferring its support plane from a nearly zero-length rounded edge was wrong. The implementation now clips solid facets against the original power support equations. No geometric tolerance was relaxed. Failed runs `initial200` and `initial200_v2` are preserved; `initial200_v3` passed with 7,912 cells and uniform density/pressure on the actual cut volumes.

Native evolution `pilot200` (t=0 to .01) and `early200` (t=.01 to .1) completed. Retained cell masses are bitwise identical; U and V remain positive. Wall work is exactly zero at every reported stage. At t=.1, the physical energy change is -1.8249855e-5 relative to initial energy, from the pre-existing contact velocity cancellation. The energy accounting residual is <=7.11e-15; it must not be presented as physical energy conservation at that smaller level.

At that early checkpoint, the endpoint comparison was pending. The later coarse t=3.50025 comparison is recorded below; the high-resolution qualification is still running. This document does not declare wake blotches resolved or the method validated.

## Source integration and ongoing batch

The 14 boundary-related files were integrated into canonical `Exam/` with prior contents saved under `canonical_before/` and an exact delta in `canonical_fixed_wall.patch`. Existing canonical inverse policy and unrelated optimizations were preserved. The isolated canonical build passed C/CUDA compilation; a restart step from t=.25025 to .25225 produced all 13 stored state columns bitwise equal to the frozen production-based candidate.

At t=.25025, 166 actual fluid-wall segments cover the same 128 SU2 facets. Maximum normal displacement is 2.21e-16, maximum fractional facet gap 2.22e-16, and overlap zero. Reconstructed obstacle area is .7850827892386891 versus SU2 .785082789238688. The early density image is `comparison_01000_v2/density.png`; it is low-resolution and does not qualify the mature wake.

`matched200_t350025` is the coarse continuation to 3.50025. `matched800_t025025` advances the 126,696-cell run. The finite batch runner `continue_fixed_wall_batch.py` waits for both checked results, then continues the high-resolution checkpoint to 3.50025 with .25 checkpoint intervals and produces a same-time SU2 comparison. It stops on failed run/qualification gates. `batch_status.json` records actual progress. This is a live process, not an assertion that model reasoning continues unattended.

## Coarse endpoint at t=3.50025

The 7,912-cell run completed to the exact SU2 snapshot time 3.50025. Its final continuation took 1,674 accepted steps and 980.06 wall seconds while the fine job shared the GPU. All retained cell masses are bitwise unchanged; U and V stayed positive. The actual 160 fluid-wall segments cover all 128 fixed SU2 facets, with no detected overlap, maximum fractional gap 2.22e-16, and maximum normal position error 2.59e-16. The reconstructed obstacle area is .7850827892386874 versus SU2 .785082789238688.

Wall total-energy work is zero. Across all four continuation segments from t=0, the existing contact operation removes 0.0019203013409 energy, or 3.4630947e-05 relative to initial energy. The summed energy accounting residual is -3.79471e-14; it must not be substituted for the physical loss. See `coarse_t350025_cumulative_audit.json` for segment provenance.

The density comparison is in `Docs/figures/cylinder_fixed_wall_20261109/coarse_vs_su2_t350025.png` and `.pdf`. It uses unsmoothed native cell density and the same colour scale for both panels. The native coarse spacing is h=0.1, versus h=0.025 for the SU2 reference. The bow shock is present, but native wake patches are visibly cell-sized, less symmetric, and do not reproduce the smooth SU2 wake structure. Therefore the stationary-boundary defect is corrected, while the wake morphology is not yet qualified. This result does not identify those patches as physical vortices or attribute all of them to a single operator.

The display-grid density relative L1 difference is 8.6542% in the displayed window. This is not an exact cell-overlap error norm and cannot be compared directly to earlier fine-resolution overlap metrics. The 126,696-cell run is still executing; its same-time comparison is required before assessing the effect at comparable resolution.

The finite continuation runner had a self-copy error when invoked from its own archived path. The copy is now conditional on distinct paths. An actual archived-path restart from .25025 to .25225 completed with exit 0 and all run gates passed (`runner_selfpath_check`). This is a runner correction; the hydrodynamic executable and the active fine calculation are unchanged.

The idle batch supervisor was replaced with a version that detects an absent native process without a completion record and validates the target time, fixed-solid gap, and restart files before continuing. The native fine calculation was not signalled or restarted. Supervisor state and replacement provenance are in `batch_process.json`; completed status still means execution/geometric checks only, not scientific morphology acceptance.


## Fine endpoint completed and visually compared at t=3.50025

The126,696-cell continuation completed exit0; the real supervisor/native PIDs
4131872/4145694 are absent and the final completion marker/checkpoint agree on
3.50025. Runtime for this continuation was32,777.74s, not a matched GIZMO speed
measurement. No calculation was restarted to create the comparison.

`comparison_14000_fine/density.png` and PDF compare actual power-cell density
against SU2 flow_14000 at the same actual time3.50025, common scale0..6, no
smoothing. The128 fixed polygon facets and obstacle area agree; max gap1.249e-14,
overlap7.438e-15, and endpoint relative volume error6.871e-16. All14,296 RHS
wall-power records in the continuation are zero. Retained material masses are
exactly fixed, and U/V stay positive. Open-boundary total mass is not constant.

The physical relative energy drift is-7.126e-7, with contact kinetic loss
-3.95138e-5. The much smaller ledger residual-8.501e-14 includes that measured
loss and must not be described as exact physical energy conservation.

At nominal spacing.025 both methods use different discretizations: native
126,696 moving material cells versus SU2148,656 vertices/295,264 triangles.
Both have the same stationary slip polygon, periodic y and supersonic x
inflow/outflow; native material-column injection/exterior stream and SU2 fixed
mesh boundary implementations differ. Native uses its adaptive SSP update;
SU2 uses BDF2 dt.00025 and HLLC/MUSCL/Barth.

The shared display-grid density relative L1 is2.839%, compared with8.654% for
the coarse native result at the same time. This is a sampled diagnostic at.005
pixel spacing, not exact cell-overlap integration or a measured convergence
order. The bow shock and broad wake agree visually, but patchy downstream
material-cell densities remain visible while SU2 is smoother. Remaining patches
are not certified as physical vortices by this comparison. Wall kinematics is
resolved; downstream thermodynamic/shape accuracy remains a separate criterion.

See endpoint `fixed_wall_audit.json`, comparison `completed.json` and
`review.json`. The plot was inspected. No claim of <=2x GIZMO, complete2D
qualification or permission to bypass2D before3D follows.
