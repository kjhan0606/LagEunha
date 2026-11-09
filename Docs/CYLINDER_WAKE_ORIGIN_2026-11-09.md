# Cylinder wake density texture: boundary and operator audit, 2026-11-09 KST

## Finding and scope

The resolved wake contains shear in both native and SU2. Its native density
blotches cannot all be identified as physical vortices: selected deep density
troughs lack a rotation signal in the measured boxes, and several have
entropy-proxy differences with much smaller pressure than density differences
against SU2. This is not proof of local neighbour pressure equilibrium. Numerical contributions
are plausible, but their shares relative to physical shock entropy and advected
inhomogeneity have NOT been isolated. The current discrete stationary-wall
representation has a measured physical defect and must remain in the diagnosis.

No physical source, limiter, viscosity, floor, entropy repair, default method,
or production trajectory was changed. This request used saved-state diagnostics
and three short controls with the existing native binary. Original checkpoints
and both restart companions are unchanged. The bounded controls all completed;
no new production continuation beyond t=3.50025 was launched by this audit.

## Matched data and boundary conditions

The native endpoint and both SU2 fields have t=3.50025. Native has 126,696
owned cells, initial spacing .025; SU2 has nominal h=.025 and .05. Native local
cell size changes substantially. These are not uniformly equal-resolution runs.
The actual native solver is fixed-mass material direct-FV, reconstructed exact
Riemann response and SSPRK2, with endpoint volume-coordinate inverse active
only since t=2.5415. Old comments mentioning KDK/HLLC do not describe this path.
SU2 uses Euler HLLC, MUSCL, Barth–Jespersen and BDF2, dt=.00025.

Both have domain 20 by 4, circle centre (2.5,2), radius .5, upstream rho=1,
u=(1,0), P=.07999999821186066, gamma=1.399999976158142 and Mach about2.988.
Both are inviscid: native Re=0 means viscosity disabled in this implementation,
not physical Reynolds number zero. There is no matched finite physical Reynolds
number. SU2 smoothness cannot be attributed to a lower prescribed Reynolds number.
Numerical dissipation and perturbation seeds are not matched.

| Boundary | Native | SU2 | Assessment |
|---|---|---|---|
| Upstream | Constant supplied upstream lattice | Supersonic inlet | Same intended state; distinct discrete treatment |
| Downstream | Outgoing samples retained ballistically | Supersonic outlet | Same outward-flow intent; distinct algorithms |
| Top/bottom | Periodic owner aliases/wrap | Periodic markers, shift (0,4) | Same intended periodic condition |
| Cylinder | Radial-image power-cell faces and original circle-contact reflection | Fixed polygonal Euler slip wall | Material geometric difference |

No explicit perturbation was added in either of these production lineages.
Cartesian native IC and triangular SU2 mesh can seed different perturbations.
Separate earlier seeded runs are not ancestors of this result. Matching periodic
conditions does not prove that periodic images have no effect in either solution.

At t=3.50025 native has convex cells, zero circle intrusion and no boundary gap,
but the excluded obstacle area is .7903427952 versus pi/4=.7853981634
(excess .6296%). Maximum boundary radius is .5246762 instead of .5.
The maximum geometric normal wall speed is .86036 and instantaneous wall
pressure power is -1.73125e-4. SU2 has a fixed polygon; its fine obstacle area
is .785082789. Thus no-penetration and convexity do not establish a stationary
circular wall. The native wall-work plus contact-kinetic change over the last
production interval is -8.64721e-5. This known defect has not been causally
assigned a fraction of the wake texture. Simply zeroing wall work would break
the existing geometry/work relationship and is not an approved remedy.

## Density versus rotation

The diagnostic integrates the ORIGINAL piecewise-constant velocity exactly along
square boundaries of side L=.1,.2,.4, splitting at actual native faces and SU2
median-dual boundaries. It does not differentiate a smoothed density image.
For CCW boundaries, n ds=(dy,-dx), G=integral(u tensor n ds)/L^2 and
omega=G_yx-G_xy=circulation/L^2. The plotted rotation diagnostic is
lambda_ci=.5 sqrt(max(0,-[(G_xx-G_yy)^2+4 G_xy G_yx])).
This is an invariant of the finite-area MEAN gradient, not the average of a
pointwise invariant. Affine translation, shear, rotation and expansion checks
pass. Unresolved roundoff-length intersections are recorded; their errors have
not been fully propagated into lambda_ci. Positive-loop counts are not vortex
counts. Nonzero vorticity also occurs in shear, and nonzero lambda_ci can be
numerical rotation. Zero on these scales does not exclude all smaller vortices.

On 341 common locations, rotation-diagnostic RMS is:

| L | Native | SU2 h=.025 | SU2 h=.05 |
|---|---:|---:|---:|
| .1 | .342316 | .149919 | .101734 |
| .2 | .212517 | .120218 | .090426 |
| .4 | .064437 | .048208 | .048975 |

Native excess is strongest at small scales. SU2 itself changes with resolution,
so this does not certify either native small vortices or SU2 as the exact flow.
Upper and lower shear layers have the same vorticity signs in all three runs.
Increasing rotation measures from t=3.00025 to3.50025 in a fixed region are not
an advected-vortex history or a measured KH growth rate.

The following descriptive density troughs were selected from the central region
x=3.4..5.6, |y-2|<=.3; no accuracy threshold was tuned to them. Differences are
against fine SU2 integrated over the SAME actual native cells. K=P/rho^gamma
is a cell-average entropy proxy, not a measurement of mass-averaged entropy.
Activation means t=2.5415, when endpoint volume inversion began.

| ID | Final (x,y) | rho difference | P difference | K_activation/K_IC | K_final/K_activation | lambda_ci at L=.1,.2,.4 |
|---|---|---:|---:|---:|---:|---|
| 54480 | (4.487, 1.842) | -61.30% | -1.68% | 10.180 | 1.048 | 0.000, 0.000, 0.000 |
| 48093 | (5.475, 1.837) | -64.08% | +0.37% | 7.705 | 1.008 | 0.728, 0.000, 0.000 |
| 72880 | (4.521, 2.151) | -58.22% | +3.18% | 9.801 | 1.039 | 0.000, 0.000, 0.000 |
| 64079 | (3.440, 2.142) | -27.21% | +280.28% | 2.842 | 3.960 | 0.000, 0.000, 0.000 |

At the centres of density holes54480 and72880, lambda_ci is zero in all
three measured boxes. At48093 it is nonzero only for L=.1. These averages
do not prove the absence of a rotating core or that all texture is numerical;
pressure agreement with SU2 is not a local pressure-equilibrium test. ID64079 explicitly has a large pressure discrepancy and ongoing
growth, unlike the mostly inherited examples.

All7,062 selected material IDs were traced through20 actual-ancestor snapshots
from IC to final state, with masses bit-exact and the final extracted state
identical to the geometry export. The old unswept wake probes are a DIFFERENT
trajectory and were not substituted. Some representative trough cells already
had K/K_IC around8--10 before t=2.5415; their subsequent increases are only
about1--5%. K rises during expansion as well as recompression. The coarse
snapshot cadence cannot locate each shock/contact or prove that this rise is
unphysical. Typical selected cells reached sqrt(V)=.125--.166 during expansion,
versus initial spacing .025, making local spatial resolution a serious issue.
sqrt(V) is an area-equivalent scale, not a directional resolution certificate.

On790 complete neighbour stencils in the descriptive central region, RMS of
log(rho_i) minus the face-length-weighted neighbour mean is .22991 native versus
.083879 fine SU2. For log(P), it is .066438 versus .132863; for log(K), .317511
versus .047304. The same native-cell averages and neighbour graph are used.
These measure texture, not L1 errors. log(K)=log(P)-gamma log(rho), so the third
metric is not independent causal evidence. Other wake cells have large density
AND pressure differences; small pressure differences relative to SU2 apply only to the selected examples.

## Existing pressure/work operator and short control

For each directed material face, the existing operator uses
F=-A p_star n, Vdot=A u_star, Edot=-A p_star u_star, with opposite exchanges
on its neighbour, no mass transfer, and Udot=Edot-v_i dot F. The weight-rate
solve matches instantaneous geometric volume changes. The endpoint inverse
changes weights to match the trapezoidal step volume; it does not repair U,
mass, position or momentum. In endpoint mode the predictor volume is not
inverted; the existing stages mode also inverts that predictor.

The logger quantity H=Udot+P Vdot satisfies d(log K)/dt=H/U for fixed mass
in the semidiscrete equations. It includes signed subcell pressure/velocity
work, not only irreversible heating. With closed faces, uniform face pressure
equal to cell pressure cancels this contribution. Nonuniform reconstructed
pressure, face response and cell-average velocity need to be examined together;
total energy conservation alone does not certify accurate thermodynamics.

Three controls restarted the EXACT same t=3.00025 checkpoint and companions,
using the same native binary and params, ending at3.00525. Environment differences
were only timestep cap and endpoint/stages inverse choice. No restart IC or
new physics was used. All retained masses remain exact; all U,V remain positive.
The half-step run has20 steps plus one3.1e-15 endpoint roundoff remainder.

| Case | Accepted steps | K increase, ID64079 | delta(log K) minus baseline |
|---|---:|---:|---:|
| endpoint_coarse | 10 | 1.15198548% | 0 |
| endpoint_halfdt | 21 | 1.15199503% | 9.44e-08 |
| stages_coarse | 10 | 1.15199064% | 5.1e-08 |

Caps are .0005/.00025/.0005. First RHS values are identical in all three:
Udot=1.74542549e-4, Vdot=.04571512, P Vdot=7.84871149e-5,
H=2.53029664e-4, H/U=2.14258839. This cell is expanding yet its U increases.
Its kinetic-energy rate is -4.20698383e-4 and total face-energy rate is
-2.46155834e-4; energy accounting closes. Direct circle-contact displacement
of this selected cell is zero in these controls. That does not remove indirect
boundary influences carried by its state or neighbours.

Halving dt changes this selected delta(log K) by only8.24e-6 RELATIVE to its
increase; stage-volume inversion changes it by4.45e-6. Thus this local increase
is already present in the spatial RHS and is not explained mainly by the new
time step or predictor-volume mismatch. This does NOT exclude historical time
errors, other intervals, or physical shock entropy. No formal time order is
claimed from these two caps and this short common inherited state.

Current finite-step production volume closure is about9.6e-12, and short-control
RHS volume closure is about2.2e-12. These do not establish spatial thermodynamic
accuracy or a correct stationary wall. The earlier inverse-off5.14% volume
defect and limiter-jump records belong to older trajectories; their causal
role in THIS saved state has not been re-established.

## Decision before any patch

Do not add viscosity, density smoothing, entropy floors or a new limiter to
hide the texture. The evidence supports a shared physical shear layer plus
additional cell-scale density/entropy-proxy structure, with some small-scale
rotation. It does not support calling all density blotches resolved KH vortices.

The next causal comparison should separate (1) the existing nonstationary
discrete wall pressure work from (2) reconstructed intercell pressure/velocity
work during resolved expansion/recompression, using the same saved material
IDs and face-level budgets. A candidate stationary-wall treatment must satisfy
the geometry/volume/work relationship as well as convexity and no penetration.
Earlier failed clipping, ghost fitting and wall-work suppression must not be
repeated without addressing their recorded failure. Restarting a late checkpoint retains previously formed entropy-proxy
contrasts; that restart alone cannot isolate how they originated.
Matched local spatial resolution and both SU2 meshes remain required before
classifying the residual small-scale structures as physical or numerical.

## Evidence and reproducibility

- [Density / pressure / circulation / rotation figure](/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1/Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1/contact_failure_t2838_v1/positive_gap_v1/entropy_trace_v1/prevention_t265_v1/continuation_t350025_v1/wake_origin_v1/fields_t350025/density_pressure_rotation.png)
- [Final field diagnostics and source hashes](/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1/Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1/contact_failure_t2838_v1/positive_gap_v1/entropy_trace_v1/prevention_t265_v1/continuation_t350025_v1/wake_origin_v1/fields_t350025/comparison.json)
- [Earlier matched field diagnostics](/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1/Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1/contact_failure_t2838_v1/positive_gap_v1/entropy_trace_v1/prevention_t265_v1/continuation_t350025_v1/wake_origin_v1/fields_t300025/comparison.json)
- [Material history and neighbour texture](/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1/Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1/contact_failure_t2838_v1/positive_gap_v1/entropy_trace_v1/prevention_t265_v1/continuation_t350025_v1/wake_origin_v1/material_history_v2/summary.json)
- [Selected density-trough diagnostics](/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1/Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1/contact_failure_t2838_v1/positive_gap_v1/entropy_trace_v1/prevention_t265_v1/continuation_t350025_v1/wake_origin_v1/spot_comparison_v1/summary.json)
- [Short-control endpoint comparison](/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1/Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1/contact_failure_t2838_v1/positive_gap_v1/entropy_trace_v1/prevention_t265_v1/continuation_t350025_v1/wake_origin_v1/short_operator_control_v1/comparison.json)
- [Actual boundary geometry](/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1/Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1/contact_failure_t2838_v1/positive_gap_v1/entropy_trace_v1/prevention_t265_v1/continuation_t350025_v1/assessment/endpoint_geometry.json)
- [Original production endpoint audit](/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1/Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1/contact_failure_t2838_v1/positive_gap_v1/entropy_trace_v1/prevention_t265_v1/continuation_t350025_v1/endpoint/audit.json)
- [Diagnostic and frozen-source SHA256 manifest](/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1/Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1/contact_failure_t2838_v1/positive_gap_v1/entropy_trace_v1/prevention_t265_v1/continuation_t350025_v1/wake_origin_v1/evidence_manifest.json)

Frozen implementation / reference configurations:

- [Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1/contact_failure_t2838_v1/positive_gap_v1/entropy_trace_v1/volume_coordinate_v1/precision_stop_v1/symbolic_reuse_v1/native/source/Exam/NativeU64/gfs_material_cylinder_rhs2d.h](/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1/Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1/contact_failure_t2838_v1/positive_gap_v1/entropy_trace_v1/volume_coordinate_v1/precision_stop_v1/symbolic_reuse_v1/native/source/Exam/NativeU64/gfs_material_cylinder_rhs2d.h)
- [Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1/contact_failure_t2838_v1/positive_gap_v1/entropy_trace_v1/volume_coordinate_v1/precision_stop_v1/symbolic_reuse_v1/native/source/Exam/gfs_material_direct_flux.h](/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1/Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1/contact_failure_t2838_v1/positive_gap_v1/entropy_trace_v1/volume_coordinate_v1/precision_stop_v1/symbolic_reuse_v1/native/source/Exam/gfs_material_direct_flux.h)
- [Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1/contact_failure_t2838_v1/positive_gap_v1/entropy_trace_v1/volume_coordinate_v1/precision_stop_v1/symbolic_reuse_v1/native/source/Exam/gfs_material_image_rhs2d.h](/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1/Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1/contact_failure_t2838_v1/positive_gap_v1/entropy_trace_v1/volume_coordinate_v1/precision_stop_v1/symbolic_reuse_v1/native/source/Exam/gfs_material_image_rhs2d.h)
- [Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1/contact_failure_t2838_v1/positive_gap_v1/entropy_trace_v1/volume_coordinate_v1/precision_stop_v1/symbolic_reuse_v1/native/source/Exam/Cylinder/material_boundary.inc](/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1/Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_endpoint_compare_v1/wall_support_v1/transition_v1/full_near_v1/production_v1/spatial800_v1/contact_failure_t2838_v1/positive_gap_v1/entropy_trace_v1/volume_coordinate_v1/precision_stop_v1/symbolic_reuse_v1/native/source/Exam/Cylinder/material_boundary.inc)
- [Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_reference_su2_v1/reference_h0025_bdf2_barth_t5/cylinder.cfg](/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1/Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_reference_su2_v1/reference_h0025_bdf2_barth_t5/cylinder.cfg)
- [Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_reference_su2_v1/reference_h005_bdf2_barth_halfdt_t5/cylinder.cfg](/gpfs/kjhan/LagForce/Git_publish_native20pi_contact_20261104_v1/Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_reference_su2_v1/reference_h005_bdf2_barth_halfdt_t5/cylinder.cfg)

External diagnostic context (not evidence of our particular root cause):

- [SU2 official boundary documentation](https://su2code.github.io/docs_v7/Markers-and-BC/) defines the Euler slip wall.
- [Canivete Cuissa & Steiner, 2020](https://arxiv.org/abs/2005.12871) discusses why vorticity alone is sensitive to shear and why an additional rotation diagnostic is useful. Our finite-loop version has the limitations stated above.

The first material extraction attempt is preserved under material_history/failed.json:
the reader refused a reused exclusive-create temporary filename. material_history_v2
uses unique temporary outputs and extracted all20 states. No failed extraction was
silently replaced by data from another trajectory. Analysis scripts accompany each
artifact. The large pre-existing working-tree changes were not staged or committed.
