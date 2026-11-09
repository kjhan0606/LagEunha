# Fine fixed-wall cylinder: transient swept-CFL event

## Status and scope

The preceding goal turn made progress: archived Sedov IC reproduction, exact source-moment audit, and a verified finite waiting runner. This turn preserves all live trajectories and adds read-only evidence on a newly observed cylinder slowdown. Native PID 4145694 remained live. No timestep coefficient, solver tolerance, force, weight, or mesh velocity was changed. No restart was launched.

## Observed event

The 126,696-cell fixed-wall run slowed near t=1.496. In [1.496,1.497), 63 accepted steps contain 55 retries, minimum dt=1.961963227863163e-6, median dt=6.61868017972858e-6, and approximately 321 seconds of recorded accepted-step wall time. The two geometric limiting IDs are 64121 and 63321. The common geometric coefficient remains 0.1. Its bound uses cell volume divided by the face-integrated absolute normal sweep relative to the generator, including affine face rotation. It is not an acoustic signal-speed cap.

The trajectory passed the event without intervention. At the captured t=1.5016547846295663 state, dt returned to .0005 with zero retries; wall pressure work remains zero. The contact kinetic-energy term is approximately -1.915e-5, so a small ledger residual must not be called exact total energy conservation.

## Saved-state geometry

The new accepted checkpoint cylmat.003217.dat has SHA256 74d8a26d884a70a570dba36be2c988c2bedd02bced9ada7b9404f1185f331e55. Its critical generators lie near (4.0991,1.9996), downstream of the radius-.5 cylinder centered at (2.5,2). Their distance is 1.2923858451700332e-4, or .00516954338 of nominal h=.025. Their shared face remains .03447109548 long and both native volumes are about 6.15e-4.

Independent power half-plane reconstruction of these two interior cells certifies omitted owned sites using periodic-y KD candidates and matches native volumes within 5.1e-14 relative. Both polygons are convex. Generator 64121 is outside its own power cell; 63321 is inside. The relative weight ratio |Wi-Wj|/d^2 is 28.13685. This is consistent with the earlier cell-survival audit: positive Laguerre volume does not require generator containment. It is not by itself a geometry failure, nor grounds for silently reinstating a weight clip.

At the saved state the pair is separating (radial separation rate .0234429); the bisector orientation rate from generator velocities is 1.89648 per unit time. These values are AFTER the minimum-dt event. The existing log has no facewise weight-rate/sweep decomposition, so this does not prove that pair proximity alone caused the minimum timestep. The minimal wall gap is not the apparent trigger here: these are interior downstream cells, not obstacle cut cells.

## Existing implementation and history

The frozen cylinder adapter sets mesh[i]=state[i], using each cell mean velocity for its generator. It does not locate the physical cell containing that generator. Docs/UNIFIED_HYDRO_KEP_FILTER_POINT_2026-11-05.md already records a different experimental FV path that evaluates point velocities in their containing cell. Its full strong-shock tests did NOT all pass. This is prior evidence to review, not a qualified cure to transplant or an excuse to discard the successful Kepler lineage. The current source also has common swept-geometry safeguards; simply disabling them would remove the measurement that detected this event.

Next diagnosis should isolate weight-rate and relative generator motion contributions using retained states and the existing operators, and assess whether a single physical velocity field can define generator motion without steering or new pressure terms. Any candidate must preserve the force/work/volume relation and be compared across the already required tests. The live t3.50025 run remains unchanged; its final SU2 density comparison is still pending.

Artifacts: Research/MaterialVolume_20261105/material_pressure_acceleration_20261107_v1/cylinder_stationary_wall_20261109_v1/cfl_event_t1497_v1/{report.json,pair_geometry.json,audit.py,pair_geometry.py,observed_log.txt,state_ids.bin}.

## Power ownership and cell-mean versus point velocity

The previous goal turn made progress by completing the KH1.2 and RT.3 assessments and qualifying synchronized GIZMO comparison. This continuation adds read-only transport diagnostics, not another physical model.

All126,696 owned generators were checked against periodic-y owned sites with certified omitted-candidate bounds. At t=.25025,209 generators are demonstrably outside their own power cells; at t=1.5016547846 this count is1,721 (1.35837%). Within the fixed downstream ROI x in(3.2,8),y in(1,3), the counts are4/15,320 and87/13,334 respectively. The best competing OWNED site suffices to disprove own-cell containment; unqueried prescribed exterior sites could still change the global containing owner near open boundaries. No global exterior-complete owner claim is made. Candidate bounds remain positive, and no ties lie in the roundoff ambiguity band. This CPU diagnostic took approximately .61s and3.34s; it is not the cost of a production stencil walk or a hydrodynamic speedup.

For the retained downstream pair, both generator positions lie in power cell63321. Their actually assigned cell-mean velocity difference has magnitude .0234441895. A read-only reconstruction uses the frozen native face-length/centroid-distance least-squares metric, independently rebuilt actual neighboring cell polygons, and the same centroid offsets. In the containing cell's UNLIMITED affine velocity field, the difference between those same two positions is2.2764265e-5, a ratio .000970998. The two local normal matrices are well conditioned (1.47 and1.24), and constant/affine reproduction checks pass at rounding precision. This narrows a candidate mechanism: label-attached cell mean velocities need not approximate a single pointwise velocity field when generators nearly coincide.

This is NOT the actual limited Riemann reconstruction, a proof of the event's cause, or a successful new time integration. The native method prepares a face-specific normal limiter (shared ordered-jump factor per face), whereas the older FV point-advection experiment used a cellwise limited velocity field. Thus its point locator cannot simply be transplanted by replacing mean velocities with a containing-cell mean or an unrestricted extrapolation. Such a replacement would conflate geometric generator motion with the conserved cell momentum. A candidate must define generator transport separately, keep cell masses/momenta/energies intact, and use that same geometric transport in volume and pressure-work coupling. It must retain the existing no-steering constraint and pass affine, cold-shear, contact and shock checks before any adoption.

No pair repulsion, weight cap, new viscosity, mean-velocity replacement, or live binary change was made. Original0.25 and1.50 checkpoints remain intact. The first local-fit postprocessor had a report-loop indentation error and exited before writing a result; its artifacts and failure.txt were retained. local_velocity_fit_v2 is the corrected completed diagnostic. Evidence: cfl_event_t1497_v1/point_ownership/report.json and local_velocity_fit_v2/report.json, with exact C helper, build command, hashes and source scripts.
