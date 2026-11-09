# Actual-cell primitive positivity diagnostic — 2026-11-09

This is a saved-state operator experiment, not an adopted hydro method. Canonical C/CUDA sources and running RT evolution are unchanged. The successful regular/glass Sedov comparisons use the original frozen native operator; their figures must not be attributed to this experiment.

## Failure mechanism isolated

The archived `common_primitive_cell_geometry_20261109_v1` replays the primitive-v3 QR reconstruction and the original scalar pair projection. At the early glass failure, face354–4018 has owner pressures near1e-6, but side4018 reconstructs pressure−1.8750563 at the actual face and−0.13811974 at its tangential projection. Both tangential sample points are inside their actual cells, so outside-cell extrapolation is not required for this failure. The pair projection clips its reversed jump to zero and clamps the negative common mean to9.094947e-19. Positive cell averages therefore do not imply useful positive face traces.

Tangential points can nevertheless lie outside other deformed cells:60 such faces in the frozen control-glass endpoint,123 in the rejected log-candidate endpoint and92 in the regular endpoint. The raw face and tangential samples use the same physical centroid coordinates and the same gradient. The pressure-offset LS displacement is not a replacement sampling coordinate.

## Isolated candidate

After the unchanged primitive QR fit, multiply each cell scalar gradient by one factor alpha. Determine alpha from the minimum linear increment over the actual convex polygon vertices; when necessary, take the representable value immediately inside the mathematical positivity bound. Then check all actual floating-point vertex evaluations. A finite64-ULP inward search aborts if strict positivity cannot be verified; this is a roundoff check, not a chosen physical pressure floor. The existing pair projection, velocity map, exact Riemann response and force/work formulas are unchanged.

A linear field reaches its minimum at a convex polygon vertex. Scaling around the volume centroid preserves its scalar spatial mean. A positive affine field requiring no restriction is unchanged. This does not prove positivity after a time step, guarantee a non-vacuum Riemann solution, or preserve a common pressure acceleration when different cells require different factors. No new physical viscosity, pressure floor or sensor is introduced.

The edge endpoints are recovered from the physical face midpoint plus/minus half the face length times the tangent to the original power-plane normal. Reconstructed-centroid displacement is used for sampling, not to redefine that normal. This code is only exercised on the existing certified periodic2D fixtures.

## Completed frozen-state checks

`common_primitive_vertex_positive_20261109_v1/{plan.json,run.py,headers,report.json}` records exact build commands, input hashes, scalar traces and RHS arrays. With the new switch disabled, all five original primitive-v3 face and RHS arrays match bit-for-bit. These are operator evaluations, not evolved endpoints of the new candidate.

|Frozen input|Nonpositive vertex cells before|Riemann failures before|After|Minimum face pressure after|
|---|---:|---:|---:|---:|
|glass_initial|30|0|0|9.4278555e-07|
|glass_early_failed|42|2|0|9.4323358e-07|
|control_glass_t01|463|4|0|8.2740458e-07|
|candidate_glass_t01|517|3|0|8.9690634e-07|
|control_regular_t01|452|0|0|9.9531237e-07|

All tested raw face pressures and vertex pressures become positive, and the nearly zero final face pressures disappear in these five fixtures. Nine failing faces across three frozen inputs now have material-contact solutions. The candidate initial glass acoustic dissipation rate is95.02576 rather than the original primitive-v3 value40.04240. Such a substantial change must be judged by evolved shock accuracy; it is not evidence of a physical heating improvement on its own. The summed force and energy rates remain at roundoff in these periodic operator tests.

## Disappearing-face check

Eight predeclared fixed random positive pressure fields on the existing certified shear geometries were evaluated at both signs of1e-4,1e-6,1e-8 around a vanishing face. The largest relative force changes for the candidate are2.703893e-4,2.704138e-6,2.704141e-8, respectively. Thus the measured force difference decreases with the geometric perturbation; the finite jump from hard-neighbor extrema is not reproduced here. This finite test does not prove all topology transitions or rank-degenerate stencils continuous.

## Limits and next gate

Do not adopt the candidate from these results alone. Positive traces may approach vacuum; cell-dependent scaling can spoil common acceleration balance. The2D local formula has not yet passed a common1D implementation/evolution, affine/free-stream regression, C/CUDA parity, boundary or full2D evolved accuracy check. These are the next gates before any native benchmark. All3D evolution remains deferred until2D is complete. The original working method and all rejected candidates remain preserved.

Independent RT monitoring during this diagnostic verified the original control reachedt=3.3 with fixed masses, positive state, gravity-energy relative change−9.60417e-15 and maximum volume mismatch9.61867e-12. Its existing supervisor advanced to the3.6 segment without intervention. This is separate from the new reconstruction experiment.


## Subsequent full-method decision

The common helper was implemented in the existing1D and native C/CUDA paths, tested through all9 original1D endpoints and both64-square Sedov endpoints. It is rejected for worsened2D density/pressure error and morphology; see `COMMON_VERTEX_PRIMITIVE_DECISION_2026-11-09.md`. The earlier frozen-state successes do not authorize adoption.
