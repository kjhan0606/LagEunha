# Stable reference input check before remaining native integration

The objective remains long Kepler evolution beyond t6--10 in the original
CUDA/MPI/OpenMP solver. This check is a prerequisite for the pending transfer,
not a demonstration of physical stability or a new physical model.

Use the fixed v7 remaining panel plans and already saved regular/oblique real
Gauss branch values. Reuse all completed whole-target results and original
target68526. No new event roots, arrangements, full fields or physical moment
queries are evaluated. The computation runs only on lageunha with OMP16/MPI1.

The input-only sample includes the first/last Gauss nodes on each axis of
every remaining whole panel and the first/last eta panels and eta Gauss nodes
at every already saved implicit free node. It is explicitly not exhaustive
coverage of the remaining21629184 physical queries. Unselected clipped
panels with no implicit free group are not represented by this sample.

`material_shared_runtime_reference.h` expresses the original reference point,
Jacobian and barycentric coordinates from the retained Duffy/affine chart.
Its Cartesian coordinate expression matches the original sd::basis graph.
Regular uncertain roots and oblique branch corrections remain intervals.

At the same actual point compare:

1. Existing binary80 Cartesian inverse for the three barycentric centers.
2. Direct conversion of the original chart's barycentric centers to binary80.

Both use `r80::cache(center, original_chart_lambda)` to retain the input
uncertainty. Thus the old path is given the more favorable original chart
enclosure; this comparison does not claim the original inverse enclosure
would be tighter. Check each original triangle endpoint pair sum has a
strictly positive lower bound before its edge quotient. Save every failed
old or direct input, including original points, lambda balls and pair sums.
Keep field guards, floating precision, order, profile and tau unchanged.

This checks a necessary input condition only. Positive pair sums do not prove
native field accuracy or all12 forward-error acceptance. A later optional
direct-lambda evaluator must use the same canonical function and gradient
caches, preserve the legacy default, and carry original point/input errors.

No state mutation, heat reset, output replacement, accepted RK4 step or
restart of previous physical calculations is authorized by this input check.

The first runner (v1) stopped before compilation because an input basename
was `real_boundary_values` instead of the actual
`real_regular_boundary_values`. Its immutable envelope and partial input
manifest are retained. The corrected v2 uses the verified producer manifest;
no successful numerical input checks were repeated.
