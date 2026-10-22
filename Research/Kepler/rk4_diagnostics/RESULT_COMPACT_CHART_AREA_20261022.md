# Curved target-chart area and positive quadrature — 2026-10-22

## Implemented

`material_compact_chart_area.h` now constructs positive physical-area nodes
over the complete Duffy chart, evaluates their positions from the same source
exterior curve, and checks the chart area against an independent Green-boundary
integral. It also provides positive 8×8 and 12×12 Gauss reference/physical
weights for later compact mass and field integration. Both constructions bind
to the same candidate-set, Jacobian, source epoch and weight seals.

The curve correction was rewritten relative to its endpoint chord. The cubic
part is factored as `du² t(t−1)` times an affine function of `t`; source-edge
analytic residuals and endpoint enclosures remain explicit. This is an
algebraically equivalent evaluation of the accepted curve and reduced
interval wrapping in the whole-chart Jacobian check. The native rezone object
compiles with an entrypoint that requires actual degree-one source coverage
and the same immutable fixed-mass law before constructing area rules.

## Synthetic result

A10 Slurm job `434702` completed, exit 0, in 10 s (fixture elapsed 8.510 s).
The curved-square fixture has 48 global chart patches, 18 curved patches and
analytic area `25/6 = 4.166666666666667`. Its certified minimum whole-chart
reference Jacobian lower bound is `0.266521756657919`.

The interval area nodes are positive and the independent interior/Green
enclosures both contain the analytic area. At the explicitly coarse
`0.02` per-patch interval budget, those sums are respectively
`[4.0966786307286034, 4.2339168148228294]` and
`[4.0496774463192073, 4.2886421381932598]`. These widths cannot support
an actual conservation tolerance. They remain geometry consistency checks.

The 8×8 rule used 3,072 positive points across all patches; the 12×12 rule
used 6,912. Both gave area center `4.1666666666666661`, agreeing with the
analytic synthetic polynomial boundary within displayed binary64 rounding.
This order agreement does not certify integration error for the actual
nonpolynomial source fields or analytic rank-1/rank-2 edges.

The fixture binary, native object, source hashes and Slurm receipt are in
`material_compact_chart_area_v4/build_contract_v4.json`. Earlier v1–v3
experiments are preserved; the v2 absolute-coordinate residual weakened the
interval Jacobian bound, while the final factored residual improved it.

## Physical-state boundary and next work

Accepted state remains step5 / `t=0.00011263377801805954`. No actual source
was imported, no source/target overlap packet was constructed, and no
density, heat, momentum, angular momentum or kinetic-energy projection or
physical step was run. The `t=6–10` shear-heating wall remains untested by this
candidate.

Next connect actual source-triangle/target-chart overlap routing to one common
positive packet stream. The true source field and reference-measure rho/e,
shared p/L/K projection, curved quadrature error budget, rounded native chart
bridge, commit/restart and amortized cost gate still need implementation.
Only the completed bounded transaction should use the accepted checkpoint.
