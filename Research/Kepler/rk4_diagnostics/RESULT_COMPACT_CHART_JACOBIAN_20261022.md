# Curved chart disk, Jacobian and degree gates — 2026-10-22

Later same-day work refactored the interval Jacobian evaluation and added
positive chart-area quadrature. The v23 numbers below remain the frozen v23
historical result; see `RESULT_COMPACT_CHART_AREA_20261022.md` and its v4
contract for the current synthetic bound and native area component.

## Why this stage exists

This is a preprocessing stage for a conservative remap from the accepted
checkpoint. It constructs target-cell geometry and the maps needed to transfer
mass and state fields; it does not advance physical time or restart RK4 yet.
The purpose is to prevent a folded or incomplete chart from contaminating the
Kepler evolution that later tests the artificial shear heating after `t=6`.

## Work completed

- Added a Q128 outward-interval bound for the Jacobian on the entire curved
  collar reference triangle. For
  `G(s,t)=A+s[(1-t)B+tC]+s^4[Gamma(t)-((1-t)B+tC)]`, the reference area factor
  `s` is removed symbolically before interval evaluation. The verifier
  subdivides the full `[0,1]^2` parameter square and accepts only boxes with a
  strictly positive Jacobian lower bound. No grid sampling can pass the gate.
- Straight patches have the exact affine value `J=1`.
- Bound the proof to the actual source-curve descriptor and root-parameter
  enclosures retained by the curved-domain chart.
- Added a candidate-set seal to the global connectivity certificate. The
  Jacobian wrapper checks that the node enclosures, patch edges, source epoch,
  weight epoch and fixed-M law are still the same as those certified earlier.
  The fixture mutates a node after connectivity and expects rejection.
- Added a global manifold guard: every interior vertex link must be one cycle,
  every boundary link one path with two exterior edges, and the connected
  complex must have Euler characteristic 1. The fixture has 33 vertices,
  80 edges, 48 global patches and 16 boundary edges; all 33 vertex stars pass.
- Added a degree-one gate that combines that disk certificate, the positive
  chart Jacobians, and the globally injective positive-J accepted-source
  geometry certificate. A fixture with no accepted-source certificate is
  rejected. It does not turn synthetic geometry into an actual coverage claim.
- Compiled both the synthetic fixture and native rezone component from the
  final frozen v23 source set. The frozen set contains 35 source files and 45
  project dependency hashes; its fixture binary and native object reproduce
  byte-for-byte.

## Runtime evidence

H200 job `434680` could not start because of `QOSGrpGRES` and was cancelled.
A10 job `434682` exposed a missing `intel/tbb` module prerequisite before the
fixture launched. A10 job `434683` ran an earlier Jacobian-only version; A10
job `434687` ran an earlier disk/degree-gate version. Job `434686` was a
script path error and did not execute the fixture. Two exploratory identity
square fixtures (`434688`, `434691`) were rejected on unresolved coincident or
endpoint roots, as designed. The final frozen v23 regression, Slurm job
`434692`, completed on `syn01` in 5 s (Bash elapsed 2.423 s, max RSS 3772 KiB).
It retained the fixture dimensions: 5 partition cells, 24 arcs, 8 curved arcs,
8 shared straight arcs and 48 global chart patches.

The continuous interval verifier covered the full charts with 144 boxes and
64 subdivisions. Every reference Jacobian lower bound was positive; the
minimum certified lower bound was `0.2132873241228527`. The chart has
`V-E+F=33-80+48=1`, and all 33 vertex links passed. The altered chart set was
rejected after connectivity certification; the unbound source-geometry
degree gate was also rejected. This result is synthetic only: the fixture
reports `degree_one_physical_coverage=false`, with no actual source import or
physical step.

## Physical-state boundary

The accepted state remains step5 / `t=0.00011263377801805954`. This work has
not imported the accepted source, proven its actual target partition, run the
projection, or advanced a physical step. It therefore provides no evidence
that the `t=6–10` Kepler wall is solved.

## Next gates

1. Complete the missing conservative source/target overlap integration,
   positive rho/e transfer, coupled shared p/L/K projection, field/entropy
   budgets, compact commit/restart wiring, and full transaction cost estimate.
2. Once that executable and cost gate are complete, run one bounded transaction
   on the actual accepted source. Bind the disk/Jacobian/degree proofs to that
   source epoch and exterior registry inside the transaction, then perform
   overlap, projection and conditional commit. Do not run a geometry-only
   source import or replay old RK stages as a separate job.

## Frozen build

`material_compact_rezone_build_v23/build_contract_v23.json` records 35 frozen
source files, 45 project dependency hashes, compile/runtime receipts,
binary/object hashes, and the original repository HEAD.
