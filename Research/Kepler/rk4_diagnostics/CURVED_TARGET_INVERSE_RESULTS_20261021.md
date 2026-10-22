# Target-local inverse integration: same tiny target passes

2026-10-21. Driver follows `ASTRA_TARGET_INVERSE_PLAN_20261021.md` without
another agent review. The completed original-error bundle remains in
`CURVED_ORIGINAL_ERROR_RESULTS_20261021.md`; none of its failures is erased.

## Result and scope

**The SAME lens target101 passes all12 original moment criteria** under the
stated IEEE / eight-ulp-libm engineering model. Its worst total error bound
is.004075452596 times the UNCHANGED tolerance. The previous polar/ray error
model for this target was211.5989 times tolerance. This completes the FIRST
tiny-target gate of the selected inverse-map bundle.

This is a native C++/MPI contract run on lageunha: binary80 FMA field
evaluation with fixed binary128 geometry/inverse predicates. It is not a
Python simulation, CUDA kernel test, full source/all-target transfer or
physical time evolution. The physical225/137 sources, shared-edge gluing,
source baselines and eventual CUDA/MPI t6/8/10 / GIZMO longest duration are
still outstanding. No full24, weight solve or physical timestep released.

## New expression, unchanged physical reference

The actual physical power polygon is constructed from the original stored
site/weight binary80 constants. Normals and constants are computed in a
target-centered expression, without subtracting two source-scale squared
positions. Binary128 ball arithmetic encloses exact plane construction and
vertex intersections. The three vertices satisfy the other plane strictly;
the resulting positively oriented bounded triangle is the entire tiny cell.

Coordinates are x=c+ell*xi, ell=2^-12. At each physical point, invert the
ORIGINAL stored P3 map. The lens actually has zero residual coefficients,
so its original mathematical map is a cubic; original-basis algebraic
coefficient expansion is exact up to enclosed arithmetic, not a fit or
new physical state. Its coefficients are expanded around the inverse center
and evaluated by binary80 FMA Horner. Conversion of each original Q128
coefficient interval to the native binary80 midpoint is charged explicitly.
Independent original nodal-basis evaluation checks every new moment.

No degree6 source-ray conversion, polar root or128-epsilon angular guard is
called on this new path. Those error terms disappear because the operations
are absent. New plane, inverse, native field, quadrature-node/weight,
geometry and summation errors are included instead.

## Inverse, positivity and source coverage

One complex infinity-norm box of reference radius.01 gives a uniform
preconditioned contraction for the whole physical target and the analytic
quadrature remainder domain:

- contraction upper:.00100000000000000359;
- invariant-box image radius upper:.00190525588832576529 <.01;
- chi-I Lipschitz upper on the WHOLE convex reference source:.10000000000000155 <1;
- local J lower:.999999999999999803735;
- inverse box's reference-boundary margin:.2600000000000000133;
- maximum point inverse radius:5.543561331830813e-33.

The first two prove inverse existence and uniqueness in the box; the whole
source Lipschitz bound excludes a second source preimage. The entire box is
strictly in the exact original unit source triangle. Consequently this
target encounters no shared material edge and uses no source-minus-outgoing
baseline. Source-boundary/gluing ambiguity and source-baseline errors are
zero HERE for this proved reason. This does not certify them for225/137.

## Actual integral remainder, not rule-difference agreement

Use ONE16x16 positive Duffy triangle rule,256 points. All16 Gauss roots are
certified in disjoint Q128 sign brackets with nonzero derivative; positive
weights and their uncertainties are propagated. No order/cap/tolerance scan.

For x=A+u(B-A)+(1-u)v(C-A), radius1 Cauchy discs about every real u,v in
[0,1] are covered by |u|,|v|,|1-u|<=2. The same uniform complex inverse
contraction applies there. Polynomial J is separated from0 and its negative
real branch; the potential argument x^2+y^2+eps^2 has a positive real-part
lower bound4.169801124845304. Absolute holomorphic integrand bounds include
the Duffy Jacobian, rho/J, velocity, pressure, Udot and potential.

The Gauss remainder formula and Legendre norm are given in
[NIST DLMF3.5.19/21](https://dlmf.nist.gov/3.5#E19). Scaling to[0,1] and
applying the Cauchy derivative bound gives C16=(16!)^4/[33*(32!)^2], with
C16<=8.387276125592904e-20. The two tensor directions are bounded by
2*C16*M for the certified complex integrand bound M. Node/weight arithmetic
uncertainty is charged separately. Agreement of two quadrature rules is
not used as an integration proof.

These analytic remainder bounds retain the stated finite-arithmetic/libm
model; they are not a formal proof of the CPU or GPU elementary library.

## Target moments and unchanged normalizer

Existing U_lower=6.9848668650464756e-13 and M_lower=4.6566128730773926e-10
are unchanged. No overlapping proof was added or normalizer enlarged.
Energy tolerance=6.9848668650464756e-23. Other moments use the existing
1e-12 mass/length/velocity form with this conservative target mass lower
bound; Udot uses1e-9 U_lower.

| Moment | Native integral | Total error bound / tolerance |
| --- | ---: | ---: |
| mass | 5.1961524227066334e-8 | .00112399 |
| angular momentum | 1.8565332991088532e-9 | .00050048 |
| K | 1.2990381056766585e-9 | .00019540 |
| U | 7.7942286340599504e-11 | .00001308 |
| potential | -2.5441784210458727e-8 | .00407545 |

All12 pass, including linear momentum, volume/position moments and Udot.
For potential, native field/geometry/node/sum error is2.6759094338922211e-25;
analytic integral remainder1.7073994574088427e-26; total2.8466493796331054e-25.
Original-basis same-rule oracle excesses:0. Positive cell area is
5.1961524227066324e-8 with Q128 area radius2.34875e-38.

## Preserved old calculation versus the continuous reference

`CURVED_TARGET_INVERSE_COMPARISON_20261021.json` compares exact decoded
binary80 old/new moments. The NEW integral enclosure now bounds the old
nominal calculation's deviation from the same continuous target integral.
The maximum old error upper / tolerance is.07471392291, for potential.
Its old-minus-new native difference is-4.934003106908099e-24.

Thus the old lens nominal happens to be accurate enough on this fixture,
while its old error MODEL could not establish that. No observed error is
substituted for a bound; this conclusion uses the new integral enclosure.
It says nothing about long physical Kepler evolution or the broken shared
traces of the actual225/137 state.

## Reproducibility / exact next work

- Final run:`curved_target_inverse_v2_lens`, build/native exit0, .8204s;
  frozen envelope:`curved_launch_RYdY3Rf5`.
- First run v1 is retained: same native moments and exact same positive
  binary packet SHA; v2 adds outward treatment of three proof-margin
  scalar expressions. No quadrature parameter or physical state changes.
- New packet SHA256:
  `53d8aad29e2d3643a74f64d31ce28c6aff777f89b70e819908c844e48b6ab4c9`.
  Header6uint64, magic0x474654494e563031/version1/target101/count256/
  sizeof(SavedPoint)/Q width16. Records are native Q128 physical/ref
  coordinates, positive physical weight and input radii, then12 native
  binary80 moments/radii. Host ABI is recorded by frozen source; CSV is
  also retained. This is separate from the old56-byte quadrature packet.
- Source, original packet and old positive/moment/normalizer SHA checks pass.
  Both inverse runs are terminal. No monitor/background process exists.
- Files:`material_target_inverse.h`, `test_curved_target_inverse_v1.cpp`,
  immutable launcher `run_curved_target_inverse_v1.sh`.
- Contract, inverse certificate, all target moments and new positive packet
  are in the final run directory. The original-basis helper gained a
  Q128-coordinate entry; the old double-coordinate wrapper is preserved.

Continue the SAME Astra bundle with physical225/137 affected targets:
original nonzero-residual P3+initial-v0 map, certified inverse and reference
source membership, local target boundary intersections, full source/gluing
ledger. The current zero-residual/unit-source helper must reject those
inputs until that actual general path is implemented. Do not apply its
polynomial proof to a nonzero-residual source or call a tiny-target pass a
full physical transfer. No additional user approval or Astra microreview.
