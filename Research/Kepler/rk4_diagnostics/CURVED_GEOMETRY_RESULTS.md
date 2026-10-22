# Actual curved geometry: interval bounds and complete angular coverage

2026-10-21. Follow-on implementation under ASTRA_CURVED_REMAP_PLAN.md. This
is progress toward the actual curved remap, not a completed remap, weight
solve or time evolution. The preceding opinion-only turn changed no runtime
state; this turn adds native geometry code and authoritative execution data.

## Preserved inputs and finite state

Use `material_p3_finite_v2/mpi1_cpu/curved_map_control.bin` and
`curved_map_algebra_v2/advected_sites.csv`. The control is the SAME actual
current P3/h2 chi,v at h=2e-4, not a new initial condition or future exact
orbit. Inputs and each frozen executable are hashed. Original native runtime
25 source hashes still match `enriched_v5_build/source.sha256` as checked
this turn. No historical restart/output was changed.

## Local map bounds and conditional boundary geometry

`material_curved_bounds.h` encloses the analytic map and derivative over
each whole reference triangle using outward long-double interval arithmetic.
Complete affine reproduction is removed algebraically before evaluating
the higher modes; this reduces dependency overestimation without changing
the map. Current coefficients, stored geometry and residual nodal data enter
the enclosures. An additional explicit engineering envelope for native
double evaluation is recorded separately (position about9.1e-13, gradient
about3.69e-12). These engineering factors are not a formal proof of every
operation in the CPU/CUDA implementation or of libm correctness.

The successful `curved_bounds_v3` resumes the complete384-row interval
ledger from v2 rather than rerunning that calculation. Packet hash is
rechecked and parent-ledger hash recorded. It then checks native samples,
boundary geometry and the inverted-map fixture. Results:

- All384 triangle enclosures complete without subdivision at the final
  affine-separated formulation. Interval J lies in[0.8916673,1.1104773].
- Maximum norm bound on grad(chi-X) is0.130586.
-9600 native point/gradient queries lie inside their explicit envelopes.
- Reference exterior has5 simple loops,80 edges;536 internal shared edges.
-3080 nonadjacent exterior-edge pairs have positive separation margin,
  minimum0.0041382.80 adjacent pairs have corner contraction factor at
  most0.184676<1. Disconnected components remain separated by these bounds.
- Native shared positions differ by at most1.78e-15 in the sampled trace
  checks. The inverted-map fixture is rejected.

**Global injectivity is NOT proved.** The boundary argument is conditional
on exact shared-trace gluing; the current continuity check is sampled.
An analytic trace mismatch bound/continuous representative with quantified
geometry error remains required before calling this a complete geometric
certificate or promoting a Newton solve.

### Retained failures and recovery

v1 and v2 terminated with exit2 at the final inverted-map assertion, after
their actual interval ledgers were written. The fixture result aggregator
initialized max_J=0, so it could not report an entirely negative J range.
It now starts at negative infinity. v1 also used a looser polynomial
interval expression; v2 extracted the exact affine part. All outputs remain
unchanged. v3 resumes v2's ledger, recomputes only the small independent
sample/boundary/fixture checks, and terminates exit0/AUDIT_COMPLETE.

## Interval ray polynomials and angular partition

`material_curved_roots.h` constructs interval coefficients of degree<=6
power-plane numerators. Two charts centered at0 andpi with u in[-1,1]
cover the reference circle using exactly0,+1,-1 trigonometric constants.
Interval radius is supported for the forthcoming radial-panel enclosure.
Coefficient interval arithmetic is retained through every polynomial
operation. Physical position uncertainty widens the plane by its actual
normal times the explicit native evaluation envelope.

Bernstein subdivision certifies sign where possible and retains ALL
remaining root bands. Endpoints, even roots, nearby roots and flat regions
are never removed by a face-length cutoff. At a fixed cap/budget an unresolved
interval remains in the coverage set and the run reports failure. Root bands
can include several roots; this implementation does not falsely equate the
number of bands to the number of roots or claim simple-root uniqueness.

Native `curved_roots_v1` uses OpenMP2 and MPI1/4 on authorized lageunha:

- Five bounded fixtures cover six distinct roots, close roots, an even
  tangent root, chart endpoints, and a small positive-width lens.
- One centroid reference radius per each of384 actual triangles, with
  W=0 and a small prescribed nonzero weight pattern. No weight Newton yet.
-768 source angular intervals;90 candidate pair polynomials requiring
  isolation;90 retained target-root bands;858 classified angular fragments,
  including90 whose target differs from their source cell.
-167494 total root panels, zero budget failures, zero native target-label
  mismatches and zero native power-plane enclosure excess.
- MPI1/4 angular-partition CSV files bitmatch exactly. Warm wall times are
  0.848s and0.944s respectively; no strong-scaling improvement is claimed.
- Maximum retained uncertain angular fraction bound is4.74e-8. This is
  an ANGULAR coverage bound at the sampled radii, not a mass/K/U accuracy
  certificate. In particular it has not passed cold K-error/U requirements.

The current working source adds a FUNNELED-safe fatal path after this frozen
successful run: OpenMP workers terminate their process on an error instead
of invoking MPI_Abort from a worker. No successful arithmetic is changed.

## Next required implementation

1. Replace sampled shared-trace continuity by an analytic bound and account
   for any representation/gluing uncertainty in transported moments.
2. Use the interval-radius coefficients to detect or bound missing radial
   topology events. A radial quadrature comparison alone is insufficient.
3. Integrate actual common-refinement m,p,L,K,U,V and report missing/uncertain
   contributions, including absolute K error normalized by cold U. Record
   the actual M_i(0)-m_i0 witness; W0 fixed-mass reduction is still unresolved.
4. Then finite screened weight solve/coarea Jacobian, fragment transfer,
   one admissible cold compression and one reverse remap, per Astra's plan.

After that comes sparse/matrix-free native composition, density-wave/sliding
and shock gates, and full CUDA/MPI Kepler evolution through t6/8/10 and the
longest GIZMO paper duration. No short geometry control substitutes for those
requirements. No new Kepler time evolution or Slurm job occurred this turn.
