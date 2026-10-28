# H28 polar run: domain expansion, not inner-disk refinement

## Finding

I compared the accepted weighted 55,351-cell initial snapshot with the full
384² polar initial snapshot, using the native snapshot reader. The disk is
defined here by the same active annulus, (0.5\le r<2). Both inputs contain
exactly 49,038 generators in that annulus, with identical counts in every
sub-annulus below. Median cell areas also agree to the shown precision.

| Initial radial band | 55,351-cell input | 384² polar input | Median area, 55k | Median area, 384² |
|---|---:|---:|---:|---:|
| 0.5–0.75 | 3,931 | 3,931 | 0.000240378 | 0.000240437 |
| 0.75–1.0 | 6,050 | 6,050 | 0.000240379 | 0.000240370 |
| 1.0–1.5 | 15,865 | 15,865 | 0.000240386 | 0.000240377 |
| 1.5–2.0 | 23,192 | 23,192 | 0.000240364 | 0.000240360 |
| **0.5–2.0 total** | **49,038** | **49,038** | — | — |

The 384² input has 147,456 generators total and tiles a 6×6 domain. Its extra
generators are predominantly inside the central cavity and outside the active
disk: for example, the 2.35–2.75 band has 26,296 cells in the polar input and
257 in the sparse 55k input. This matches the prior H28 plan: retain the
GIZMO-effective interior resolution and widen the periodic halo/domain.

## Consequence for the running experiments

Jobs 437076, 437250, 437251, 437639, and 437640 form a W0/Voronoi full-domain
and halo-width chain. They remain useful for checking periodic-boundary
effects through (t=20), but they do not test whether more particles in the
inner disk suppress shear error. The staged t=10 and t=20 manifests have been
corrected to record this scope. They remain gated by each preceding 15-check
native orbit/energy/entropy/area audit and the zero-injection ledger.

At the latest verified Slurm read, job 437076 remained RUNNING toward
`t=2.221`; its log had reached `t=1.95547` (global step 26248). A nearby
complete ledger at `t=1.95290` reports `dEtot/|Etot0|=-2.523e-10`, zero
nonpositive internal-energy cells, and zero floor injection. The local/global
step offset remains 10,856, the accepted restart step. An independently
audited native `t=1.500026676` checkpoint passes the W0, fixed-mass, area,
positivity, EOS, and conservation checks; see the
[segmented-chain record](KEPLER_Q2_HALO28_SEGMENTED_T65_2026-10-27.md). This
does not yet provide a completed `t=2.221` endpoint receipt.

## Remaining approved test

The approved inner-resolution hypothesis still needs its own initial
condition. The next run must increase the generator count specifically inside
the active inner disk while keeping the outer domain, halo width, gravity,
EOS, and flux operator fixed. It must preserve the intended radial mass and
pressure profiles, keep each cell's mass fixed after initialization, begin
with (W=0) for a Voronoi control, and pass the same complete-snapshot,
positive-state, conservation, orbit, and density gates at (t=6.5). A global
384²-to-larger-grid change would also raise inner counts, but it is a different
and more expensive comparison than local refinement.

The existing polar ring-spacing and hybrid-polar-glass generators rearrange
or polarize a fixed total point count; they do not by themselves establish a
larger inner-disk count. No refined IC has been generated or submitted yet.
The proposed matched local-refinement profile and its runtime gate are
recorded in the [inner-resolution test design](KEPLER_INNER_RESOLUTION_TEST_DESIGN_2026-10-27.md).
