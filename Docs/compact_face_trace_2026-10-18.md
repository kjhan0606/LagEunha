# Saved donor face traces have small segment counts

Reuse the previously retained one-sided faces and exact donor/polygon
intersections at t2 for weight fractions .25 and1. Count intervals with constant
donor velocity, and compare density line quadrature orders2,4,8 against12.

Inner faces: 4768 retained,7876 segments, mean1.652 segments/face, maximum3.
Edge faces: 680/671 retained,773/766 segments, mean1.137/1.142, maximum3.
These counts apply only to retained unambiguous faces. About24% of eligible
one-sided faces are excluded due to coincident donor boundaries or coverage
ambiguity; no full-mesh memory estimate or universal three-segment bound follows.

With constant donor velocity v, fixed diagnostic face speed w, and normal n,
let I0=integral_segment rho ds. Then the segment mass flux is

    Fm = ((v-w).n) I0,
    Fp = v Fm,
    FK = |v|^2 Fm/2.

The same identity holds in a consistently transformed reference frame. Thus
for this field, one density integral and the donor velocity per segment encode
the advective mass/momentum/kinetic flux without truncating velocity moments.
The diagnostic evaluates the equivalent quadrature sums; it does not yet
implement a persistent compact face data structure or performance benchmark.

Worst retained-face results across fractions/regions/components:

| Density line quadrature | Max absolute flux difference vs12 | Max relative L2 difference |
|---|---:|---:|
| 2 | 1.216e-9 | 4.515e-7 |
| 4 | 6.574e-13 | 2.284e-10 |
| 8 | 6.632e-15 | 2.113e-12 |

The fitted density includes the smooth external-potential feature; it is not
generally polynomial. These are empirical comparisons, not exactness orders.
Using2-point quadrature merely because velocity is constant is not justified
at tighter relative flux tolerances. Adaptive density integration is an option.

## Next gate

For an affine normal face speed w.n=w0+w1*s and constant donor velocity, add
I1=integral rho*s ds. The mass flux then uses (v.n-w0)I0-w1 I1; momentum and
kinetic follow from the same shared Fm. This targets actual moving-face
kinematics without introducing a cell-volume third-moment closure. It still
requires a common face convention, donor-side selection, thermal/pressure
work, and treatment of coincident boundaries.

The donor fragments represent existing spatial information, not a law for
its physical evolution. Segment discovery, maintenance across topology
changes, growth over repeated remaps, and MPI ownership are unresolved costs.
Do not claim a fast complete solver or original-model reduction from this
static segment census. Existing long Kepler calculations remain unchanged.

Evidence: `canonical_compact_face_trace_2026-10-18.json`.
Script: `audit_compact_face_trace.py`.
Remote: checkpoint remap directory, `compact_face_trace_audit.json`.
