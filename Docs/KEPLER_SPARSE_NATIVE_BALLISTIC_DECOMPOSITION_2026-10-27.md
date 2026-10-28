# Sparse Kepler: native density versus exact-orbit Voronoi geometry

The native `t=10` sparse gate reports disk density change but does not
separate the change produced by shearing fixed-mass Voronoi cells on exact
softened circular orbits. The independent
[diagnostic](audit_kepler_sparse_native_ballistic_20261027.py) makes that
comparison by cell ID. It reads complete native `t=0` and endpoint
snapshots, checks fixed mass, W=0, positive state, `rho=m/V`, and periodic box area,
and constructs unweighted Voronoi polygons from **all initial generators
moved analytically** to the endpoint time. Every selected disk polygon
must be bounded and pass `distance to box edge > 2 × circumradius`, so
omitted periodic images cannot affect its area.

The diagnostic reports three different density quantities: native change
from the initial state, exact-orbit geometric change from the initial native
cell area, and the **pointwise** native-versus-exact-orbit density residual.
The difference of the first two absolute means is not used as a residual.
It also reports position error against the analytic orbit. This is a
geometric diagnostic; entropy, energy, heating, and long-time stability
still require the separate native physics gates.

## Check on accepted native t=.825 sparse snapshots

| Initial disk, mass weighted | Wide, 55,931 sites | Narrow, 54,312 sites |
|---|---:|---:|
| Native absolute density change | `1.1652064e-4` | `1.1659595e-4` |
| Exact-orbit geometric density change | `1.1654484e-4` | `1.1654484e-4` |
| Native-versus-orbit absolute density residual | `8.47206e-8` | `2.31479e-7` |
| Native-versus-orbit position error | `2.44020e-9` | `4.12776e-9` |

Both [wide](KEPLER_SPARSE_WIDE_NATIVE_BALLISTIC_T0825_STATECHECK_2026-10-27.json)
and [narrow](KEPLER_SPARSE_NARROW_NATIVE_BALLISTIC_T0825_STATECHECK_2026-10-27.json)
read-only receipts passed their geometry prerequisites and periodic-image
safety check. Their maximum relative `rho=m/V` error is `5.90e-8`, and the
minimum omitted-image safety ratio is `41.08` in each case. The geometric
values agree with the independently computed
[long-time exact-orbit table](KEPLER_SPARSE_FULL_LONG_GEOMETRY_2026-10-27.json)
at `.825`; the native density norms agree with the existing
[same-ID comparison](KEPLER_SPARSE_NATIVE_T0825_SPEED_ACCURACY_2026-10-27.md).
The early native residual is small, but no `t>6` result follows.

After an A10 or H200 `t=10` sparse endpoint is complete and its physics and
state audits pass, run this diagnostic on that exact snapshot. Compare the
pointwise density residual, its radial bins, and position error with the
early baseline. A density rise that is explained by exact-orbit geometry
has a different interpretation from a growing native-versus-orbit residual;
neither quantity alone measures irreversible heating.

The first ongoing-A10 output at actual `t≈1.650242` has now been checked
with both the state and geometry auditors; see the
[wide/narrow intermediate comparison](KEPLER_SPARSE_A10_T165024_INTERMEDIATE_2026-10-27.md).
