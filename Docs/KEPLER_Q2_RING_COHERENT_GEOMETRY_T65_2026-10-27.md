# Bounded q2 ring-coherent IC: exact-orbit geometry at t=6.5

The [IC manifest](KEPLER_Q2_BOUNDED_RINGCOHERENT005_R275_N1024_2026-10-27.json)
defines one random angular offset, at most 5% of that ring's site spacing,
shared by every site in each of 70 polar rings. Radii, ring populations,
exterior sites, pressure, and orbital prescription match the accepted q2
bounded controls. Its independently calculated initial Voronoi areas set the
fixed cell masses. The minimum site separation remains positive.

The [read-only exact-orbit audit](KEPLER_Q2_BOUNDED_RINGCOHERENT005_BALLISTIC_T65_2026-10-27.json)
advects all 55,351 sites on softened circular orbits to `t=6.5` and rebuilds
the 49,038 disk Voronoi cells. It predicts the following mass-weighted mean
absolute `A(0)/A(6.5)-1` density changes:

| IC | Exact-orbit disk density change |
|---|---:|
| aligned | 1.85727e-4 |
| 5% ring-coherent phase offset | 2.05525e-4 |
| independent 5% site-angle offset | 1.47879e-3 |

The ring-coherent value is only 1.11 times aligned and is 7.20 times smaller
than independent site jitter. The comparison uses the same softened orbit
and Voronoi-area method for all three controls; the ring result's periodic
image safety ratio is 41.07. It supports the hypothesis that breaking the
regular spacing *within* a ring, rather than merely shifting adjacent rings
relative to each other, excites most of the current geometric density noise.
It does not prove the native pressure-force response or long-time stability.
The queued native ring-coherent preflight and `t=6.5` continuation provide
that separate test.

The independent [t=10 exact-orbit check](KEPLER_Q2_BOUNDED_RINGCOHERENT005_BALLISTIC_T10_2026-10-27.json)
gives `2.05637e-4` for ring-coherent, against `1.94290e-4` aligned and
`1.47192e-3` independent site jitter. The geometric contrast therefore
persists through the old `t=6–10` interval, although it is still only an
orbit/area prediction.
