# Same-IC fresh-weighted Kepler comparison at t≈7.425

The saved 5%-jitter, bounded native W0 and fresh-weighted trajectories have
the same initial hash and times `7.4252993734264185` and
`7.425299506861853`, respectively: a difference of `1.33435e-7`. CPU job
**437593** independently verified the actual binary hashes, exact saved
snapshots, origin, entropy and each branch's native physics gates; all 19
W0 checks and all 21 weighted checks pass. The checks are branch-specific:
the W0 density safety criterion differs from the weighted `1e-3` equilibrium
criterion, so 'pass' does not mean identical density thresholds.

The [comparison receipt](comparison.json) records disk mass-weighted mean
absolute density changes `1.83783e-3` W0 and `2.62810e-5` weighted, and mean
absolute log-entropy changes `1.64065e-5` and `7.55488e-6`. Total-energy
relative errors are `5.99393e-9` and `5.99229e-9`. The signed entropy state
changes are not automatically irreversible heating.

The [density-field figure](density.png), CPU job **437594**, shows the actual
Voronoi/power winner per pixel using common color scales. It states the small
time difference explicitly. Lower panels show each cell's density relative
to that same cell's initial density. The two trajectories use different
binaries and W0 versus evolving weights; both retain fixed individual mass.

This is progress beyond `t=6` for the bounded test, not a demonstrated local
moving-face mass/energy flux, full 384² periodic solution, evolved cell spin,
or `t=600` GIZMO reproduction. The independent
[fresh t=6.5 local-work and weight-swept-mass analysis](../KEPLER_FRESH_WEIGHT_LOCAL_T6501_2026-10-27/README.md)
sets out a concrete reason to test local transport before asserting physical
accuracy.
