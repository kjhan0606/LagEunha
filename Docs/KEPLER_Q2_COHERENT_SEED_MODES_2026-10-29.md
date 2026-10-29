# Existing Q2 ring-phase control at t=6.5

The preserved aligned and 5% coherent ring-phase-shifted Q2 W0 runs were
compared at the same t=6.5. The phase shift rotates all particles within a
ring together by a small ring-dependent offset; it leaves the ring particle
counts unchanged. The two ICs represent the same continuum profile but do
not have identical discrete masses. Hash-verified initial and endpoint states,
fixed individual mass, W0, and uniform ring spacing were checked. See the
[receipt](KEPLER_Q2_COHERENT_SEED_MODES_2026-10-29.json) and
[reader](audit_kepler_q2_coherent_seed_modes_20261029.py).

On each of the same 11 rings initially in `0.5<=r<0.75`, the largest
initial-cell-area Fourier mode is `m=13` in both ICs. The largest mode of
the t=6.5 specific-angular-momentum change is also `m=13` for **all 11**
rings in both ICs. The relative phase of the latter to the initial-area
mode ranges from about -1.59 to -1.42 radians; amplitudes differ by ring and
between ICs. At the later aligned-only t=30 and t=47.5 snapshots the dominant
L modes instead include `m=65,78,91,94,115`; see the
[late-mode receipt](KEPLER_Q2_INITIAL_GEOMETRY_MODES_2026-10-29.json).

The existing small ring-phase control thus does not remove the `m=13`
pattern and cannot test whether changing the ring-population pattern changes
the late high-mode selection. These are post-hoc descriptive comparisons,
not causal evidence or a measurement of the late mode's growth rate. The
phase comparisons alone do not distinguish direct forcing, geometric
symmetry, mode coupling or another correlated mechanism. The t=6.5 control
does not demonstrate that the coherent case remains stable to t=30.

Past independent site jitter worsened t=6.5 errors; it is not a proven
replacement IC. A discriminating future population/geometry comparison must
match radial profile, local resolution, pressure and runtime controls, and
verify geometry and conservation before attributing differences to ring
counts. It should not alter the ongoing Q1 or pressure sensitivity branches.
The angular Fourier transform is on ordered ring IDs; no periodic spatial
boundary was introduced.
