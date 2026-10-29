# Initial ring geometry and late Q2 angular-momentum modes

The aligned bounded Q2 initial condition consists of azimuthally uniform
rings. The stored initial cell areas and angular momentum of the same IDs at
t=30 and t=47.50031653561342 were analyzed along each of the 11 original
rings in `0.5<=r0<0.75`. This is an angular coordinate around each ring,
not a periodic spatial boundary. Input hashes and each ring's values are in
the [receipt](KEPLER_Q2_INITIAL_GEOMETRY_MODES_2026-10-29.json); the
[reader](audit_kepler_q2_initial_geometry_modes_20261029.py) uses the existing
native snapshots without running a new simulation.

All 11 rings have initial stored-area spectral peak `m=13`. Adjacent selected
disk rings generally differ by 12 or 13 sites; for the innermost selected ring
the inward background is not represented by the disk-ring neighbor count.
At both late times, 8 of the 11 late angular-momentum peak modes are multiples
of 13: `m=65,78,91` (the others are `m=94,94,115`). For those eight rings,
the phase of late angular momentum relative to the *initial area* at the
corresponding late peak is close to `-pi/2`. At t47.5 their initial-area
power ranks at these peaks are between 5 and 13; their individual power
fractions are small, about 0.00017–0.00135. The strongest initial area mode
is not the strongest late angular-momentum mode on any of the 11 rings.

These exploratory, after-the-fact comparisons suggest a possible connection
between aligned-ring geometry and which perturbations are seeded or amplified.
They do **not** demonstrate that the initial area pattern causes the late
instability: small initial area harmonics, phase correlations, nonlinear
coupling, shearing dynamics and growing numerical modes need controlled
comparison. The 11 rings are not independent statistical samples and no
significance claim is made. Existing jitter-versus-aligned t10 data show
that removing coherence alone does not establish accuracy. A meaningful next
test should vary ring phase or population independently while preserving the
radial profile, fixed-mass/no-steering conditions and the same pressure,
then inspect both initial spectral seeds and equal-time native evolution.
It must use the ongoing pressure comparison and preserved checkpoints rather
than silently reinterpreting this post-processing as a new physical run.
