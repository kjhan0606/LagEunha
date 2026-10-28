# Prepared coupled position/velocity seeds with fixed entropy

Prepared reference and +/-0.001 times actual t2 displacement and velocity
departure from the exact initial circular trajectories, restricted to current
r<.6. Removed mass-weighted translation from both seed fields. This is an
isolated diagnostic IC perturbation, not steering during evolution.

Independent Voronoi areas for r<2.4 match original areas within3.25e-14.
Rebuilt plus/minus geometry changes areas by at most3.75e-5 relative.
Let R=V_new/V_old from ratios of independently reconstructed areas. Retained
mass and entropy using V*=R, rho/=R, P*=R^-gamma, U*=R^(1-gamma),
cs*=R^((1-gamma)/2), leaving stress.K unchanged. Entropy relative error is
4.45e-16 or less; no volume change measured beyond r1. Reference snapshot
is byte-identical to the preserved t2 source.

Actual solver terminal geometry/force probes then reproduce positions and
velocities exactly, masses exactly, and interior(r<1.9) volumes/pressures
within2.40e-14 relative. No RK evolution has yet been performed.

The first analyzer failed because raw snapshot IDs include particle-type
bits while exported IDs do not. Corrected comparison masks IDs with2^56-1.
Reanalyzed existing outputs; no probes were repeated. Reference exit0 was
collected. Plus/minus return codes were not collected before the original
runner assertion terminated; both terminal probe completion markers and
complete numeric arrays were verified. Do not invent missing return codes.
Runner now collects every exit code before any numerical assertion.

Next: evolve the three prepared case directories using the existing orbit
work binary, with OMP4 and dt=.0025 over .1, and compare central responses
including displacement, circular-departure velocity and acoustic pressure.
This seed holds current entropy fixed; it is not a full unstable eigenvector.
No physical solver cure or long-term pass is claimed.

Remote root: /gpfs/kjhan/LagForce/Codex_coupled_onset_seed_20261018.
Evidence: canonical_coupled_onset_seed_2026-10-18.json,
canonical_coupled_onset_initialization_2026-10-18.json.
Scripts: prepare_coupled_onset_seed.py, check_coupled_onset_initialization.py.
