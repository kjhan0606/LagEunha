# Fixed-mass material-boundary 1D paper evidence

Authoritative suite: `/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/material_full_original_nine_20261105_v1`.

The manuscript reports the exact Riemann / primitive PLM-minmod / SSPRK2 native C comparison. The method evolves shared material boundaries with immutable masses. No RK4, KDK, variable-weight mass exchange, or spin evolution is used in this suite. Multidimensional Laguerre equivalence and production CUDA integration are not established.

`accuracy.json` and `assessment.json` preserve measured cell-mean errors, endpoints, boundary-work conservation balances, and original-data hashes. `profiles.npz` preserves physical-domain cell samples and reference profiles for all nine cases. `make_figures.py` regenerates the three publication PDF/PNG panels from the archived profiles without running a solver. The plotting script retains raw cell segments and center points, including narrow spikes.

Order in the profile archive: Sod, blast, Noh, Shu–Osher, Lax, double rarefaction, collision, stationary contact, moving contact. Data columns: left face, right face, mass, density, velocity, pressure, total energy. Reference columns: density, velocity, pressure.

Analytic profiles are represented on 40,000 intervals for reference integration. Blast and Shu use independent finite-resolution Athena++ references at 6400 and 12800 cells. Such numerical references are not exact solutions.

`resolution_analysis.json` and `acoustic_convergence.json` are separate reconstruction comparisons. The four-case resolution sequence predates the full original nine suite and has different contact and boundary specifications. The manuscript only quotes matched Sod/blast reconstruction comparisons and smooth-wave convergence, not an all-nine resolution qualification.

`prior_gfs_1d_section.tex.txt` preserves the historical 1D GFS prose and figure references replaced in the main results. Original figure files remain untouched. Historical GFS RK4, 2D KH/RT KDK, and Kepler RK4/KDK results are separate methods, not tests of the new material-cell update.
