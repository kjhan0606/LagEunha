# Net-mass ALE first-RHS representation audit

See ../../Docs/UNIFIED_HYDRO_ALE_REPRESENTATION_2026-11-05.md. The original
native small geometry and positive degree-three-exact face quadrature supply
ten analytic-field controls in 1D/2D/3D. Small reference mass-balance solves
never subtract the RHS mean. Every physical face mass/momentum/energy flux
is retained. Periodic cuts have explicit external boundary ledgers.

Transport, centroid, analytic Euler, Galilean, native geometry finite-difference,
and variance-rate identities pass. Net cell mass zero does not eliminate
surface mass circulation. The first-RHS proxy Edot-ubar dot pdot contains
Kvar_dot even if the initial discrete energy omits variance. The measured term
has mixed signs and local h^2 scaling. This is not a positive dissipation
measurement, a whole-PDE convergence test, or a no-go against ordinary FV.
No physical cold-flow accuracy or production scalar-energy closure is claimed.

The old Hopkins initialization experiment is not repeated. No evolved moment,
reservoir, energy reset, large global solve, native128 spectrum, cost campaign,
or PDE is introduced. The compile-time geometry hook defaults off; old/new
unhooked diagnostic binaries and output are bit-identical. Native geometry
and all source, commands, FD directions, arrays, and checks are archived.

Frozen legacy voro.h headers are stored as lossless .gz files to retain their
original bytes, including whitespace. Decompress to the original name before
rebuilding a frozen source tree. Uncompressed SHA256 hashes are recorded.
