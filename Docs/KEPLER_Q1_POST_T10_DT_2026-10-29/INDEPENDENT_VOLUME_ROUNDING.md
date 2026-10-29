# Endpoint volume representation and entropy

For each of the three certified t=10.1 controls, an independent all-site clipped Voronoi reconstruction was evaluated at two guard radii. Both reconstructions agree within 2.22e-10 relative and cover the fixed polygon to 1.19e-12 absolute area. The native stored float volumes agree with independent areas to better than 5.96e-8 relative in every case. The t=10 reference uses the separately audited double-volume initialization at exactly the same generator positions.

Replacing only the initial and final stored float volumes with independent double geometry changes the disk mass-mean absolute per-cell incremental log-entropy diagnostic by 1.34–1.35e-8. The disk absolute-mean result with native versus reconstructed volumes is 3.01252e-6 versus 3.01248e-6 for coarse, 6.02580e-5 versus 6.02582e-5 for half, and 1.076029e-4 versus 1.076030e-4 for quarter. The endpoint volume readout rounding cannot by itself explain the increase across the timestep controls.

This comparison leaves the internal energy, actual trajectory and all intermediate operations untouched. It does not rule out accumulated volume or energy rounding during the integration, physical pressure work, or other numerical effects. Full results and input hashes are in `INDEPENDENT_VOLUME_ROUNDING.json`; the separate H200 thermal64 run tests a different arithmetic path.
