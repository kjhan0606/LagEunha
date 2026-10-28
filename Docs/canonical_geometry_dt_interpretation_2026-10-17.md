# Canonical geometry timestep investigation

The existing t=10.0020667693911 pre-HG CPU RHS dump was analyzed on lageunha,
without restarting or modifying any of the four live calculations.

Among emitted active TCF faces, the smallest reconstructed geometry timestep
is 0.0023311118 at IDs 808/688, radii 0.449713/0.455500. Separation is
0.006157793, face sound speed 0.3192802, and normal relative speed 0.01091592.
The sound term contributes 96.69% of the denominator. Separation is only
0.2184/0.1859 times the respective square-root cell volume. Thus this face's
restriction combines small separation and elevated face sound speed; normal
approach speed is a minor contribution at this instant.

Cell sound speeds reconstructed from pressure, volume and mass agree with
the printed face values within 4.78e-7 relative. No active emitted face has
a clamped factor, so clamp ambiguity does not affect this dataset.

Limitations: the dump includes emitted TCF faces, not necessarily every
neighbor examined by the timestep limiter. It samples one pre-HG RHS, not
all RK stages. Its minimum must not be equated with the actual accepted
timestep, or used alone to explain the smaller timestep at t>10.3.
The next diagnostic must identify the actual minimizing neighbor and stage.

Live process check confirmed PIDs 2286374, 2295268, 2352123 and 2432921.
Subsequent accepted log endpoints were baseline t=10.43485817,
HG-boundary-corrected t=10.41156364, cap .00125 t=5.52553815,
and cap .0025 t=4.40553815. All latest entries had finite=1 and dM0=0.
These are execution observations, not physical-accuracy passes.
