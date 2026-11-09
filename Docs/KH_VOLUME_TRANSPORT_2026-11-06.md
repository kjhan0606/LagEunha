# Volume-null disturbances are not transported consistently by native W=0 shear

## Scope and native checks

Construct exact unperturbed row-shear configurations from the existing
32x64/64x128 zero-seed IC: x(t)=x(0)+U(y)t, y(t)=y(0), modulo the periodic
box. These are diagnostic snapshots, not a new integrator. Evaluate actual
C/CUDA K1 and export faces at t=1 and t=1+-1e-4,+-5e-5. No W, force, mass,
or energy correction is applied to a production trajectory.

Reconstruct B=dV/dx from these native W=0 faces, ordered by particle ID.
The pressure force B^T P reproduces CUDA acceleration within1.11e-13 on
the exact shear configurations and1.63e-13 on the earlier actual disturbed
t=1 checkpoint. Uniform pressure force vanishes, and B times the base shear
velocity is below1.15e-16. All face neighbor matches are checked against
BOTH endpoints of the Voronoi bisector.

An initial diagnostic assertion failed on a1.19e-10-long edge: its endpoint
normal predicted the neighbor location with1.34e-10 absolute error. The
diagnostic now uses that normal only to locate the neighbor, verifies the
bisector endpoints, and uses the actual generator separation for the normal.
The maximum location prediction discrepancy in the disturbed geometry is
1.15e-9; force agreement above is the independent final check. This changed
postprocessing only, not native tessellation or the physical force.

## Material transport test

Let xi be a small generator displacement and A the base velocity gradient:
A maps each cell's y displacement to U'(y) times its x displacement.
A^2=0 for this shear. The Eulerian velocity disturbance is

    u = dot(xi) - A xi.

At uniform background P and fixed entropy per particle, linearized pressure is

    p = -gamma P V^-1 B xi.

For a material divergence operator transported by this shear, a relabeling
with u=0 and p=0 should not create pressure. The discrete diagnostic is

    R = dot(B) + B A.

Start with B xi=0 and dot(xi)=A xi. Then the native geometric pressure
derivative is -gamma P V^-1 R xi. A nonzero result exposes failure of this
discrete material-transport property, even though base shear volume and
constant-pressure force are correct.

Generate a high-k transverse displacement localized near one shear layer,
project it into ker(B) using sparse least squares, and normalize ||xi||2=1.
The projection residual ||B xi|| is1.12e-13 at32 and1.14e-13 at64.
Nevertheless ||R xi|| is0.62608 and0.097882 respectively. The induced
divergence RMS ||V^-1 R xi||/sqrt(N) is14.1667 and4.42965. These absolute
values use a unit discrete displacement norm; they are not observed physical
growth rates and should not be used as a spatial convergence estimate.

Halving the geometry differentiation interval changes R xi by6.42e-7 and
5.32e-6 relatively. Thus this defect is not explained by differentiation
roundoff or the failed tiny-edge matching assertion.

## Why total energy conservation does not settle it

The linearized native force about uniform pressure is M ddot(xi)=B^T p.
Define the Eulerian perturbation energy using native cell masses and volumes:

    E = (1/2) u^T M u + (1/2) p^T [V/(gamma P)] p.

Direct differentiation gives

    dot(E) = -u^T M A u - p^T R xi.

The first term is the continuum shear-energy exchange. The second is an
additional discrete geometry-transport contribution. This identity is checked
numerically with independent random displacement/velocity vectors in
analyze.py. It does not violate total native energy conservation: energy
can be exchanged with the background shear by a different discrete mechanism.

This is a concrete spatial mechanism to investigate, not proof that it
accounts for the full measured instability or that an arbitrary null-space
vector is a faithfully reconstructed smooth continuum disturbance. The
finite-resolution volume/velocity representation is precisely under test.

## Consequence for dynamic weights

Weight differences need not move all faces equally. Their volume operator L
can, in principle, compensate a zero-sum volume-transport defect. Any candidate
must transport the intended material volume/pressure response and include
its corresponding momentum and energy terms; simply adding pressure-driven
weight inertia does not establish this property. No weight solve or correction
has yet been implemented by this audit, and no extra damping term is adopted.
Next evaluate that kinematic role against admissibility, work consistency,
and solve cost, rather than changing a pressure coefficient to hide growth.

Artifacts: Research/MaterialVolume_20261105/kh_volume_transport_20261106_v1/.
Native input/output, results.json, metrics.json, analysis.log, and the saved
null-mode arrays are retained. Full2D, later3D, and <=2x GIZMO remain open.
