# Coupled massless-weight acoustic screen: conservative but shear-heating

This candidate explicitly solves neighboring weight rates together; it does
NOT prescribe equal velocities for all faces. Existing native1024-cell
Gresho geometries are reused, with the kinematic operators D,N,K,B=DN from
the previous diagnostic. Ti,Tj reconstruct the two side normal velocities.
The impedance is the ideal-gas linear acoustic value Z=rho*c, not a fitted
coefficient. No auxiliary weight momentum or kinetic energy is introduced.

Let w=Nv+Kq, ri=w-Ti v, rj=w-Tj v and
R=.5 sum_edges A[Zi ri^2+Zj rj^2]. At fixed state this is quadratic in v,q.
Massless weight force balance gives

    K^T diag(A(Zi+Zj)) K q
      = K^T {A Delta P - A Zi(N-Ti)v - A Zj(N-Tj)v}.

A mean-zero gauge fixes q. The paired physical force and internal-energy rate:

    F=B^T P -(N-Ti)^T A Zi ri -(N-Tj)^T A Zj rj,
    Udot_i=-P_i Vdot_i + sum_incident_faces A Zi ri^2,
    Vdot=Dw.

Since the q equation implies L^T P=partial_q R, Euler's identity for the
quadratic R gives v.F+sum Udot=0. Cell entropy production numerator is
nonnegative. These are semidiscrete identities, not a nonlinear timestep or
shock accuracy guarantee. All side velocities differ as required by their
neighbor relations. This is linear-acoustic, not the exact nonlinear piston
solver already used for the heated-pair reference.

Numerical results on initial, deformed-frozen and deformed-dynamic geometries:
- Uniform P=1 rest: q=0, heat=0, force <4.6e-18.
- Saved Gresho fields: energy-rate residual <4e-16; positive dissipation.
- Cold smooth periodic shear vx=sin(2pi y), vy=0, P=1e-6:
  regular initial geometry heat/U ~2.31e-29 per unit time;
  deformed-frozen heat/U=.46829 per unit time;
  deformed-dynamic heat/U=.51581 per unit time.
These are instantaneous rates, not measured finite-time pressure increases.
The shear case retains each saved cell's mass/volume density; it is a
manufactured operator input, not a new exact nonlinear solution on the evolved
Gresho density. Ideal inviscid flow nevertheless has no physical viscous
heating from this shear. Reconstruction/geometry incompatibility is being
converted into positive numerical heat by this closure.

Decision: do not implement this unchanged as the common PDE solver or add
another sensor to rescue it. Energy conservation plus positive entropy is
insufficient. The geometrically unresolved part of smooth material motion
must not automatically be treated as acoustic dissipation. This closes this
specific simple impedance closure, not all coupled weight laws.

Next physically motivated route: examine whether the needed internal motion
can be represented by the existing cell velocity moments and their actual
kinetic metric. Arbitrary diagonal weight inertia lacks that mapping, while
instantaneous relaxation thermalizes the missing motion. A physical metric
must be derived from represented fluid velocity, and screened for cost before
any production port. Do not automatically add new state variables or invoke
an expensive per-step quadrature; existing moment/rotation work is preserved.

Artifact: Research/MaterialVolume_20261105/coupled_weight_acoustic_20261106_v1.
audit.py and results.json, nine frozen-state cases, all complete. No native
PDE changes, no new simulation, no3D or GIZMO speed claim.
