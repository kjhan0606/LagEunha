# Reversible-pressure contribution to the coupled response at t2.8

Used saved +/- states at checkpoint503, reconstructed independent finite
Voronoi faces, and evaluated the linear pressure operator L(x)P. Interior
volumes reproduce solver values within2.85e-14; masses match initial
values and all weights are zero. No simulation was launched or modified.

For central response delta_F=(L+P+ - L-P-)/(2 epsilon), split exactly:

    delta_F_pressure = (L+ + L-)(P+ - P-)/(4 epsilon)
    delta_F_geometry = (L+ - L-)(P+ + P-)/(4 epsilon)

Separately use P_iso=P_seed*(V_seed/V)^gamma, with seed at t2.0029,
to divide each branch pressure force into seed-entropy isentropic and
subsequent entropy-excess components. Both force splits close to relative
L2 below9.5e-12 on the initial interior region.

For initial r<.5, instantaneous sum delta_u dot delta_F is:

| Component | Departure-response power |
|---|---:|
| Total reversible pressure | +1.7070519e-8 |
| Pressure-value variation | +1.6622748e-8 |
| Geometry-operator variation | +4.4777082e-10 |
| Seed-entropy isentropic force | +1.7300101e-8 |
| Entropy excess since seed | -2.2958283e-10 |

The first split contributes about97.38% through pressure values and2.62%
through explicit operator geometry. These are two decompositions of the
same pressure force; do not add rows across the two decompositions.
Pressure values themselves depend on geometry/volume, so small explicit
operator contribution does not exclude a geometric origin of pressure
variation. Seed entropy already includes everything accumulated before t2;
the small negative excess is not a statement about all heating since t0.

The initial .3-.5 band supplies1.57695e-8 of this pressure response power.
The directional quantity -sum delta_x dot delta_F is positive overall
(3.03645e-8 inner), despite a smaller negative explicit-geometry component.
Do not interpret positive work on velocity as a negative restoring stiffness;
nor is this finite directional quantity a full Hessian spectrum.

This measures pressure's contribution to the derivative of the central
circular-departure kinetic norm, not its complete derivative. TCF, HG,
gravity and the reference-velocity position derivative are omitted.
Transient compressible/shearing response remains possible; no instability
or numerical cure is established. Continue the running refined comparison
before using the quantitative response to choose a physics modification.

Evidence: canonical_coupled_pressure_response_2026-10-18.json,
canonical_coupled_pressure_geometry_2026-10-18.json.
Script: audit_coupled_pressure_response.py --step503.
