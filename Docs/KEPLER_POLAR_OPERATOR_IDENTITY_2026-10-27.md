# Active polar Kepler pressure operator: verified 2026-10-27

The q1/q2 H200 binary SHA-256
`7584426920df56413fa05bdc950b82f890412cfaa793269f4edcc66cd67c0284`
prints `[VARPRESS_W0] experimental smooth fixed-mass variational pressure
active` in its actual run log. The isolated source and binary contain
`GFS_VARPRESS_W0`; the run enables phase 1, AV mode 5, W=0, fixed mass,
and geometric-CFL coefficient 0.1. It uses the user's CUDA+MPI+OpenMP RK4
solver with one MPI rank and two OpenMP threads in these jobs.

## What the force and work actually use

In the active `Exam/exam_gpu.cu` production force kernel, the geometric
Voronoi face velocity is computed from both generator velocities and the
face-centroid correction. For a directed face with area vector `A`, length
`L`, generator separation `d`, and centroid offset
`o=c-(x_i+x_j)/2`, the variational branch forms

`b_i = A/2 + L o/d`, `F_i = (P_i-P_j) b_i`,

`Udot_i(face) = -P_i [(w_g-v_i)·A]`.

It also forms the total-energy face contribution as `v_i·F_i + Udot_i`.
The sum over the complete closed cell faces recovers geometric `-P_i dV_i/dt`
because `Σ_face A=0`. This is the intended adjoint geometric pressure
operator. The actual log and the source-level compile branch establish the
operator identity; they do not prove long-time accuracy.

The AV5 HLLC calculation still computes a star pressure and contact speed,
and the local short-pair code still evaluates a `p_pair` and increments the
`n_pair` counter. Under `GFS_VARPRESS_W0`, the force and phase-1 total-energy
work use `F_i` and `Udot_i` above instead of that star/pair pressure.
Therefore **`n_pair>0` shows a close-pair criterion was evaluated, not that
its restoring pressure acted on the dynamics**. Likewise this active run
must be identified as variational W0 pressure, not an HLLC pressure-flux
Kepler test. HLLC is present in the executable but does not supply the
pressure force in this smooth branch. The branch's shock behavior requires
its own validation; no GIZMO HLLC equivalence follows from the Kepler run.

## Consequence for the Kepler program

The [material/face-motion audit](KEPLER_FACE_MOTION_GAP_INTERPRETATION_2026-10-27.md)
measures geometric non-materiality and fixed-mass density variation. It is
not a missing pressure-work term in this variational branch. At `t>6`,
the tests must separate geometric density drift, entropy and positivity,
variational pressure torque, and outer periodic-seam forcing. If short
pairs become a failure mechanism, an actual conservative pair-force/work
extension needs an explicit design and native comparison; changing the
counter or relying on the currently bypassed `p_pair` is insufficient.

The old documentation called these runs `RK4/HLLC` and treated `n_pair` as
an active restoring force. Those labels are superseded by this source and
runtime audit. Existing binary, jobs, and checkpoints were not changed.

## Independent interior adjoint check

The [small-mesh reproducer](audit_kepler_varpress_adjoint_20261027.py) and
[machine-readable result](KEPLER_VARPRESS_ADJOINT_SMALL_MESH_2026-10-27.json)
check the active W=0 force and face-work formula against independently rebuilt
SciPy Voronoi polygons. Three perturbed 7-by-7 meshes, two pressure fields per
mesh (uniform and nonuniform), and five interior sites per field give 30 site
checks. A central finite difference of the fixed-entropy polygon internal
energy supplies the reference force; a separate finite difference of polygon
area under generator drift supplies the reference volume rate. All five
predeclared checks pass. Across the 30 sites, the largest absolute errors are
`4.55e-9` for force versus negative energy gradient, `1.24e-9` for the area
rate, and `1.77e-9` for thermal work versus `-P dV/dt`; closed-face area-vector
error is at most `7.83e-16`. The JSON records both source and reproducer hashes.

This establishes the smooth **interior geometric identity** for the tested
source formula. It does not check native GPU face extraction/pairing, periodic
seams, free/ghost boundaries, shocks, or stability at late Kepler times. Those
remain separate acceptance gates; in particular the finite-difference result
cannot promote a short polar run into a `t>6` or GIZMO-length validation.

A separate [periodic 5-by-5 reproducer](audit_kepler_varpress_periodic_adjoint_20261027.py)
and [result](KEPLER_VARPRESS_PERIODIC_ADJOINT_SMALL_MESH_2026-10-27.json)
extend the mathematical check across periodic image faces. For three perturbed
meshes with uniform and nonuniform pressure (150 central-cell evaluations),
each generator perturbation moves all nine periodic images before the reference
energy is rebuilt. The maximum force-versus-energy-gradient error is
`4.43e-9`; face area-rate and thermal-work errors are `1.69e-9` and
`2.47e-9`. The box-area partition error is `1.78e-14`, total force residual
`1.23e-15`, and total kinetic-plus-thermal power residual `8.19e-16`.
All nine declared checks pass, including zero force at uniform pressure.
This supports the periodic **continuum geometry formula**, but still does not
exercise native CUDA image pairing or the discontinuous Kepler velocity field
at the outer seam. A live GPU checkpoint and boundary-control comparison are
still required before assigning late heating to the pressure formula.

## What the polar IC tests especially well

The Hopkins-profile IC sets nearly uniform `P=1e-6`. In the active
variational operator `F_i=(P_i-P_j)b_i`, its initial pressure force is
therefore zero up to storage and reconstruction roundoff; the softened
point gravity and assigned circular velocity determine the early motion.
The q2 regular rings also keep exact-orbit Voronoi areas unusually steady.
These two facts make q2 a legitimate but favorable equilibrium-preservation
control for this particular operator. This is an inference from the force
formula and the independent geometry audit, not evidence that generic
shearing ICs or the long-time native run are stable. A small site-jitter
control with the same physics and boundary is needed to test how strongly
the result depends on ring regularity.

The [native q1 `t=.825` strict audit](/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q1/halo12_from061876_to0825/audit_long_t0825_strict.json)
measures initial-disk orbital angular-momentum drift `-9.31e-13` and signed
mass-weighted entropy-log change `-8.45e-6`; the absolute entropy-log error
is `1.53e-4`. At this early checkpoint there is no evidence of a substantial
disk torque or net irreversible heating for cell spin to cure. This does
not exclude a late shear defect. The large apparent outer-atmosphere
angular-momentum change in a periodic-box, center-referenced diagnostic
must not be treated as a conserved global angular momentum: positions
jump under periodic wrapping, and the outer velocity field itself fails
the torus matching condition.
