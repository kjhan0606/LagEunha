# Gresho physical error comparison at t=0.1

Artifact: Research/MaterialVolume_20261105/gresho_weight_accuracy_20261106_v1.
Existing native C/MPI/CUDA KDK now also accepts a frozen-feature checkpoint
for controlled W=pi=0 comparison; weight kicks/drift are disabled only when
the checkpoint explicitly carries the frozen feature. No existing Kepler
entry point was changed. One process, two OpenMP threads, A10.

## Metric definition

rho L1 = sum V_i |m_i/V_i - 1| / sum V_i.
v L1 = sum V_i |v_i-v_exact(x_i)| / sum V_i.
radial RMS = sqrt(sum V_i v_radial,i^2 / sum V_i).
The standard Gresho maximum velocity is one and density is one. Velocity
reference is point-sampled at current generator positions, consistent with
the IC velocity sampling, not an exact integral over each deformed cell.
All intervals below are t=0.1, not a full-period vortex test.

| Grid | W mode | density L1 | vector velocity L1 | vector velocity L2 | radial RMS |
|---|---|---:|---:|---:|---:|
|16x16|frozen|.00394572|.01700297|.02107934|.01414167|
|16x16|dynamic|.00226004|.01514424|.02578867|.02396571|
|32x32|frozen|.00147269|.00717092|.00979194|.00690223|
|32x32|dynamic|.00061976|.00501130|.01057047|.01016898|

New frozen16/frozen32/dynamic32 calculations use the same dt cap 1e-4.
The dynamic16 result is reused from its previous adaptive-dt checkpoint to
avoid needless replay; it is not an exactly matched-time-grid comparison.
The 32-grid pair is matched in timestep cap and restart endpoints.

## Interpretation

Both modes improve with spatial refinement on these metrics. At 32^2 dynamic
weights lower density L1 by about 58% and velocity L1 by about 30%, but raise
velocity L2 by about 8% and radial RMS by about 47%. Thus a mean density or
energy plot alone would conceal a worse rotational-flow error.
Weight kinetic energy / translational kinetic energy falls from 1.59% at16
to .543% at32; this supports further investigation but does not establish a
physical interpretation of auxiliary inertia or an asymptotic order.

Do not select dynamic weights as the universal method on these results.
Keep both paths as controlled candidates. Need localization of radial error,
longer vortex/instability tests, consistent 1D shock/entropy behavior, and
end-to-end GIZMO timing. 2D still incomplete; 3D remains deferred.
