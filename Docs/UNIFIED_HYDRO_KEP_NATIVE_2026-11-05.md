# Kinetic-compatible native coupling: the remaining cold heat is stress

The selected cell-mean KEP/EC core is now coupled to the same actual
native 1D/2D/3D power geometry, physical generator reconstruction and
birth-mu mass--weight law. It removes the previous cold birth-RHS heat
source. Eleven actual evolution endpoints complete, including cold64 to
t=0.25. However the registered pressure criterion still fails, so this
core plus the selected compression stress is not yet adopted.

## Fixed implementation and evidence

`GFS_NATIVE_KEP` explicitly selects the new core in
`Exam/gfs_fv_native.c`; the existing default HLLC/HLL path is retained.
The new experimental checkpoint version is 2 and rejects version-1
continuations, preventing silent use of late HLLC states under a changed
operator. Every new trajectory starts from its analytic birth condition.
The initial cold64 and Gresho64 tables are byte-identical to their HLLC
counterparts. No energy conversion, floor or reset is applied.

The mass flux is linear in mesh normal speed. Each RHS uses one scalar
SPD solve with the existing immutable mu, PCG cap 64 and true residual
criterion. Newton is absent. The stress and its explicit timestep bound
are exactly those in `UNIFIED_HYDRO_KEP_ALE_PLAN_2026-11-05.md`; coefficients
are not reduced to fit the cold threshold.

`kep_native_build_v2` is the tested build, binary SHA256
`0495742cdf9b0846a87e4cac7afa355f71e39ad98aff30b0782db33d7f87541f`.
Its source hashes match this implementation. v1 preserves the first
birth-RHS check; v2 adds actual-geometry entropy and boundary ledgers
before any PDE run. Executables remain in GPFS; compact evidence is in
`Research/KEPNative_20261105`.

## Mathematical checks are also evaluated on the actual mesh

The optimized C and C++ pure-flux audits pass 240 dimensional fixtures.
On actual native cold64/128 birth geometry, thermal rate, relative
pressure rate, kinetic pair defect, compression heat and the entropy
balance defect are all exactly zero in the recorded arithmetic.

Every actual RHS additionally computes the extensive mathematical entropy
rate with the changing cell volume included:

    eta_dot_i = v_entropy_i dot Qdot_i - rho_i Vdot_i.

The numeric boundary entropy flux and the signed stress dissipation are
recorded separately. The sum closes to the declared roundoff tolerance;
the volume term is not omitted. On accepted PDE states, maximum kinetic
pair defects are about 1.8e-15 for Sod, 1.1e-15 for Gresho, and below
4e-19 for cold shear. Semidiscrete entropy defects remain below 3.51e-13
in the cold tests and below 3e-16 in Sod/Gresho. These do not establish
exact finite-step entropy conservation or positivity for arbitrary data.

The driver separately logs

    deltaK_RK = K_final-K_initial-dt*(Kdot_stage0+Kdot_stage1)/2,

including the accumulated absolute defect. This is distinct from GCL.
`kinetic_ledger.json` identifies the core, checkpoint version and segment
scope, so a resumed segment is not mistaken for a full-history ledger.

## Actual endpoints and reference norms

`kep_native_group_v1` completes all eleven endpoints: stationary and
periodic translating contacts in 1D200, 2D8^2 and 3D4^3 to t=0.25;
1D Sod200/400 to t=0.2; cold64 with the registered step and quarter-step
to t=0.25; and Gresho64 to t=1. Six contact noise checks pass 1e-9.
Maximum conservation ledger defect is 1.40e-15 and the residual-corrected
mass--weight invariant defect 1.39e-12. No floor or repair is used.

The exact Riemann solution is integrated on actual Sod intervals:

| N | Density L1 | Mass-mean velocity L1 | Physical mean pressure L1 |
|---|---:|---:|---:|
| 200 | 0.00619038 | 0.0139286 | 0.00724473 |
| 400 | 0.00336976 | 0.00735139 | 0.00393889 |

All three errors decrease. They are larger than the retired reconstructed
HLLC results at the same resolutions. This cell-mean EC core has not
retained the full high-order face reconstruction of that warm-flow method.

Actual-polygon analytic means and physical node velocities give:

| Run | Physical node velocity L1 | Cell-mean velocity L1 | Physical pressure relative L1 |
|---|---:|---:|---:|
| Cold64, t=0.25 | 0.000249533 | 0.0000268182 | 0.0158379 |
| Cold64 quarter-step, t=0.25 | 0.000249554 | 0.0000267983 | 0.0185049 |
| Gresho64, t=1 | 0.0369119 | 0.0368108 | 0.00636790 |

Both cold pressure norms exceed 0.01. The quarter-step control does not
pass, and no smaller cold-only CFL is selected. Positive cold pressure
is maintained; the failure is physical accuracy rather than endpoint
collapse. Gresho improves over the older W0 velocity error near 0.073,
but worsens from the HLLC value 0.0126674. A 128^2 KEP vortex refinement
or matched GIZMO timing is not claimed.

## The residual cold heating has a measured source

For cold64 the mean-state thermal increase is 2.37568862e-8. Integrated
compression-stress heating is 2.37568997e-8, while accumulated absolute
RK kinetic defect is only 1.09e-14. The quarter-step values are
2.77573129e-8, 2.77573151e-8 and 3.33e-14, respectively. Thus almost all
remaining global cold heating is the registered stress, rather than a
failure of the advective kinetic identity.

The stress is exactly off at Cartesian birth, but responds as the smooth
sine-shear mesh deforms. The current law vanishes for affine reconstructed
velocity, which does not imply it vanishes for nonlinear incompressible
shear. The integrated stress changes with step sampling; a small kinetic
time defect alone does not rule out every temporal contribution.

Gresho has thermal increase 0.0138540 and stress production 0.0143537;
pressure work accounts for the remaining balance. Its RK kinetic defect
is 2.62e-8. Stress-induced dissipation and spatial accuracy remain serious
limitations, and the failed cold gate is not relaxed because it is close.

## Cost, regression and remaining requirements

Cold64 takes 2 accepted steps and 0.0251 evolution seconds; its quarter-step
control takes 8 steps and 0.107 seconds. Gresho64 takes 4182 accepted steps,
2676 rejected attempts and 77.39 seconds, compared with 54.26 seconds for
the retired HLLC warm run. The one-SPD solve reduces weight-solve cost, but
does not make the whole calculation faster. Cold timings are not evidence
of accuracy qualification, and none of these is a matched GIZMO benchmark.

A separately compiled default HLLC path still reproduces both contact200
and Sod200 final tables AND checkpoints byte-for-byte. Original production
OpenMP/MPI/CUDA paths, accepted Kepler histories and manuscript results are
not replaced by this prototype.

Original blast/collision and remaining 1D cases, explicit affine/alternating
native mode controls, 3D shocks, higher-resolution vortex/cold refinement,
the full 2D suite, full native GPU/MPI evolution and matched <=2x GIZMO
runtime remain unfinished. The next grouped evaluation must address the
compression-stress response to resolved shear without coefficient fitting,
while preserving the now verified mass/momentum/energy kinetic closure.
