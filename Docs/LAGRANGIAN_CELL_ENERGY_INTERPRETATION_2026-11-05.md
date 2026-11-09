# Lagrangian cell interpretation and coherent kinetic energy

The user's current direction is to reason from moving, deforming fluid cells (a balloon analogy), including flattening, pressure evolution and pressure transmission to neighbors. The next derivation must use material-volume kinematics and mechanical work. The fixed spatial-field projection below is a diagnostic of representation, not a replacement physical model or proof that Laguerre cells follow material boundaries.

For a sealed material cell, mass is fixed. Its smooth adiabatic pressure obeys dU=-P dV and P V^gamma=constant. Shape deformation at constant volume is distinct from compression; inviscid smooth shear and rotation must not generate irreversible heat. Smooth compression can increase temperature without increasing entropy. A shock must produce the appropriate entropy and satisfy Rankine–Hugoniot relations. Pressure transmission is finite-speed acoustic dynamics, not instantaneous pressure equalization.

A shared face must transmit equal-and-opposite momentum and consistent mechanical energy. Generator motion with fluid velocity alone does not establish that the face is a material surface: normal face speed, fluid normal speed, volume rate and mass flux must be audited together. Retain the user's no-steering requirement and the preference for fixed mass. If the previously authorized weight-linked mass exchange is used, explicitly account for transported momentum and energy; do not silently treat it as a sealed balloon or convert coherent transport into heat.

The exact continuum-cell energy decomposition is

E_continuum = |p_cell|^2/(2 m_cell) + K_within_cell + U_thermal,

K_within_cell = (1/2) integral_cell rho |u-u_mass_mean|^2 dV.

K_within_cell can describe coherent shear/rotation; it is not inherently turbulence or thermal energy. The current native cold initialization explicitly uses U_thermal plus mean-motion kinetic energy and separately records the continuum projection deficit. This is an approximate discrete initial representation, not proof that a new evolved physical field is required.

Read-only native-geometry results for the same analytic steady sine shear:

| Quantity | Initial geometry | Saved final geometry |
|---|---:|---:|
| Analytic total kinetic energy | 0.25 | 0.25 |
| Exact cell-mean kinetic energy | 0.24979926688481 | 0.24978966794708 |
| Within-cell variance energy | 0.00020073311519221 | 0.00021033205291836 |
| Variance / prescribed thermal energy | 133.8221 | 140.2214 |
| Saved simulation mean kinetic energy | 0.24979926688481 | 0.24979926397588 |
| Saved velocity L1 vs exact cell means | 1.725725e-15 | 2.683334e-05 |

The analytic total is unchanged, but its mean/within-cell partition changes with native cell geometry. The evolved candidate nearly preserves mean kinetic energy and has small but nonzero velocity error relative to exact cell means on the new geometry. Neither thermal accuracy alone nor preservation of the mean-motion kinetic invariant is a complete accuracy certificate. The partition identity does not establish that all previous artificial heating was caused by omitted internal kinetic variance.

Using an affine velocity field, the within-cell contribution is approximately m trace(G C G^T)/2, with the velocity gradient G and geometric covariance C of the cell. Native face quadrature gives C via the divergence theorem, without a new tessellation. This affine estimate matches the exact analytic variance with relative cellwise L1 about0.37% initially and0.44% on the saved final geometry. The measured local moment/gradient pass takes milliseconds here; this is not an end-to-end solver cost guarantee or a proof of a valid closure.

Next work: derive the moving-cell mass, volume, momentum, angular-momentum and energy budgets first. Track coherent deformation/rotation, reversible pressure work, and irreversible shock heating separately in the diagnostic. Determine whether compatible reconstruction and existing geometric moments suffice before introducing any independent spin/elastic/subcell-energy variable. A candidate must also retain the successful pressure-mode response, improve density/wall peaks, pass the same 1D/2D/3D shock tests and satisfy the GIZMO speed requirement. Do not add an assumed rubber membrane shear modulus from the analogy.

References for the interpretation, not validation of this candidate: [Aluie, coarse-grained compressible flow](https://arxiv.org/abs/1012.5877) and [Zier & Springel, cold shear on a moving mesh](https://doi.org/10.1093/mnras/stac1783). The former distinguishes coarse-grained momentum/kinetic transport; the latter demonstrates that accurate face integration matters for moving-mesh shear noise.

All audit inputs are read-only. No pressure reset, new state variable, or production-code change was made for this physical interpretation.

Current computation status: physical11's original nine 1D endpoints, cold64 t=.25, and existing dimensional contact controls have completed. Cold pressure L1 is0.193928%, within the unchanged1% gate. Collision density peak remains8.56% above the exact high state. Independent strict-no-floor Athena PLM Blast references at N1600/3200/6400 reached.038; density mean self-difference from3200 to6400 is1.1267%, so this is a numerical reference with measured resolution sensitivity. The new dimensional fixture harness has been built; its RH/wall test matrix has not yet been launched. Those tests and full speed qualification remain outstanding.

Audit artifact directory: /gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/cold_shear_energy_decomposition_20261105_v2


# Material-face kinematics on existing physical11 cold-shear snapshots

This is a read-only audit of the actual common native RHS, not a new model or evolution. Original checkpoint hashes are recorded before/after. The existing physical11 method and all prior results remain intact.

## Question

Does nearly fixed cell mass imply an impermeable material face, and does face slip itself explain artificial heating? Inspect the existing 64² periodic sine-shear test, rho=1 and P=1e-6, using saved initial and t=0.25 checkpoints. Evaluate geometric face speed including solved weight rates, exact sine-shear normal velocity at each quadrature point, reconstructed normal velocity, actual EC+acoustic mass flux, and the thermal-work ledger.

## Evidence

| Quantity | Initial | t=0.25 |
|---|---:|---:|
| Area-mean absolute exact material-normal slip | 0.013495 | 0.008372 |
| Absolute unique-face mass transfer rate | 0.859281 | 0.516809 |
| Sum absolute cell net mass rate | 1.26e-15 | 4.76e-7 |
| Sum absolute cell volume rate | 0 | 2.56e-6 |
| Total thermal energy rate | 6.77e-19 | 1.57243410e-8 |
| Sealed-cell sum -P_i dV_i/dt | 0 | -3.56e-15 |
| Acoustic thermal rate | 0 | 1.57243445e-8 |

Mass transfer counts each face once and integrates the absolute quadrature flux; it is not a net domain mass change. Cell net mass rates sum to zero to roundoff. Reassembled mass flux agrees with the native edge flux to 1.30e-18 or better. At the final checkpoint the exact analytic sine field is a reference, not an assertion that the numerical field is exactly analytic; reconstructed-trace slip is 0.008375, close to the analytic 0.008372.

## Interpretation and limitations

1. Fixed or nearly fixed cell mass does not imply that every face follows the same fluid parcels. Transfers can cancel around a cell.
2. Significant transfer already exists at the initial state with negligible heating. Slip or exchange alone is therefore not a sufficient explanation of the observed spurious heating.
3. At the final state, pressure-volume work is negligible compared with the acoustic thermal rate. The outstanding heat is associated with the numerical acoustic transfer, not an imposed physical viscosity (Re remains infinite).
4. The closed material-cell law dU=-P dV is a useful baseline but cannot be imposed on the present open discrete cells without including transported energy and momentum. This audit does not invalidate ALE or establish that a new internal velocity/rotation state is necessary.
5. The user's deforming-fluid-balloon viewpoint should explain the existing method and identify the simplest remaining force/work inconsistency. It does not authorize discarding previous successes or adding an elastic membrane law.
6. Instantaneous source diagnostics at two snapshots are not a causal trajectory decomposition, nor a benchmark against GIZMO. The prior integrated acoustic ledger is independent supporting evidence.

## Next discriminating step

For the same saved cells, compare the compatible kinetic-energy work with the reconstructed face motion and the analytic coherent-shear transport. Separate reversible pressure work, transport across moving faces, and irreversible acoustic dissipation. Any alternative must retain pressure-checkerboard damping, constant-pressure density contacts, original shock endpoints, and full reflected-wall validation in 1D/2D/3D; no new sensor or fitted coefficient is justified by this audit alone.


## User clarification: material faces, not merely fixed cell totals

The user explicitly requested fixed cell masses and then clarified that matter must not cross individual faces. This supersedes the current experimental use of weight-linked mass transport as the selected development direction. Preserve previous variable-mass successes and failures as comparison evidence; do not discard earlier work or reinterpret accepted production trajectories.

The newly tested physical12 branch at `/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/acoustic_fixed_cell_mass_stage_v1` constrains only the net flux of each cell and therefore DOES NOT satisfy the clarified no-face-transfer requirement. Sod200 and Blast200 reached original targets with bit-identical birth masses, and stationary contact in 1D/2D/3D completed. Cold64 failed the entropy ledger at an unaccepted first Euler trial (defect 1.56e-9); it is not a successful shear result. Do not continue presenting this net-only branch as the requested material-cell model.

For an outward face normal n_if and shared face pressure P_f, the material-cell baseline is

    dm_i/dt = 0,
    d(momentum_i)/dt = -sum_f A_f P_f n_if,
    d(total_energy_i)/dt = -sum_f A_f P_f u_f dot n_if,
    dV_i/dt = sum_f A_f w_geometry,f dot n_if.

A material face requires its physical fluid normal velocity u_f dot n_if to equal its actual geometric normal velocity. A shared pressure and face velocity make internal pair momentum and total-energy exchanges antisymmetric. The implied thermal update follows from total energy minus cell mean kinetic energy; it must not be independently patched. Smooth closed adiabatic motion and shock entropy generation must be checked separately. Numerical shock dissipation is permitted when derived consistently with this material Riemann/pressure-work problem; no physical shear viscosity at Re=infinity is implied.

Use existing acoustic information as a possible pressure/face-velocity closure, not an automatic reason for mass advection. Merely zeroing mass flux while retaining unrelated advective momentum/energy terms is not the requested derivation. Equally, assigning independently chosen Riemann velocities to all faces does not establish that one Voronoi/Laguerre tessellation with shared generator and weight rates can realize them. Audit this compatibility explicitly before selecting a closed model; do not silently introduce steering, elastic membrane physics, or heat reservoirs.

Immediate design task: distinguish (a) material pressure and work, (b) actual shared-cell geometry, and (c) compatibility between their face motions. Preserve the existing pressure-mode/contact/shock diagnostics as acceptance tests for the resulting common dimensional operator.
