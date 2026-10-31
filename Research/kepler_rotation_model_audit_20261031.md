# Kepler rotation model audit — 2026-10-31

## Scope

This is a source audit and a conservative-model checklist. It is not a claim that the current Kepler run contains cell spin, and it does not change the running binary.

## What the current rotation code does

- `Exam/exam.c:75-83` defines `GFS_LAGUERRE_ROTATION` as an extension of the existing face-normal-velocity correction from Voronoi faces to weighted Laguerre faces. The W=0 midpoint/centroid correction is a separate, already-active branch.
- `Exam/exam.c:4729-4754` and `Exam/exam_gpu.cu:600-626` use
  `u_n(c_f) = u_n(a) - (v_j-v_i)·(c_f-a)/d`,
  where `a` is the midpoint for W=0 or the weighted bisector anchor
  `x_i + [(1 + (w2_i-w2_j)/d²)/2](x_j-x_i)` for a Laguerre face. This corrects the normal speed used by the face flux when the face centroid is offset from its anchor.
- `Exam/exam.c:5336` also enables the Laguerre-face correction while dynamic W2 is active.
- `ParticleSoA` in `Exam/exam_gpu.h:74-132` stores translational position, velocity, mass, thermodynamic fields, gradients and stress. `GPUPhysicsParams.laguerre_rot` is a switch for the face correction; it is not a cell angular-velocity state. No per-cell spin momentum, moment of inertia, or spin energy is present in these structures.

Therefore, the present correction accounts for face kinematics. It does not add a persistent cell rotation degree of freedom, a torque evolution equation, or rotational kinetic energy. The active W=0 run tests the baseline face correction, fixed mass, RK4, and the disk evolution; it does not test spin dynamics.

An earlier run-data audit (`spin_balance_audit_20261030.md`) also found that reconstructed-field spin energy is a diagnostic component of the resolved kinetic energy, not a new reservoir to add unchanged. Any independent spin state must either remove that rigid-rotation component from the resolved velocity field or define a separate subgrid degree of freedom. That audit further derives the active `GFS_VARPRESS_W0` force and work ledger; its restrictions are included below.

## Constraints for a conservative spin model

1. Define the cell reference point before defining spin. A Voronoi/Laguerre generator is not generally the mass centroid. A generator-centered split can contain offset and cross terms; a centroid-centered split requires the centroid's motion to be consistent with the conserved cell momentum.
2. In 2D, store an extensive internal angular momentum `L_spin` and a clearly defined inertia `I`; if a rigid-spin closure is used, `omega = L_spin/I` and `E_spin = L_spin²/(2I)`. For a uniform-density polygon, its centroid and second area moment can be computed from its vertices without subcell quadrature. Any nonuniform subcell closure needs a separate cost and accuracy case.
3. Derive torque from the actual discrete force operator. The familiar torques `(x_f-x_i)×F` and `(x_f-x_j)×(-F)` apply only to a canonical shared-face traction that is equal and opposite. The active Kepler `GFS_VARPRESS_W0` branch is a directed-cell `BᵀP` operator; its two cell-owned contributions are not generally `F` and `-F`. The canonical pair formula is therefore not a production spin update for this branch. First establish the complete global torque identity, including cell-face closure and the polygon-wall term, or derive a new couple-stress operator together with its energy ledger.
4. Couple the active work flux to translation, spin and internal energy in one budget. The current `GFS_VARPRESS_W0` branch uses cell-centered `v_i·F_i - P_i Vdot_i`; it does not already use a single shared `F·u_f` work expression. A face-velocity correction and a spin-induced contact velocity must be derived with that exact work path so face work is neither omitted nor counted twice.
5. If inertia changes with cell shape, include `dE_spin/dt = omega·dL_spin/dt - 0.5*omega²*dI/dt`. A torque-work term alone does not close the spin-energy ledger for a deforming cell.
6. The user's `Re → ∞` condition means zero limiting shear viscosity. It supplies no viscous spin-down mechanism. Spin can exchange angular momentum and energy through conservative pressure/stress fluxes, but any dissipative relaxation must come from a deliberately specified finite-viscosity or subgrid closure, not an arbitrary damping factor.

## Gates before a Kepler spin run

- A rigid-rotation patch test checks preservation of uniform spin and the expected face velocity.
- A canonical two-cell off-axis traction test checks pairwise momentum, total angular momentum (orbital plus spin), and total energy to roundoff; it is only a unit test of that canonical operator.
- Before using the production W=0 operator, compare its global cell torque sum against the closed-cell face-closure identity and measured polygon-wall torque, and against a fixed-wall virtual-rotation finite difference of the discrete internal energy. Do not infer global angular-momentum failure from non-antisymmetric individual directed faces.
- A shearing patch test checks the resolved stress work and the partition among orbital kinetic, spin kinetic, and internal energy.
- A Kepler comparison records cell-centroid motion, `L_spin`, `E_spin`, total angular momentum, total energy, entropy/floor injection, and density-field error against the same t=0 reference. Compare the face-kinematics-only baseline with spin enabled at matched resolution and timestep criteria.
- Retain the H200 GPU path and measure added wall time per RK4 step; reject a formulation that needs expensive subcell quadrature or materially slows the AREPO-scale evolution.

The current long-time W=0 continuation is the baseline for these comparisons. Its result must not be relabeled as a spin validation.
