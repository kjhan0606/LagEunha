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

## Constraints for a conservative spin model

1. Define the cell reference point before defining spin. A Voronoi/Laguerre generator is not generally the mass centroid. A generator-centered split can contain offset and cross terms; a centroid-centered split requires the centroid's motion to be consistent with the conserved cell momentum.
2. In 2D, store an extensive internal angular momentum `L_spin` and a clearly defined inertia `I`; if a rigid-spin closure is used, `omega = L_spin/I` and `E_spin = L_spin²/(2I)`. For a uniform-density polygon, its centroid and second area moment can be computed from its vertices without subcell quadrature. Any nonuniform subcell closure needs a separate cost and accuracy case.
3. Apply every face force at one shared face point. For force `F_ij` at `x_f`, the spin torques about the two chosen cell reference points are `(x_f-x_i)×F_ij` and `(x_f-x_j)×(-F_ij)`. Use the same force and point on both sides. This makes the spin torque cancel the orbital torque in the global angular-momentum budget, subject to the chosen reference-point and moving-volume bookkeeping.
4. Couple the work flux to translation, spin and internal energy in one budget. A face-velocity correction and a spin-induced contact velocity must be derived together so the face work is counted once. The current geometric correction must not be mistaken for `E_spin` or added a second time without a derivation.
5. The user's `Re → ∞` condition means zero limiting shear viscosity. It supplies no viscous spin-down mechanism. Spin can exchange angular momentum and energy through conservative pressure/stress fluxes, but any dissipative relaxation must come from a deliberately specified finite-viscosity or subgrid closure, not an arbitrary damping factor.

## Gates before a Kepler spin run

- A rigid-rotation patch test checks preservation of uniform spin and the expected face velocity.
- A two-cell off-axis force test checks pairwise momentum, total angular momentum (orbital plus spin), and total energy to roundoff.
- A shearing patch test checks the resolved stress work and the partition among orbital kinetic, spin kinetic, and internal energy.
- A Kepler comparison records cell-centroid motion, `L_spin`, `E_spin`, total angular momentum, total energy, entropy/floor injection, and density-field error against the same t=0 reference. Compare the face-kinematics-only baseline with spin enabled at matched resolution and timestep criteria.
- Retain the H200 GPU path and measure added wall time per RK4 step; reject a formulation that needs expensive subcell quadrature or materially slows the AREPO-scale evolution.

The current long-time W=0 continuation is the baseline for these comparisons. Its result must not be relabeled as a spin validation.
