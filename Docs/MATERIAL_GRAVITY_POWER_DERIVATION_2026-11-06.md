# Candidate gravity coupling derived from the existing material power identity

Derivation only; not adopted or implemented in evolution. Preserve the successful Kepler baseline. Constant gravitational field is the initial scope; current RT source treatment is not qualified.

Let Epot(x,W)=sum_i m_i Phi(c_i), with fixed masses and actual volume centroids, uniform density within each cell, and linear Phi. Let gx=dEpot/dx and gW=dEpot/dW. The current material geometry imposes

    Vdot = B v + L q = C v + D P,

where q=Wdot and L is the symmetric weight graph. Solve L eta=gW, modulo a constant (sum gW=0 from weight-gauge invariance). Then

    Epotdot = [gx+(C-B)^T eta] dot v + eta^T D P.

A candidate compatible coupling replaces the dissipative pressure drive with Pe=P-eta and uses

    Vdot = C v + D Pe,
    F = C^T Pe - K v + B^T eta - gx,
    Udot_i = -P_i Vdot_i + H_i,
    sum H_i = Pe^T D Pe + v^T K v.

These equations cancel total kinetic+internal+potential power algebraically. The acoustic impedance must still come from physical rho/P, not from an arbitrary potential-gauge-dependent Pe. Positive local heat requires the same squared-jump distribution with the effective pressure drive. Boundary pressure offsets cannot simply be layered on top. Actual geometric Wdot must realize Vdot, and active constraints require their own rederivation with these source terms; the present projector is not automatically qualified.

At uniform rho and fixed domain, the virtual-work identity yields gx=-B^T(rho Phi(c)) and gW=-L(rho Phi(c)); hence eta=-rho Phi(c)+constant. The force becomes C^T(P+rho Phi(c))-K v. Therefore P+rho Phi(c)=constant has zero force/heat, provided C^T1=0 and the corresponding trace jump of constants vanishes. This is derived from potential work, not subtraction of a measured hydrostatic residual. It does not yet prove accuracy for density jumps or nonlinear central gravity.

## Exact face derivatives in 2D

For a unique interior face oriented i->j, let n=(xj-xi)/d, area A (segment length), face centroid f, midpoint m=(xi+xj)/2, tangent t perpendicular to n, and g=grad Phi. Define

    psi = rho_i [Phi(f)-Phi(c_i)] - rho_j [Phi(f)-Phi(c_j)],
    k = A^3 (rho_i-rho_j) (g dot t) t / (12d).

The exact linear-potential contribution for rotating straight faces is

    gx_i += A psi [n/2+(f-m)/d] + k,
    gx_j += A psi [n/2-(f-m)/d] - k,
    gW_i += A psi/(2d),
    gW_j -= A psi/(2d).

The k term is the covariance of affine potential and affine face velocity, not an empirical correction. Fixed impermeable walls contribute zero. These formulas need independent finite-difference checks of actual Laguerre cells, including unequal density and face rotation. Products of centroid values alone miss k. In3D, the analogous term needs the face second-moment tensor; no3D implementation is authorized before2D qualification.

Cost: one additional graph solve for eta per RHS unless an equivalent cheaper reduction is found. Reuse of an existing matrix does not prove the user's runtime requirement. Benchmark before adoption; no per-step volume quadrature or new physical parameter is proposed. This derivation is a candidate within the existing fixed-mass/weight framework, not a claim that the method is complete.

Native geometric derivative follow-up: Exam/gfs_material_potential2d.h implements gx/gW only. Independent half-plane reconstruction on4irregular weighted16-cell partitions and16virtual-motion directions matches native gradients to4.56e-12 atdifference step1e-5; omitted face covariance causes error up to.0032293. Equal/unequal density, face rotation and weight changes included. No PDE integration or runtime gravity adoption. Evidence material_potential_gradient_20261106_v1.

Native graph-power reduction also verified on the same4irregular graphs: existing PCG solves for eta and q give direct/reduced potential power agreement <=9.76e-18. This leaves the physical-impedance hydro coupling, constrained variant, time integrator, native RT adapter and cost unqualified; no evolution adoption.

## Active weight constraints with compatible gravity

The generalized pressure drive is Pe=P-eta, with the force source B^T eta-gx. Physical pressure still sets the impedance and internal pressure work. The existing projector now accepts this optional drive/source, with its old API delegating to NULL inputs. Potential power is evaluated from the *projected* Wdot, not the unconstrained rate.

Native C check_active.c tests an 8x8 nonuniform weighted partition, periodic x and stationary reflected y walls. For each of two perturbed states (uniform and varying density), every right/top interior face is tried with both constraint signs. Of 480 single-constraint trials, 241 outward cases accept and 239 inward cases reject consistently with independently calculated unconstrained rates. Maximum normalized total kinetic/internal/potential power residual is 6.024e-14; maximum constraint-rate residual is 2.711e-13. Accepted cases have nonnegative local heat and net volume residual below 1e-9. Unconstrained uniform and two-layer hydrostatic cases, and zero-gravity equivalence, also pass.

The old eight-state wall projection regression reproduces every saved diagnostic exactly, including 11 accepted orientations, 5 rejected orientations and zero same-RHS cache discrepancy. Evidence: material_gravity_rhs_20261106_v1/active_results.jsonl and wall_drive_regression.json.

Scope: instantaneous single-active-constraint algebra only. This is not RT evolution, a multiple-active-set proof, CUDA/MPI gravity qualification, time-discrete energy conservation, or a GIZMO speed comparison. Time integration and the native RT adapter remain pending.

## SSPRK2 planar time integration check

`material_gravity_rhs_20261106_v1/check_time.c` calls the existing gravity RHS and existing SSPRK2 stage helper, rebuilding exact rectangular power cells at both stages. Native C evolves 8x8 cells to t=0.01 at dt=0.001, 0.0005, 0.00025, 0.000125. The fixture is uniform along periodic x, nonuniform and weighted along bounded y, with smooth density/pressure/velocity perturbations. Mass is copied unchanged; no extra gravity kick or thermal repair is applied. Geometry positivity, positive internal energy and planar symmetry are checked. Adjacent-resolution velocity differences give temporal orders 2.00087 and 2.00044. Total kinetic/internal/centroid-potential energy residual is at roundoff (maximum absolute 2.665e-15); this special planar result does not prove exact energy preservation for general moving polygons.

An initial nonuniform-x variant failed the planar-fixture assumption at the first predictor: transverse forces up to about 1.3e-8 appeared. Its source is retained as `check_time_nonuniform_rejected.c`; it is not a passing evolution result. Do not silently constrain those forces or use the rectangular geometry for that trajectory. General 2D evolution requires the actual polygon geometry adapter. The first monotonic-energy convergence assertion was also inappropriate once errors reached roundoff; qualification instead uses velocity self-convergence plus an energy bound.

This is a planar reduction of the 2D operator, not KH/RT completion, native production integration, GPU/MPI qualification or the runtime benchmark.
