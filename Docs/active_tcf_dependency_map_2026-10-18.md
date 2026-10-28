# Actual active TCF dependency map and FACE_Z override

Source: `Exam/exam.c` TCF assembly around9016--9138, with caller setup around
8447--8472. Existing active flags include QBP, CHI, EDGE_INNER=3 and FACE_Z.

| Stage | Inputs | Actual output / constraint |
|---|---|---|
| Common face geometry | area A, separation d,e, face centroid | tangent offset t projected perpendicular to e |
| QBP reconstruction | unlimited velocity gradients; Hessians/third derivatives at arc interfaces; validity flags | delta=du_cell+d e.(Gi+Gj).e/2; invalid gradients give qin=0 and count failure |
| Compression blend | symmetric gradients, divergence | chi=max of two cell compression shares; qin=(1-chi)delta+chi*du_cell |
| Limiter | qin, du_cell | mm=minmod(qin,du_cell) |
| FACE_Z | cell rho/P, density-centroid projections, d, gamma, mm | fi/fj floored at1e-12; rho_face=rho*f, P_face=P*f^gamma; Zq=max(ZD,ZH) |
| Dissipative pressure | applied Zq, mm, Q0/QLAW flags | default Q=Zq*mm |
| Force / work | cell P,v and A,d,e,t,Q | compatible pair force, total-energy work, heat and geometric volume rate |

With FACE_Z enabled, Zq is overwritten by the face-state impedance. The earlier
tcf_pstar difference quotient is calculated and available diagnostically but
does not directly supply the applied Zq in this branch. Therefore generic-HLLC
versus-TCF pressure-helper differences are not direct active-force differences.

The force/work expressions are

    F_i=-A[(Pbar+Q)e+(Pj-Pi)t/d],
    Edot_i=-A[(Pj*ui+Pi*uj)/2+(Pj*vi-Pi*vj).t/d+Q(ui+uj)/2],
    heat_i=A*Q*(ui-uj)/2,
    Vdot_i=A[(ui+uj)/2+(vi-vj).t/d].

For paired common geometry, the dissipative kinetic work is balanced by the
two nonnegative heat contributions. HG and subsequent RK-stage handling remain
separate parts of the full solver; these expressions are not the entire RHS.

## Saved-face readback

Read existing t2/t4/t6/t10-named diagnostic directories,221159 directed faces
each. In83965/84215/82779/84125 nonzero-Q faces, applied Zq differs from the
legacy stored Z by more than relative1e-5. Actual Q agrees with
Zq*minmod(qin,du_cell) to relative L2 about1.88--1.90e-7, consistent with the
roughly7-significant-digit text dump. Directed heat is nonnegative.
Replacing Zq by legacy Z in a frozen diagnostic would change summed heating
to about.9495--.9498 of recorded heat. This is NOT an evolved alternative or
recommended change. Directory names alone are not used here as fresh timing
verification of each snapshot.

## GPU/MPI implementation boundary

A faithful GPU path needs the above reconstruction/geometry caches, ownership,
pair symmetry and force/work accumulation, not merely a CUDA tcf_pstar helper.
The existing GPU unsupported guard remains. Next isolate and compare the
active algebra using recorded face inputs, then separately verify extraction,
cache freshness and distributed ownership before any full-kernel enablement.
Existing reconstruction/heat tests should be reused; do not rerun long Kepler
trajectories to rediscover these scalar identities.

No production code or running calculation changed. CPU physical stability,
full CUDA/OpenMP/MPI parity and time-dependent Laguerre reduction remain open.

Evidence: `canonical_active_face_z_dependency_2026-10-18.json`.
Script: `audit_active_face_z_dependency.py`.
