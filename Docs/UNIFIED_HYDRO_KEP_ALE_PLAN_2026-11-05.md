# Selected replacement: cell-mean kinetic-compatible ALE flux

One grouped Astra evaluation of the actual failed cold-flow bundle selects
this candidate. It does not reject all FV or weight-linked mass methods.
There is no additional timestep sweep of the retired HLLC candidate.
The driver implements and runs this plan without per-test agent calls.

## A common cell-pair closure

Use the existing stored m,p,E, declared mean-state EOS, actual power cells,
limited physical generator velocity, immutable birth mu and SSPRK2.
Let averages be arithmetic unless marked logarithmic, and set

    beta_i = rho_i/(2 P_i).
    rho_log = logmean(rho_i,rho_j).
    e_log = 1/[2 (gamma-1) logmean(beta_i,beta_j)].
    P_hat = rho_bar/(2 beta_bar).
    u_bar = (u_i+u_j)/2; kappa = dot(u_i,u_j)/2.

Apply the core to the STORED CELL MEANS, not independently reconstructed
face states:

    Fm = rho_log (dot(u_bar,n)-wn).
    Fp = Fm u_bar + P_hat n.
    FE = Fm (e_log+kappa) + P_hat dot(u_bar,n).

The spatial KEP/EC construction follows
[Chandrashekar](https://arxiv.org/abs/1209.4994); the ALE temporal entropy
state must be included in the identity, not approximated by a physical
Riemann star state. These algebraic identities do not prove positivity or
finite-step stability. See also the distinction between flux identities
and robust evolution in [Ranocha](https://arxiv.org/abs/1701.02264).

## One sign-compatible compression stress

Use the existing velocity reconstruction/limiter to define the actual
cell-pair normal jump delta and reconstructed face jump delta_rec.

    s = max(0,-minmod(delta,delta_rec)).
    shat = [(h_i+h_j)/2]/|q_j-q_i| * s, h_i=V_i^(1/d).
    Pi = rho_bar (c_bar shat + 2 shat^2).
    Fp += Pi n; FE += Pi dot(u_bar,n).

Do not change mass flux or add another divergence switch. This vanishes
for affine velocity reconstruction and detects an alternating mode with
zero slopes. Positive Pi implies delta<0. Pair kinetic loss is
`A Pi delta<=0`, thermal production is nonnegative, and contraction with
the entropy-variable jump is `2 beta_bar Pi delta<=0`.
This is numerical shock entropy production at physical Re=infinity,
not a physical viscosity parameter or a thermal reset.

## Simpler weight solve and explicit limits

Because the core mass flux is linear in wn, the unchanged birth-mu law
requires one scalar SPD solve per RHS:

    [diag(mu)+L_rho_log] z = b, edge coefficient A rho_log/(2 r).

Keep PCG cap 64, the true residual criterion and the integrated
mass--weight invariant. No Newton, mu scan, RHS mean removal or dynamic
weight recentering is needed. Keep wave CFL 0.4, cone safety 0.9 and
whole-step rollback cap 8. Add the stress stability bound

    D_f=A Pi/(-delta) for delta<0; otherwise D_f=0.
    dt_nu <= 0.4 min_i m_i/sum_faces D_f.

No implicit viscosity solver is added. Boundary and residual ledgers stay
explicit. Reduced solver structure is not a measured GIZMO speed result.

## Separate spatial kinetic work from RK error

On every actual RHS, record the pair kinetic work

    Kdot_ij=A dot(delta_u,Fp-u_bar Fm).

On every accepted full step, record

    deltaK_time=K_final-K_initial-dt*(Kdot_stage0+Kdot_stage1)/2.

This is separate from the existing GCL defect. SSPRK2 is not exactly KEP
in time. Keep all rejected work in the measured cost.

## Registered execution and stopping criteria

1. Native dimensional kinetic/entropy/Galilean/contact algebra, actual
   geometry and boundary entropy ledgers, affine shear and alternating
   modes. Separate core/stress on the existing analytic cold birth fields
   at 64/128. The Cartesian sine-shear birth RHS must eliminate the old
   kinetic-scale thermal source. Algebra-only success is not PDE success.
2. Same-settings actual 1D Sod200/400, blast200, Mach8 collision200; cold64
   to t=0.25 and Gresho64 to t=1. Register one cold64 step-factor-1/4
   control. Preserve original shock ICs and endpoints.
3. If these pass, cold128/Gresho128, the remaining original 1D tests, and
   small 3D contact/Sod. Every new trajectory begins at its analytic birth
   IC. Do not continue failed late HLLC cold states under the new core.
4. Only after physical qualification, run isolated matched GIZMO cold/
   Gresho timing. Same analytic IC, EOS, resources, endpoint and output
   policy; report initialization+evolution and evolution-only wall time.

The cold physical pressure relative L1 limit remains 1%. Positivity-driven
time collapse, PCG/rollback cap, or failed shock refinement rejects this
candidate. If only the quarter-step control passes, time accuracy remains
unresolved; do not select a cold-only CFL. Gresho must improve over W0 and
with refinement; disclose differences from the retired HLLC warm-flow
result as well. No floor, hidden reservoir, kinetic cap or energy repair.

## Implementation checkpoint

`Exam/gfs_kep_ale.h` implements the pure fixed-five-slot cell-pair core.
`Exam/HydroCore/kep_flux_audit.c` passes 240 mixed native dimensional
fixtures in optimized C and C++: ALE entropy discrepancy 3.05e-16,
kinetic-chain discrepancy 5.81e-16, boost discrepancy 4.36e-15, and
compression entropy discrepancy 1.17e-14. Evidence is `kep_flux_v1`.
The subsequent native coupling, registered stress, one-SPD solve,
kinetic/entropy ledgers and eleven initial PDE endpoints are documented
in `UNIFIED_HYDRO_KEP_NATIVE_2026-11-05.md`. Both cold64 pressure gates
still fail (1.58%/1.85%), so the candidate is not adopted. Original blast/
collision, higher-resolution qualification and the remaining suite are
unfinished; the selected kinetic core and the remaining stress defect
are distinguished rather than discarded together.
The previous primitive CUDA audit belongs to the HLLC/HLL core; it is
not reused as a CUDA qualification of this new flux.
