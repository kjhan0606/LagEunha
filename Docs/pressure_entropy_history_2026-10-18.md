# Early versus late pressure work: reuse of actual saved geometry

Reused actual t2 and t6 RHS/common-face outputs, and the completed t10
analysis. No new trajectory was launched or modified. The generalized
audit_pressure_entropy_force.py accepts --root and --label, retains
exclusive creation of the output, and checks the reconstructed force
against saved hydro-minus-Q on r<1.9.

For fixed mass, P_iso=P0*(V0/V)^gamma and P_excess=P-P_iso. Linearity of
the reversible pressure operator gives an exact algebraic force split.
Initial-ID material masks are identical across the three analyses.

| Time | Inner actual pressure work | Isentropic | Entropy excess |
|---|---:|---:|---:|
| t2 | +1.40147e-9 | +1.44419e-9 | -4.27131e-11 |
| t6 | +4.44146e-6 | +2.61790e-6 | +1.82356e-6 |
| t10 | +2.80264e-5 | +3.16280e-6 | +2.48636e-5 |

Work is instantaneous sum F_pressure dot(v-v_circular), not an integrated
energy increment. Entropy-excess fractions of inner pressure work are
about -3.05%,41.06%,88.72%. These are signed algebraic fractions; they are
not cumulative heating fractions or fractions of the total growth rate.
The entropy-excess component includes numerical entropy error, so it is
not independently identified as AV heat in this calculation.

Reconstructed pressure-force relative L2 errors are4.51e-8 at t2 and
1.09e-7 at t6, consistent with rounded saved Q. Linear-force splitting
closes to1.07e-21 and1.91e-21 absolute respectively.

At t2 the pressure excitation is predominantly the isentropic component;
entropy excess slightly damps inner departure at that instant. By t6 the
excess contributes positively and by t10 it dominates the pressure work.
This separates early excitation from later thermodynamic feedback. It
contradicts treating late entropy dominance as proof that initial AV
heating alone initiates growth. It does not prove a unique instability
mechanism or distinguish geometric forcing from a coupled growing mode.

The recent accepted t2.105-to2.205 budget further shows direct TCF/HG
damping while other hydro work and shear supply departure. Next causal
priority is a coupled position/pressure response at early t2, because
the earlier velocity-only response damped. Do not erase entropy, remove
pressure work, or steer circular trajectories as a stabilization fix.
Any candidate must still pass the historical t6-10 failure interval and
long evolution with the approved weight/mass reduction constraints.

Evidence: canonical_pressure_entropy_history_2026-10-18.json plus linked
per-time source JSON files. Long-term physics stability remains unresolved.
