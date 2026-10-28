# Original and inherited entropy decomposition at matched t6

Read-only analysis of complete legacy/adaptive checkpoint1783 at identical
time6.000405917763673. Initial and source IDs/masses match, weights are zero;
the t2 source path and SHA are verified against the adaptive run manifest.
No trajectory, pressure or entropy was changed.

At each actual final volume V, define P0iso=P0(V0/V)^gamma and
P2iso=P2(V2/V)^gamma. Then

    P-P0 = (P0iso-P0) + (P2iso-P0iso) + (P-P2iso).

These are initial-entropy compression, inherited entropy before the t2
restart, and entropy gained after the restart. Normalized pointwise closure
is at most1.78e-15. The exact source time is2.0029059177634028.

For each component c, its signed projection is <c,e>/<e,e>, where
e=(P-P0)/P0 and brackets are initial-mass-weighted over804 inner-background
cells. Projections sum to1; RMS values do not add.

| Component | Legacy RMS | Adaptive RMS | Legacy projection | Adaptive projection |
|---|---:|---:|---:|---:|
| Initial-entropy compression | .613674 | .671282 | .362463 | .382855 |
| Inherited pre-t2 entropy | .00783990 | .00808263 | .00481753 | .00477507 |
| Post-t2 entropy | .917956 | .939181 | .632720 | .612370 |
| Total pressure error | 1.317010 | 1.374947 | 1 | 1 |

The directly retained pre-t2 entropy component is small in the final pressure
error projection. This does NOT imply that the early seed or inherited
position/velocity errors were dynamically irrelevant. Post-t2 entropy and
geometric compression both contribute substantially on the actual evolved
geometry. Adaptive reconstruction has not reduced either component's RMS
at this endpoint. A smaller fractional projection is not a smaller absolute
error: adaptive post-t2 RMS is higher despite its lower projection fraction.

These are pressure-state decompositions, not force-work budgets, heating
source attribution, independent causal fractions, or a justification for
removing entropy. Next apply the actual reversible pressure operator to
these components and include TCF/HG, gravity and reference-motion work.
Preserve the ongoing adaptive t10 run and completed legacy reference.

Evidence: canonical_t6_entropy_pressure_history_2026-10-18.json.
Script: audit_t6_entropy_pressure_history.py.

## Actual pressure force and work

Two terminal MPI1/OMP1 probes from the matched t6 snapshots completed with
exit0 and no accepted evolution. Source hashes were unchanged. The adaptive
binary is the existing validated build, with adaptive/smooth flags enabled
for adaptive and disabled for legacy. Exported x,v,P,V,m match snapshots.
The pressure force reconstructed from common-face geometry matches hydro
RHS minus rounded TCF faces to relative L2 <=1.034e-7. Three-component force
closure is <=2.12e-21.

Inner-background sum((v-v_circular) dot F_pressure), in simulation units:

| Pressure source | Legacy | Adaptive |
|---|---:|---:|
| Initial entropy at actual volume | 2.63342e-6 | 3.06219e-6 |
| Inherited entropy before t2 | 1.68949e-8 | 1.94849e-8 |
| Entropy added after t2 | 1.77295e-6 | 1.94691e-6 |
| Total | 4.42327e-6 | 5.02858e-6 |

Both main pressure components supply departure kinetic energy at this
instant. The initial-entropy component accounts for about61% of adaptive
pressure work, even though post-t2 entropy dominates the pressure-error
projection. Pressure-error proportions therefore cannot substitute for work
budgets. Different trajectories/geometry/velocity contribute to differences;
this is not a same-state operator counterfactual. TCF/HG damping and
gravity/reference-motion contributions remain necessary before interpreting
the net growth. No pressure or entropy reset is justified by these results.

Evidence: canonical_t6_entropy_pressure_work_2026-10-18.json.
Scripts: run_t6_entropy_force_probes.py, analyze_t6_entropy_force_probes.py.
Remote root: /gpfs/kjhan/LagForce/Codex_t6_entropy_force_probes_20261018.

## Complete instantaneous departure generator

Actual HG mode assembly was reconstructed with frozen-support work/heat
accounted for. Its kinetic-plus-thermal balance residual is <=5.79e-21.
Analytic circular-reference directional derivatives agree with independent
central differences to2.71e-10 at the smaller increment.

Inner-background derivative of half sum(m*|v-vcirc|^2):

| Term | Legacy | Adaptive |
|---|---:|---:|
| Reversible pressure | +4.42327e-6 | +5.02858e-6 |
| Hydro RHS minus pressure | -3.58436e-6 | -4.02727e-6 |
| HG generator | -1.50795e-6 | -1.58154e-6 |
| Gravity plus moving circular reference | +7.83036e-7 | +9.93271e-7 |
| Total | +1.13995e-7 | +4.13037e-7 |

Both states have positive instantaneous departure growth. Hydro-minus-pressure
and HG damp it; pressure and the combined orbital term supply it. The larger
adaptive total is a difference of large terms on different evolved states,
not an instability eigenvalue or a same-state operator effect. Do not quote
the gravity term without its largely cancelling circular-reference term.
Hydro remainder includes TCF and any other assembled hydro terms; rounded
TCF face reproduction was checked separately but is not exact term purity.
This is a continuous instantaneous generator, not a validated finite-step
HG/RK4/HG energy budget. The unresolved prior generator/map discrepancy is
not closed by this analysis. Long physical tests remain necessary.

Evidence: canonical_t6_complete_departure_work_2026-10-18.json.
Script: analyze_t6_complete_departure_work.py.
