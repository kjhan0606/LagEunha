# Hard reconstruction support crossing identified

Production Exam/exam.c selects candidates within d<=2.5sqrt(V), fits
cubic weighted least squares (nine columns), and expands distance shells
only if rank is insufficient. The row multiplier is sqrt(V)/d. Thus an
entering/leaving candidate has nonzero weight at the support boundary.

Independent fits from actual saved RK call2/3 states reproduce all nine
selected cells' production gradients within4.93e-7, consistent with their
rounded gradient output. Both stages have rank9, not a rank fallback.
For the minus branch:

| Cell | Lost candidate | d/sqrt(V), before -> after | Count | Condition number |
|---|---|---|---|---|
|155|214|2.49996590 -> 2.50005640|21 -> 20|5.89 -> 7.01|
|848|626|2.49990837 -> 2.50006203|20 -> 19|6.26 -> 6.07|

Face-neighbor sets remain unchanged. Reconstruction support differs from
face topology. The plus branch loses neither candidate over this stage.
Production candidate lists were not exported; the full-point-cloud fits
and their matching gradients independently reproduce the transition.

## Local remedy candidate, not yet adopted

Apply an additional row multiplier
T=(1-(d/(2.5sqrt(V)))^2)^2 for d<2.5sqrt(V), and zero outside.
The least-squares objective therefore uses T^2 as the additional weight.
This makes a sample's influence vanish continuously at the existing
support radius. These reconstruction weights are NOT Laguerre cell weights.

Gradient change norms between actual RK call2 and3:

| Branch/cell | Existing hard support | Fixed union control | Smooth taper |
|---|---:|---:|---:|
|minus155|.02102880|.00094109|.00094308|
|plus155|.00094109|.00094109|.00094308|
|minus848|.00454002|.00042295|.00041939|
|plus848|.00042295|.00042295|.00041938|

The fixed-union control isolates the effect of changing sample membership;
it is not a proposed frozen-stencil evolution. Smooth fits remain rank9
and recover cubic polynomial coefficients to<1e-10 in these local stencils.
No steering, physical force removal, or cell mass change was used.

This explains a localized Q/reconstruction transition and its contribution
to the short-step diagnostic discrepancy. It does not establish the cause
of historical t6-10 physical failure or prove a stable long-time solution.
The next implementation should be opt-in, retain legacy behavior by
default, explicitly handle shell expansion/fallback, and verify candidate
search coverage. Test polynomial/rigid-motion consistency, actual single
steps and accepted budgets, then matched saved-state trajectories through
historical failure times before claiming a remedy. Verify backend scope;
the current production hydro path remains CPU, CUDA is used for neighbors.

Evidence: canonical_coupled_shell_transition_2026-10-18.json,
canonical_coupled_smooth_reconstruction_2026-10-18.json.
Scripts: audit_coupled_reconstruction_shell.py,
audit_coupled_smooth_reconstruction.py.
