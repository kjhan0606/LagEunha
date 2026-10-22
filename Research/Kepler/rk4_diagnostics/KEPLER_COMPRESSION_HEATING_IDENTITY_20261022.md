# Separating bulk compression from subcell compression energy

2026-10-22. Derived from the actual native generating-context thermal
functional, not an additional physical closure or a proposed heat correction.

## Exact discrete identity for one unchanged positive rule

The evolved CUDA kernels use, at each reference point q:

    j0_q = generating_J_q + generating_delta_cell > 0
    j_q  = current_J_q    + current_delta_cell    > 0
    lambda_q = j_q/j0_q
    U_q = P0*w_q*j0_q/(gamma-1) * lambda_q^(1-gamma)

For this fixed rule define

    V0 = sum_q w_q*j0_q
    V  = sum_q w_q*j_q
    alpha_q = w_q*j0_q/V0,  sum alpha_q = 1
    m = V/V0 = sum alpha_q*lambda_q
    U0 = P0*V0/(gamma-1)

Since gamma>1, f(lambda)=lambda^(1-gamma) is strictly convex. Jensen's
inequality, directly applied to this positive finite sum, gives

    U/U0 = sum alpha_q*f(lambda_q) >= f(m)
    U_uniform = U0*(V/V0)^(1-gamma)
    U - U_uniform >= 0

Equality requires uniform relative compression lambda_q=m at every positive
node. This is about relative j/j0, not the variance of raw J alone. Initial
j0 need not equal1. The native fresh constant-P0 generating law makes these
weights applicable even if the initial density varies between reference nodes.

At m=1 and small fluctuations e=lambda-1, with weighted mean e=0:

    (U-U0)/U0 = gamma*(gamma-1)/2 * <e^2>_alpha + O(<abs(e)^3>)

For gamma=5/3 the quadratic coefficient is5/9. The native gamma is the exact
double conversion of its recorded float value; the saved-data analysis reads
that value rather than substituting an ideal5/3. This leading expansion does
not measure pointwise variance when point values were not saved.

## Physical interpretation and limits

Exact circular inviscid Kepler flow has j/j0=1 and no compression heating.
An approximate material map can preserve the whole-cell volume while creating
nonuniform subcell compression. Its positive U increase is then consistent
with the same fixed-entropy thermal law and pressure work. Energy/torque
conservation alone does not exclude such an inaccurate flow. It need not be
irreversible entropy production or artificial viscosity; the current smooth
operator keeps the generating entropy fixed. Whether this effect dominates
the historical failure requires actual trajectory evidence, not this identity.

This diagnostic separates the uniform bulk-compression energy from the
remaining convex compression energy. It must NEVER be used to reset U to
U_uniform: doing so would erase resolved thermal physics and change the
conserved state/entropy. No correction or damping is implemented.

## Saved-data implementation and precision scope

`analyze_kepler_compression_heating_v1.py` uses all71080cells of the already
accepted original-rule step4. It obtains U0 from the existing audited initial
cell ledger and reconstructs V0=(gamma-1)*U0/P0 using the exact recorded binary
law. Current V and U come from the existing endpoint cell CSV. It computes

    bulk_fraction = (V/V0)^(1-gamma)-1
    total_fraction = U/U0-1
    signed_Jensen_residual = total_fraction-bulk_fraction

Only CSV/binary-header reads and55-digit Decimal arithmetic are performed.
No new initial state, point evaluation, mass solve or RK4 step is executed.
All source signs, including negative near-roundoff residuals, are retained.
Canonical V and summed thermal measure agree only to their finite numerical
accuracy. Hence tiny negative residuals are not clamped, and tiny positive
residuals cannot automatically be identified as physical or artificial heat.
The same saved summaries do not recover the actual pointwise j variance.

This first diagnostic is explicitly restricted to the unchanged ORIGINAL
positive rule. A future compressed-rule state requires its own initial
positive-measure context and the separate numerical quadrature discrepancy
ledger; the old U0 cannot silently stand in for a changed initial measure.

See material_full_compression_heating_audit_v1/contract.json and cells.csv for
actual values, precision limits and the exact source hashes. t6/8/10 and
longest-GIZMO stability remain unverified.
