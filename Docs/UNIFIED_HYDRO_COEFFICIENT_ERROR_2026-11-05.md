# Constrained-weight coefficient error and frozen physical response

This finite stage continues `fd91a39`. The overall dimensional-hydrodynamics
objective remains incomplete. The constrained-weight candidate is not released
for PDE evolution. The previous CPU coefficient cost failure remains in force.
The completed work distinguishes coefficient-solve error from the spatial
response of the proposed operator itself.

## Unchanged native capture limits

The C/OpenMP coefficient probe now has an optional diagnostic capture path,
`MATERIAL_VOLUME_COEFFICIENT_CAPTURE`. The default probe still stops on the first
failed coefficient. Capture mode can observe the second RHS after nonconvergence
only within the same shared time budget. A first failure is never cleared.
Capture I/O occurs after the candidate timer. No production fluid code changed.

One native capture was performed on the original 128² Gresho t=0.5 graph.
The numerical limits remained 256 iterations per RHS, relative tolerance 1e-10,
and a shared 0.076127-second budget. Both fields reached 256 iterations with
status 3. True relative residuals were 2.6413754e-7 and 5.6032240e-7. The timed
preparation and two incomplete solves took 0.043036 seconds. R/R-transpose,
weighted geometry, and PDE evolution were not timed or released. No repeat was
used to select a favorable capture. This result is not a passed cost screen.

The independent sparse reference from the previous snapshot was reused without
another global solve. Saved-reference residuals remain about 4.6e-14 and 5.9e-14.
Mean-zero coefficient differences relative to the reference have L2 fractions
2.6172e-5 and 4.7944e-5. A residual of a few times 1e-7 therefore does not mean a
coefficient error of the same magnitude.

## Residual-to-operator identities

For each traceless symmetric strain basis, define

    e^b = psi_approx^b - psi_exact^b.
    r^b = L e^b.
    s_i^b = G_i(v):S^b.
    (Q_b v)_i = s_i^b.

At the frozen W=0 state, coefficient errors produce

    delta(Wdot) = sum_b e^b s^b.
    delta(Cv) = sum_b [diag(s^b) r^b + [L,diag(s^b)] e^b].
    delta(F) = sum_b Q_b^T diag(e^b) L P.

The commutator vanishes for a constant affine strain component. Spatially
varying strain retains the coefficient error itself. The paired force and work
remain adjoints even with an inaccurate coefficient field. Conservation alone
therefore does not certify coefficient accuracy.

Fourteen 2D/3D small-native-geometry checks inject lowest-positive and highest
eigenmode errors with the same residual L2 norm of one. Every strain basis is
tested. For the first basis the results are:

| Dimension | Mode | Coefficient-error L2 | Volume-rate-error L2 | Force-error L2 |
|---|---|---:|---:|---:|
| 2 | lowest positive | 2.03320 | 0.39690 | 0.99653 |
| 2 | highest | 0.26400 | 0.03588 | 0.11616 |
| 3 | lowest positive | 12.10015 | 0.48843 | 0.98285 |
| 3 | highest | 1.13535 | 0.03057 | 0.08377 |

The experiments demonstrate different physical responses at equal residual.
The measurements do not prove low-mode dominance for every possible state.
The reference force/work and affine-constant-strain identities pass.

## Analytical upper estimates without another elliptic solve

The positive reciprocal native conductances are averaged for an undirected
analysis graph. No face is discarded. The maximum matrix-entry change from the
raw runtime graph is 2.66e-15. Runtime conductances are unchanged. Both approximate
and reference coefficients are compared to the same compatible mean-zero target.
Reference uncertainty is retained separately.

One deterministic spanning tree routes each residual into a feasible edge flux.
For a mean-zero coefficient error, the residual-flow energy gives

    e^T L e <= sum_tree flux_edge^2 / conductance_edge = E_upper.
    |e_i|^2 <= E_upper [R_path(i,root) + mean_j R_path(j,root)].

The pointwise estimates propagate through the commutator, transpose force, and
relative weight-rate difference across each face. Calculations use long double
with 63 mantissa bits and monitor residual closure and tree-flux divergence.
The tree divergence defects are below 1.2e-28 for the native coefficient errors.
The estimates are analytical upper estimates evaluated with controlled roundoff,
not interval-arithmetic certificates. Ritz estimates are not used as certified
spectral-gap bounds. Conservative estimates can be too loose to select a useful
physical tolerance.

## Frozen native physical-error scales

The original checkpoint verifies equal masses 1/16384, W=0, and matching
positions, volumes, pressure, and velocity. The captured common fit reproduces
the independent native fitted divergence to 5.3e-13. Actual and exact nodal
Gresho velocity fields are tested separately. Exact nodal pressure is a generator
sample rather than a cell-average projection.

The original native step at the same checkpoint supplies the reporting scale
dt=0.00070072076050564647. The new candidate has no qualified timestep bound.
First-order frozen one-step scales are defined as

    |delta(U)|/U ~ (gamma-1) dt |delta(Cv)|/V.
    |delta(v)| ~ dt |delta(F)|/m.
    |delta(face_normal_position)| ~ dt |delta(Wdot_i)-delta(Wdot_j)|/(2 d_ij).

For actual velocity and pressure, the observed coefficient differences and
the approximate-solution upper estimates give:

| Quantity | Observed weighted mean | Observed maximum | Upper estimated mean | Upper estimated maximum |
|---|---:|---:|---:|---:|
| relative thermal increment | 3.87e-8 | 2.11e-5 | 4.41e-5 | 1.55e-2 |
| velocity increment | 2.41e-7 | 1.07e-4 | 3.72e-4 | 9.43e-2 |
| face displacement / h | 1.73e-8 | 8.51e-6 | 1.67e-5 | 6.69e-3 |

The thermal and velocity means are volume weighted. The face means are area
weighted and h=1/128. Observations use the sparse reference; reference uncertainty
has separate bounds in the archived results. The observed volume-rate L1 error
is 8.2790e-5, while the residual-only term is 3.5889e-7 and the commutator term is
8.2754e-5. The commutator dominates the residual-only estimate by about 231. The
decomposition defect is below 1.2e-21. Force/work differences remain below 1e-22
despite the nonzero operator errors.

Exact nodal Gresho velocity and pressure give smaller coefficient-error scales.
The mean relative-volume-rate error is 1.3965e-6 and the mean acceleration error
is 5.2874e-7. These frozen results do not establish a long-time PDE error bound.
No numerical scientific error budget has been declared. No new solve tolerance
is selected from the observed residuals. The loose tree bounds mean that a
tight tolerance is not certified; loose bounds do not establish large actual
errors.

## Spatial-response warning independent of coefficient accuracy

Using the accurate reference coefficients reveals a separate issue. On the
same saved geometry, the candidate reduces the relative volume-rate L1 for
exact nodal Gresho velocity from 0.54168 to 0.06168. The candidate increases the
response for the actually evolved velocity from 1.34192 to 6.54532.

Linear decomposition separates the deviation from the nodal reference:

| Frozen input deviation | Original B operator | Reference candidate C operator |
|---|---:|---:|
| velocity deviation, relative volume-rate L1 | 1.68097 | 6.55577 |
| pressure deviation, force-density L1 | 4.92320 | 21.88994 |
| exact nodal pressure, force-density error from analytic gradient | 0.01569 | 0.03244 |

The norm components need not add. The original state was evolved with B rather
than C. The measured amplification is a frozen response and is not proof of a
candidate instability. The amplification persists with accurate coefficients,
so tighter coefficient solves alone cannot remove the measured response.
The next decision must assess the spatial operator before investing in solver
tuning or a long candidate trajectory.

The bundled Astra evaluation selects a finite linear-mode and frozen-native
response stage. The checks must distinguish known longitudinal compression
from transverse, divergence-free motion. Larger response than B alone is not
a rejection condition. The next stage reuses the reference coefficients and
records pressure null modes, long-wave sound speed, false transverse
compression, and positive-acoustic stiffness. No PDE, AV change, new global
reference solve, or coefficient-solver tuning belongs to the next bundle.

## Provenance and open scope

Run root is `/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/material_volume_error_propagation_v1`.
Native capture occurred once. Initial analysis is preserved. The authoritative
analysis is `analysis_v2/propagation.json`, with canonical graph edges used
consistently in the commutator and bounds and an added exact/deviation response
decomposition. No native capture, reference solve, or PDE run was repeated for
the analysis revision. Source snapshots distinguish both analysis versions.

Compact evidence is archived under `Research/CoefficientError_20261105`.
Large arrays and native input remain on GPFS with SHA256 hashes. The manuscript
keeps the previously verified trajectory. New-method vortex/shock accuracy,
weighted time integration, dimensional PDE regressions, the remaining 2D suite,
and matched GIZMO throughput are still open.
