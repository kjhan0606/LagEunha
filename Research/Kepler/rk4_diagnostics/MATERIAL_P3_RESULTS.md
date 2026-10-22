# Fixed P3 material space: bounded initial accuracy gate passed

2026-10-21. Implements `ASTRA_MATERIAL_P3_PLAN.md` after the independent
hybrid/source-curvature bundle. The native initial spatial engineering gate
passes. **No Kepler timestep, curved remap, shock, t6-10 or paper-duration
success is claimed.** The long goal remains active.

## Preserved physical problem and actual implementation

Use the same24 physical polygons and original continuous rho/v/P, potential,
softening and initial cell masses. No generator or Laguerre weight changes.
Replace the internal kinematic P1 space with complete continuous P3 plus the
same two initial-velocity interpolation residuals, now v0-I3(v0). Coefficients
remain free dynamic kinematic unknowns; no gravity/future-orbit function is
added. Parent-cell mesh vertices/edge nodes are shared topologically. Split
each reference triangle into4 for the one h/2 comparison, preserving original
physical boundaries and24 physical cells. This is NOT a new graded-generator
IC or a repeat of the old8x inner-generator refinement.

Native sources: `test_material_p3_hybrid_v1.cpp`, `material_p3_mesh.h`,
`material_p3_trace.h`. New C/CUDA headers in the canonical repository:
`Exam/gfs_material_p3.h`, `Exam/gfs_material_p3_radial_input.h`. The latter
stabilizes evaluation of the ORIGINAL initial velocity residual only. The
user's production native runtime has not been replaced; its25 audited source
hashes still match `enriched_v5_build/source.sha256`.

| Internal representation | P3 | P3, h/2 |
|---|---:|---:|
| reference triangles | 96 | 384 |
| physical cells | 24 | 24 |
| shared polynomial nodes | 497 | 1853 |
| broken scalar DOFs | 648 | 2088 |
| independent physical trace rows | 139 | 223 |
| actual conforming scalar DOFs | 509 | 1865 |
| cell mass condition maximum | 863.03 | 1924.10 |
| Schur condition | 6412.08 | 12241.84 |
| trace Gram condition maximum | 19221.76 | 69251.57 |

## General trace and hybrid operator

Use actual four P3 edge nodal functions and two input-residual traces.
Positive48/96-point face Gram rules, projection against the polynomial trace,
and a residual Gram eigensystem determine the rank. A small candidate mode
is omitted only if its actual sampled function is also at the explicitly
estimated arithmetic noise level; uncertain rank fails. There are no
coordinate-axis incidence exceptions. Rank1 residual faces count8/16;
rank2 faces20/40 at the two representations.

Polynomial node rows use deterministic incidence independence. General
residual constraint rows use load-independent reorthogonalized rank selection.
The chosen canonical row gauge precedes and does not optimize work/torque or
pressure error. The same block-cell saddle equations determine acceleration
and multiplier together. F contains volume stress/gravity and only EXTERNAL
known traction. Actual dual-basis face quadrature independently recovers
B^T lambda, work and angular transfer. Internal triangulation traces are also
checked, not merely physical cell interfaces.

In v4,312/1464 physical-cell polynomial bubbles have no multiplier load.
Their measured virtual-force residuals are <=2.1e-19 of load scale. Direct
traction quadrature on bubbles is <=4.5e-15 of load scale across all controls.

## Important numerical correction, retained failed v1

`material_p3_hybrid_v1` exposed an arithmetic failure in OUTER zero-gravity
controls: direct subtraction of order-one velocities produced residuals as
small as2.5e-13. Values and derivatives lost their compatible precision,
giving chi up to2.26 even with constant physical pressure and no gravity.
This was not accepted as a physical P3 result. The frozen run is preserved.

Evaluate the mathematically identical residual using

    q=v0-T3(v0),       q-I3(q)=v0-I3(v0).

T3 is a degree3 polynomial of the INITIAL velocity, so I3(T3)=T3. No source
or heat is modified. For the radial initial input, write L=2 a.dot(delta)/z,
Q=delta^2/z, s=L+Q, z=a^2+eps^2. A binomial tail of (1+s)^(-3/4) starting
at n=4 evaluates the velocity remainder and its analytic first derivatives
without subtracting order-one values. The tail has explicit geometric bounds,
convergence checks and a validity rejection at abs(s)>=.75. These are basis
evaluation checks, not physical cutoffs.

An independent80-digit DIRECT v0-I3(v0) calculation (not the tail formula)
and direct derivatives validate31 native queries: maximum normalized errors
1.09e-14 and6.25e-14 for P3/h2. See `check_p3_remainder.py` and
`material_p3_hybrid_v3/remainder_oracle.json`. It performs small algebraic
analysis only; all hydro/operator assembly is native C++/MPI.

Initial velocity differs from the original native continuous-input evaluation
by <=8.8e-14, consistent with coordinate/rounding precision. Relative initial
mass and physical U differences against the retained SAME24 native material
reference are <=1.20e-14 and3.69e-14. Maximum cell K difference/U is4.37e-12.
Physical thermal energy was never obtained by subtracting polynomial K from E.

## Curvature convergence and predeclared gate

Use C_h, the error in (D^2 P/Dt^2)/P from the actual reconstructed acceleration
and its gradient. Define chi=abs(C_h)/Omega^2. Include all positive volume
quadrature points plus zero-measure vertex/P3-node/edge diagnostic queries;
these extra samples do not alter any integral. Reported maxima are sampled
engineering diagnostics, not rigorous continuous-space supremum proofs.

Finer64x28 integration and original32x20 integration both pass:

| Actual Kepler diagnostic | P3 | P3, h/2 |
|---|---:|---:|
| max abs C_h, per time^2 | 25.07476 | 2.949230 |
| max cell volume L1 abs C_h, per time^2 | 1.751831 | .2255238 |
| max chi | .005376232 | **.0006291199** |
| max cell volume mean chi | .0003101900 | **.00003960264** |

Final bounds were fixed beforehand: max chi<=.002 and EVERY cell's mean
chi<=.001. Maximum and L1 curvature errors improve8.50x and7.77x. Even the
minimum individual-cell improvements are4.84x and4.79x, so no numerical-floor
exception is needed for the actual Kepler cases. The uniform boost case
passes the same gates. Two levels do not establish asymptotic order.

The old P1 source curvature maximum1216.89/time^2 has therefore fallen to
2.94923/time^2 on this initial fixture. This is not a measured evolution error
or a t10 extrapolation. Initial pressure/source acceleration derivative error
was improved, not just a traction ledger renamed.

## Numerical uncertainty and native parallel/device evidence

- 32x20 to64x28: final max chi changes8.5e-14; max cell mean chi changes
  1.73e-8. Face acceleration changes <=3.1e-12. These are measured integration
  sensitivities, not formal error bounds.
- Read-only saddle residual-correction probe: actual-gravity chi sensitivity
  <=1.08e-14. The probe is NOT applied to acceleration, lambda or thermal
  state. All v4 native state binaries match v3 bit for bit.
- Final-space finite-map check at h=5e-6: normalized curvature difference
  <=2.01e-7, well below the engineering bounds. It is an adiabatic geometric
  derivative check, not an evolved pressure reset.
- Native MPI1/4: all8 state binaries per accepted version bitmatch, as do
  sampled physical-face accelerations and tractions.
- Actual-gravity cell numerical energy residual/U <=4.7e-11/time and angular
  residual/2K <=3e-15. Physical-pressure/source defects remain separately
  reported; their final maxima are about3.23e-5 U/time and2.12e-10 of2K.
- CUDA check on actual31 queries, RTX5000Ada: invalid0; normalized residual
  CPU/GPU difference3.94e-14; polynomial scaled difference2.83e-15. This is
  device arithmetic for the new basis/input helper. No GPU mass/Schur solve,
  production CUDA speedup or complete CUDA hydro evolution is implied.

Accepted frozen runs: `material_p3_hybrid_v2` (32x20), `v3` (64x28), `v4`
(same64x28 equations plus bubble/linear-error probes). The original basis
prerequisite is `material_p3_basis_v1`. CUDA records: `material_p3_cuda_v1`.
Machine evidence: `MATERIAL_P3_FINAL_BUNDLE_20261021.json` and
`MATERIAL_P3_ACCEPTANCE_20261021.json`.

## Cost and current limitation

On lageunha MPI1 the finer reference assembly takes about.43/1.28 seconds
for P3/h2 and trace/Schur construction about.14/1.24 seconds. The fine
quadrature reference has peak rank0 RSS about.97 GiB, including cached points.
These are small-patch reference costs, not a production parallel benchmark.
The dense reference Schur and full point cache must not be used as the
full71080-cell scalable design.

`gfs_material_map.h` STILL has maximum5 basis entries. P3 has12. This initial
reference assembles at F=I and does not call that helper with12 entries.
Before finite P3 material geometry, implement a compatible12-entry C/CUDA map
and symmetric stress-work kernel without changing the preserved5-entry API.
Then connect actual curved map -> Laguerre overlap/remap -> cold re-expression.

The radial helper describes only the INITIAL basis function. After transport
and remap, do not reconstruct a new state by calling the original Kepler IC
again. The transported density, momentum, physical thermal state, resolved
kinetic content and angular moments must survive re-expression. No orbit
reset, steering, heat refill or unresolved kinetic reservoir is authorized.

## Next dependency and goal status

Proceed under the existing Astra plan to the actual nonuniform curved Kepler
material-map/remap and cold re-expression bundle. The prior constant-density
affine overlap test is insufficient. Then composed density-wave/sliding,
shock entropy and native scalable CUDA evolution are required. Historical
t6/t8/t10 and longest GIZMO-paper duration remain the acceptance target in
`../GIZMO_COMPLETION_CONTRACT.md`. The goal is active; no new time evolution
was started and all native sessions in this bundle are terminal.
