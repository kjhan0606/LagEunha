# Native ALE shadow: mass, work and affine moment evidence

2026-10-20. This is a completed read-only operator bundle. No ALE mass flux,
spin energy, moment correction or new physical timestep was applied to the
evolving state. The long Kepler goal remains incomplete.

## Native implementation and observation check

The confirmed source is `/gpfs/kjhan/LagForce/Kepler_stencil5x5_Nx256_A100_call33720_20261020/src`.
The driver added `Exam/gfs_ale_flux.h` and `Exam/gfs_ale_shadow.inc` to the user's
CUDA/MPI body. `GFS_ALE_SHADOW_CALL=1` captures one RHS. `GFS_ALE_SHADOW_ONLY=1`
in the opt-in SSP path exits collectively after current-stage geometry -> EOS
-> halo -> gradient/force preparation, before any accepted evolution step.
Defaults remain off.

One rank, owning the lower global ID, constructs each shared pair packet.
One CUDA thread per pair computes a full conservative Euler HLLC flux with
Davis bounds; an inadmissible HLLC state uses the complete HLL flux. No
pressure-ratio100 switch is applied to this shadow solver. The same moving
face frame and two Gauss points produce mass, two momenta and total energy.
Native reconstruction clamps are retained and counted, not hidden.

The packet records full flux, pressure work, relative advection and an explicit
Riemann residual. In an HLL fallback, its numerical flux need not equal the
physical flux of its HLL-averaged state; that difference is retained. Differences
from the production operator also contain changed Riemann closure and face
quadrature, so they are not attributed entirely to omitted advection.

Exam and KH libraries were explicitly built on user-authorized lageunha.
Frozen binary SHA256:
`fbdb5d2bea0354b30c34d82ee60198c101f51ece1d5b2d22a0888edb53217d48`.
Build sources/hashes/logs are in `ale_shadow_v1_build/`.

As the single observation check, the already specified three-step fine interval
was replayed with only the shadow capture added. Every recorded field and all
metadata match the prior SSP v2 endpoint bitwise. It used5.35seconds; no new
physical parameter case was introduced. See `ale_shadow_observer_v1/manifest.json`.
All other captures stop after the first RHS with no time evolution. Original
initial/checkpoint files retain their SHA256 hashes.

## Distorted 64^2 patches: exact integration and mass-volume contract

Five captures: uniform, rigid rotation, trace-free linear shear, uniform-boosted
shear, and the same shear under MPI4. The position jitter is22% of the spacing;
all cases use the same periodic geometry and initial polygon-integrated masses.
Each has4096cells and12288 unique canonical pairs. Native execution took about
2.4seconds per MPI1 capture and4.68seconds for MPI4, including launch and I/O.
These times are not a scaling benchmark.

The manufactured exact affine reconstruction is stored separately from the
actual native slope-limited reconstruction. The central1024cells are away
from the nonperiodic affine velocity discontinuity at the periodic seam.

- Two-point CUDA quadrature agrees with an independent three-point analytic
  polynomial integration to ~1.5e-15 of the natural flux scale.
- Shear reference max|mdot-rho*Vdot|=1.13570e-17, relative to its peak mass
  rate1.09209e-13.
- Shear max|Vdot|=1.03994e-4. With rho=1 and fixed mass, its instantaneous
  max|rho_dot| would be .410738 per unit time. The full mass flux cancels
  that geometric density change in the uniform-density reference.
- Actual native reconstructed shear mass-GCL error is6.92398e-13. Uniform
  rho is exactly1, while EOS float P ranges9.99999884e-7 to1.00000011e-6;
  actual uniform mass-GCL residual3.54210e-12 therefore is not reported as
  exact machine-zero uniform preservation.
- Adding velocity(3,-2) gives full mass/momentum/energy boost-transform errors
  <=5.07e-15 relative to peak flux. MPI1 versus MPI4 canonical full fluxes
  are bitwise identical.
- Pair antisymmetry is exact by the single-packet construction. Global cell
  sums have ordinary floating accumulation residuals, which are recorded.

The first analysis assertion used a peak-relative normalization even for a
vanishing uniform mass flux. It failed on ~2.75e-19 absolute roundoff. Only
the analysis was corrected: the gate now uses1e-9*peak +128*machine-epsilon*
natural-flux scale and retains both absolute and relative values. No native
patch was rerun or changed to make that analysis pass.

**Meaning:** at this finite mesh, a geometrically changing cell in an exactly
uniform-density incompressible affine flow requires nonzero mass transfer.
Keeping each mass fixed while treating the same actual polygon as the exact
Euler control volume cannot reproduce that reference exactly. This is not a
general impossibility theorem for fixed-mass particle methods or their
continuum convergence.

Data: `ale_shadow_patches_v1/manifest.json` and `analysis.json`.

## Two actual cold states: full ALE is not a demonstrated cure

Both use71080 cells/213240 canonical pairs. Only their first RHS was captured;
no additional physical evolution was performed. Their executions took3.02 and
2.93seconds. Global conservative flux sums are at roundoff, and surface closure
is <=8.34e-16. There are8/16 clamped reconstructed pairs and525/545 HLL-fallback
pairs; these interventions are explicitly visible.

Representative cell68526:

| Rate | Captured normal state t=.10196265 | SSP final normal state t=.13509104 |
|---|---:|---:|
| Existing native Udot | -2.36902e-8 | +3.06994e-9 |
| Hypothetical full ALE Udot | +3.90306e-7 | -8.37755e-7 |
| Hypothetical mdot | +2.08998e-6 | -2.04909e-6 |
| Actual geometric Vdot | +1.19475e-3 | -2.88789e-3 |
| Advection part of Udot difference | +3.96834e-7 | -8.32805e-7 |
| Pressure closure/quadrature part | +1.82277e-8 | -1.47288e-9 |
| Native work-speed part | ~0 | -6.01634e-9 |
| Explicit approximate-Riemann residual | -1.06557e-9 | -5.31512e-10 |

The thermal conversion always includes
`Udot=Edot-v·pdot+.5*v²*mdot`. The four cells show the same qualitative result.
The full-minus-native signed identity closes at~1e-21; canonical-vs-actual
native Udot differences ~1e-15 include float-stored native rates.

The large negative second ALE thermal rate is **not** concealed. A conservative
flux on the existing point-value/internal-U representation is not thereby a
validated thermal or positivity model, and the old fixed-mass state is not
silently reinterpreted as a consistent finite-volume mean state. No mass or
energy was changed using this result.

Data: `ale_shadow_real_v1/analysis.json`.

## Internal affine velocity, angular momentum and pressure work

Additional read-only analysis uses the same packets; no new simulation/capture
was launched. Exact polygon centroids/covariances are integrated from the
actual edges. A separate, unlimited, face-area/distance weighted affine fit
is obtained from the current point velocities. It is **not** the actual
native limited gradient, an existing spin state, or a proof that new energy
must simply be added to the particle model.

The kinetic estimate uses uniform density within each polygon, so its area
centroid/covariance also serve as mass moments in that proposed local model.
They are not the true density-weighted mass moments of an unknown subcell
density reconstruction. In particular, the quoted fitted cell-mean velocity
offset has this uniform-density affine-model assumption.

For cell68526 the affine within-cell kinetic variance estimate is1.22794e-8
and1.36881e-8, versus stored U7.23612e-12 and3.71630e-12 (ratios1697 and3683).
The fitted rotation dominates this scale, but strain and an anisotropic
rotation-strain cross term are separately nonzero. The estimated difference
between generator velocity and fitted cell-mean velocity is(.1599,-.01317)
and(-.1578,.2109). Replacing one interpretation by the other without a
momentum/energy contract is not valid.

For the actual native pressure force, its power against this fitted internal
affine field decomposes as:

| Cell68526 signed power | Captured state | SSP final state |
|---|---:|---:|
| Native Udot, canonical faces | -2.36902e-8 | +3.06994e-9 |
| Rotation pressure power about generator | -2.74967e-8 | -6.51533e-9 |
| Symmetric-strain pressure power | +4.52008e-9 | +4.02221e-10 |
| Remaining work after affine decomposition | -7.13543e-10 | +9.18305e-9 |

Thus at the first state most of the signed native pressure work is compatible
with coherent internal affine motion. This strengthens the case for examining
a rotation/strain moment representation. It is **not** permission to subtract
these values from thermal Udot: dKvariance/dt also requires moment evolution,
moving-origin/centroid terms, gravity and transport. About the true polygon
centroid the rotation/mean-motion allocation changes substantially; at the
second state even the rotation-power sign changes. A generator-torque scalar
alone would therefore leave an unresolved origin/energy contract.

Full advective momentum flux also has noncentral pair force. For any common
face point, pair orbital torque `(xj-xi) cross Fp` is cancelled by the sum of
the two relative-point intrinsic torque rates. The pair identity closes at
~2.12e-22. This algebra requires an intrinsic angular-momentum account if
generator orbital momentum alone is used. Our per-cell face-midpoint torque
allocation is approximate when flux varies along the face; it is not an exact
angular-flux quadrature or a spin-energy evolution law.

Data: `ale_shadow_real_v1/moments.json`, `affine_work.json`,
`angular_balance.json`. Analysis used one CPU thread on lageunha.

## Next decision

This bundle closes the requested read-only contract checks. It does not close
the fixed-mass/material-volume/thermal/angular-moment evolution model. These
results are sent together to Astra to choose the next implementable physical
operator, preserving no steering, Re=infinity, the approved W0 reduction and
the prior weight-DAE/sliding-instability evidence. No graded-IC or long restart
is launched from a diagnostic pass. Final scope remains the t6-10 wall and the
GIZMO-matched longest orbital interval.

### Bundled Astra decision

Astra reviewed the completed bundle and selected an opt-in affine-moment ALE
model together with a gauge-fixed weight–mass relaxation, preserving the old
W0 fixed-individual-mass path. It requires positive density moment realization,
consistent total/affine/thermal energy and true angular flux, rather than a
torque-power correction. The full specification, preconditions and stop gates
are recorded in `ASTRA_AFFINE_MOMENT_PLAN.md`. That new evolving model has not
yet been implemented or claimed stable.
