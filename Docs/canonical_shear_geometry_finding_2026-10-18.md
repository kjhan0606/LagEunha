# Initial shear and geometric volume mismatch

The local circular field is v(x)=Omega(r) J x, with J the 90-degree rotation.
At generator i its Taylor gradient is

G_i = Omega_i J + (Omega'_i/r_i) (J x_i) x_i^T.

Both terms have zero trace. The second term represents differential rotation.
Applying the existing geometric volume row to velocities
v_j = v_i + G_i (x_j-x_i) therefore tests response to a locally divergence-free
affine field. Every row has its own Taylor field; no globally modified velocity
field or new pressure operator is introduced.

The common-face geometry gives the existing row exactly:

Vdot_i = sum_j A_ij [ (v_i+v_j).e_ij/2
                    + (v_i-v_j).(c_ij-d_ij/2)/|d_ij| ].

Here c_ij is the face centroid relative to generator i. Reconstructed actual
Vdot/V agrees with the original C diagnostic within 3.90e-14 for r<1.9.
All inputs are saved t0 outputs; no new simulation or production edit occurred.

At r=.2-.3, mass RMS rates are:

| Field | Vdot/V RMS |
|---|---:|
| Actual circular velocity | .16515613 |
| Local affine approximation | .15679606 |
| Actual minus local affine | .01536817 |

The local rigid-rotation response is below 8.11e-14 throughout r<1.9.
The affine term and remainder are correlated; their RMS values do not sum,
and .1568/.1652 is not a causal percentage of the nonlinear instability.

This refines the earlier h/r explanation: most of the instantaneous signal
in this radial band survives without higher-order velocity curvature. A
Euclidean Voronoi partition does not generally transform as material volumes
under volume-preserving shear. The current volume derivative can therefore be
geometrically correct while not maintaining the material equilibrium.

No direct replacement Vdot=0, pressure suppression, or reduced AV is justified.
Any correction must retain compatible pressure work, mass/momentum/energy
accounting, no steering, and the approved equal-weight fixed-mass reduction.
Prior material-volume DAE secular failure and nonzero constant-pressure forces
of the raw gradient/transpose candidate remain relevant rejection evidence.

Data: `canonical_volume_shear_split_2026-10-18.json`.
Script: `audit_canonical_volume_shear_split.py`.

## Uniform affine-field follow-up

`canonical_global_affine_volume_2026-10-18.json` applies one common velocity
gradient to all generators, using image-correct neighbor coordinates. Results
are assessed only for r<1.9, away from domain boundaries. No boundary-condition
compatibility or evolved affine flow is claimed.

Rotation and isotropic dilation reproduce material divergence to about
1.3e-13 maximum error. Unit simple shear v=(y,0) has zero material divergence,
but the geometric volume-rate mass RMS is .0307867 in inner background,
.00939208 at the disk inner edge, and .00526369 in the disk core. Thus the
effect is not introduced solely by using different Taylor fields in each row.
These are differences from material divergence, not failures to compute the
derivative of the Voronoi volume itself.

Historical constraints on correction remain explicit:

- `dae_periodic_verification_2026-10-09.md`: second-order flux DAE initially
  suppresses growth, but weights grow secularly and suppression degrades near
  t10 even without frozen/arc boundaries. Reintroducing it unchanged is not a
  justified next production candidate.
- `canonical_shared_flux_weight_response_2026-10-17.json`: scalar weight rates
  can fit cell-volume rates at one state while failing to match individual
  material face speeds. Nonzero interior cycle sums certify that distinction
  for the tested reconstruction.
- Consequently, matching divergence or heating alone is insufficient. A
  weight-coupled mass relaxation proposal must also account for the residual
  shared-face transport of momentum and energy, and explicitly prove the
  approved reduction for persistently equal weights. No such coupled cure is
  established by these diagnostics.
