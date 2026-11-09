# Acoustic contact response and cold-shear rejection, 2026-11-05

Status: **experimental, not adopted; the user goal remains incomplete**. This work uses the existing native C/OpenMP geometry and common 1D/2D/3D finite-volume operator. The production OpenMP/MPI/CUDA solver is preserved. No simulation output was overwritten, no trial checkpoint was used as a restart, and no work used `/scratch`.

## What has been established

The earlier stationary-contact initialization defect and the shock-generated contact oscillation are different problems. The former is documented in CONTACT_1D_CAUSAL_AUDIT_2026-11-04.md. For the latter, the cell-mean KEP core has an exact alternating-pressure null mode on regular periodic meshes: nonuniform alternating pressure has zero acceleration and zero pressure rate in 1D, 2D and 3D. Increasing a compression stress cannot act on this zero-velocity state. HLLC damps it while retaining a true constant-pressure density contact.

A new complete acoustic correction also damps this pressure mode in all three dimensions. On the same original 1D target times, it strongly reduces pressure oscillations near the shock-generated contacts. This supports missing acoustic response as a mechanism sustaining the previous oscillation. It is not a proof of every source of the oscillation: flux reconstruction and boundary treatment also differ, and density errors do not improve uniformly.

## Why the first cell-mean acoustic replacement was rejected

Read-only audit: `/gpfs/kjhan/LagForce/UnifiedHydro_20261104_v1/acoustic_exact_shear_audit_20261105_v1`.

The saved physical9 cold64 geometry was supplied with analytic cell means of u=(sin(2 pi y),0), rho=1 and P=1e-6. The original saved state was not modified. Exact cell integration uses native faces, the divergence theorem and a sinc integral. The inherited gamma is the original single-precision representation of 5/3.

At constant rho and P, the unreconstructed rank-two acoustic correction has

D_m = 0,

D_p = -rho lambda delta(u_n) n / 2,

D_E = ubar dot D_p.

Its global discrete cell-mean kinetic-to-thermal transfer is sum_face area rho lambda delta(u_n)^2 / 2. This is positive for smooth divergence-free shear on oblique faces. On the saved 4096-cell mesh it gives heat rate 5.9933695e-4, or 399.55794 initial thermal energies per unit time. This is an instantaneous diagnostic, not a fitted trajectory. The same analytic shear on a regular grid is effectively zero; translation is exactly zero. The independent cell and face energy sums agree to 2.1e-18, and the analytic formula agrees to 4.5e-18.

The original diagnostic's absolute mass-correction gate narrowly fails: 1.04095e-12 versus the registered 1e-12 threshold. This cancellation residue and failed result are preserved, not hidden by relaxing the threshold. It is far too small to explain the heat, which agrees with the formula having exactly zero mass correction. A crosscheck on this same geometry reproduces the previous filtered-stress heat rate 1.3735108e-8.

## Reconstruction diagnostics and experimental common operator

The EC core follows [Chandrashekar](https://arxiv.org/abs/1209.4994). Sign-compatible entropy reconstruction is motivated by [Fjordholm, Mishra and Tadmor](https://epubs.siam.org/doi/abs/10.1137/110836961). The trace-local acoustic amplitude construction below is an experimental construction here; those references do not qualify this particular solver.

Literal projection of reconstructed entropy variables onto the cell-mean acoustic basis creates unwanted mass transfer from the squared tangential velocity. The native frozen-geometry diagnostics are saved in `acoustic_reconstruction_audit_20261105_v1`. They retain raw and limited gradients, with and without a sign constraint. The experiment therefore reconstructs the two acoustic amplitudes in the trace pair's own mean frame, retaining the cell-mean basis for the conservative correction.

For beta=rho/(2P), beta_log the logarithmic mean, a^2=gamma/(2 beta_log), and the two cell means, the exact entropy projections simplify to

A = delta(log P) - delta(beta) |delta u|^2 / 4,

B = a (beta_L+beta_R) delta(u_n),

z_minus = A-B,

z_plus = A+B.

This form eliminates the large cancelling kinetic terms in a literal entropy-variable dot product. Compute reconstructed amplitudes in the same way from the two face traces. Bound each reconstructed amplitude to the signed interval between zero and its cell-mean amplitude. This is a numerical sign constraint with no fitted coefficient. For c=rho_log lambda/(2 gamma), the complete correction is

D = -c [r_minus zrec_minus + r_plus zrec_plus] / 2.

The entropy contraction is exactly -c (zmean_minus zrec_minus + zmean_plus zrec_plus)/2 <= 0. The correction is shared by mass, momentum and energy. It is independent of actual face velocity, so the existing affine mass-flux derivative and scalar screened weight solve remain unchanged. The mass correction is included in that solve, not merely added after it.

This replaces the compression/LSQ-residual stress; the filter and centroid-strain macros are disabled. No new evolved physical variable or adjustable coefficient is introduced. Physical node velocities and the SSPRK2 harness remain as before. Reflecting boundaries use mirrored cell-mean Riemann states, rather than forcing both incident hydrodynamic traces to zero normal velocity before solving the wall problem. Boundary reconstruction is explicitly first order. This boundary treatment is separate from keeping the physical wall face stationary.

Source snapshots, compiler commands and hashes:

- `acoustic_reconstruction_stage_v1`: first experimental physical10 implementation and complete first attempt.
- `acoustic_reconstruction_stage_v2`: same flux and reconstruction, numerical time-step correction only.
- `acoustic_reconstruction_stage_v2/whole_physical_path.json`: accepted checkpoint segment provenance.

The frozen-face audit has 12000 state pairs across dimensions. It verifies original-flux recovery when traces equal means (scaled difference 5.25e-14), entropy contraction (9.60e-16), orientation reversal (2.78e-15), Galilean transformation (4.60e-14), and exactly zero constant-pressure/constant-velocity contact correction. This is not a proof of positivity or fully discrete entropy stability.

The actual coupled RHS damps the alternating pressure amplitude at approximately -0.19025847 in every dimension for initial amplitude 0.01. The constant-pressure density checkerboard remains exactly stationary. These are initial semidiscrete rates, not measured full-evolution decay constants.

## Time-step implementation error and correction

The first implementation incorrectly estimated a viscous restriction by dividing total acoustic traction by the normal velocity jump. Acoustic traction includes pressure-driven acceleration; this quotient is not a valid diffusion coefficient and becomes singular as the velocity jump vanishes. At the saved cold failure, the quotient restriction was 3.41e-16 while the acoustic CFL permitted 0.210 and the admissible-state interval permitted 0.0968. The corresponding normal jump was 6.21e-31.

The second implementation removes that quotient and includes the actual acoustic dissipation speed in the wave CFL. The conservative flux and its thermal/entropy transfers are unchanged. The three interrupted cases continue from their accepted physical10 checkpoints; they do not restart from birth or reset their thermal history. The diffs and numerical_change.json preserve this correction. The error is separate from the physically unacceptable shear heating, which persists after correction.

## Results at original target times

All nine original 1D cases were attempted with the same physical10 operator. Shu uses N800 and the other eight use N200. Eight reach their original target times, including double rarefaction to 0.15 and Noh to 0.6. Blast stops at 0.0083890513 before its 0.038 target, after an unaccepted trial's screened solve reaches the PCG cap with true scaled residual 2.08e-12 versus its unchanged 1e-12 requirement. Its preceding extremely small steps and low-density cell also require diagnosis; this is not evidence of a qualified strong-shock solution.

Stationary and moving 1D contacts and the 2D/3D stationary contact controls preserve pressure/velocity to roundoff. This does not certify multidimensional shocks.

Contact-window comparison uses actual cell overlap with x_contact +/- 0.03, exact constant star pressure and matched target times:

| Test | Previous pressure relative L1 | Acoustic pressure relative L1 | Previous density cell-average L1 | Acoustic density cell-average L1 |
|---|---:|---:|---:|---:|
| Sod | 1.4931% | 0.01298% | 2.2188% | 2.3786% |
| Lax | 0.9904% | 0.005071% | 1.6219% | 1.0014% |
| Collision | 0.3572% | 0.001301% | 1.4986% | 2.3462% |

The pressure improvement is substantial. **Contact accuracy as a whole is not solved:** collision has a 6.84% density overshoot above the exact high-density state. Density field L1 is distinct from the cell-average norm and is also retained. Scripts and metrics: `acoustic_reconstruction_stage_v1/compare_contacts.py`, `contact_comparison.json`, `contact_density_accuracy.py`, `contact_density_accuracy.json`. Figure: `contact_comparison.png` and PDF, with old/new density and pressure shown together.

Noh at 0.6 has domain-integrated relative field L1 approximately 0.997% density, 0.908% pressure and 0.191% velocity against the finite-pressure mirrored Riemann solution. However, its maximum density is 8.358 versus the exact postshock 3.99999, and maximum pressure is 1.63081 versus 1.33334. The wall-cell peak remains unacceptable as a claim of uniformly accurate wall capture. Whole-domain averages must not conceal it.

The cold64 shear reaches its original 0.25 target after the numerical correction, but pressure relative L1 is **17.9693%**, failing the registered 1% gate. Thermal energy increases by 2.69540e-7. The integrated acoustic thermal transfer is 2.69622e-7, accounting for 100.030% of the increase; summed absolute time-integration kinetic defect is only 8.19e-11. The failure is dominated by the explicit acoustic dissipation, not the time-step error, a floor, or global energy nonconservation. No floors or repairs were used.

## Geometric reconstruction followup and next decision

`acoustic_midpoint_reconstruction_audit_20261105_v1` evaluates a centroid-midpoint dissipative jump with the same gradients, without changing the EC core or introducing extra derivatives. On the exact saved shear it reduces the instantaneous heat only from 1.7794e-6 to 1.5461e-6 with limited gradients (raw: 1.7621e-6 to 1.5284e-6). It is therefore insufficient evidence for another full solver branch; no midpoint evolution was launched.

The next design question is the work pairing of the dissipative reconstruction. The current face force depends on a reconstructed residual, but its global kinetic work multiplies the full cell-mean velocity jump. It therefore need not scale as the square of the unresolved residual. A mathematically compatible adjoint residual operator is a candidate to investigate: on a frozen stencil, a force of the form -R^T D R u has nonpositive work -|D^(1/2) R u|^2. Conservation, the constant/affine null spaces, finite-temperature entropy, pressure-equilibrium preservation, boundary treatment and actual cost all require verification before adoption. This is a research direction, not an implemented cure or permission to drop any gate.

The pressure response should be retained as a demonstrated requirement. The next change must improve density and cold shear together, not merely suppress the plotted pressure curve. Do not add fitted sensors or tune separate coefficients by test. Do not relabel the present candidate as accepted, restart failed trials, weaken solver/physics tolerances, or substitute short diagnostic times for the requested endpoints.

## Remaining user requirements

- A single accepted contact/shock/shear/wall operator: **not achieved**.
- 1D contact pressure mechanism: supported by null-mode and evolved intervention evidence; density and full causal decomposition still incomplete.
- 1D full suite: 8/9 endpoints only; wall peaks and density errors remain.
- Common 1D/2D/3D code: exercised for RHS/contact controls, **not** full multidimensional shock accuracy.
- Supersonic upstream/subsonic downstream RH shocks and hard walls in all dimensions: previous HLLC evidence cannot be assigned to this new candidate; still required.
- End-to-end <=2x matched GIZMO runtime: **not established**. No matched throughput claim is made. The native 3D N8 control still enumerates up to 511 geometry candidates; this experimental path is not a scalable substitute for the production parallel code.

No test or helper process from this batch remains running at the time of this record. All accepted states, failure states, commands and source snapshots are retained under `/gpfs`.
