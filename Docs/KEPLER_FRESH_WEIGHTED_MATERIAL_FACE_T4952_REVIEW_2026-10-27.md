# Fresh weighted Kepler power faces remain nonmaterial at t≈4.952

The [kinematic face audit](KEPLER_FRESH_WEIGHTED_MATERIAL_FACE_T4952_2026-10-27.json)
uses the same accepted fresh-weighted one-step native snapshots and hashes
as the [local-work study](KEPLER_Q2_WEIGHTED_FRESH_LOCALWORK_T4952_REVIEW_2026-10-27.md).
It reconstructs the midpoint generators and power diagram independently,
uses the secant `Wdot=(W_end-W_start)/0.001`, and compares each disk-disk
face's geometric normal speed with the analytic softened Kepler velocity
at its midpoint. There are 145,671 disk-disk faces, including 49,038
same-initial-ring and 96,633 cross-ring faces. These are geometric face
speeds, not measured native Riemann mass fluxes.

| Same generator positions and velocities | Integrated `ΣL|v_face−v_K,n|` | Length-weighted slip RMS |
|---|---:|---:|
| W=0 Voronoi | `0.192565` | `2.01376e-4` |
| Actual power weights held fixed | `0.194148` | `2.02936e-4` |
| Actual power weights and measured `Wdot` | `0.174979` | `1.81059e-4` |

Moving weights reduce integrated absolute face slip by **9.13%** and RMS
slip by **10.09%** relative to W=0 at these same generator states. They
leave about 90% of the original slip, even while almost freezing cell
volumes. Same-ring integrated slip falls from `0.114857` to `0.100999`
(12.1%); cross-ring slip falls from `0.077708` to `0.073980` (4.80%).
This directly distinguishes a volume constraint from the stronger condition
that every face move with the fluid.

The separate [endpoint total-PdV check](KEPLER_WEIGHTED_ONE_STEP_PDV_CLOSURE_2026-10-27.json)
finds local internal-energy/PdV consistency to p99 `8.36e-8` relative
over the same smooth step. Thus a nonmaterial face does not mean a missing
term in the code's current adiabatic cell-energy update. It does mean that
fixed mass attached to these moving power cells cannot be interpreted as
exact material mass on every face. The analytic Kepler face speed and
midpoint/secant reconstruction are proxies; stage-exact native transport
and general-flow convergence remain unverified. The result argues against
relying on **weight rates alone** to make the tessellation fully material
through the historical t>6 interval.
