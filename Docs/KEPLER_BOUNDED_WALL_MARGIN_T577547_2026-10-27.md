# Exact polygon-wall generator clearance through t=5.775473

The [read-only wall-margin auditor](audit_kepler_bounded_wall_margin_20261027.py)
uses the native fixed 1,024-gon: side normal angle
`(k+1/2)2π/1024`, apothem `2.75 cos(π/1024)`, centered at `(3,3)`.
For each accepted generator it computes the minimum clearance over wall
supporting planes and the velocity along the nearest side normal. A direct
evaluation of all 1,024 planes for the ten closest late sites agrees with
the fast angular formula exactly to printed precision. The
[machine-readable receipt](KEPLER_BOUNDED_WALL_MARGIN_T577547_2026-10-27.json)
hashes each snapshot and checks fixed mass, W=0, and positive wall gap.

| Accepted time | Minimum site-to-wall gap | Sites with gap < 0.2 | Largest outward normal speed among those sites |
|---:|---:|---:|---:|
| 0 | 0.10312633 | 72 | 0.00189868 |
| 3.30024072 | 0.10299799 | 72 | 0.00190830 |
| 5.77547349 | 0.10189275 | 72 | 0.00136295 |

At the latest accepted checkpoint, the nearest site is ID 55293, with
gap `0.10189275`; its instantaneous wall-normal velocity is only
`+6.57e-5` although its radial velocity is `+7.00e-4`. The two quantities
differ because a straight polygon side is not radial. All generators
remain inside the fixed polygon, and the minimum gap has decreased by
about `0.00123` from initialization. These observations support the
short-run no-impact gate through `t=5.775`; they do **not** predict the
time of a future collision or show that wall waves cannot affect the
disk at `t=120` or `600`. A tested collision/reflection rule and
time-integrated wall impulse/torque ledger are still required before
promoting the closed-wall branch as a long-time boundary method.
