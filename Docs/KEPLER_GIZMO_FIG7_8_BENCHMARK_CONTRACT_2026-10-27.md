# GIZMO Kepler benchmark contract for the LagForce campaign

The primary source is [Hopkins, *A New Class of Accurate, Mesh-Free
Hydrodynamic Simulation Methods*, §4.2.4 and Figures 6–8](https://arxiv.org/pdf/1409.7395).
Figure 6 shows surface density at `t=120`; Figure 7 explicitly shows the
MFM density field at **`t=600`**. The Figure 7 caption says the inner
`r≈0.5` disk has executed more than 250 orbits. Since a Kepler period at
`r=.5` is about 2.22 in these units, `t=600` corresponds to about 270
inner orbits. The prose elsewhere gives a different approximate orbit
count; the fixed physical time `t=600` is the reproducible endpoint.

The paper initializes a 2D `γ=5/3` disk with

```
Σ(r) = 0.01 + (r/0.5)^3                 for r < 0.5
     = 0.01 + 1                         for 0.5 ≤ r ≤ 2
     = 0.01 + [1+(r-2)/0.1]^-3         for r > 2,
P = 10^-6,  Φ = -(r²+ε²)^-1/2,
v_c(r) = r(r²+ε²)^-3/4.
```

The authors describe effective `256²` resolution: initial points are
evenly spaced over `−2<x,y<2`, then particles are removed inside
`r<0.5` and outside `r>2` to match the taper profile. Our current
bounded q2 control has a comparable
49,038 disk sites, but uses deliberately polar ring placement, a sparse
core/exterior, a fixed polygonal wall, and `ε=0.01175`. It is a useful
mechanism test, **not the same IC or boundary problem**. The three-IC
aligned/ring-coherent/jitter comparison measures sensitivity to this
layout, not reproduction of the paper's setup.

The [accepted bounded-q2 glass spacing audit](KEPLER_Q2_BOUNDED_RING_SPACING_2026-10-27.json)
uses the original SHA-256-verified positions, not a nominal spacing formula.
Across its 68 disk rings, the radial gap is `0.02209708` and the same-ring
neighbor chord is about `0.01088`: their ratio is `2.029–2.034`. For all
49,038 disk generators, the nearest generator on an adjacent disk ring is
farther than a same-ring neighbor (sitewise cross/same ratio `2.029–2.093`).
Thus the conjecture that neighboring rings have closer generator centers
than a given ring's neighbors does not apply to this q2 IC. This distance
ordering does not eliminate cross-ring Voronoi faces or their shear-induced
geometry changes. The [initial face-motion audit](KEPLER_FACE_MOTION_GAP_INTERPRETATION_2026-10-27.md)
found 96,593 disk-disk cross-ring faces and 49,038 same-ring faces; their
integrated material-versus-geometric normal-speed gaps are `0.03617` and
`0.06934`, respectively. Neither distance nor face-speed gap alone measures
irreversible heating.

For Figure 8, the paper defines the arithmetic particle means in the
**current-position** annulus `0.8<R<1.2`:

```
L1(v) = mean |v_phi,i − v_K(R_i)| / v_K(R_i)
L1(ρ) = mean |Σ_i − 1|.
```

The fixed-time velocity comparison is at `t=20π` (10 orbits at `R=1`).
The paper also charts first crossing of `L1(v)=0.01` and `L1(ρ)=0.25`.
These are **different from our mass-weighted full-disk density and
exact-orbit-position gates**. Figure 7 has no single published scalar
pass threshold; its density image, inner/outer ring retention, orbit and
angular momentum diagnostics must be shown together at `t=600`.

The [hash-verified metric baseline](KEPLER_GIZMO_FIG8_METRIC_BASELINE_2026-10-27.json)
applies the published annular formulas to available native snapshots and
also reports an initial-ID annulus and `mean |Σ/Σ₀−1|`. For the independent
jitter IC at `t=0.825`, paper-style `L1(v)` is `9.13e-8` W=0 versus
`4.61e-11` fresh weighted. The literal `L1(ρ)` is near `0.01` for both,
because the IC itself has `Σ=1.01` in this annulus: while all density
values remain above one, `mean |Σ−1|` mostly measures the constant
`0.01` floor and can hide signed small-scale noise. The own-t0 relative
metric reveals `1.460e-3` W=0 versus `1.285e-5` weighted. This additional
metric supplements rather than silently replaces the paper's definition.

At `t=30`, the aligned q2 case has annular paper-style `L1(v)=1.82e-8`
but a much larger **inner** `0.50≤r₀<0.65` mean orbit error
`2.42e-5`; the paper's `0.8<R<1.2` annulus alone can miss early damage
to the inner ring. Future release gates should retain both radial and
full-disk measures. Neither `t=30` nor the current weighted `t≈0.825`
state tests Figure 6 at `t=120` or Figure 7 at `t=600`.

## Full-q2 W=0 early metric sample (2026-10-27)

The paper-style annular metrics were also evaluated read-only on the complete
384² W=0 q2 checkpoint at `t=2.0000159421`. In `0.8<R<1.2`, both the current
position and initial-ID annuli contain 10,453 sites. The result is
`L1(v)=7.08e-9`, literal `L1(Σ)=0.01000004`, and own-initial-density-relative
`L1(Σ/Σ₀)=1.745e-4`. The literal density value is almost entirely the
prescribed `0.01` floor. At `t=2`, the run has completed only `0.318` orbit at
`R=1`; passing the paper's `L1(v)=0.01` threshold this early is not evidence
of long-time stability. This remains a mechanism/control IC with polar ring
placement and a periodic halo, not a direct reproduction of GIZMO's uniform
effective-256² particle layout. Input and endpoint hashes and the metric
definitions are recorded in
[the machine-readable result](KEPLER_GIZMO_FIG8_Q2_FULL_W0_T2_METRICS_2026-10-27.json).

## Full-q2 exact `t=2.221` gate (2026-10-28)

The completed 384² full-domain periodic q2 W=0 checkpoint
`khout.030070.dat` has a passing 15-check global receipt and a passing
five-check radial receipt. The global receipt verifies exact fixed mass and
W=0, positive state, `rho=m/V`, periodic-area tiling, zero floor/dual-energy
injection, and the paper-style accuracy thresholds. In the current-position
`0.8<R<1.2` annulus, the published-formula diagnostics are
`L1(v)=7.8992e-9`, literal `L1(Σ)=0.01000004`, and plateau-corrected
`L1(Σ)=1.6166e-4`. The literal density norm is dominated by the specified
`0.01` density floor. The initial-disk mass-weighted density change is
`1.4081e-4`, entropy-log change `3.5605e-4`, orbit-position error
`2.4628e-8`, angular-momentum drift `-6.3937e-13`, and total-energy drift
`-3.7761e-10`.

At `epsilon=0.01175`, `t=2.221` is `0.9994` orbit at `r=0.5` and `0.3534`
orbit at `R=1`. The Fig. 7 reference is `t=600`, whose caption reports more
than 250 inner orbits at `r≈0.5` (about 270 with this softened period).
Therefore the excellent q2 checkpoint is an early equilibrium check; it
does not approach the paper's long-time result.

The radial receipt also shows where the current periodic setup departs from
that inner-disk result. Initial-radius bins `2.35≤r₀<2.75` and
`2.75≤r₀<5` have mass-weighted density changes `0.1747` and `0.3343`;
the latter has absolute radial-velocity change `0.1826`. In the current
outer-edge selection, the angular-momentum change is `-0.6030`. These bins
include the low-density exterior and periodic seam; they do not by themselves
show that the scientific `0.5≤r₀<2` disk has failed. They make the matched
radial diagnostics and density fields at the already queued `t=4.5` and
`t=6.5` endpoints essential. The polar 384², 6×6 periodic setup remains a
distinct global initial-value problem from the paper's effective-256²
particle setup and must retain that boundary-condition caveat.

Machine-readable evidence:
[global t=2.221 receipt](/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/halo12_after_failure_to0825/to2p221_h28/audit_t2p221.json),
[radial t=2.221 receipt](/gpfs/kjhan/LagForce/Kepler_PolarIC384_H200_20261026/q2/halo12_after_failure_to0825/to2p221_h28/audit_radial_t2p221.json),
and [q2 density image](KEPLER_Q2_H200_T2.221_DENSITY_2026-10-27.png).
