# Kepler comparison units and error measures: saved-input audit

2026-10-22. The full t>6 and longest-GIZMO goal remains unachieved. This audit
does not advance a trajectory, select a shorter endpoint, modify an initial
condition, or release long production. The existing positive-rule producer
and H200-first paired-step dispatcher remain unchanged.

## Primary paper evidence

[Hopkins, arXiv:1409.7395, section 4.2.4 and Figures 6–8](https://arxiv.org/pdf/1409.7395)
was read again. Figure times are 120 and 600. The text additionally reports
approximately 1000 orbits for angular-momentum conservation, without a uniquely
specified reference radius in that statement; no 1000-orbit density image was
identified in this section. Figure 8 uses the annulus 0.8<R<1.2, the mean
absolute fractional azimuthal-speed error, and mean absolute density departure
from unity. Its failure thresholds are respectively 0.01 and 0.25. The velocity
convergence comparison is at ten periods at R=1. Equation 34 has a 0.01 density
floor plus a plateau contribution of 1, although the density-error definition
uses unity. Preserve that convention when reproducing the figure.

The main-text estimate of 160 inner orbits at t=600 conflicts with Figure 7's
count above 250 near R=0.5. With the caption's unsoftened period formula, the
calculated counts are 95.493 at R=1 and 270.095 at R=0.5. These are arithmetic
checks, not a determination of the paper's actual numeric softening.

## Actual retained original-code configuration

The frozen native driver sets rin=2, rout=8 and a box of side 24. Its reference
field has

    rho0(r) = 0.001 + 2/(r+1e-6)
                       * sigmoid(8*(r-2)) * sigmoid(8*(8-r)),
    Omega(r) = (r*r + eps*eps)^(-3/4), GM=1.

The saved generating-law binary header independently gives eps=0.047,
P0=1e-6 and gamma=1.6666669845581055. The latter is the retained native
value, not an instruction to alter an evolved checkpoint to exact 5/3.
The paper's plateau profile and the current tapered 1/r profile differ.
Aligning their nominal inner/outer radii by a factor of four does not make
the initial conditions equivalent. With GM fixed and all lengths scaled
together by four, Keplerian times would scale by eight; such a complete
similarity of these actual configurations has not been established.

For the actual current potential,

    T(r) = 2*pi*(r*r+0.047^2)^(3/4), N(t,r)=t/T(r).

| Reference radius | Current-code period | Current-code duration of 1000 periods |
| --- | ---: | ---: |
| 2, nominal inner disk edge | 17.7788919909 | 17778.8919909 |
| 4, scaled counterpart of paper R=1 | 50.2706872012 | 50270.6872012 |
| 8, nominal outer disk edge | 142.1759343783 | 142175.9343783 |

At the current inner reference radius, historical t=6/8/10 are only
0.33748/0.44997/0.56246 periods. They remain essential historical failure
checkpoints, and passing them would not establish even one inner orbit.
Neither of the 1000-period numbers above is selected as the final endpoint
by this audit. The longest comparison still requires a declared reference
radius and matched-setup evidence.

## Consequences for later result analysis

For the historical configuration, report departures from its own rho0(r),
circular-speed law and original accepted history. In the nominally scaled
paper annulus 3.2<r<4.8, rho0 is approximately 0.626 to 0.418. Applying the
paper's departure-from-unity threshold directly would already declare a large
error in this correct initial state. Such a result would measure an IC
difference rather than accumulated numerical error.

For a matched paper reproduction, retain the paper's annulus, definitions,
times and thresholds. Distinguish the paper-style sample mean from mass- or
volume-weighted diagnostics. Report radial velocity, phase/orbit departure,
spatial density and heating as well: a small global angular-momentum error
alone cannot demonstrate preservation of the disk. Density images must share
stated scales and units with the comparison images.

## Reproducible evidence and current execution

The adjacent `audit_kepler_paper_time_units_20261022.py` reads only an 80-byte
saved header and two frozen source files, records hashes, and evaluates the
analytic expressions. Its exclusive-create output is
`KEPLER_PAPER_TIME_UNITS_AUDIT_20261022.json`. It ran successfully; this is
saved-input analysis, not another small evolution test. No GPU kernel or RK
stage was executed by this audit.

At 02:25:12 UTC the original producer launcher, four MPI ranks, and dispatcher
were confirmed live with their retained process start identities. Preprocessing
had completed 787896/1138432 triangles (69.21%). No terminal exit code or
physical-transition claim existed. The latest accepted physical state remained
step 4, t=8.2355880701376864e-5. The paired accuracy/cost decision still comes
before any new physical continuation or conditional map-guard work.
