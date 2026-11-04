# Kepler acceptance relative to GIZMO, 2026-11-04

The user's latest instruction supersedes the earlier absolute cohort-error
interpretation: compare errors with GIZMO at the corresponding physical time;
errors need not be driven toward zero if they are not materially larger.
Do not add complexity or costly tests solely to satisfy an arbitrary tighter
internal threshold. Prefer the simpler and faster adequately accurate method.

## Primary comparison

Use the GIZMO MFM result at the same time and effective disk resolution,
with the same error definition and radial selection. Report the error ratio
where a reference measurement exists. No new arbitrary ratio cutoff is set.
MFV can be shown as a secondary reference, labeled separately.

Hopkins (2015), section4.2.4 and Figure8, defines equal-element mean errors
in the **current-radius annulus .8<R<1.2**:

    L1(v) = mean(abs(v_phi - v_K(R))/v_K(R))

    L1(rho) = mean(abs(rho - 1))

The fixed-time velocity comparison is at t=20*pi (ten R=1 orbital times),
with effective N1D=256 for the standard case. Figure8's MFM/MFV endpoint
velocity errors at this resolution are approximately6e-3, read from the
graph rather than a published numerical table. Treat this as approximate.

The .01 velocity and .25 density values in Figure8 are **thresholds used
to plot first-crossing times**, not measured GIZMO errors at t=20*pi.
Being below those thresholds alone does not establish parity with GIZMO.
The figure does not tabulate the density L1 value at that exact time.

The original profile has a1.01 central plateau including its.01 background.
Retain literal density error relative to1 for the published definition and
also report change from the actual initial profile. Do not silently subtract
the initial offset in only one code.

Source: https://arxiv.org/pdf/1409.7395, PDFpages20 and22 (Figure8 and metric
definitions). The downloaded PDF/hash and page rendering are preserved in
/gpfs/kjhan/LagForce/Kepler_Q1_BoundedGlobal_20261029/gizmo_error_comparison_20261104_v1.

## Internal diagnostics and essential integrity

The old original-radius cohort density-change1% screen remains useful as a
diagnostic of where errors arise. Exceeding it **does not by itself fail
the user's GIZMO-relative accuracy criterion**. Retain historical results
and flags as measured; qualify their meaning instead of deleting or changing
the measurements. Original-radius mass-weighted cohort errors are not the
same metric as the published current-radius equal-element annular L1.

Continue checking state validity, mass, energy/angular-momentum accounting,
and disk morphology. Record external-wall impulses and boundary/IC differences
when interpreting errors. A custom wall-impulse budget is not itself a
published GIZMO-relative accuracy measurement. Keep integrity defects distinct
from the size of a valid numerical solution's error.

## Available evidence and next endpoint

Completed native RK4/circle branch448969 is at exact20*pi. Its published-
definition velocity L1 is3.39535e-7, far below the approximate6e-3 graph
value. Its annular literal density L1 is.0112173, but no exact-time published
GIZMO density value has been extracted, so no direct density-error ratio is
claimed. Its disk-cohort density change1.4386% is a different diagnostic,
and is not grounds for rejection by the old1% screen under this instruction.

These numerical magnitudes support adequate velocity accuracy in that
realization, not an assertion of intrinsically superior hydrodynamics:
particle placement and boundary conditions differ from the paper reference.

Active449508 is an RK4-prefix/global-KDK/block2-KDK continuation, with a
different post-t40 boundary from448969. Finish its already running exact20*pi
endpoint, then apply this comparison to its own output. Do not substitute
an intermediate t~60 result for the fixed-time endpoint.

The separately running GIZMO449138 has the same Q1 t0 positions/masses/
velocities and softened force but a periodic box8 rather than a reflective
polygon. Its existing first snapshot is att120, not20*pi. It is valuable
for a future same-time t120 comparison; do not compare its t120 errors to
our t20*pi errors as if they were matched, or repeat its evolution merely
to manufacture a missing earlier snapshot.
