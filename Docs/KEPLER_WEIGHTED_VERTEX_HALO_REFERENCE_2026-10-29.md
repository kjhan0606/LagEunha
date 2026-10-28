# Weighted vertex halo certificate: proof and isolated host reference

The vertex-disk halo completeness condition extends from W=0 to weighted
Laguerre planes, including the native clamp to [0.05,0.95]. Astra reviewed
the proof against the actual native clamp; the driver implemented a separate
host reference and adversarial tests. No running executable, native source,
frozen W0 candidate, queued batch or checkpoint dispatcher was changed.

## Sufficient condition

Let p own a finite convex candidate polygon, let v=p+u be one of its vertices,
and let D be the rectangle actually covered by complete generator data. Let
Wmax be an upper bound on every generator weight at the same RK stage,
including generators outside D and periodic images. W denotes the stored
power weight (native w2), with units of length squared. Set

    delta = Wmax - Wp >= 0
    R(v)^2 = |u|^2 + delta.

If the closed disk centered at every v with radius R(v) lies strictly inside
D, no outside generator can cut the polygon. This is sufficient, not
necessary; rejection means completeness remains unproven.

For q outside D, write r=q-p and d2=|r|^2>0. Disk containment implies

    |q-v|^2 > |u|^2 + delta,
    u dot r < B = (d2-delta)/2.

The ordinary Laguerre plane has right side A=(d2+Wp-Wq)/2. Since Wq<=Wmax,
A>=B. The actual native clamped right side is

    Ac = max(0.05*d2, min(0.95*d2, A)).

Also B<=0.5*d2<=0.95*d2, so both arguments of the inner minimum are at least
B, and therefore Ac>=B. The upper clamp can reduce A, but cannot reduce it
below the required lower bound B. The outside plane thus preserves every
vertex; its affine inequality preserves the entire convex hull. The proof
does not use a radial search approximation or alter any cell force/work.

## Preconditions that the certificate cannot establish

- Every relevant generator inside D must have been searched completely.
  A missing internal generator is not protected by this proof. Rounded
  linked-cell grid limits must not substitute for actual halo coverage.
- The polygon must be finite and convex, with no active artificial-box
  support or capacity truncation. Check every vertex in the same unwrapped
  coordinates as the periodic images used for clipping.
- Wmax must bound the weights actually used for this stage. A fresh global
  owner reduction is sufficient only if halo/image weights agree with their
  owners. Time-dependent weights require a new reduction or another proven
  bound after each change; the previous accepted step's maximum is not enough.
- Negative weights are allowed. Initializing a maximum reduction to zero
  is a safe overbound when all weights are negative, but can reject otherwise
  sufficient halos. An underestimated rank-local maximum is unsafe.
- Coincident generators, nonfinite input, and native floating clipping
  errors need their own checks. This is an exact-arithmetic proof for the
  represented candidate polygon, not a numerical error bound for its creation.

The reference header rejects Wmax<Wp, nonfinite values and overflow. It
bounds delta and radius squared upward, and clearance and its square
downward, before a strict comparison. If Wmax==Wp it calls the existing W0
function directly, preserving its exact decision, including at Wp=Wmax=0.
Equal negative weights with an overbound Wmax=0 remain safe but need not give
the same decision as a tighter equal-weight bound.

## Validation and counterexamples

The isolated C++ test uses long-double clearance inequalities and directly
evaluates outside-site standard and clamped planes, rather than using the
disk certificate as its plane oracle. It includes 100000 exact W0/equal-weight
decision comparisons, 20000 near-boundary weighted cases, explicit lower and
upper clamp regions and transition-adjacent inputs, tangent rejection,
large-coordinate cancellation, underflow/overflow and invalid inputs.

It also demonstrates why a current global weight bound is mandatory. With
p=(0,0), v=(0.1,0), D=[-0.5,0.5]^2, q=(0.6,0), Wp=0 and Wq=0.36, a stale
Wmax=0 passes the disk check even though q cuts v under both standard and
clamped planes. The correct Wmax=0.36 rejects the halo. A second example uses
a whole finite triangle and an unclamped cutting plane. Rejection does not
prove that a missing generator exists; these are explicit constructed cutters.

Reproduce the small host check from this repository:

```sh
g++ -std=c++17 -O2 -fno-fast-math -ffp-contract=off \
  Docs/test_gfs_weighted_vertex_halo_reference_20261029.cpp \
  -o /tmp/test_gfs_weighted_vertex_halo_reference_20261029
/tmp/test_gfs_weighted_vertex_halo_reference_20261029
```

The [receipt](KEPLER_WEIGHTED_VERTEX_HALO_REFERENCE_2026-10-29.json) records
the exact test output, source hashes and compiler command. The reference is
not connected to the solver, and has not been built or executed on CUDA.
This establishes neither native clipping/face/force parity, MPI completeness,
speed improvement nor long-time stability. Native integration remains deferred
until the existing W0 comparison supplies its own geometry and physics evidence.
The original t>6 and GIZMO long-time objectives remain unverified.
