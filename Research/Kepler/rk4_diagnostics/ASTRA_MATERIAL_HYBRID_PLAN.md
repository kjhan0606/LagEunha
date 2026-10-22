# Material hybrid interface bundle

2026-10-21. Records the completed bundled Astra planning/evaluation decision
from the shared-H1/source-control review. Coding and execution remain with the
driver; no new micro-review or production restart is requested.

The preceding opinion turn reread existing refinement evidence but changed no
authoritative implementation or next decision: classify it as no progress.
The available safe next action is the independent hybrid operator below.

## Equations and representation

Use broken cell-local P1 plus the two initial-velocity residual functions.
Derive B from actual physical face traces and compare ker(B) to the old shared
space before comparing accelerations. The old component-wide residual sharing
may impose stronger conditions than H1 conformity. Cell bubbles must not be
constrained. Disconnected components must remain disconnected.

Solve independently, without fitting a previously computed force residual:

    M a = F + B^T lambda,   B a = 0,
    (B M^-1 B^T) lambda = -B M^-1 F.

M is block diagonal by cell. F contains volume stress/gravity and only known
EXTERNAL traction. Internal pressure traction is not counted twice. Fixed
material-reference B has Bdot=0. Choose a deterministic independent-row gauge
before looking at forces, energy, or torque. No gauge optimization toward a
physical-pressure target is allowed.

For this axis-aligned fixture, endpoint trace constraints and analytically
independent velocity-residual traces give incidence rows. On axis-radial faces
one residual is identically zero and must not be constrained. Verify positive
Gram rank of all remaining actual trace functions by positive quadrature.
Record this fixture restriction; it is not a generic curved-face rank solver.

## Physical interface and gates

Represent each multiplier as an actual face traction in the L2-dual trace
basis: G eta=mu, G=int(phi phi^T)dS0. Check its direct quadrature load against
B^T lambda. Shared position and velocity traces must give opposite momentum,
work and angular transfer. Retain smooth physical Udot=-int P div(v) and the
original source. No heat/source redistribution, independent kinetic reservoir,
viscosity, or damping is introduced.

Report separately (a) cell energy/angular closure using numerical traction and
(b) numerical-minus-physical pressure traction, including its work and torque.
The prior15.7 U/time defect is not actual heating and may reappear in (b).
Old residual representability requires vanishing work on ALL ker(B) modes,
including bubbles, not just vanishing global summed force.

Bounded cases: gravity off, representable linear gravity, actual softened
Kepler gravity, and actual gravity plus uniform velocity perturbation. Run
native MPI1/4; report rank, Schur positivity/conditioning, acceleration and
traction agreement, and quadrature sensitivity. A finite-map traction measure
check is a geometry control, not a time integrator.

Stop promotion if physical-trace representation, pair work/torque, bubble
annihilation, or discrete ledgers fail, or if pressure/source approximation is
hidden by renaming a residual as traction. A passing initial operator does not
establish stability. Curved Kepler remap/re-expression, density-wave/sliding,
shock stress, native CUDA runtime and historical t6-10 remain dependencies.
The final target remains the longest GIZMO paper Kepler evolution.
