# Avoiding double-counted spin in a reconstructed velocity field

The native solver already stores a velocity-gradient tensor `A`. Adding
`omega × (x-C)` directly on top of that full gradient would generally count
its rigid-rotation component twice. The isolated
[host/CUDA-compatible operator](gfs_affine_spin_reconstruction_20261027.h)
instead separates the rigid part using the mass central second-moment tensor
`Q`, about the mass centroid `C`:

```
I = Qxx + Qyy
S_A = Ayx*Qxx - Axy*Qyy + (Ayy-Axx)*Qxy
omega_A = S_A/I
B = A - omega_A*J,   J = [[0,-1],[1,0]]
A_effective = B + (S_state/I)*J.
```

`B` has zero integrated intrinsic angular momentum under the assumed
moment tensor. Its non-rigid affine kinetic energy is retained separately,
and is Q-orthogonal to the rigid part, so
`K(A_effective) = K(B) + S_state^2/(2I)`. Setting `S_state=S_A` reconstructs
the original affine gradient algebraically. This is a kinematic closure
building block, not a physical torque or a completed hydrodynamic update.

The C++ [host test](test_gfs_affine_spin_reconstruction_20261027.cpp)
passes anisotropic-Q projection, the zero-state-spin and unchanged
gradient cases, the energy split, and invalid-tensor rollback. The separate
[CUDA device test](test_gfs_affine_spin_reconstruction_cuda_20261027.cu)
completed on A10 as Slurm **437559**. Its output reports `PASS: GPU affine
spin projection, energy split and invalid rollback`; the staged source and
executable checksums passed. The earlier test **437549** was canceled before
allocation. These tests do not establish the runtime cost on a full disk,
face-force/torque physics, RK4 staging, local flux conservation, MPI exchange,
or t>6 stability.

An additional 100,000 deterministic randomized positive-definite, often
anisotropic moment tensors and gradients passed the host test's independent
checks that the effective gradient has the requested intrinsic spin and that
its integrated affine kinetic energy equals the spin plus non-rigid terms.
This is an algebraic invariance test; it does not establish physical flux
accuracy or a stable numerical evolution.

An important reduction distinction remains: when this spin reconstruction
is **disabled**, the native solver must use its original `A` unchanged. If
the reconstruction is enabled and `S_state=0` while the original gradient
has `S_A!=0`, the operator intentionally removes that gradient's rigid
component and therefore does not automatically reproduce the old solver.
Any native integration must state its W0/S=0 reduction explicitly, satisfy
the user's W0 fixed-mass condition, and avoid claiming equivalence merely
because the code compiles. Mass-centroid/geometry closure, initialization of
the extended total energy, and a physically specified traction remain open.

### W=0 and fixed-mass meaning

In this project, `W=0` means every generator's Laguerre weight is zero
(stored as `w2_i=0` in the geometry path), so the power bisectors reduce to
ordinary Voronoi bisectors. It does not mean zero cell mass and does not
mathematically imply fixed mass. The current q1 checkpoint's `params.dat` has
`KH Kappa=0`; its `GAS w2 power=0.2` does not make a nonzero weight when
`Kappa=0`. The tested RK4 build also omits the `LAGMFV` compile flag, which
excludes the optional mass-flux accumulation and half-step mass update. Thus
the accepted W=0 reduction here is fixed-mass by the selected solver path,
not by the value W=0 alone. A future time-varying-weight path must preserve
that mass invariant unless mass transfer is explicitly introduced; audit
per-cell mass by particle identity as well as total mass. This stipulated
reduction is not the exact geometric ALE mass flux: a Voronoi face is
generally not material, so that flux changes individual cell masses under
shear. The measured counterexample and its modeling implications are in the
[W=0 ALE fixed-mass contract](KEPLER_W0_ALE_FIXED_MASS_CONTRACT_2026-10-27.md).

An additional algebraic option stores **excess** spin
`delta_S = S_physical - S_A`. Its local face gradient is
`A_effective = A + (delta_S/I)*J`; setting `delta_S=0` leaves the native
gradient bitwise unchanged in the host test, while the physical intrinsic
spin is still `S_A`, not zero. This could help an exact no-extra-spin
reconstruction reduction, but it does not solve energy initialization or
evolution: because `S_A` changes with the gradient and cell shape, the
correct state equation is `d(delta_S)/dt = physical_torque - d(S_A)/dt`.
The corresponding reference-energy variation must also enter the energy
ledger. Treating `delta_S` itself as the physical spin or omitting those
time derivatives would be incorrect. This excess-spin option is a local
prototype, not a proposed production closure.

The separate [two-cell affine-spin ledger](gfs_spin_affine_face_exchange_20261027.h)
extends the existing pair exchange with an explicit non-rigid kinetic-energy
reservoir. Its residual thermal energy is total minus bulk, spin and that
reservoir; it rejects an update transactionally if either cell's thermal or
reservoir energy would become invalid. The caller supplies the physical face
transfers and non-rigid-energy changes: the ledger does not invent a stress,
friction or dissipation law. Its [host test](test_gfs_spin_affine_face_exchange_20261027.cpp)
passes pair mass, momentum, first mass moment, angular momentum, total
energy and positivity rollback. The separate actual-device test completed
on A10 as **437550**. Its output reports `PASS: GPU nonrigid/spin/thermal
pair conservation and rollback`; the staged source and executable checksums
passed. This remains accounting infrastructure, not a native Kepler solver
or a demonstration of non-rigid energy evolution.
The host check also exercises 50,000 deterministic randomized pairs and
verifies mass, linear momentum, first mass moment, total angular momentum,
total energy and the explicitly exchanged non-rigid reservoir. These are
accounting invariants for supplied exchanges, not proof that the supplied
fluxes model inviscid fluid physics.
The same randomized pair test also applies a common Galilean velocity boost,
transforming momentum and energy transfers consistently; centroid, intrinsic
spin and thermal energy agree with the unboosted result. This checks the
instantaneous ledger algebra, not Galilean invariance of the complete RK4,
mesh-motion, gravity and MPI implementation.

## Consequence of the requested `Re -> infinity` limit

`Re -> infinity` removes physical viscous shear stress; it does not prescribe
a spin-relaxation time. A smooth, shock-free Kepler flow therefore must not
convert cell spin into heat through an added drag or damping term. Entropy
generation at a captured shock is a separate Euler weak-solution effect and
must be diagnosed separately from smooth-orbit heating.

A conservative spin update can be derived from the angular-momentum balance
of the same moving control volume used by the native finite-volume scheme.
For a cell with mass centroid `C`, bulk momentum `P`, intrinsic angular
momentum `S`, face normal `n`, and mesh velocity `w`, its total angular
momentum is `L = C x P + S`. The Euler angular-momentum flux through a moving
face is the moment of the existing momentum flux:

```
dL/dt = - integral_boundary x × [ rho v ((v-w) · n) + p n ] dA
         + integral_volume x × (rho g) dV.
```

This says how the cell angular-momentum budget must be closed; it does not
invent a tangential viscous force. In a discrete update, mass, momentum,
centroid motion, and this moment flux must use the same RK stage and face
quadrature. The spin increment is then the change in `L - C x P`. Pairwise
face contributions cancel in the global budget, while each cell can exchange
angular momentum with its resolved bulk motion. If the code stores excess
spin over the angular momentum already represented by its native velocity
gradient, the corresponding gradient-spin change and reference-energy
change must remain in the same ledger.

The rotational energy `S^2/(2I)`, with `I` the mass second moment about `C`,
is only one part of the kinetic-energy partition. The native total-energy
flux must remain authoritative; the RK4 update must repartition its change
between bulk, affine/non-rigid, spin, and thermal reservoirs without adding
energy or counting the native gradient's rigid component twice. No independent
spin-to-thermal relaxation is admissible in the stated inviscid smooth-flow
limit.

This derivation narrows the next implementation review but is not a native
solver update or a validated closure. Before a Kepler production run, the
driver should demonstrate: (1) exact original-solver reduction when the spin
extension is disabled; (2) cell and global angular-momentum closure against
an independently integrated moving-face flux; (3) total-energy and
positivity closure through every RK4 stage; and (4) no secular thermal rise
in a smooth affine/shear control before testing the Kepler disk beyond `t=6`.

## Native integration contract for the driver

The next step is to connect the tested ledger to the existing RK4 face flux,
not to add a relaxation or friction term. For each cell, form the mass-central
moment tensor Q from the same density model and polygon used by the native
finite-volume update, with 2D moment of inertia I=Qxx+Qyy. Let S_A be the
intrinsic angular momentum already implied by the native affine gradient A.
Store excess spin delta_S=S_physical-S_A and reconstruct

    A_effective = A + (delta_S/I) J .

This keeps the current face reconstruction unchanged at delta_S=0. The state
equation must nevertheless evolve

    d(delta_S)/dt = torque_from_native_face_fluxes - d(S_A)/dt .

For a momentum impulse delta_p transferred from cell i to j through a face
with centroid x_f, the pair torque increments about their mass centroids are

    delta_S_i = -(x_f-C_i) cross delta_p
    delta_S_j = +(x_f-C_j) cross delta_p .

These are the signs and moment arms used by the tested two-cell ledger. With
the matching equal-and-opposite bulk momentum update, orbital plus intrinsic
angular momentum is conserved. For fixed-mass Kepler runs there is no mass
transfer term; if the operator is later used with mass flux, the centroid
shift terms in the ledger must also be retained.

The energy state must include the subcell affine kinetic part
K_sub=K_nonrigid+S_physical^2/(2I), with
K_nonrigid=0.5 tr(B Q B^T) and B=A-(S_A/I)J. The native ALE energy flux
must then update total energy consistently with the same stage face fluxes;
thermal energy is the residual after bulk and subcell kinetic energies.
Initializing this extra reservoir by subtracting it from the old thermal
energy while holding the old total-energy reference fixed is not valid. The
initial extended energy must instead be derived consistently from the
continuous Kepler IC and checked for positive thermal energy.

For the requested Kepler limit Re=infinity, the spin equation has no
phenomenological damping or relaxation. Spin may exchange with resolved
orbital motion through physical face torques, while total energy remains
conservative; any decay must be measured as numerical transfer, not imposed
friction. The isolated tests above do not supply the face torque, the
non-rigid energy transfer, Q/I stage evolution, gravity coupling, or MPI
halo exchange.

Native release gates for a spin-enabled branch are therefore: the disabled
branch reproduces the current W=0 fixed-mass trajectory and ledgers; rigid
rotation advects without spurious heating; an off-centre face impulse
conserves pair orbital-plus-intrinsic angular momentum and total energy;
each RK4 stage keeps Q positive definite and thermal energy positive; and a
short native Kepler continuation passes the independent density, orbit,
angular-momentum, energy, and power-geometry audits before any t>6 spin run.

## Low-cost geometry path

The current GPU tessellator already visits every polygon corner to compute
the shoelace area and already stores the clipped vertices in PolygonOut.
For a uniform-within-cell density closure, the mass centroid and Qxx, Qxy,
Qyy follow from the standard signed polygon-moment sums in that same
per-cell pass. This adds work proportional to the existing corner count and
needs no extra neighbor search, global solve, or subcell-force integration.
Center the raw moments about the resulting mass centroid before using
I=Qxx+Qyy.

Keep the raw clipped polygon area separate from any numerical volume floor.
If the geometry is hidden or the floor is invoked, reject or flag its moment
tensor instead of treating the artificial floor as physical inertia. Compare
the GPU centroid and Q tensor to an independent polygon-moment audit first.
The current postprocessing assumes uniform density; any later density-gradient
moment correction should be justified by a measured change in the Kepler
density and spin-energy gates.

The native particle position is the mesh generator; the two-cell ledger's
position is the mass centroid. Those roles must remain distinct unless a
measured gate shows they coincide closely enough. Fixed cell mass does not
fix the mass centroid: polygon deformation or a nonuniform density
reconstruction can move it. Each RK4 stage must therefore use a consistent
centroid, central moment tensor, face centroid and bulk-angular-momentum
reference, while keeping the generator coordinates for tessellation.

### Signed polygon formula for the first native geometry prototype

The CUDA clipping kernel stores vertex coordinates relative to the generator,
so the moments can be accumulated without subtracting large global
coordinates. For vertex `i`, let `j=(i+1) mod n` and
`c_i=x_i*y_j-x_j*y_i`. With the current counter-clockwise orientation:

```
A  = 1/2 Σ c_i
cx = Σ (x_i+x_j)c_i / (6A)
cy = Σ (y_i+y_j)c_i / (6A)

Jxx = Σ (x_i²+x_i*x_j+x_j²)c_i / 12
Jyy = Σ (y_i²+y_i*y_j+y_j²)c_i / 12
Jxy = Σ (2*x_i*y_i+x_i*y_j+x_j*y_i+2*x_j*y_j)c_i / 24
```

These `J` values are geometric area moments about the generator. Centralize
them with `Jcxx=Jxx-A*cx²`, `Jcyy=Jyy-A*cy²`, and
`Jcxy=Jxy-A*cx*cy`; for uniform cell density, the mass centroid is the
generator plus `(cx,cy)`, `Q=(m/A)*Jc`, and `I=Qxx+Qyy`. Keep signed sums
until positive area and polygon orientation pass. A cell whose clipped area
is invalid or whose effective volume came from the numerical floor must have
invalid moments and be excluded from spin evolution.

The first implementation gate is an independent host `long double` audit of
all returned GPU centroid/Q components, plus analytic rectangle and triangle
controls. Only after that passes should the state arrays, RK4 staging, restart
format, or face torque be changed. Measure transfer and full-step wall-time
overhead because the GPU kernel runs at every RK stage.

The host reference is now implemented in
[`gfs_polygon_moments_reference_20261028.h`](gfs_polygon_moments_reference_20261028.h)
with rectangle, asymmetric triangle, translated-triangle, reversed-orientation
and degenerate-polygon controls in
[`test_gfs_polygon_moments_reference_20261028.cpp`](test_gfs_polygon_moments_reference_20261028.cpp).
It compiled with `g++ -std=c++17 -O2 -Wall -Wextra -pedantic` and reported
`PASS: signed polygon area, centroid and central Q reference`. This validates
the independent formula implementation only; no GPU moments or spin dynamics
have been added to the native solver yet.
