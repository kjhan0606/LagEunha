# Two-dimensional screen of static weight condensation

This follows WEIGHT_STATIC_CONDENSATION_2026-11-06.md. It tests a proposed
extension, not the existing native solver and not a KH/RT evolution.
No production physics, checkpoints, or successful Kepler results changed.

## Representation tested

On a unit square lattice each displacement component is quadratic only in
its own local coordinate. Its mean and the two normal face values determine
the polynomial. A single scalar weight per cell determines different face
offsets through neighbor differences. Equal speeds on all faces are neither
required nor imposed. This separable representation is not asserted to
describe general Laguerre polygons or their rotating faces.

For Fourier amplitudes z=(Xx,Xy,W), define rows

    Da z = i sin(theta_a) Xa + (1-cos(theta_a)) W,
    Sa z = (cos(theta_a)-1) Xa - i sin(theta_a) W.

The component field is Xa + (Da z) ra + 3(Sa z)(ra^2-1/12).
Exact cell integration gives kinetic matrix

    M = diag(1,1,0) + sum_a [Da^H Da/12 + Sa^H Sa/20]

and linear compressive energy matrix

    A = (Dx+Dy)^H(Dx+Dy) + 3 sum_a Sa^H Sa.

Four-point tensor Gauss integration independently checks both polynomial
integrals for 100 random fields, maximum absolute error 1.43e-14.
Condense W by W=-A_WW^-1 A_WX X, and lift BOTH matrices. The grid-axis
case reproduces the prior 1D result and has a zero transverse eigenvalue.

## Failure exposed by an oblique transverse mode

At theta_x=theta_y=t, a transverse displacement Xx=-Xy has condensed W=0.
Although its mean divergence is zero, its reconstructed local divergence
is not. The energy therefore gives an artificial elastic restoring force:

    omega_T^2 = 3(1-cos t)^2 /
                [1 + sin(t)^2/12 + (1-cos t)^2/20].

This analytic expression agrees with the generalized eigenvalue calculation.
For a wave oriented at 45 degrees, its artificial transverse speed divided
by sound speed is 0.04249, 0.16842, 0.32748, 0.58927 at wavelengths of
64, 16, 8, 4 cells respectively. It tends to zero under refinement; this
is a finite-resolution consistency error, not a proof of a nonzero continuum
shear modulus. It is elastic, not dissipative: do not report it as heating.
Nevertheless the anisotropic response is unacceptable as an unqualified
common closure for the intended low-pressure shear tests.

This identifies a defect of the separable internal reconstruction, NOT a
proof that scalar weights, static condensation, or the native solver must
have the same defect. A divergence-compatible vector reconstruction would
need independent evidence before adoption. Adding a shear sensor or subtracting
this energy after a step would not repair the underlying representation.

## Complexity check

Let qa=1-cos(theta_a), Q=qx+qy. The weight coefficient is

    A_WW = Q^2 + 3[sin(theta_x)^2+sin(theta_y)^2].

After factoring Q, its symbol is 6+Q-3(qx^2+qy^2)/Q, between 2 and 6.
In 2D this quotient is direction dependent; the 1D local cyclic tridiagonal
equation cannot simply be reused. A Laplacian preconditioner is plausible,
but no native solve or end-to-end timing was measured. The sampled maximum
reduced frequency is sqrt(10); this is not a proof of the global maximum
and is not a GIZMO runtime comparison.

## Decision

Do not port this separable extension into the C/CUDA evolution. The next
closure must preserve the divergence-free subspace as well as pressure work,
before investing in KH/RT runs. Existing native trajectories and their known
weight-admissibility failure remain the empirical reference. Full 2D testing,
3D validation after 2D, and <=2x GIZMO runtime are still incomplete.

Reproducer and full output:
Research/MaterialVolume_20261105/weight_condensation_2d_screen_20261106_v1/
audit.py and results.json. This was a local Python algebraic diagnostic,
not a replacement simulation implementation.
