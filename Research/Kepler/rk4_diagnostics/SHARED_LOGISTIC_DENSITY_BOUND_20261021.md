# Same density, two-sided analytic exponential representation

The historical complex_density_bound requires both |exp(z_inner)|<1 and
|exp(z_outer)|<1. This is a sufficient geometric-series condition only on
one side of each transition. In particular, a real positive exponent is not
itself a logistic pole. The original profile is left unchanged:

rho(r)=rho_bg+rho_amp/(r+eps_rho) L(z_inner) [1-L(z_outer)],
L(z)=exp(z)/(1+exp(z)), z_j=2(r-R_j)/taper.

For the existing complex disk z=c+w, |w|<=R and real c, choose E=exp(-z)
when c>0, E=exp(z) otherwise. With the original holomorphic_exp/native libm
premise this gives a disk E=a+delta, |delta|<=b. The stored strictly positive
gap 1+a-b excludes denominator zeros everywhere in the disk. Use:

- c<=0: L=E/(1+E), 1-L=1/(1+E).
- c>0: L=1/(1+E), 1-L=E/(1+E).

These are exact holomorphic identities, with the same actual poles
z=(2k+1) i pi, original taper/eps/profile and original disk. They neither
enlarge a graph neighborhood nor relax an accuracy criterion. If the disk
bound for 1+E is not positive the result remains unproved. No automatic
refinement or fallback is authorized by this identity alone.

The first validation consumes existing failed density-radius records, with
zero new sd::Point or physical quadrature queries. Passing it certifies only
that old failed radius disk under the new equivalent expression. It cannot
certify the unattempted suffix of that slot or any all12/time-evolution result.
