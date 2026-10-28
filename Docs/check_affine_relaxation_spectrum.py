"""Frozen shared-face affine damping spectrum and a conservative RK4 bound.

gamma_face=alpha*A*harmonic(rho*c_s) is an experimental numerical-viscosity
scale, not an Euler tangential traction law. Here alpha=rho=c_s=1.
"""
import json,sys
import numpy as np
sys.path.insert(0,'/gpfs/kjhan/LagForce/Codex_moving_cell_spin_20261018')
from affine_moment_patch import bounded_mesh,integrate
rows=[]
for n in [4,8]:
    xx,yy=np.meshgrid((np.arange(n)+.5)/n,(np.arange(n)+.5)/n)
    x=np.column_stack((xx.ravel(),yy.ravel()))+np.random.default_rng(731).uniform(-.15,.15,(n*n,2))/n
    polys,edges=bounded_mesh(x,np.zeros(n*n))
    mom=[integrate(p,lambda z:np.zeros_like(z)) for p in polys]
    roots=[];inverse=[]
    for a in mom:
        values,vectors=np.linalg.eigh(a['S']);assert values.min()>0
        roots.append((vectors*np.sqrt(values))@vectors.T)
        inverse.append((vectors/np.sqrt(values))@vectors.T)
    dimension=6*n*n;matrix=np.zeros((dimension,dimension));max_pair=0.;basis_rows=[]
    nodes,weights=np.polynomial.legendre.leggauss(3)
    for i,j,a,b,normal,length in edges:
        if j<0:continue
        tangent=(b-a)/length
        for z,w in zip(nodes,weights):
            point=.5*(a+b)+.5*z*(b-a);row=np.zeros(dimension)
            for owner,sign in [(i,1),(j,-1)]:
                arm=point-mom[owner]['center'];dual=inverse[owner]@arm
                row[6*owner:6*owner+6]=sign*np.r_[tangent/np.sqrt(mom[owner]['area']),np.outer(tangent,dual).ravel()]
            gamma=.5*length*w
            matrix+=gamma*np.outer(row,row);max_pair=max(max_pair,gamma*np.dot(row,row))
            basis_rows.append(row)
    eigen=np.linalg.eigvalsh(matrix);maximum=float(eigen[-1]);bound=float(np.max(np.sum(abs(matrix),axis=1)))
    assert eigen[0]>-1e-11 and maximum<=bound*(1+1e-12)
    dt=.25/bound;z=-dt*eigen
    amp=1+z+z*z/2+z**3/6+z**4/24
    assert np.max(abs(amp))<1+1e-12
    null_errors=[]
    # Six global affine velocity degrees of freedom, including rigid rotation.
    for k in range(6):
        v0=np.zeros(2);G=np.zeros((2,2))
        if k<2:v0[k]=1
        else:G.ravel()[k-2]=1
        state=np.concatenate([np.r_[np.sqrt(s['area'])*(v0+G@s['center']),(G@root).ravel()] for s,root in zip(mom,roots)])
        null_errors.append(float(np.max(abs(np.array(basis_rows)@state))))
    assert max(null_errors)<1e-12
    cs=np.sqrt(float(np.float32(1.666667))*((float(np.float32(1.666667))-1)*(5*float(np.float32(1.666667))/2)*1e-5)/.01)
    rows.append(dict(n=n,degrees_of_freedom=dimension,lambda_max=maximum,
       absolute_row_sum_bound=bound,max_individual_contact_rate=max_pair,
       assembled_over_max_contact=maximum/max_pair,minimum_eigenvalue=float(eigen[0]),
       count_eigenvalues_abs_below_1e_10=int(np.sum(abs(eigen)<1e-10)),
       max_affine_null_slip=max(null_errors),conservative_frozen_dt=dt,
       max_RK4_amplification_at_bound=float(np.max(abs(amp))),
       homogeneous_cold_background_cs=cs,
       homogeneous_cold_background_scaled_lambda_max=cs*maximum,
       homogeneous_cold_background_scaled_dt=dt/cs))
print(json.dumps(dict(scope=__doc__,cases=rows,
    limits=['Frozen linear mechanical block only; not nonlinear hydro/gravity/geometry stability or heat positivity.',
    'Bound uses kinetic-normalized affine coordinates; raw-moment row sums are not equivalent.',
    'Uniform background scaling is analytic on these small meshes, not the actual 192-square disk.',
    'No chosen alpha or physical calibration, no production timestep changed.']),indent=2))
