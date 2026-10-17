import argparse,json
from pathlib import Path
import numpy as np
from power_cells_3d import Geometry3D
from check_3d_coupled_rhs import initial
from kinetic_moment_rhs_3d import rhs
from moment_state_3d import unpack


def run(library):
    geo=Geometry3D(library);rng=np.random.default_rng(271);sites=np.array([[i,j,k] for i in (.25,.75) for j in (.25,.75) for k in (.25,.75)])+rng.uniform(-.035,.035,(8,3));weights=np.zeros(8)
    raw,L=initial(geo,sites,weights,lambda p:np.zeros_like(p));gradient=np.array([.1,-.2,.15])
    for r in raw:r[34]=(r[0]+r[1:4]@gradient)/(5/3-1)
    dr,gen,z,diag=rhs(geo,sites,weights,raw,L,face_order=3,pressure_reconstruction='linear');exact=np.zeros_like(raw)
    for i,r in enumerate(raw):
        m,q,Q,p,D,Pi,E=unpack(r);exact[i,13:16]=-m*gradient;exact[i,16:25]=-np.outer(q,gradient).ravel()
    error=float(np.max(abs(dr-exact)));assert error<2e-12
    return dict(passed=True,max_raw_rhs_error=error,mass_flux=diag['mass_flux'],pressure_model=diag['pressure_model'],scope='Affine pressure with exact 3D Euler moment acceleration')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('library');a=p.parse_args();result=run(a.library);Path('linear_pressure_3d_report.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
