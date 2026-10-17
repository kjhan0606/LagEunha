"""Actual C geometry plus positive-factor 3D midpoint evolution gates."""
import argparse,json
from pathlib import Path
import numpy as np
from power_cells_3d import Geometry3D
from check_3d_coupled_rhs import initial
from moment_evolution_3d import step
from moment_state_3d import unpack


def run(library):
    geo=Geometry3D(library);rng=np.random.default_rng(271)
    sites0=np.array([[i,j,k] for i in (.25,.75) for j in (.25,.75) for k in (.25,.75)])+rng.uniform(-.035,.035,(8,3));velocity=np.array([.2,-.1,.3]);rows=[]
    for dt in (.05,.025,.0125):
        sites=sites0.copy();weights=np.zeros(8);raw,L=initial(geo,sites,weights,lambda p:np.tile(velocity,(len(p),1)));original=raw.copy();ledger=np.zeros(7);peak=0.;mass=raw[:,0].copy()
        bounds=lambda t:np.array([[0.,0.,0.],[1.,1.,1.]])+t*velocity
        for k in range(round(.1/dt)):
            sites,weights,raw,L,d=step(geo,k*dt,sites,weights,raw,L,dt,bounds,lambda p:velocity)
            ledger+=d['ledger_increment'];peak=max(peak,d['stage_mass_flux'])
        exact=original.copy();shift=.1*velocity
        for i,r in enumerate(original):
            m,q,Q,p,D,Pi,E=unpack(r);exact[i,1:4]=q+m*shift;exact[i,4:13]=(Q+np.outer(q,shift)+np.outer(shift,q)+m*np.outer(shift,shift)).ravel();exact[i,16:25]=(D+np.outer(shift,p)).ravel()
        error=float(np.max(abs(raw-exact)));position_error=float(np.max(abs(sites-sites0-shift)))
        assert error<3e-11 and position_error<3e-11 and np.array_equal(raw[:,0],mass)
        rows.append(dict(dt=dt,t=.1,max_raw_error=error,position_error=position_error,mass_exact=True,stage_mass_flux=peak,ledger_max=float(max(abs(ledger)))))
    return dict(passed=True,rows=rows,scope='3D actual moving cells and positive midpoint uniform-flow evolution; not a nonlinear accuracy or dense-path certification')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('library');a=p.parse_args();result=run(a.library)
    Path('time_evolution_3d_report.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
