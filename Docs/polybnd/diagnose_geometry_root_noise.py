"""Read-only conditioning and root-tolerance check on a saved lifted path."""
import os,json
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic,restore
from simple_wave_evolution import install_domain
from geometry_mass_lift import lifted_state
import geometry_mass_lift as geometry
import vector_kinetic_root as vector


def run():
    install_domain();cfg=json.loads(Path('root_noise_config.json').read_text());source=json.loads(Path(cfg['initial']).read_text());mass=np.asarray(source['mass'])
    record=restore(json.loads(Path(cfg['trial']).read_text())['record']);rows=[];original=geometry.solve_arrays
    for s in [.5,.5+1e-8,.5-1e-8,.5+1e-10,.5-1e-10]:
        box={}
        def capture(data,N):
            z,values,norm,iterations=original(data,N)
            residual,L,_=vector.evaluate(data,N,z);eigenvalues=np.linalg.eigvalsh(L.toarray()[:-1,:-1])
            tight,tv,tn,ti=original(data,N,tolerance=1e-16)
            box.update(default_norm=norm,default_iterations=iterations,tight_norm=tn,tight_iterations=ti,
                maximum_rate_change=float(np.max(abs(tight-z))),minimum_gauge_eigenvalue=float(eigenvalues.min()),maximum_gauge_eigenvalue=float(eigenvalues.max()),
                weight_rate=z,tight_weight_rate=tight)
            return z,values,norm,iterations
        geometry.solve_arrays=capture
        state,L,_=lifted_state(record,mass,s);rates,diag=geometry.rates(record['t']+record['dt']*s,state,L,mass)
        box.update(s=s);rows.append(box)
    geometry.solve_arrays=original
    result=dict(rows=rows,dt=record['dt'],source_time=record['t'],note='Same saved numerical fluid/geometry path, no time integration. Compare the original mass root tolerance with1e-16 and inspect its actual gauge Jacobian.')
    atomic('root_noise_report.json',result)
    for row in rows:print({k:v for k,v in row.items() if k not in ['weight_rate','tight_weight_rate']},flush=True)


if __name__=='__main__':run()
