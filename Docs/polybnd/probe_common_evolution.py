"""Reproduce rejected common-contact path nodes without changing their rates."""
import json
from pathlib import Path
import numpy as np
from shapely.geometry import Polygon
from adaptive_coupled_mass import atomic,restore
from coupled_factor_path_audit import dense_state
from predictive_evolution_batch import initialize_case,density_options
from predictive_boundary_density import BoundaryDensity
from predictive_mass_path import evaluate_path
import moment_ale_evolve as model


def run():
    config=json.loads(Path('probe_config.json').read_text());source=Path(config['source']);rows=[]
    for case,rejection,nodes in config['probes']:
        cp=json.loads((source/'inputs'/case/'checkpoint.json').read_text());cfg=cp['config'];mass,_=initialize_case(cfg)
        rejected=json.loads((source/'evolution'/case/'short'/('rejection_'+str(rejection)+'.json')).read_text())
        record=restore(rejected['record'])
        for node in nodes:
            state,L,dy=dense_state(record['y'],record['L'],record['d0'],record['dm'],record['diag'],mass,record['dt'],node/64)
            t=record['t']+record['dt']*node/64;a=state[:-4].reshape(-1,20);rates=dy[:-4].reshape(-1,20)
            polys,_=model.geometry(a[:,:2],a[:,2],t);means=np.array([m/Polygon(v).area for m,v in zip(mass,polys)])
            tau=32*np.finfo(float).eps*max(1.,float(np.max(abs(means))))
            gaps=np.minimum(means-means.min(),means.max()-means)
            variants=config.get('variants',[dict(mode='moment'),dict(mode='neighbor')])
            for variant in variants:
                options=density_options(cfg,True);options.update(variant);mode=options['mode'];provider=BoundaryDensity(**options)
                try:
                    residual,diag=evaluate_path(t,state,L,mass,(rates[:,:2],rates[:,2]),provider)
                    row=dict(case=case,rejection=rejection,node=node,t=t,mode=mode,options=options,means=means,gaps_over_tau=gaps/tau,
                        normalized_mass_rates=residual/mass,closure=diag,priors=provider.last_priors,fits=provider.last_fits)
                    print(case,rejection,node,mode,'massrate',float(np.max(abs(residual)/mass)),flush=True)
                except Exception as exc:row=dict(case=case,rejection=rejection,node=node,t=t,mode=mode,options=options,error=repr(exc))
                rows.append(row);atomic('probe_report.json',dict(completed=False,rows=rows))
    atomic('probe_report.json',dict(completed=True,rows=rows,
        note='Rejected states/rates preserved. Alternative prior is read-only diagnosis, not accepted evolution.'))


if __name__=='__main__':run()
