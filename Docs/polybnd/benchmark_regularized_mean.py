"""Check a covariance-rank regularity hypothesis on a stored fluid trial."""
import os,json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import numpy as np
from adaptive_coupled_mass import atomic,restore
from geometry_mass_lift import raw_path,integrate
from smooth_mass_path import sample_record_file,initialize_worker
from smooth_path_quadrature import node_indices,assess_values
from predictive_evolution_batch import record_digest
from simple_wave_evolution import install_domain
from quadratic_velocity_reconstruction import prepare,correction
import moment_ale_evolve as model


def run():
    install_domain();cfg=json.loads(Path('regularized_mean_config.json').read_text());source=json.loads(Path(cfg['initial']).read_text());mass=np.asarray(source['mass'])
    record=restore(json.loads(Path(cfg['trial']).read_text())['record']);old=record.pop('geometry_lift',None)
    quiet=np.asarray([np.linalg.norm(v)==0 for v in record['L']]);rows=[]
    for mode in ['projected','projected_regularized']:
        os.environ['LAGFORCE_QUADRATIC_BUDGET']=mode;values=[]
        for s in [1e-8,1e-6,1e-4,.01]:
            state,L,_=raw_path(record,mass,s);a=state[:-4].reshape(-1,20);x=np.asarray(a[:,:2],float);weights=np.asarray(a[:,2],float)
            cells,edges=model.geometry(x,weights,record['t']+record['dt']*s);states=[model.decode(p,np.r_[m,r[3:]]) for m,r,p in zip(mass,a,cells)]
            prepare(states,cells,edges,mass,L);added=np.asarray([np.linalg.norm(correction(st,p)[0]) for st,p in zip(states,x)])
            values.append(dict(s=s,maximum_quiet_correction=float(added[quiet].max()) if np.any(quiet) else 0.,maximum_correction=float(added.max())))
        rows.append(dict(mode=mode,samples=values));atomic('mean_rank_samples.json',rows);print(rows[-1],flush=True)
    os.environ['LAGFORCE_QUADRATIC_BUDGET']='projected_regularized';started=time.monotonic()
    record['geometry_lift']=integrate(record,mass,rtol=1e-12,atol=1e-14,root_tolerance=1e-16);seconds=time.monotonic()-started
    path=Path('regularized_trial.json');atomic(path,dict(record=record));digest=record_digest(record)
    budgets=json.loads(Path(cfg['config']).read_text());cache={};history=[]
    with ProcessPoolExecutor(max_workers=8,mp_context=mp.get_context('spawn'),initializer=initialize_worker) as executor:
        for count in [9,17,33,65,129,257]:
            ids=node_indices(count,1024);missing=[j for j in ids if j not in cache]
            result=list(executor.map(sample_record_file,[str(path.resolve())]*len(missing),[digest]*len(missing),[mass]*len(missing),[j/1024 for j in missing]))
            cache.update(zip(missing,result));atomic('samples.json',{str(k):v for k,v in cache.items()})
            checked=assess_values([cache[j] for j in ids],mass,record['dt'],budgets['end']-source['t'],budgets,True)
            row={k:checked[k] for k in ['ratio','physical_ratio','mass_ratio','generator','quadrature_ratio']};row['nodes']=count;history.append(row);print(row,flush=True)
            if checked['ratio']<=1:break
    report=dict(completed=True,rank_samples=rows,old_geometry_evaluations=old['evaluations'] if old else None,new_geometry_evaluations=record['geometry_lift']['evaluations'],seconds=seconds,history=history,
        note='Conditional comparison only: same saved primary moment/factor trial, different fixed64 projected-gradient normalization. No approved coupled evolution repeated.')
    atomic('regularized_mean_report.json',report);print(report,flush=True)


if __name__=='__main__':run()
