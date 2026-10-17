"""Separate moment-rounding noise from root accuracy on one saved path."""
import json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import numpy as np
from adaptive_coupled_mass import atomic,restore
from geometry_mass_lift import integrate,lifted_state,rates
from smooth_mass_path import sample_record_file,initialize_worker
from smooth_path_quadrature import node_indices,assess_values
from predictive_evolution_batch import record_digest
from simple_wave_evolution import install_domain


def run():
    install_domain();cfg=json.loads(Path('smoothed_geometry_config.json').read_text());source=json.loads(Path(cfg['initial']).read_text());mass=np.asarray(source['mass'])
    record=restore(json.loads(Path(cfg['trial']).read_text())['record']);old_evaluations=record['geometry_lift']['evaluations'];noise=[]
    for mode in ['float64','longdouble']:
        current=dict(record)
        if mode=='longdouble':current['dense_moment_arithmetic']='longdouble'
        values=[]
        for s in [.5,.5+1e-8,.5-1e-8,.5+1e-10,.5-1e-10]:
            state,L,_=lifted_state(current,mass,s);g,d=rates(current['t']+current['dt']*s,state,L,mass);values.append(g[:,2])
        noise.append(dict(mode=mode,second_difference_1e8=float(np.max(abs(values[1]+values[2]-2*values[0]))),second_difference_1e10=float(np.max(abs(values[3]+values[4]-2*values[0])))))
    atomic('preliminary_noise.json',noise)
    current=dict(record);current.pop('geometry_lift');current['dense_moment_arithmetic']='longdouble';started=time.monotonic()
    current['geometry_lift']=integrate(current,mass,rtol=cfg.get('geometry_rtol',3e-14),atol=cfg.get('geometry_atol',1e-15),root_tolerance=cfg.get('geometry_mass_root_tolerance',2e-13));seconds=time.monotonic()-started
    path=Path('extended_trial.json');atomic(path,dict(record=current));digest=record_digest(current);cache={};history=[]
    budgets=json.loads(Path(cfg['config']).read_text())
    with ProcessPoolExecutor(max_workers=4,mp_context=mp.get_context('spawn'),initializer=initialize_worker) as executor:
        for count in [9,17,33,65,129,257]:
            indices=node_indices(count,1024);missing=[j for j in indices if j not in cache]
            result=list(executor.map(sample_record_file,[str(path.resolve())]*len(missing),[digest]*len(missing),[mass]*len(missing),[j/1024 for j in missing]))
            cache.update(zip(missing,result));atomic('samples.json',{str(k):v for k,v in cache.items()})
            a=assess_values([cache[j] for j in indices],mass,current['dt'],.05,budgets,True)
            row={k:a[k] for k in ['ratio','physical_ratio','mass_ratio','generator','quadrature_ratio']};row['nodes']=count;history.append(row);print(row,flush=True)
            if a['ratio']<=1:break
    report=dict(completed=True,noise=noise,original_geometry_evaluations=old_evaluations,new_geometry_evaluations=current['geometry_lift']['evaluations'],new_geometry_seconds=seconds,history=history,
        note='Same stored source, RK2 moment polynomial and factor coefficients. Extended intermediate arithmetic only; accepted states and stored source precision are not changed. Original approved coupled evolution not repeated.')
    from kinetic_speed_smoothing import parameters
    report['kinetic_smoothing_coefficient'],report['kinetic_smoothing_length']=parameters()
    report['note']='Conditional geometry comparison on the same stored fluid/factor path, with explicitly changed smooth kinetic splitting. Not a full PDE accuracy result. All budgets and independent polynomial-derivative audits retained.'
    atomic('smoothed_geometry_report.json',report);print(report,flush=True)


if __name__=='__main__':run()
