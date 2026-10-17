"""Choose auxiliary ODE tolerances using unchanged posterior path budgets."""
import json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import multiprocessing as mp
import numpy as np
from adaptive_coupled_mass import atomic,restore
from geometry_mass_lift import integrate
from smooth_mass_path import sample_record_file,initialize_worker
from smooth_path_quadrature import assess_values
from predictive_evolution_batch import record_digest
from simple_wave_evolution import install_domain


def run():
    install_domain();cfg=json.loads(Path('budget_geometry_config.json').read_text());source=json.loads(Path(cfg['initial']).read_text());mass=np.asarray(source['mass']);budgets=json.loads(Path(cfg['config']).read_text())
    original=restore(json.loads(Path(cfg['trial']).read_text())['record']);baseline=original.pop('geometry_lift',None);dt=original['dt'];duration=budgets['end']-source['t'];N=len(mass);rows=[]
    for divisor in [4,16,64,256]:
        record=dict(original);position_tol=budgets['generator_path_budget']*dt/divisor;weight_tol=2/np.sqrt(N)*position_tol;atol=np.tile([position_tol,position_tol,weight_tol],N)
        started=time.monotonic();record['geometry_lift']=integrate(record,mass,rtol=1e-12,atol=atol,root_tolerance=1e-16);seconds=time.monotonic()-started
        path=Path('trial_'+str(divisor)+'.json');atomic(path,dict(record=record));digest=record_digest(record);cache={};history=[]
        extra=set()
        for segment in record['geometry_lift']['segments']:
            extra.update([float(segment['t_old']/dt),float(segment['t']/dt),float((segment['t_old']+segment['t'])/(2*dt))])
        with ProcessPoolExecutor(max_workers=8,mp_context=mp.get_context('spawn'),initializer=initialize_worker) as executor:
            for count in [9,17,33,65,129,257]:
                phases=[float(j/(count-1)) for j in range(count)];missing=sorted((set(phases)|extra)-cache.keys())
                result=list(executor.map(sample_record_file,[str(path.resolve())]*len(missing),[digest]*len(missing),[mass]*len(missing),missing));cache.update(zip(missing,result))
                atomic('samples_'+str(divisor)+'.json',{float(k).hex():v for k,v in cache.items()})
                a=assess_values([cache[s] for s in phases],mass,dt,duration,budgets,True)
                peakgen=max(v['diag']['generator_defect'] for v in cache.values());peakrate=duration*max(float(np.max(abs(np.asarray(v['flux']))/mass)) for v in cache.values())
                a['generator']=peakgen;a['rate']=peakrate;a['physical_ratio']=max(a['mass_ratio'],peakgen/budgets['generator_path_budget'],peakrate/budgets['path_rate_budget']);a['ratio']=max(a['physical_ratio'],a['quadrature_ratio'])
                row={k:a[k] for k in ['ratio','physical_ratio','mass_ratio','generator','rate','quadrature_ratio']};row.update(nodes=count,actual_samples=len(cache));history.append(row);print(divisor,row,flush=True)
                if a['ratio']<=1 or peakgen>budgets['generator_path_budget'] or peakrate>budgets['path_rate_budget']:break
        result=dict(divisor=divisor,position_atol=position_tol,weight_atol=weight_tol,evaluations=record['geometry_lift']['evaluations'],seconds=seconds,history=history,passed=bool(a['ratio']<=1))
        rows.append(result);atomic('budget_geometry_progress.json',dict(rows=rows))
        if result['passed']:break
    atomic('budget_geometry_report.json',dict(completed=True,rows=rows,old_geometry_evaluations=baseline['evaluations'] if baseline else None,
        note='Identical saved moment/factor coefficients and unchanged mass/generator/rate/quadrature budgets. Only auxiliary ODE tolerances differ. Extra DOP segment endpoints/midpoints checked. No approved coupled evolution repeated.'))


if __name__=='__main__':run()
