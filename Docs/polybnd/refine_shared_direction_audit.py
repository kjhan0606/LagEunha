"""Independent257-node plus focused65-node audit of completed trial paths."""
import json,time
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import multiprocessing as mp
import numpy as np
from adaptive_coupled_mass import atomic,simpson,budget_ratio
from predictive_evolution_batch import initialize_case,initialize_path_worker,path_sample
from predictive_boundary_density import BoundaryDensity


def inspect(mode,root):
    root=Path(root);source=json.loads((root/'initial_config_checkpoint.json').read_text());cfg=source['config'];mass,_=initialize_case(cfg)
    trial=json.loads((root/(mode+'_trial.json')).read_text());record=trial['record'];dt=record['dt'];options=trial['options']
    saved=json.loads((root/(mode+'_nodes.json')).read_text());samples={i*4:saved[str(i)] for i in range(65)};started=time.monotonic()
    with ProcessPoolExecutor(max_workers=16,mp_context=mp.get_context('spawn'),initializer=initialize_path_worker,
            initargs=(cfg['n'],cfg['amplitude'],True,True,options)) as executor:
        futures={executor.submit(path_sample,record,mass,i/4):i for i in range(257) if i not in samples}
        for f in as_completed(futures):samples[futures[f]]=f.result()
        atomic(mode+'_nodes257.json',samples)
        # Known stalled time is in the first1/64 of this larger trial.
        focus={i*16:samples[i] for i in range(5)}
        futures={executor.submit(path_sample,record,mass,i/64):i for i in range(65) if i not in focus}
        for f in as_completed(futures):focus[futures[f]]=f.result()
        atomic(mode+'_focus65.json',focus)
    flux=np.array([samples[i]['mass_rate'] for i in range(257)]);fflux=np.array([focus[i]['mass_rate'] for i in range(65)])
    signed=simpson(flux,dt);absolute=simpson(abs(flux),dt)
    first_signed=simpson(flux[:5],dt/64);first_absolute=simpson(abs(flux[:5]),dt/64)
    focused_signed=simpson(fflux,dt/64);focused_absolute=simpson(abs(fflux),dt/64)
    corrected_signed=signed-first_signed+focused_signed;corrected_absolute=absolute-first_absolute+focused_absolute
    # Conservative sum of nested-grid signed changes, including replacing the
    # initial segment by focused nodes. This retains the original1e-12 allocation.
    qdelta=abs(signed-simpson(flux[::2],dt))+abs(signed-simpson(flux[::4],dt))
    qdelta+=abs(focused_signed-first_signed)+abs(focused_signed-simpson(fflux[::2],dt/64))
    all_samples=list(samples.values())+list(focus.values())
    peak=1e-4*max(float(np.max(abs(np.asarray(s['mass_rate']))/mass)) for s in all_samples)
    generator=max(s['generator'] for s in all_samples);pressure=min(s['pressure'] for s in all_samples)
    factor_R=min(s['factor_R'] for s in all_samples)
    ratio=budget_ratio(peak,float(np.max(corrected_absolute/mass)),float(np.max(qdelta/mass)),generator,dt,1e-4,5e-11)
    row=dict(mode=mode,passed=bool(ratio<=1 and pressure>0 and factor_R>=-1e-20),source_time=record['t'],dt=dt,
        whole_nodes=257,focus_nodes=65,focus_time_width=dt/64,ratio=ratio,normalized_peak_mass_rate=peak,
        integrated_absolute_mass_defect=float(np.max(corrected_absolute/mass)),quadrature_error_screen=float(np.max(qdelta/mass)),
        signed=corrected_signed,absolute=corrected_absolute,quad_delta=qdelta,seconds=time.monotonic()-started)
    atomic(mode+'_refined_report.json',row);print('refined path',mode,row['passed'],ratio,flush=True);return row


def run():
    cfg=json.loads(Path('refined_config.json').read_text());root=Path(cfg['source'])
    original=json.loads((root/'shared_direction_report.json').read_text())
    if not original['completed']:raise RuntimeError('original trials incomplete')
    with ProcessPoolExecutor(max_workers=2,mp_context=mp.get_context('spawn')) as executor:
        futures=[executor.submit(inspect,mode,str(root)) for mode in ('moment','global_moment')]
        rows=[f.result() for f in futures]
    report=dict(completed=True,passed=all(r['passed'] for r in rows),results=rows,
        note='Prescribed rates on completed saved trials, no weight re-solve. Same physical gates plus nested257/129/65 and targeted65/33 replacement in the known stalled-time segment. Sampled evidence, not a continuous-path proof.')
    atomic('refined_path_report.json',report);raise SystemExit(0 if report['passed'] else 2)


if __name__=='__main__':run()
