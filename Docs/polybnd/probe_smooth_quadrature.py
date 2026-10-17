"""Refine a saved rejected path; reuse its stages and all nine old samples."""
import json,os,multiprocessing as mp
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor,as_completed
import numpy as np
from adaptive_coupled_mass import atomic,simpson
from saved_curved_contact import checked_json
from predictive_evolution_batch import record_digest
from smooth_mass_path import sample_record_file,initialize_worker


def run():
    cfg=json.loads(Path('smooth_quadrature_config.json').read_text())
    source=checked_json(cfg['initial'],cfg['initial_sha256']);mass=np.asarray(source['mass'])
    trial=checked_json(cfg['trial'],cfg['trial_sha256']);old=checked_json(cfg['samples'],cfg['samples_sha256'])
    cache_path=Path('refined_samples.json')
    samples=json.loads(cache_path.read_text()) if cache_path.exists() else {str(int(j)*16):v for j,v in old.items()}
    record=trial['record'];digest=record_digest(record);dt=record['dt'];rows=[]
    with ProcessPoolExecutor(max_workers=8,mp_context=mp.get_context('spawn'),initializer=initialize_worker) as executor:
        for count in (9,17,33,65,129):
            indices=list(range(0,129,128//(count-1)));missing=[j for j in indices if str(j) not in samples]
            futures={executor.submit(sample_record_file,cfg['trial'],digest,mass,j/128):j for j in missing}
            for f in as_completed(futures):samples[str(futures[f])]=f.result();atomic(cache_path,samples)
            flux=np.asarray([samples[str(j)]['flux'] for j in indices]);signed=simpson(flux,dt);absolute=simpson(abs(flux),dt)
            delta=signed-simpson(flux[::2],dt);absdelta=absolute-simpson(abs(flux[::2]),dt)
            row=dict(nodes=count,absolute_ratio=float(np.max(absolute/mass))/(1e-7*dt/.05),
                signed_quadrature_ratio=float(np.max(abs(delta)/mass))/(1e-9*dt/.05),
                absolute_quadrature_ratio=float(np.max(abs(absdelta)/mass))/(1e-9*dt/.05),
                signed=signed,absolute=absolute,signed_delta=delta,absolute_delta=absdelta)
            rows.append(row);atomic('smooth_quadrature_report.json',dict(completed=count==129,rows=rows,
                source_time=record['t'],dt=dt,reused_source_stages=True,reused_nine_samples=True))
            print({k:row[k] for k in ('nodes','absolute_ratio','signed_quadrature_ratio','absolute_quadrature_ratio')},flush=True)


if __name__=='__main__':run()
