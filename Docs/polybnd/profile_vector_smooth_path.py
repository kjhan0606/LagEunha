"""Read-only comparison of identical scalar/cell-batched dense path fluxes."""
import hashlib,json,os,time
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic
from physical_field_accuracy_audit import environment
from simple_wave_evolution import install_domain
from smooth_mass_path import sample_record


def run():
    install_domain();cfg=json.loads(Path('vector_smooth_path_profile_config.json').read_text());rows=[]
    for entry in cfg['cases']:
        model=json.loads(Path(entry['config']).read_text());os.environ.update(environment(model));os.environ.update(LAGFORCE_KINETIC_SMOOTHING=str(model['kinetic_smoothing_coefficient']),LAGFORCE_KINETIC_SMOOTHING_LENGTH=str(model['kinetic_smoothing_length']))
        mass=np.asarray(json.loads(Path(entry['initial']).read_text())['mass']);record=json.loads(Path(entry['trial']).read_text())['record'];data=[]
        for s in [0.,.37,1.]:
            values=[];seconds=[]
            for flag in ['0','1']:
                os.environ['LAGFORCE_VECTOR_SMOOTH_PATH']=flag;start=time.monotonic();v=sample_record(record,mass,s);seconds.append(time.monotonic()-start);values.append(v)
            difference=float(np.max(abs(np.asarray(values[0]['flux'])-np.asarray(values[1]['flux']))))
            if difference>2e-12 or values[0]['diag']!=values[1]['diag']:raise RuntimeError('vector path physical equivalence failed')
            data.append(dict(s=s,flux_maximum_absolute_difference=difference,scalar_seconds=seconds[0],cell_batched_seconds=seconds[1],speedup=seconds[0]/seconds[1]))
        rows.append(dict(label=entry['label'],source_time=record['t'],source_dt=record['dt'],n=model['n'],samples=data,trial_sha256=hashlib.sha256(Path(entry['trial']).read_bytes()).hexdigest()))
        print(entry['label'],data,flush=True)
    atomic('vector_smooth_path_profile_report.json',dict(completed=True,rows=rows,note='Same continuous moment density, velocity covariance allocation, positive samples and face quadrature. Only group evaluations and summations by cell/edge. Actual saved dense path input; source time/primary arrays unchanged. One observed timing per sample, not a general benchmark. Default remains scalar.'))


if __name__=='__main__':run()
