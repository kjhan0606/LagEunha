"""Read-only equivalence and timings at actual saved physical states."""
import hashlib,json,os,time
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic,restore
from physical_field_accuracy_audit import environment
from simple_wave_evolution import install_domain
from kinetic_factor_rhs import rhs


def run():
    install_domain();cfg=json.loads(Path('vector_smooth_path_profile_config.json').read_text());rows=[]
    for entry in cfg['cases']:
        model=json.loads(Path(entry['config']).read_text());os.environ.update(environment(model));os.environ.update(LAGFORCE_VECTOR_KINETIC_ROOT='1',LAGFORCE_KINETIC_SMOOTHING=str(model['kinetic_smoothing_coefficient']),LAGFORCE_KINETIC_SMOOTHING_LENGTH=str(model['kinetic_smoothing_length']))
        mass=np.asarray(json.loads(Path(entry['initial']).read_text())['mass']);record=restore(json.loads(Path(entry['trial']).read_text())['record']);data=[]
        for flag in ['0','1']:
            os.environ['LAGFORCE_VECTOR_MOMENT_FLUX']=flag;start=time.monotonic();dy,diag=rhs(record['t'],record['y'],mass,covariance_factors=record['L']);data.append((dy,diag,time.monotonic()-start))
        old,do,oldtime=data[0];new,dn,newtime=data[1];difference=float(np.max(abs(old-new)));cov={}
        for key in ['cov_A','cov_K','cov_Qfactor']:
            cov[key]=max(float(np.max(abs(a-b))) for a,b in zip(do[key],dn[key]))
        if max([difference]+list(cov.values()))>2e-12:raise RuntimeError('batched moment flux equivalence failed')
        np.testing.assert_array_equal(old[:-4].reshape(-1,20)[:,:3],new[:-4].reshape(-1,20)[:,:3])
        if do['mass_flux']!=dn['mass_flux']:raise RuntimeError('moment batching changed mass root')
        row=dict(label=entry['label'],n=model['n'],t=record['t'],scalar_seconds=oldtime,batched_seconds=newtime,speedup=oldtime/newtime,rhs_maximum_absolute_difference=difference,covariance_factor_rate_maximum_differences=cov,mass_root_and_generator_rates_exact=True,trial_sha256=hashlib.sha256(Path(entry['trial']).read_bytes()).hexdigest());rows.append(row);print(row,flush=True)
    atomic('vector_moment_flux_profile_report.json',dict(completed=True,rows=rows,note='Face states/HLLC/mass root unchanged. Shared moment flux products and positive incoming covariance blocks batched across the same records, preserving record/receiver column order. No physical interval advanced or repeated. One observed timing per case; default remains scalar.'))


if __name__=='__main__':run()
