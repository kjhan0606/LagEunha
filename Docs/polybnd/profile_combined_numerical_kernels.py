"""Validate optional equivalent kernels on immutable actual physical paths."""
import cProfile,hashlib,io,json,os,pstats,time
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic,restore
from physical_field_accuracy_audit import environment
from simple_wave_evolution import install_domain
from smooth_mass_path import sample_record
from kinetic_factor_rhs import rhs


def flags(fast):
    value='1' if fast else '0'
    return dict(LAGFORCE_VECTOR_MOMENT_FLUX=value,LAGFORCE_VECTOR_SMOOTH_PATH=value,LAGFORCE_SEEDED_POWER_CELLS=value,LAGFORCE_FAST_PRESSURE_POLYNOMIAL=value)


def covariance_rates(diag,factors):
    rows=[]
    for A,K,F,L in zip(diag['cov_A'],diag['cov_K'],diag['cov_Qfactor'],factors):
        R=L@L.T;rows.append(A@R+R@A.T-L@K@L.T+F@F.T)
    return np.asarray(rows)


def run():
    install_domain();cfg=json.loads(Path('vector_smooth_path_profile_config.json').read_text());rows=[]
    for entry in cfg['cases']:
        model=json.loads(Path(entry['config']).read_text());os.environ.update(environment(model));os.environ.update(LAGFORCE_VECTOR_KINETIC_ROOT='1',LAGFORCE_KINETIC_SMOOTHING=str(model['kinetic_smoothing_coefficient']),LAGFORCE_KINETIC_SMOOTHING_LENGTH=str(model['kinetic_smoothing_length']))
        mass=np.asarray(json.loads(Path(entry['initial']).read_text())['mass']);stored=json.loads(Path(entry['trial']).read_text())['record'];record=restore(stored);data=[]
        for fast in [False,True]:
            os.environ.update(flags(fast));start=time.monotonic();dy,diag=rhs(record['t'],record['y'],mass,covariance_factors=record['L']);data.append((dy,diag,time.monotonic()-start))
        old,do,oldtime=data[0];new,dn,newtime=data[1];difference=float(np.max(abs(old-new)));cov=float(np.max(abs(covariance_rates(do,record['L'])-covariance_rates(dn,record['L']))));samples=[]
        for s in [0.,.37,1.]:
            values=[];seconds=[]
            for fast in [False,True]:
                os.environ.update(flags(fast));start=time.monotonic();v=sample_record(stored,mass,s);seconds.append(time.monotonic()-start);values.append(v)
            diff=float(np.max(abs(np.asarray(values[0]['flux'])-np.asarray(values[1]['flux']))));samples.append(dict(s=s,flux_maximum_absolute_difference=diff,scalar_seconds=seconds[0],optimized_seconds=seconds[1],speedup=seconds[0]/seconds[1]))
        if max([difference,cov]+[v['flux_maximum_absolute_difference'] for v in samples])>=1e-12:raise RuntimeError('combined kernel physical equivalence failed')
        row=dict(label=entry['label'],n=model['n'],t=record['t'],scalar_rhs_seconds=oldtime,optimized_rhs_seconds=newtime,rhs_speedup=oldtime/newtime,rhs_maximum_absolute_difference=difference,physical_covariance_rate_maximum_absolute_difference=cov,path_samples=samples,trial_sha256=hashlib.sha256(Path(entry['trial']).read_bytes()).hexdigest());rows.append(row);print(row,flush=True)
        if model['n']==64:
            os.environ.update(flags(False));profile=cProfile.Profile();profile.runcall(sample_record,stored,mass,.37);profile.dump_stats('scalar_wave64_path.prof');stream=io.StringIO();pstats.Stats(profile,stream=stream).strip_dirs().sort_stats('cumtime').print_stats(35);Path('scalar_wave64_path_profile.txt').write_text(stream.getvalue())
    result=dict(completed=True,rows=rows,rhs_maximum_absolute_difference=max(v['rhs_maximum_absolute_difference'] for v in rows),independent_mass_flux_maximum_absolute_difference=max(s['flux_maximum_absolute_difference'] for v in rows for s in v['path_samples']),physical_covariance_rate_maximum_absolute_difference=max(v['physical_covariance_rate_maximum_absolute_difference'] for v in rows),note='Same mathematical model/budgets. Initial clipping intersection seeded by nearby sites; every site remains in global violation search and all-constraint verification. Same pressure polynomial/Horner, shared moment flux accumulation and cell-grouped independent path. Covariance compared in physical units, allowing incoming factor column permutations from equivalent polygon vertex/edge order. No physical interval advanced or repeated; timing samples are not a general benchmark.')
    atomic('combined_numerical_kernel_validation_report.json',result)


if __name__=='__main__':run()
