"""Validate certified reuse against fresh LPs along an actual accepted path."""
import json,time,cProfile,pstats
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic,restore
from adaptive_moment_batch import setup
from coupled_factor_path_audit import dense_state
from predictive_boundary_density import BoundaryDensity
from predictive_mass_path import evaluate_path


def run():
    record=restore(json.loads(Path('reference_step.json').read_text())['record'])
    source=json.loads(Path('inputs/n4_r100/checkpoint.json').read_text());cfg=source['config']
    setup('shear_contact',cfg['n'],cfg['amplitude']);mass=np.asarray(source['mass']);N=len(mass)
    import sys
    vector='--vector' in sys.argv
    fresh=BoundaryDensity();warm=BoundaryDensity(reuse_witness=not vector,vector_geometry=vector);rows=[]
    for s in (0.,.25,.5,.75,1.,.0625,.125):
        y,L,dy=dense_state(record['y'],record['L'],record['d0'],record['dm'],record['diag'],mass,record['dt'],s)
        rates=dy[:-4].reshape(N,20);args=(record['t']+s*record['dt'],y,L,mass,(rates[:,:2],rates[:,2]))
        started=time.monotonic();ff,df=evaluate_path(*args,fresh);cold_seconds=time.monotonic()-started
        started=time.monotonic();fw,dw=evaluate_path(*args,warm);warm_seconds=time.monotonic()-started
        flux_error=float(np.max(abs(ff-fw)));trace_error=max(float(np.max(abs(fresh.last_profiles[k]-warm.last_profiles[k]))) for k in fresh.last_profiles)
        bound_error=max(max(0.,w['baseline_bound']-f['baseline_bound']) for f,w in zip(fresh.last_fits,warm.last_fits))
        excess=max(max(0.,w['objective']-f['baseline_bound']-1e-5) for f,w in zip(fresh.last_fits,warm.last_fits))
        passed=flux_error<1e-10 and trace_error<1e-10 and bound_error<1e-8 and excess<1e-8 and dw['closure_moment_error']<1e-8 and dw['closure_min_density']>0
        row=dict(s=s,passed=passed,fresh_seconds=cold_seconds,warm_seconds=warm_seconds,lp_calls=dw['closure_lp_calls'],witness_reuses=dw['closure_witness_reuses'],
            flux_difference=flux_error,trace_difference=trace_error,dual_lower_bound_excess=bound_error,actual_volume_cap_excess=excess)
        rows.append(row);atomic('reuse_progress.json',dict(rows=rows));print(row,flush=True)
    report=dict(passed=all(r['passed'] for r in rows) and (vector or sum(r['witness_reuses'] for r in rows)>0),rows=rows,vector_geometry=vector,
        fresh_seconds=sum(r['fresh_seconds'] for r in rows),warm_seconds=sum(r['warm_seconds'] for r in rows))
    atomic('reuse_report.json',report)
    return report


if __name__=='__main__':
    profile=cProfile.Profile();profile.enable();result=run();profile.disable()
    profile.dump_stats('reuse_profile.stats')
    with open('reuse_profile.txt','w') as out:pstats.Stats(profile,stream=out).strip_dirs().sort_stats('cumulative').print_stats(25)
    raise SystemExit(0 if result['passed'] else 2)
