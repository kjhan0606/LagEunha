"""Compare one saved physical source with dense and vector/sparse roots."""
import cProfile,json,time,os,pstats
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic
from saved_curved_contact import checked_json
from simple_wave_evolution import install_domain
import kinetic_factor_rhs as kinetic


def run():
    cfg=json.loads(Path('vector_root_config.json').read_text());install_domain()
    source=checked_json(cfg['initial'],cfg['initial_sha256']);saved=checked_json(cfg['rhs'],cfg['rhs_sha256'])
    mass=np.asarray(source['mass']);state=np.asarray(source['state']);L=[np.asarray(v) for v in source['factors']]
    rows=[];results=[]
    for mode in ('0','1'):
        os.environ['LAGFORCE_VECTOR_KINETIC_ROOT']=mode
        profiler=cProfile.Profile();started=time.monotonic()
        dy,diag=profiler.runcall(kinetic.rhs,0.,state,mass,covariance_factors=L,face_order=4)
        elapsed=time.monotonic()-started;stats=pstats.Stats(profiler).stats
        top=sorted([dict(file=key[0],function=key[2],calls=v[1],self_seconds=v[2],cumulative_seconds=v[3]) for key,v in stats.items()],key=lambda v:v['cumulative_seconds'],reverse=True)[:25]
        rows.append(dict(backend=diag['mass_root_backend'],profiled_seconds=elapsed,mass_flux=diag['mass_flux'],top=top));results.append((dy,diag))
        print('root benchmark',rows[-1]['backend'],elapsed,diag['mass_flux'],flush=True)
    old,do=results[0];new,dn=results[1];N=len(mass)
    exact=float(np.max(abs(old-np.asarray(saved['dy']))))
    difference=abs(new-old);a=difference[:-4].reshape(N,20)
    physical=np.delete(a,2,axis=1);gram=max(float(np.max(abs(x@x.T-y@y.T))) for x,y in zip(do['cov_Qfactor'],dn['cov_Qfactor']))
    passed=exact==0. and np.max(physical)<1e-10 and np.max(a[:,2])<1e-9 and gram<1e-10 and max(do['mass_flux'],dn['mass_flux'])<2e-13
    report=dict(passed=bool(passed),legacy_saved_rhs_difference=exact,physical_rhs_difference=float(np.max(physical)),
        weight_rate_difference=float(np.max(a[:,2])),incoming_covariance_gram_difference=gram,
        ledger_rate_difference=float(np.max(difference[-4:])),rows=rows,profiled_speed_ratio=rows[0]['profiled_seconds']/rows[1]['profiled_seconds'],
        note='Same numerical density, quadrature, pressure force, positive factors, mass equations and root tolerance. Linear algebra/summation changes are measured; no bitwise solver equivalence claimed. Profiled single-source cost, not whole-trajectory timing.')
    atomic('vector_root_benchmark.json',report);print({k:report[k] for k in ('passed','legacy_saved_rhs_difference','physical_rhs_difference','weight_rate_difference','profiled_speed_ratio')},flush=True)
    raise SystemExit(0 if passed else 2)


if __name__=='__main__':run()
