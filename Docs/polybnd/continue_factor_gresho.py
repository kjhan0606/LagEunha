"""Resume verified corrected-loss checkpoints; keep cumulative boundary ledgers."""
import json,time,hashlib
from pathlib import Path
import numpy as np
from compatible_2d_suite import initial,compare
from adaptive_moment_batch import diagnostic
from kinetic_factor_rhs import rhs
from factor_moment_evolve import rebuild,advance_factors_sylvester as advance
SOURCES={4:Path('/gpfs/kjhan/LagForce/Codex_B_sylvester_2d_20261010_004842/gresho_n4_dt0.001.checkpoint.json'),6:Path('/gpfs/kjhan/LagForce/Codex_D_finiteR_20261010_010542/corrected_loss/gresho_n6_dt0.001.checkpoint.json')}

def run(n,dt,end=.1):
    path=SOURCES[n];raw=path.read_bytes();z=json.loads(raw);start=z['t'];mass=np.array(z['mass']);y=np.array(z['state']);L=[np.array(x) for x in z['factors']]
    m0,y0=initial('gresho',n) # Set the same domain and construct the original ledger reference only.
    np.testing.assert_allclose(mass,m0,rtol=0,atol=0)
    tag=f'gresho_n{n}_dt{dt:g}';history=[];wall=time.monotonic();t=start;max_flux=0.;calls=0
    def save(suffix):
        Path(tag+suffix).write_text(json.dumps(dict(t=t,mass=mass.tolist(),state=y.tolist(),factors=[a.tolist() for a in L],source=str(path),source_sha256=hashlib.sha256(raw).hexdigest()),indent=2))
    try:
        for k in range(round((end-start)/dt)):
            t=start+k*dt
            dy,d=rhs(t,y,mass,covariance_factors=L);calls+=1
            halfL=advance(L,d,dt/2);half=rebuild(y+dt*dy/2,mass,halfL)
            dh,dd=rhs(t+dt/2,half,mass,covariance_factors=halfL);calls+=1
            nextL=advance(L,dd,dt);ny=rebuild(y+dt*dh,mass,nextL)
            mid=diagnostic(t+dt/2,half,mass,y0);last=diagnostic(t+dt,ny,mass,y0)
            minR=min(float(np.linalg.eigvalsh(a@a.T).min()/m) for a,m in zip(halfL+nextL,list(mass)*2))
            max_flux=max(max_flux,d['mass_flux'],dd['mass_flux'])
            if minR < -1e-20 or min(mid['min_pressure'],last['min_pressure'])<=0 or max_flux>1e-10 or last['ledger_error']>1e-10:raise RuntimeError('continuation physics gate failed')
            history.append(dict(t=t+dt,min_factor_R=minR,min_pressure=min(mid['min_pressure'],last['min_pressure']),ledger_error=last['ledger_error']))
            y,L=ny,nextL;t=start+(k+1)*dt
        save('.checkpoint.json');result=dict(status='completed',errors=compare('gresho',t,y,mass),final=diagnostic(t,y,mass,y0),min_factor_R=min(a['min_factor_R'] for a in history),min_pressure=min(a['min_pressure'] for a in history))
    except Exception as e:
        # y,L still refer to the last accepted state; time is the attempted step start.
        save('.last_accepted.json');result=dict(status='failed',error=repr(e),last_accepted_time=t)
    Path(tag+'.history.json').write_text(json.dumps(history,indent=2))
    result.update(n=n,dt=dt,start=start,end=end,source=str(path),source_sha256=hashlib.sha256(raw).hexdigest(),max_mass_flux=max_flux,calls=calls,seconds=time.monotonic()-wall)
    return result
if __name__=='__main__':
    for n in (4,6):
        for dt in (.001,.0005):print(json.dumps(run(n,dt)),flush=True)
