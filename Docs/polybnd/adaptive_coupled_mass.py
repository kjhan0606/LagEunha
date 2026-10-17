"""Adaptive coupled midpoint/factor steps with a sampled mass-defect budget.

Evolution always uses material sites and the current-state mass solve. The
diagnostic quadratic path is accepted only if both its rate and integrated
absolute mass defect meet local allocations. A separate 65-node pass audits
all saved steps. Neither sample grid supplies a rigorous continuous bound.
"""
import argparse
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from adaptive_moment_batch import setup,diagnostic
from factor_moment_evolve import rebuild,advance_factors_sylvester as advance
from kinetic_factor_rhs import rhs
from coupled_factor_path_audit import dense_state


def plain(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,dict):return {k:plain(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [plain(v) for v in value]
    if isinstance(value,np.generic):return value.item()
    return value


def atomic(path,data):
    path=Path(path);temp=path.with_suffix(path.suffix+'.tmp')
    temp.write_text(json.dumps(plain(data)));temp.replace(path)


def simpson(values,width):
    if len(values)<3 or len(values)%2!=1:raise ValueError('odd node count required')
    a=np.asarray(values)
    return width/(3*(len(a)-1))*(a[0]+a[-1]+4*a[1:-1:2].sum(0)+2*a[2:-1:2].sum(0))


def budget_ratio(rate,absolute_defect,quad_difference,generator_defect,dt,end,
                 mass_budget=5e-11):
    if not 0<dt<=end or mass_budget<=0:raise ValueError('invalid budget allocation')
    return max(rate/5e-9,absolute_defect/(mass_budget*dt/end),
               quad_difference/(1e-12*dt/end),generator_defect/5e-9)


def restore(record):
    r=dict(record)
    for k in ('y','d0','dm'):r[k]=np.asarray(r[k])
    r['L']=[np.asarray(x) for x in r['L']]
    r['diag']={k:[np.asarray(x) for x in v] for k,v in r['diag'].items()}
    return r


def assess(record,mass,end,nodes):
    r=restore(record);N=len(mass);rates=[];peak=0.;gen=0.;pressure=float('inf');gram=0.
    worst_s=0.
    for s in np.linspace(0.,1.,nodes):
        state,L,dy=dense_state(r['y'],r['L'],r['d0'],r['dm'],r['diag'],mass,r['dt'],s)
        velocity=dy[:-4].reshape(N,20)
        _,d=rhs(r['t']+s*r['dt'],state,mass,covariance_factors=L,
                audit_geometry_rates=(velocity[:,:2],velocity[:,2]))
        f=d['mass_rates']
        if not np.all(np.isfinite(f)):raise RuntimeError('nonfinite mass residual')
        rate=end*float(np.max(abs(f)/mass))
        if rate>peak:peak=rate;worst_s=(r['t']+s*r['dt'])/end
        gen=max(gen,d['generator_velocity_defect']);pressure=min(pressure,d['min_pressure'])
        gram=min(gram,min(float(np.linalg.eigvalsh(a@a.T).min()/m) for a,m in zip(L,mass)))
        rates.append(f)
    signed=simpson(rates,r['dt']);coarse=simpson(rates[::2],r['dt'])
    absolute=simpson(np.abs(rates),r['dt'])
    return dict(rate=peak,generator_defect=gen,min_pressure=pressure,min_factor_R=gram,
                signed=signed,absolute=absolute,quad_delta=signed-coarse,worst_s=worst_s)


def trial(t,y,L,mass,dt):
    d0,diag0=rhs(t,y,mass,covariance_factors=L)
    halfL=advance(L,diag0,dt/2);half=rebuild(y+dt*d0/2,mass,halfL)
    dm,diag=rhs(t+dt/2,half,mass,covariance_factors=halfL)
    nextL=advance(L,diag,dt);nexty=rebuild(y+dt*dm,mass,nextL)
    keys=('cov_linearization_factors','cov_A','cov_K','cov_Qfactor')
    record=dict(t=t,dt=dt,y=y,L=L,d0=d0,dm=dm,diag={k:diag[k] for k in keys})
    return record,nexty,nextL


def configuration(kind,n,amplitude,end,budget):
    if kind not in ('wave','shear','contact','moving_contact','shear_contact') or not isinstance(n,int) or n<3:
        raise ValueError('unsupported case or integer n<3')
    if not np.all(np.isfinite([amplitude,end,budget])) or end<=0 or budget<=0:
        raise ValueError('finite parameters and positive duration/budget required')
    if amplitude<0 or (kind=='wave' and amplitude>=.6):
        raise ValueError('wave amplitude must keep initial pressure positive')
    if kind=='shear_contact' and amplitude<1:
        raise ValueError('shear contact density ratio must be >=1')
    return dict(kind=kind,end=end,n=n,amplitude=amplitude,budget=budget,
                controller_nodes=17,audit_nodes=65,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())


def run(kind,reference=None,*,n=3,amplitude=.01,end=.05,budget=5e-11):
    config=configuration(kind,n,amplitude,end,budget)
    root=Path(kind);root.mkdir(exist_ok=True);stepsdir=root/'steps';stepsdir.mkdir(exist_ok=True)
    mass,y0=setup(kind,n,amplitude);L0=[np.zeros((2,2)) for _ in mass];y0=rebuild(y0,mass,L0)
    cp=root/'checkpoint.json';start=time.monotonic()
    if cp.exists():
        state=json.loads(cp.read_text())
        if state['config']!=config:raise RuntimeError('checkpoint configuration mismatch')
        y=np.asarray(state['state']);L=[np.asarray(a) for a in state['factors']]
        t=state['t'];dt=state['next_dt'];count=state['steps'];rejected=state['rejected']
    else:
        y=y0.copy();L=L0;t=0.;dt=end/32;count=0;rejected=0
    while t<end:
        dt=min(dt,end-t)
        if dt<end*1e-9:raise RuntimeError('step size exhausted')
        record,new,newL=trial(t,y,L,mass,dt)
        d=assess(record,mass,end,17)
        metric=budget_ratio(d['rate'],float(np.max(d['absolute']/mass)),
            float(np.max(abs(d['quad_delta'])/mass)),d['generator_defect'],dt,end,budget)
        hn=diagnostic(t+dt,new,mass,y0)
        if not np.isfinite(metric) or d['min_pressure']<=0 or d['min_factor_R']< -1e-20 or hn['ledger_error']>=1e-10:
            raise RuntimeError('trial state/ledger/finite gate failed')
        if metric>1:
            rejected+=1
            with (root/'rejections.jsonl').open('a') as f:f.write(json.dumps(dict(t=t,dt=dt,ratio=metric))+'\n')
            dt/=2
            continue
        count+=1;record['accepted_ledger_error']=hn['ledger_error'];record['control']=plain(d)
        atomic(stepsdir/f'{count:05d}.json',record)
        y,L=new,newL;t=min(end,t+dt)
        factor=float(np.clip(.85/np.sqrt(max(metric,1e-12)),.5,1.5))
        dt=min(end/16,dt*factor)
        atomic(cp,dict(config=config,t=t,next_dt=dt,steps=count,rejected=rejected,
                       mass=mass,state=y,factors=L))
        if count%25==0 or t==end:
            print(json.dumps(dict(kind=kind,stage='evolving',steps=count,t=t,rejected=rejected)),flush=True)
    # Resume the independent audit by step number without repeating evolution.
    auditcp=root/'audit_checkpoint.json'
    if auditcp.exists():
        accum=json.loads(auditcp.read_text())
        if accum['config']!=config or accum['total_steps']!=count:raise RuntimeError('audit checkpoint mismatch')
    else:
        accum=dict(config=config,total_steps=count,done=0,signed=[0.]*len(mass),absolute=[0.]*len(mass),
            quad_delta=[0.]*len(mass),max_rate=0.,max_generator_defect=0.,min_pressure=1e300,
            min_factor_R=0.,max_ledger_error=0.,max_prefix_mass=0.)
    for index in range(accum['done']+1,count+1):
        record=json.loads((stepsdir/f'{index:05d}.json').read_text());d=assess(record,mass,end,65)
        for key in ('signed','absolute','quad_delta'):accum[key]=(np.asarray(accum[key])+d[key]).tolist()
        accum['max_rate']=max(accum['max_rate'],d['rate'])
        accum['max_generator_defect']=max(accum['max_generator_defect'],d['generator_defect'])
        accum['min_pressure']=min(accum['min_pressure'],d['min_pressure'])
        accum['min_factor_R']=min(accum['min_factor_R'],d['min_factor_R'])
        accum['max_ledger_error']=max(accum['max_ledger_error'],record['accepted_ledger_error'])
        accum['max_prefix_mass']=max(accum['max_prefix_mass'],float(np.max(abs(np.asarray(accum['signed']))/mass)))
        accum['done']=index
        atomic(auditcp,accum)
        if index%50==0 or index==count:
            print(json.dumps(dict(kind=kind,stage='independent_audit',done=index,total=count,max_rate=accum['max_rate'])),flush=True)
    result=dict(kind=kind,n=n,amplitude=amplitude,duration=end,mass_budget=budget,
        status='completed',steps=count,rejections=rejected,
        max_path_mass_rate=accum['max_rate'],max_generator_velocity_defect=accum['max_generator_defect'],
        integrated_relative_mass_defect=float(np.max(abs(np.asarray(accum['signed']))/mass)),
        integrated_absolute_mass_defect=float(np.max(np.asarray(accum['absolute'])/mass)),
        quadrature_mass_difference=float(np.max(abs(np.asarray(accum['quad_delta']))/mass)),
        max_prefix_mass_defect=accum['max_prefix_mass'],min_pressure=accum['min_pressure'],
        min_factor_R=accum['min_factor_R'],max_ledger_error=accum['max_ledger_error'],
        independent_samples=count*65,seconds_this_invocation=time.monotonic()-start)
    result['state_gate']=bool(result['min_pressure']>0 and result['min_factor_R']>=-1e-20 and result['max_ledger_error']<1e-10)
    result['mass_path_gate']=bool(result['max_path_mass_rate']<1e-8 and
        result['integrated_absolute_mass_defect']<1e-10 and result['quadrature_mass_difference']<1e-11)
    result['generator_path_gate']=bool(result['max_generator_velocity_defect']<1e-8)
    if kind=='wave' and reference is not None:
        if (n,amplitude,end)!=(3,.01,.05):
            raise ValueError('legacy reference is only valid for N3 amplitude .01 duration .05')
        ref=np.asarray(json.loads(reference.read_text())['state'])[:-4].reshape(n*n,20)
        actual=y[:-4].reshape(n*n,20)
        result['reference_state_max_difference']=float(np.max(abs(actual-ref)))
    else:
        result['reference_comparison']='not supplied for this case; conservation/path gates only'
    atomic(root/'result.json',result)
    print(json.dumps(result),flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('predecessor');args=parser.parse_args()
    previous=Path(args.predecessor)
    summary=json.loads((previous/'refinement_summary.json').read_text())
    passes=[r for r in summary['results'] if r.get('status')=='completed' and r.get('state_gate') and r.get('mass_path_gate')]
    if not passes:
        atomic('gate.json',dict(status='blocked_by_predecessor',reason='no passing evolving-state reference'))
        raise SystemExit(2)
    chosen=passes[-1];source=previous/f"steps_{chosen['steps']}"/'checkpoint.json'
    reference=Path('reference.checkpoint.json')
    if not reference.exists():reference.write_bytes(source.read_bytes())
    if reference.read_bytes()!=source.read_bytes():raise RuntimeError('reference changed')
    atomic('gate.json',dict(status='passed',reference_steps=chosen['steps'],source=str(source),
                          reference_sha256=hashlib.sha256(reference.read_bytes()).hexdigest()))
    results=[run(kind,reference) for kind in ('shear','wave')]
    atomic('summary.json',dict(results=results))
