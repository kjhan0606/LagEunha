"""Authorized second 2D bundle; reuse completed coarse runs, never restart them."""
import json,os
from pathlib import Path
import numpy as np
import factor_moment_evolve as short
import factor_2d_suite as spatial
from compatible_moment_rhs import rhs as donor_rhs
from kinetic_factor_rhs import rhs as kinetic_rhs

ROOT=Path.cwd()
PREVIOUS={
 'donor':Path('/gpfs/kjhan/LagForce/Codex_B_factor_moments_20261010_002319'),
 'kinetic':Path('/gpfs/kjhan/LagForce/Codex_B_kinetic_factors_20261010_003555')}

def emit(record):
    with (ROOT/'results.jsonl').open('a') as f:f.write(json.dumps(record)+'\n')
    print(json.dumps(record),flush=True)

def gate(folder):
    path=folder/'results.jsonl'
    if not path.exists():return False,'missing predecessor results'
    rows=[json.loads(s) for s in path.read_text().splitlines()]
    expected={(kind,a,dt) for kind,a in [('wave',.001),('wave',.1),('shear',0.),('contact',0.)] for dt in (.001,.0005)}
    coverage={(r.get('kind'),r.get('amplitude'),r.get('dt')) for r in rows}
    ok=coverage==expected and len(rows)==8 and all(r.get('status')=='completed' and r.get('min_factor_R',-1)>=-1e-20 and r.get('min_pressure',0)>0 and r.get('final',{}).get('ledger_error',1)<1e-10 for r in rows)
    return ok,rows

def convergence(mode,kind,a,folder):
    dt=[.001,.0005,.00025,.000125];states=[]
    for h in dt:
        parent=PREVIOUS[mode] if h>=.0005 else folder
        p=parent/f'{kind}_a{a:g}_dt{h:g}.checkpoint.json'
        if not p.exists():return dict(status='unavailable',missing=str(p))
        z=json.loads(p.read_text());assert z['t']==.02
        states.append(np.array(z['state'])[:-4].reshape(-1,20))
    groups={'sites':slice(0,2),'q':slice(3,5),'Q':slice(5,9),'p':slice(9,11),'D':slice(11,15),'Pi':slice(15,19),'E':slice(19,20)}
    result={}
    for name,sl in groups.items():
        diffs=[float(np.linalg.norm(x[:,sl]-y[:,sl])) for x,y in zip(states[:-1],states[1:])]
        order=[float(np.log2(x/y)) if min(x,y)>1e-13 else None for x,y in zip(diffs[:-1],diffs[1:])]
        result[name]={'successive_L2_differences':diffs,'observed_orders':order}
    return dict(status='evaluated',dt=dt,groups=result,note='Self convergence only; orders below roundoff cutoff omitted. No general order claim.')

def main():
    report={}
    for mode,rhs in [('donor',donor_rhs),('kinetic',kinetic_rhs)]:
        passed,detail=gate(PREVIOUS[mode]);report[mode]={'predecessor_passed':passed}
        if not passed:
            emit(dict(mode=mode,status='blocked_by_predecessor',detail=detail));continue
        folder=ROOT/mode;folder.mkdir();os.chdir(folder);short.rhs=rhs;spatial.rhs=rhs
        report[mode]['time_convergence']={};passed_new=True
        for kind,a in [('wave',.1),('shear',0.)]:
            for dt in (.00025,.000125):
                r=short.run(kind,a,dt);r['mode']=mode;r['suite']='time';emit(r)
                passed_new &= r.get('status')=='completed' and r.get('min_factor_R',-1)>=-1e-20 and r.get('min_pressure',0)>0 and r.get('max_mass_flux',1)<1e-10 and r.get('final',{}).get('ledger_error',1)<1e-10
            report[mode]['time_convergence'][kind]=convergence(mode,kind,a,folder)
        if mode=='kinetic' and passed_new:
            for kind,end in [('wave',.05),('gresho',.02)]:
                for n in (4,6,8):
                    r=spatial.run(kind,n,end);r['mode']=mode;r['suite']='spatial';emit(r)
                    if r.get('status')!='completed' or r.get('min_factor_R',-1)<-1e-20 or r.get('min_pressure',0)<=0:
                        report[mode]['spatial_stopped_at']=[kind,n];break
                if 'spatial_stopped_at' in report[mode]:break
        report[mode]['time_physics_gate']=bool(passed_new)
        os.chdir(ROOT)
        (ROOT/'convergence.json').write_text(json.dumps(report,indent=2))
    (ROOT/'convergence.json').write_text(json.dumps(report,indent=2))

if __name__=='__main__':main()
