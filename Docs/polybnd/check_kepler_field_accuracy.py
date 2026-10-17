"""Read-only Eulerian field errors of canonical Kepler fluid checkpoints.

Integrate over the fixed physical circle r<2.5 with positive split polar
quadrature. Relative norms and sampled extrema are diagnostics; this does not
certify 100-orbit evolution, density path bounds, or total potential flux.
"""
import argparse,json,os,hashlib
from pathlib import Path
import numpy as np
import moment_ale_evolve as base
from check_kepler_moment_preflight import install
from physical_field_accuracy_audit import environment
from kepler_disk_reference import polar_rule,fields,CENTER
from quadratic_velocity_reconstruction import prepare,correction
from quadratic_pressure_reconstruction import prepare as prepare_pressure,value


def inspect(path,order,progress=None,configuration_sha256=None):
    data=Path(path).read_bytes();cp=json.loads(data);mass=np.asarray(cp['mass']);a=np.asarray(cp['state'])[:-4].reshape(-1,20);L=[np.asarray(v) for v in cp['factors']]
    polys,edges=base.geometry(a[:,:2],a[:,2],cp['t'])
    density_provider=None
    if 'potential_moments' in cp:
        from potential_moment_coupled_2d import provider
        from audit_kepler_energy_step import phi
        from kepler_disk_reference import frozen
        U=np.asarray(cp['potential_moments'])
        if U.shape!=mass.shape:raise ValueError('invalid stored potential moments')
        density_provider=provider(U,phi,lambda p:polar_rule(p,16,radii=(.25,.5)),frozen(cp['t'],a[:,:2]))
    raw=[np.r_[m,row[3:]] for m,row in zip(mass,a)]
    states=[base.decode(p,r,density=None if density_provider is None else density_provider(i,p,r)) for i,(p,r) in enumerate(zip(polys,raw))]
    prepare(states,polys,edges,mass,L);pressure,_=prepare_pressure(states,polys)
    error=np.zeros(3);denom=np.zeros(3);minimum=np.full(2,np.inf);maximum=np.zeros(3);volume=0.
    identity=dict(input_sha256=hashlib.sha256(data).hexdigest(),order=order,configuration_sha256=configuration_sha256)
    if progress is not None:
        identity['source_manifest_sha256']=hashlib.sha256(Path(__file__).with_name('SHA256SUMS').read_bytes()).hexdigest()
    first=0
    if progress is not None and Path(progress).exists():
        saved=json.loads(Path(progress).read_text())
        if saved['identity']!=identity:raise ValueError('field audit progress identity mismatch')
        first=saved['next_cell'];error=np.asarray(saved['error']);denom=np.asarray(saved['denom']);minimum=np.asarray(saved['minimum']);maximum=np.asarray(saved['maximum']);volume=saved['volume']
    def accumulate_cell(i,poly,s):
        nonlocal error,denom,minimum,maximum,volume
        local=poly-CENTER;distance=np.maximum(np.maximum(local.min(axis=0),-local.max(axis=0)),0.)
        if np.linalg.norm(distance)>=2.5:return
        x,w=polar_rule(poly,order,radii=(.25,.5,2.,2.5));inside=np.linalg.norm(x-CENTER,axis=1)<2.5;x=x[inside];w=w[inside]
        if not len(w):return
        rho,v,P=fields(x);rn=s['density'].value(x);vn=s['v']+(x-s['c'])@s['B'].T+correction(s,x);pn=value(pressure[i],x)
        delta=np.column_stack([abs(rn-rho),rho*np.linalg.norm(vn-v,axis=1),abs(pn-P)])
        error+=w@delta;denom+=w@np.column_stack([rho,rho*np.linalg.norm(v,axis=1),P]);volume+=w.sum()
        minimum=np.minimum(minimum,[rn.min(),pn.min()]);maximum=np.maximum(maximum,np.max(np.column_stack([abs(rn/rho-1),np.linalg.norm(vn-v,axis=1),abs(pn/P-1)]),axis=0))
    for i,(poly,s) in enumerate(zip(polys,states)):
        if i<first:continue
        accumulate_cell(i,poly,s)
        if progress is not None and ((i+1)%32==0 or i+1==len(polys)):
            saved=dict(identity=identity,next_cell=i+1,error=error.tolist(),denom=denom.tolist(),minimum=minimum.tolist(),maximum=maximum.tolist(),volume=volume)
            temp=Path(str(progress)+'.tmp');temp.write_text(json.dumps(saved));temp.replace(progress)
            print('Kepler field audit',cp['t'],order,i+1,'/',len(polys),flush=True)
    return dict(t=cp['t'],order=order,relative_L1=dict(zip(('density','velocity_mass_weighted','pressure'),(error/denom).tolist())),volume=float(volume),volume_reference=float(np.pi*2.5**2),minimum_density=float(minimum[0]),minimum_pressure=float(minimum[1]),sampled_maximum_density_relative=float(maximum[0]),sampled_maximum_velocity_absolute=float(maximum[1]),sampled_maximum_pressure_relative=float(maximum[2]),sha256=hashlib.sha256(data).hexdigest())

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('initial');p.add_argument('checkpoint');p.add_argument('configuration');a=p.parse_args()
    cfg=json.loads(Path(a.configuration).read_text());os.environ.update(environment(cfg));install()
    rows=[];cfg_sha=hashlib.sha256(Path(a.configuration).read_bytes()).hexdigest()
    for label,path in [('initial',a.initial),('evolved',a.checkpoint)]:
        for order in (16,24):
            row=dict(state=label,**inspect(path,order,progress=f'{label}_{order}_progress.json',configuration_sha256=cfg_sha));rows.append(row)
            Path(f'{label}_{order}_result.json').write_text(json.dumps(row,indent=2))
    result=dict(scope='first-step field audit only',full_validation=False,rows=rows);Path('kepler_field_audit.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True)
