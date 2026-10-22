"""Signed native pressure work against a fitted internal affine velocity.

No energy is subtracted from U in evolution. This is an instantaneous work
decomposition, not a proof that omega*torque is the complete spin-energy RHS.
"""
import json
import numpy as np
from analyze_ale_shadow_v1 import ROOT,FOCUS,load

base=ROOT/'rk4_diagnostics/ale_shadow_real_v1'
target=base/'affine_work.json'
if target.exists():raise SystemExit('Refuse to overwrite affine work analysis')
moments=json.loads((base/'moments.json').read_text())
result={}
for name in ('fine_captured','fine_fresh_final'):
    p,c,a,b,metas=load(base/name);z=p['a']
    dx=(c['x'][b]-c['x'][a]+12)%24-12
    dy=(c['y'][b]-c['y'][a]+12)%24-12
    rx=.5*(z[:,0]+z[:,2]);ry=.5*(z[:,1]+z[:,3])
    rows=[]
    for ident in FOCUS:
        i=int(np.flatnonzero(c['id']==ident)[0]);m=next(r for r in moments[name]['focus'] if r['id']==ident)
        first=np.flatnonzero(a==i);second=np.flatnonzero(b==i)
        faces=np.r_[first,second];sign=np.r_[-np.ones(len(first)),np.ones(len(second))]
        rr=np.column_stack([np.r_[rx[first],rx[second]-dx[second]],np.r_[ry[first],ry[second]-dy[second]]])
        force=sign[:,None]*z[faces,22:24]
        B=np.array(m['B_fit']);omega=m['omega_fit'];centroid=np.array(m['generator_to_polygon_centroid'])
        affine_velocity=rr@B.T
        prot=omega*np.sum(rr[:,0]*force[:,1]-rr[:,1]*force[:,0])
        pfit=float(np.sum(force*affine_velocity));pstrain=pfit-prot
        ucanon=float(np.sum(sign*z[faces,24])-np.dot([c['vx'][i],c['vy'][i]],force.sum(axis=0)))
        uactual=float(c['Edot_native'][i]-c['vx'][i]*c['pxdot_native'][i]-c['vy'][i]*c['pydot_native'][i])
        mean_offset=float(np.dot(force.sum(axis=0),B@centroid))
        torque_centroid=float(np.sum((rr[:,0]-centroid[0])*force[:,1]-(rr[:,1]-centroid[1])*force[:,0]))
        rows.append(dict(id=ident,native_Udot=uactual,native_Udot_canonical=ucanon,
            rotation_pressure_power_about_generator=float(prot),
            symmetric_strain_pressure_power_about_generator=float(pstrain),
            affine_pressure_power_about_generator=pfit,
            remaining_pressure_work_after_affine_decomposition=float(ucanon-pfit),
            force_power_from_generator_to_fitted_mean_velocity=mean_offset,
            rotation_pressure_power_about_polygon_centroid=omega*torque_centroid,
            strain_pressure_power_about_polygon_centroid=pfit-mean_offset-omega*torque_centroid,
            identity_residual=float(ucanon-prot-pstrain-(ucanon-pfit))))
    result[name]=dict(scope='read-only signed decomposition using actual native pressure and midpoint work plus an independently fitted affine velocity',
        caution='Pressure torque power alone is not dKvariance/dt: moment evolution, mass transport, moving origin, gravity and angular flux remain. No change to U or EOS was made.',focus=rows)
target.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(output=str(target),focus68526={k:next(r for r in v['focus'] if r['id']==68526) for k,v in result.items()}),indent=2),flush=True)
