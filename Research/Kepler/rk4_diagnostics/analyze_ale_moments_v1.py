"""Cell geometric/velocity moments from existing native face packets only.

These are read-only estimates for a possible uniform-density affine/spin
representation, not stored/evolved cell rotation or an energy correction.
Area moments equal mass moments only under that uniform-density assumption.
"""
from pathlib import Path
import json
import numpy as np
from analyze_ale_shadow_v1 import ROOT,FOCUS,load


def moments(d):
    p,c,a,b,meta=load(d);z=p['a'];n=len(c)
    dx=(c['x'][b]-c['x'][a]+12)%24-12
    dy=(c['y'][b]-c['y'][a]+12)%24-12
    sums=np.zeros((n,6))
    for owners,x1,y1,x2,y2 in ((a,z[:,0],z[:,1],z[:,2],z[:,3]),
                            (b,z[:,2]-dx,z[:,3]-dy,z[:,0]-dx,z[:,1]-dy)):
        cross=x1*y2-x2*y1
        values=np.column_stack([.5*cross,(x1+x2)*cross/6,(y1+y2)*cross/6,
            (x1*x1+x1*x2+x2*x2)*cross/12,
            (2*x1*y1+x1*y2+x2*y1+2*x2*y2)*cross/24,
            (y1*y1+y1*y2+y2*y2)*cross/12])
        np.add.at(sums,owners,values)
    area=sums[:,0];assert np.all(area>0)
    mx=sums[:,1]/area;my=sums[:,2]/area
    cxx=sums[:,3]/area-mx*mx;cxy=sums[:,4]/area-mx*my;cyy=sums[:,5]/area-my*my
    assert np.all(cxx>0) and np.all(cyy>0) and np.all(cxx*cyy-cxy*cxy>0)
    # Independent face-area/distance weighted point-value affine fit. It is
    # intentionally not called the actual slope-limited native reconstruction.
    length=np.hypot(z[:,2]-z[:,0],z[:,3]-z[:,1]);dist=np.hypot(dx,dy);w=length/dist
    dux=c['vx'][b]-c['vx'][a];duy=c['vy'][b]-c['vy'][a]
    edge=np.column_stack([w*dx*dx,w*dx*dy,w*dy*dy,w*dux*dx,w*dux*dy,w*duy*dx,w*duy*dy])
    fit=np.zeros((n,7));np.add.at(fit,a,edge);np.add.at(fit,b,edge)
    xx,xy,yy,t1,t2,t3,t4=fit.T;det=xx*yy-xy*xy;assert np.all(det>0)
    bxx=(t1*yy-t2*xy)/det;bxy=(-t1*xy+t2*xx)/det
    byx=(t3*yy-t4*xy)/det;byy=(-t3*xy+t4*xx)/det
    mass=c['mass'];omega=.5*(byx-bxy);sxy=.5*(bxy+byx)
    kvar=.5*mass*((bxx*bxx+byx*byx)*cxx+2*(bxx*bxy+byx*byy)*cxy+(bxy*bxy+byy*byy)*cyy)
    krot=.5*mass*omega*omega*(cxx+cyy)
    kstrain=.5*mass*((bxx*bxx+sxy*sxy)*cxx+2*sxy*(bxx+byy)*cxy+(sxy*sxy+byy*byy)*cyy)
    kcross=kvar-krot-kstrain
    lspin=mass*(byx*cxx+(byy-bxx)*cxy-bxy*cyy)
    dvmeanx=bxx*mx+bxy*my;dvmeany=byx*mx+byy*my
    delta_point_to_affine_mean_K=mass*(c['vx']*dvmeanx+c['vy']*dvmeany)+.5*mass*(dvmeanx**2+dvmeany**2)
    result=dict(directory=str(d),scope='exact polygon moments; independent unlimited affine point-value fit; no spin state, no energy added, no evolution',
        max_polygon_volume_relative_to_stored=float(abs(area-c['V']).max()/c['V'].min()),
        max_polygon_volume_pointwise_relative=float((abs(area-c['V'])/area).max()),
        velocity_fit='sum_face (area/distance) |delta_v-B*delta_x|^2; fixed current positions/velocities',
        kinetic_definition='Kvar=.5*m*tr(B Cov(x) B^T); antisymmetric rotation, symmetric strain and their anisotropic cross term kept separately',
        focus=[])
    for ident in FOCUS:
        i=int(np.flatnonzero(c['id']==ident)[0])
        result['focus'].append(dict(id=ident,U=float(c['U'][i]),
            generator_to_polygon_centroid=[float(mx[i]),float(my[i])],
            covariance=[[float(cxx[i]),float(cxy[i])],[float(cxy[i]),float(cyy[i])]],
            B_fit=[[float(bxx[i]),float(bxy[i])],[float(byx[i]),float(byy[i])]],
            omega_fit=float(omega[i]),K_variance_estimate=float(kvar[i]),
            K_rotation_estimate=float(krot[i]),K_strain_estimate=float(kstrain[i]),
            K_rotation_strain_cross=float(kcross[i]),
            K_variance_over_native_U=float(kvar[i]/c['U'][i]),
            intrinsic_Lz_estimate=float(lspin[i]),
            affine_mean_velocity_offset=[float(dvmeanx[i]),float(dvmeany[i])],
            bulk_K_change_if_generator_velocity_were_replaced_by_affine_cell_mean=float(delta_point_to_affine_mean_K[i])))
    return result


target=ROOT/'rk4_diagnostics/ale_shadow_real_v1/moments.json'
if target.exists():raise SystemExit('Refuse to overwrite existing moment analysis')
result={name:moments(target.parent/name) for name in ('fine_captured','fine_fresh_final')}
target.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(dict(output=str(target),focus68526={key:next(r for r in value['focus'] if r['id']==68526) for key,value in result.items()}),indent=2),flush=True)
