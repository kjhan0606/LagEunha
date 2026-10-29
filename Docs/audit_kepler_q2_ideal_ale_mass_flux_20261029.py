"""Mass exchange needed by a hypothetical uniform-density ALE Kepler state."""
import hashlib
import json
from pathlib import Path
import socket
import numpy as np
from kepler_native_global_ids_20261029 import ordered_global, ID_MASK
from laguerre_nonlinear_volume_probe_t6_20261026 import power_areas
from audit_kepler_bounded_pressure_force_t10_20261027 import finite_faces_touching_active
from audit_kepler_q1_t10_pressure_timescales_20261029 import describe

docs=Path('/gpfs/kjhan/LagForce/Git_publish_20261028/Docs')
prior=json.loads((docs/'KEPLER_Q2_BOUNDED_LATE_CORE_2026-10-29.json').read_text())
out=docs/'KEPLER_Q2_IDEAL_ALE_MASS_FLUX_2026-10-29.json'
history=json.loads((docs/'KEPLER_Q2_TORQUE_FEEDBACK_2026-10-29.json').read_text())
if socket.gethostname().split('.')[0].lower() in ('grammar','syntax'):raise RuntimeError('Run geometry on compute host')
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()
initial=Path(prior['initial'])
if sha(initial)!=prior['initial_sha256']:raise ValueError('Initial changed')
a,ma=ordered_global(initial);ids=a['u4if']&ID_MASK
r0=np.hypot(a['x']-3.,a['y']-3.);active=np.flatnonzero(r0<2.)
flag=np.zeros(len(a),bool);flag[active]=True
lookup=np.full(len(a),-1,dtype=int);lookup[active]=np.arange(len(active))
cases=[]
for old in [prior['cases'][0],dict(endpoint='/gpfs/kjhan/LagForce/Kepler_BoundedPoly1024_A10_20261027/wallselect_to60_a10_v1/khout.011058.dat', endpoint_sha256=next(v['source_sha256'] for v in history['cases'] if abs(v['time']-47.50031653561342)<1e-10), time=47.50031653561342)]:
    path=Path(old['endpoint'])
    if sha(path)!=old['endpoint_sha256']:raise ValueError('Changed endpoint')
    p,meta=ordered_global(path,ids)
    if ma['t']!=0 or meta['t']!=old['time'] or not np.array_equal(a['mass'],p['mass']) or np.any(p['w2']!=0):raise ValueError('Invalid input')
    xy=np.column_stack((p['x']-3.,p['y']-3.))
    area,geom,tri,centers=power_areas(xy,np.zeros(len(a)),active,lookup,return_geometry=True)
    err=float(np.max(abs(area/p['volume'][active].astype(float)-1)))
    if geom['hidden_active'] or np.any(area<=0) or err>=1e-5:raise ValueError('Invalid independent geometry')
    i,j,c,l=finite_faces_touching_active(tri,centers,flag)
    delta=xy[j]-xy[i];dist=np.linalg.norm(delta,axis=1)
    if np.any(dist<=0):raise ValueError('Coincident generators')
    normal=delta/dist[:,None];offset=c-.5*(xy[i]+xy[j]);pressure=p['pressure'].astype(float)
    bi=l[:,None]*(.5*normal+offset/dist[:,None]);bj=l[:,None]*(.5*normal-offset/dist[:,None])
    force=np.zeros((len(p),2));dv=np.zeros(len(p));ii=flag[i];jj=flag[j];dp=pressure[i]-pressure[j]
    np.add.at(force,i[ii],dp[ii,None]*bi[ii]);np.add.at(force,j[jj],dp[jj,None]*bj[jj])
    vel=np.column_stack((p['vx'],p['vy']))
    speed=.5*np.einsum('ij,ij->i',vel[i]+vel[j],normal)-np.einsum('ij,ij->i',vel[j]-vel[i],offset)/dist
    flux=l*speed;np.add.at(dv,i[ii],flux[ii]);np.add.at(dv,j[jj],-flux[jj])
    # Same geometry and labels but the exact softened Kepler velocity at each
    # current site. Its continuum divergence is identically zero away from the center.
    rnow=np.linalg.norm(xy,axis=1);omega_ref=(rnow*rnow+.01175**2)**(-.75)
    vref=np.column_stack((-xy[:,1]*omega_ref,xy[:,0]*omega_ref))
    speed_ref=.5*np.einsum('ij,ij->i',vref[i]+vref[j],normal)-np.einsum('ij,ij->i',vref[j]-vref[i],offset)/dist
    dvr=np.zeros(len(p));flux_ref=l*speed_ref
    np.add.at(dvr,i[ii],flux_ref[ii]);np.add.at(dvr,j[jj],-flux_ref[jj])
    tangent=np.column_stack((-normal[:,1],normal[:,0]))
    def material_flux(order):
        nodes,weights=np.polynomial.legendre.leggauss(order)
        total=np.zeros(len(i))
        for node,weight in zip(nodes,weights):
            pos=c+(.5*node*l)[:,None]*tangent
            om=(np.sum(pos**2,axis=1)+.01175**2)**(-.75)
            v=np.column_stack((-pos[:,1]*om,pos[:,0]*om))
            total+=.5*l*weight*np.einsum('ij,ij->i',v,normal)
        return total
    material4=material_flux(4);material2=material_flux(2)
    # rho=1 uniformly: positive net moves mass into owner i.
    transfer=flux_ref-material4;transfer2=flux_ref-material2
    mdot=np.zeros(len(p));mdot2=np.zeros(len(p));gross=np.zeros(len(p));mflux=np.zeros(len(p))
    np.add.at(mdot,i[ii],transfer[ii]);np.add.at(mdot,j[jj],-transfer[jj])
    np.add.at(mdot2,i[ii],transfer2[ii]);np.add.at(mdot2,j[jj],-transfer2[jj])
    np.add.at(gross,i[ii],abs(transfer[ii]));np.add.at(gross,j[jj],abs(transfer[jj]))
    np.add.at(mflux,i[ii],material4[ii]);np.add.at(mflux,j[jj],-material4[jj])
    h=1e-5
    ap,_=power_areas(xy+h*vref,np.zeros(len(a)),active,lookup)
    am,_=power_areas(xy-h*vref,np.zeros(len(a)),active,lookup)
    diff=(ap-am)/(2*h);inner=(r0>=.5)&(r0<.75)
    ferr=float(np.linalg.norm((diff-dvr[active])[inner[active]])/np.linalg.norm(diff[inner[active]]))
    if ferr>1e-5:raise ValueError('geometry derivative did not match independent difference')
    graddiv=p['stress.gUxx'].astype(float)+p['stress.gUyy'].astype(float)
    geomrate=dv[active]/area;refrate=dvr[active]/area
    radius=np.hypot(xy[:,0],xy[:,1]);gravity=radius*(radius*radius+.01175**2)**(-1.5)
    omega=(radius*radius+.01175**2)**(-.75)
    accel=np.linalg.norm(force,axis=1)/p['mass'].astype(float)
    g_ratio=np.divide(accel,gravity,out=np.zeros(len(p)),where=gravity>0)
    thermal=np.abs(pressure*dv)/(p['ie'].astype(float)*omega);mass=p['mass'].astype(float)
    regions={}
    for name,lo,hi in [('background',.1,.5),('inner_disk',.5,.75),('disk',.5,2.)]:
        mask=(r0>=lo)&(r0<hi)
        if not np.all(gravity[mask]>0):raise ValueError('Undefined gravity ratio')
        regions[name]={'cells':int(mask.sum()),'pressure_acceleration':describe(accel[mask],mass[mask]),
                       'pressure_accel_over_gravity':describe(g_ratio[mask],mass[mask]),
                       'absolute_thermal_rate_over_Uomega':describe(thermal[mask],mass[mask])}
    subset=inner[active];mm=mass[active][subset];om=omega[active][subset]
    a0=area[subset];net=mdot[active][subset]/a0/om;total=gross[active][subset]/a0/om
    cancel=(mdot[active][subset]-dvr[active][subset])/a0/om
    quad=(mdot[active][subset]-mdot2[active][subset])/a0/om
    local_mass={'unit_density_net_mass_fraction_rate_over_omega':describe(abs(net),mm),
                'unit_density_gross_face_exchange_over_mass_omega':describe(total,mm),
                'post_exchange_density_rate_over_omega':describe(abs(cancel),mm),
                'quadrature_4_vs_2_abs_over_omega':describe(abs(quad),mm),
                'ratio_sum_abs_net_to_gross':float(np.sum(abs(mdot[active][subset]))/np.sum(gross[active][subset]))}
    local={'geometric_abs_over_omega':describe(abs(geomrate[subset])/om,mm),
           'stored_gradient_divergence_abs_over_omega':describe(abs(graddiv[active][subset])/om,mm),
           'geometric_minus_gradient_abs_over_omega':describe(abs(geomrate[subset]-graddiv[active][subset])/om,mm),
           'exact_kepler_velocity_geometric_abs_over_omega':describe(abs(refrate[subset])/om,mm),
           'exact_velocity_centered_difference_relative_rms':ferr}
    cases.append({'inner_mass_exchange_counterfactual':local_mass,'inner_divergence_comparison':local,'time':meta['t'],'endpoint_sha256':old['endpoint_sha256'],'active_cells':len(active),
                  'faces_touching_active':len(i),'max_relative_area_error':err,'regions':regions})
result={'initial_sha256':prior['initial_sha256'],'cases':cases,'script_sha256':sha(Path(__file__)),
        'scope':'Hypothetical full Euler ALE exchange for uniform unit density and exact Kepler velocity on saved Q2 positions; mass-rate face integrals are conservative pairwise. W0 native calculation has fixed masses, so this deliberately violates its reduction and is not an implementation or admissible W0 alternative. No weight law, nonuniform reconstruction, pressure, energy, positivity, time evolution, or Kepler stability proof.'}
with out.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
print(json.dumps(result,indent=2))
