"""Instantaneous force/heat from dumped production HG mode coefficients.

For r=L v and C=A0*(b*abs(r)+sqrt(b*b*r*r+kap*U_owner)),
F=-C*r*L and heat_owner=C*r*r. Frozen degrees of freedom are excluded
from force application with their reservoir work retained explicitly.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import struct
import re

import numpy as np
from audit_canonical_pair_history import gravity_coefficient

ap=argparse.ArgumentParser();ap.add_argument('root',type=Path);args=ap.parse_args();root=args.root
base=Path('/gpfs/kjhan/LagForce')
assert (root/'exit.code').read_text().strip()=='0'
spec=importlib.util.spec_from_file_location('reader',base/'run_canonical_av_branches_20261018.py')
reader=importlib.util.module_from_spec(spec);spec.loader.exec_module(reader)
manifest=json.loads((root/'manifest.json').read_text())
source_copy=root/manifest.get('source_copy',Path(manifest['source']).name)
initial,t=reader.read(source_copy)
assert reader.sha(source_copy)==manifest['source_sha256']
ic,_=reader.read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat')
n=len(initial)
state=np.loadtxt(root/'raw_faces.dat.hgstate.rank0000');state=state[np.argsort(state[:,0])]
assert state.shape==(n,13) and np.array_equal(state[:,0],np.arange(n)) and np.all(state[:,2]==1)
for name,col in [('x',3),('y',4),('mass',5),('vx',6),('vy',7),('ie',8)]:
    assert np.max(abs(state[:,col]-initial[name]))<1e-12,(name,float(abs(state[:,col]-initial[name]).max()))
rhs=np.loadtxt(root/'raw_faces.dat.rhs.rank0000');rhs=rhs[np.argsort(rhs[:,0])]
rhs_diff=None
if t>1:
    oldpath=(base/'Codex_TCF_rhs_identity_20261017/t6/rhs.dat' if t>5 else base/'Codex_t2_onset_rhs_20261017/t2/rhs.dat')
    oldrhs=np.loadtxt(oldpath);oldrhs=oldrhs[np.argsort(oldrhs[:,0])]
    rhs_diff=float(abs(rhs[:,:13]-oldrhs[:,:13]).max());assert rhs_diff<1e-12,rhs_diff
path=root/'raw_faces.dat.hgm.rank0000'
with path.open('rb') as f:magic,nm,ns,ncol=struct.unpack('<4q',f.read(32))
assert magic==20261017 and ncol==9 and path.stat().st_size==32+nm*72+ns*24
m=np.memmap(path,mode='r',dtype='<f8',offset=32,shape=(nm,9))
s=np.memmap(path,mode='r',dtype='<f8',offset=32+nm*72,shape=(ns,3))
assert np.isfinite(m).all() and np.isfinite(s).all()
owner=m[:,0].astype(int);start=m[:,2].astype(int);length=m[:,3].astype(int);ids=s[:,0].astype(int)
assert np.all(m[:,4:6]==0) and np.all(length>0)
assert np.array_equal(start,np.r_[0,np.cumsum(length[:-1])]) and length.sum()==ns
assert np.all((ids>=0)&(ids<n)) and np.all((owner>=0)&(owner<n))
x=state[:,3:5]-3;mass=state[:,5];U=state[:,8];v=state[:,6:8]
frozen=np.linalg.norm(x,axis=1)>=2.5
free_support=~frozen[ids]
alpha=np.add.reduceat(np.sum(s[:,1:3]**2,axis=1)/mass[ids]*free_support,start)
assert np.all(alpha>=0)
omega=np.sqrt(gravity_coefficient(np.linalg.norm(x,axis=1)))
vc=omega[:,None]*np.column_stack((-x[:,1],x[:,0]))
delta=v-vc
r0=np.hypot(ic['x']-3,ic['y']-3)
masks={'initial_inner':(r0<.5)&(ic['den']<.5),'initial_inner_edge':(r0<.6)&(ic['den']>.5),
       'initial_disk_core':(r0>=.6)&(r0<1.9)&(ic['den']>.5)}
normal=(x[804]-x[803])/np.linalg.norm(x[804]-x[803])
out={'scope':'Direct instantaneous HG force/heat reconstructed from actual single-rank C mode assembly at saved epoch; no HG time integration.',
     'time':t,'mode_count':int(nm),'support_entries':int(ns),'source_sha256':reader.sha(source_copy),
     'mode_sha256':reader.sha(path),'rhs_max_abs_difference_from_previous_dump':rhs_diff,
     'script_sha256':reader.sha(Path(__file__)),'states':{},
     'limits':['NumPy reduction order differs from C; no bitwise modal-force claim.',
               'Circular velocities reuse the same geometry, thermodynamic state and geometric mode coefficients.',
               'Instantaneous modal sum is not a finite split step or a long-time stability proof.']}
for name,velocity in [('actual',v),('circular',vc)]:
    entry_dot=np.sum(s[:,1:3]*velocity[ids],axis=1)
    residual=np.add.reduceat(entry_dot,start)
    C=m[:,6]*(m[:,7]*abs(residual)+np.sqrt((m[:,7]*residual)**2+m[:,8]*U[owner]))
    C=np.where(alpha>0,C,0.)
    assert np.isfinite(C).all() and np.all(C>=0)
    multiplier=np.repeat(C*residual,length)
    entry_force=-multiplier[:,None]*s[:,1:3]
    force=np.zeros((n,2));np.add.at(force,ids,entry_force*free_support[:,None])
    heat=np.zeros(n);np.add.at(heat,owner,np.where(frozen[owner],0,C*residual**2))
    frozen_residual=np.add.reduceat(entry_dot*(~free_support),start)
    boundary_work=float(np.sum(C*residual*frozen_residual))
    boundary_heat=float(np.sum(C[frozen[owner]]*residual[frozen[owner]]**2))
    kinetic_power=float(np.sum(velocity*force));thermal_power=float(heat.sum())
    energy_residual=kinetic_power+thermal_power-boundary_work+boundary_heat
    assert abs(energy_residual)<1e-9*max(thermal_power,abs(kinetic_power),1e-30)
    record={'heat_total':thermal_power,'kinetic_power':kinetic_power,'boundary_work':boundary_work,
            'boundary_heat':boundary_heat,'energy_balance_residual':energy_residual,
            'pair_radial_relative_acceleration':float(normal@(force[804]/mass[804]-force[803]/mass[803])),
            'cells':{},'regions':{}}
    for k in [803,804,805,933,934]:
        record['cells'][str(k)]={'acceleration':(force[k]/mass[k]).tolist(),'heat_over_U':float(heat[k]/U[k])}
    for label,mask in masks.items():
        record['regions'][label]={'heat':float(heat[mask].sum()),
                                 'work_on_actual_velocity_departure':float(np.sum(delta[mask]*force[mask])),
                                 'acceleration_mass_RMS':float(np.sqrt(np.sum(force[mask]**2/mass[mask,None])/mass[mask].sum()))}
    if name=='actual' and t>5:
        indirect=json.loads((base/'Codex_pair_history_20261018/av.json').read_text())
        record['indirect_comparison']={}
        for old in indirect['cells']:
            k=old['id'];true=force[k]/mass[k]
            coarse=np.array(old['HG_plus_unaccounted_acceleration_extrapolated_coarse'])
            fine=np.array(old['HG_plus_unaccounted_acceleration_extrapolated_fine'])
            record['indirect_comparison'][str(k)]={'direct':true.tolist(),'coarse_error_norm':float(np.linalg.norm(coarse-true)),
                                                 'fine_error_norm':float(np.linalg.norm(fine-true))}
        # Attribute force on cell804 by thermal mode owner; these are
        # nonlocal HG supports, not two-body face forces.
        mode_index=np.repeat(np.arange(nm),length)
        target=ids==804
        projected=np.sum(entry_force[target]*normal,axis=1)/mass[804]
        by_owner=np.zeros(n);np.add.at(by_owner,owner[mode_index[target]],projected)
        largest=np.argsort(abs(by_owner))[-12:][::-1]
        record['largest_mode_owner_contributions_to_804']=[{'mode_owner':int(k),'acceleration_along_pair':float(by_owner[k])} for k in largest]
    if name=='actual':
        actual_force=force.copy();actual_heat=heat.copy()
        np.savez_compressed(root/'instantaneous_actual.npz',force=force,heat=heat,residual=residual,C=C)
    out['states'][name]=record
# Full instantaneous source audit at the same saved state. At t0, if
# v=v_c and pressure is uniform, nonzero residual force is a numerical seed.
faces=np.loadtxt(root/'raw_faces.dat.rhs.faces.rank0000')
ii=faces[:,0].astype(int)
force_q=np.zeros((n,2));np.add.at(force_q,ii,-(faces[:,4]*faces[:,14])[:,None]*faces[:,6:8])
force_p=rhs[:,8:10]-force_q
gamma=float(np.float32(float(re.search(r'define Gamma_hydro\s*=\s*([\d.eE+-]+)',(root/'params.dat').read_text())[1])))
compression=-gamma*rhs[:,12]/rhs[:,6]
q_pdot=(gamma-1)*rhs[:,11]/(rhs[:,5]*rhs[:,6])
hg_pdot=(gamma-1)*actual_heat/(rhs[:,5]*rhs[:,6])
gravity=-gravity_coefficient(np.linalg.norm(x,axis=1))[:,None]*x
source_audit={}
for label,mask in masks.items():
    w=mass[mask];den=w.sum()
    ar=lambda f:float(np.sqrt(np.sum(f[mask]**2/w[:,None])/den))
    rms=lambda z:float(np.sqrt(np.average(z[mask]**2,weights=w)))
    source_audit[label]={'cells':int(mask.sum()),
        'velocity_departure_mass_RMS':float(np.sqrt(np.sum(w[:,None]*delta[mask]**2)/den)),
        'pressure_acceleration_mass_RMS':ar(force_p),'TCF_acceleration_mass_RMS':ar(force_q),
        'HG_acceleration_mass_RMS':ar(actual_force),
        'total_hydro_acceleration_mass_RMS':ar(rhs[:,8:10]+actual_force),
        'gravity_acceleration_mass_RMS':float(np.sqrt(np.average(np.sum(gravity[mask]**2,axis=1),weights=w))),
        'pressure_over_initial_min_max':[float((rhs[mask,5]/ic['pressure'][mask]).min()),float((rhs[mask,5]/ic['pressure'][mask]).max())],
        'Vdot_over_V_mass_RMS':rms(rhs[:,12]/rhs[:,6]),
        'compression_Pdot_over_P_mass_RMS':rms(compression),
        'TCF_Pdot_over_P_mass_RMS':rms(q_pdot),'HG_Pdot_over_P_mass_RMS':rms(hg_pdot),
        'total_Pdot_over_P_mass_RMS':rms(compression+q_pdot+hg_pdot)}
out['full_source_audit']=source_audit
with (root/'analysis.json').open('x') as f:json.dump(out,f,indent=2);f.write('\n')
print(json.dumps(out,indent=2))
