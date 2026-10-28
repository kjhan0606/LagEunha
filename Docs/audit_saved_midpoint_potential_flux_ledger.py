"""Saved actual midpoint potential/energy flux ledger; no fit or RHS.

Psi here is a transport rate, not an evolved state or reconstructed integral.
Endpoint reconstruction and its q24/48 uncertainty are a separate diagnostic.
"""
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path
import numpy as np

base=Path('/gpfs/kjhan/LagForce')
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args();args.output.mkdir(exist_ok=False)
potential_path=base/'Codex_weight_remap_moments_centered_20261017/audit.py'
spec=importlib.util.spec_from_file_location('numerical_potential',potential_path)
potential=importlib.util.module_from_spec(spec);spec.loader.exec_module(potential)
root=base/'Codex_compiled_midpoint_once_20261018'
manifest=json.loads((root/'manifest.json').read_text());h=float(manifest['dt'])
results={};source_paths=[Path(__file__),potential_path,root/'manifest.json']
(args.output/'manifest.json').write_text(json.dumps(dict(scope=__doc__,saved_root=str(root),
    source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),provider_calls=0,newton_fits=0),indent=2))


def sum4(values):return ((values[:,0]+values[:,1])+values[:,2])+values[:,3]


for label in ('stage0','midpoint'):
    input_path=root/label/'inputs.npz';reference_path=root/label/'reference.npz'
    source_paths.extend((input_path,reference_path))
    inputs=np.load(input_path);reference=np.load(reference_path)
    nodes=inputs['node_inputs'];owners=inputs['owners'];frozen=inputs['frozen'].astype(bool)
    active=~frozen;N=len(frozen);i=owners[:,0];j=owners[:,1];internal=j>=0
    left_active=active[i];right_active=np.zeros(len(j),dtype=bool);right_active[internal]=active[j[internal]]
    direction=left_active.astype(int)-right_active.astype(int)
    wl=nodes[:,21:25];wr=nodes[:,25:29];vl=nodes[:,5:13].reshape(-1,4,2);vr=nodes[:,13:21].reshape(-1,4,2)
    fm=nodes[:,4]*(sum4(wl)+sum4(wr))
    energyL=.5*np.sum(vl*vl,axis=2)+nodes[:,29,None]
    energyR=.5*np.sum(vr*vr,axis=2)+nodes[:,30,None]
    # Small dot kernels mirror each original node's dot products.
    energy_flux=nodes[:,4]*((wl[:,None,:]@energyL[:,:,None])[:,0,0]+(wr[:,None,:]@energyR[:,:,None])[:,0,0]+nodes[:,31]*nodes[:,32])
    Phi=potential.potential(nodes[:,:2]);Phi_flux=Phi*fm
    assert all(np.isfinite(v).all() for v in (fm,energy_flux,Phi,Phi_flux))

    def scatter(flux):
        return -np.bincount(i,weights=flux,minlength=N)+np.bincount(j[internal],weights=flux[internal],minlength=N)

    mass_rate=scatter(fm);hydro_rate=scatter(energy_flux);face_Phi_rate=scatter(Phi_flux)
    gE=reference['gravity'][:,17]
    Psi_unmasked=face_Phi_rate-gE
    Psi_evolution=Psi_unmasked.copy();Psi_evolution[frozen]=0.
    gE_evolution=gE.copy();gE_evolution[frozen]=0.
    stored_hydro=reference['original_hydro'][:,17]
    hydro_error=abs(hydro_rate[active]-stored_hydro[active])
    hydro_bound=1e-14+1e-10*abs(stored_hydro[active])
    assert np.all(hydro_error<=hydro_bound)
    assert np.max(abs(mass_rate[active]))<2e-13
    combined=stored_hydro+gE_evolution+Psi_evolution
    cancellation=gE_evolution+Psi_evolution-np.where(active,face_Phi_rate,0.)
    boundary=internal&(left_active!=right_active)
    active_exterior=(~internal)&left_active
    categories=dict(active_active=internal&left_active&right_active,
        active_frozen=boundary,frozen_frozen=internal&~left_active&~right_active,
        patch_exterior=~internal)
    transfers={}
    for name,mask in categories.items():
        transfers[name]=dict(nodes=int(mask.sum()),mass_flux_absolute_sum=float(abs(fm[mask]).sum()),
            potential_flux_absolute_sum=float(abs(Phi_flux[mask]).sum()),
            active_mass_inward_rate=float((-direction[mask]*fm[mask]).sum()),
            active_hydro_boundary_rate=float((-direction[mask]*energy_flux[mask]).sum()),
            active_potential_boundary_rate=float((-direction[mask]*Phi_flux[mask]).sum()))
    boundary_Phi=float((-direction*(boundary|active_exterior)*Phi_flux).sum())
    boundary_hydro=float((-direction*(boundary|active_exterior)*energy_flux).sum())
    combined_sum=float(combined[active].sum())
    # Conservative bound for binary64 accumulation of these stored values.
    operations=8*len(nodes)+4*N
    gamma=operations*np.finfo(float).eps/(1-operations*np.finfo(float).eps)
    bound=float(gamma*(abs(Phi_flux).sum()+abs(energy_flux).sum()+abs(gE).sum()))
    residual=combined_sum-boundary_hydro-boundary_Phi
    assert abs(residual)<=bound
    assert np.max(abs(cancellation))<=bound
    np.savez(args.output/(label+'.npz'),ids=inputs['ids'],active=active,node_Phi=Phi,node_fm=fm,
        node_energy_flux=energy_flux,mass_rate_unmasked=mass_rate,hydro_rate_unmasked=hydro_rate,
        potential_face_rate_unmasked=face_Phi_rate,gravity_energy_rate_unmasked=gE,
        Psi_transport_unmasked=Psi_unmasked,Psi_evolution_masked=Psi_evolution,combined_evolution_rate=combined)
    result=dict(time=float(inputs['time']),nodes=len(nodes),transfers=transfers,
        active_mass_residual_max=float(abs(mass_rate[active]).max()),
        hydro_flux_vs_stored_active_max=float(abs(hydro_rate[active]-stored_hydro[active]).max()),
        shared_node_pair_cancellation_max=float(abs(-Phi_flux[internal]+Phi_flux[internal]).max(initial=0.)),
        gravity_Psi_cancellation_max=float(abs(cancellation).max()),
        active_hydro_energy_rate=float(stored_hydro[active].sum()),
        active_gravity_energy_rate=float(gE[active].sum()),active_Psi_transport_rate=float(Psi_evolution[active].sum()),
        active_boundary_hydro_rate=boundary_hydro,active_boundary_Phi_rate=boundary_Phi,
        active_combined_energy_rate=combined_sum,combined_minus_boundary_residual=residual,
        accumulation_roundoff_bound=bound,frozen_Psi_unmasked_sum=float(Psi_unmasked[frozen].sum()),
        frozen_evolution_Psi_max=float(abs(Psi_evolution[frozen]).max()),
        active_patch_exterior_nodes=int(active_exterior.sum()),
        frozen_note='Frozen evolution is masked to zero; unmasked face/gravity exchange is separately retained, not claimed to evolve.')
    results[label]=result
initial=np.load(root/'initial.npz');end=np.load(root/'accepted.npz');mid=np.load(root/'midpoint/inputs.npz')
assert float(initial['time'])==0. and float(end['time'])==h and float(mid['time'])==h/2
active=~mid['frozen'].astype(bool)
a0=initial['state'][:-4].reshape(-1,20);a1=end['state'][:-4].reshape(-1,20)
delta_hydro=float((a1[active,19]-a0[active,19]).sum())
half=results['midpoint'];boundary_increment=h*(half['active_boundary_hydro_rate']+half['active_boundary_Phi_rate'])
source_paths.extend((root/'initial.npz',root/'accepted.npz'))
report=dict(status='completed',scope=__doc__,stages=results,dt=h,
    full_midpoint=dict(actual_hydro_energy_change=delta_hydro,
        expected_hydro_plus_gravity_increment=h*(half['active_hydro_energy_rate']+half['active_gravity_energy_rate']),
        hydro_update_residual=delta_hydro-h*(half['active_hydro_energy_rate']+half['active_gravity_energy_rate']),
        transported_Psi_increment=h*half['active_Psi_transport_rate'],
        hydro_plus_Phi_boundary_increment=boundary_increment,
        hydro_plus_transported_Psi_minus_boundary=delta_hydro+h*half['active_Psi_transport_rate']-boundary_increment),
    source_sha256={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_paths},
    limits='Stored moment model flux ledger only. Potential transport is not an evolved variable. Numerical reconstructed endpoint Psi/q24-q48 still required; no analytic density replacement, heat correction or long-time success claim.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
