"""Geometry and original admissibility/checkpoint audits for native candidates.

No integration or state acceptance. Reference fields are compared only; no
analytic value is substituted into the native candidate or its factors/Psi.
"""
import json,subprocess,sys
from pathlib import Path
import numpy as np
from full_moment_state_io import load as load_state,digest as state_digest,HEADER
base=Path('/gpfs/kjhan/LagForce');old=base/'Codex_full_disk_arc_roots_modules_v7_20261019';fast=base/'Codex_full_disk_fast_modules_v1_20261019'

def write(path,data):
    with Path(path).open('x') as f:json.dump(data,f,indent=2)
def command(script,arguments,log):
    with Path(log).open('x') as out:subprocess.run([sys.executable,str(script),*map(str,arguments)],check=True,stdout=out,stderr=subprocess.STDOUT)

def geometry(service,request):
    source=Path(bytes(request.state_path).decode()).resolve();root=Path(bytes(request.output_dir).decode()).resolve()
    assert source.is_dir() and str(root).startswith('/dev/shm/')
    header=HEADER.unpack((source/'header.bin').read_bytes());assert header[3]==request.epoch and header[4]==request.run_token
    state=load_state(source);root.mkdir(exist_ok=False)
    np.savez(root/'candidate.npz',**state)
    command(old/'prepare_full_curved_stage_geometry_workers.py',['--state',root/'candidate.npz','--workers',8,'--output',root/'geometry'],root/'geometry.log')
    report=json.loads((root/'geometry/report.json').read_text());assert report['status']=='passed_full_curved_stage_geometry'
    write(root/'native_geometry_manifest.json',dict(status='passed_native_candidate_geometry',epoch=int(request.epoch),run_token=int(request.run_token),source_digest=state_digest(source),state_serialization='Exact arrays; no update or replacement',accepted=False))
    return root/'geometry'

def admissibility(service,request):
    ctx=service.ctx;state=load_state(ctx['state']);initial=load_state(bytes(request.aux_path).decode())
    with np.load(Path(ctx['geometry'])/'mesh.npz') as f:mesh={k:f[k] for k in f.files}
    physics=json.loads(Path(ctx['physics_file']).read_text())['physics'];P0=physics['P0'];gamma=physics['gamma']
    n=len(initial['ids']);raw0=np.column_stack([initial['mass'],initial['state'][:-4].reshape(n,20)[:,3:]]).astype(np.longdouble)
    r0=np.linalg.norm(initial['state'][:-4].reshape(n,20)[:,:2],axis=1)
    assert np.array_equal(state['ids'],initial['ids']) and np.array_equal(state['mass'],initial['mass'])
    a=state['state'][:-4].reshape(n,20);raw=np.column_stack([state['mass'],a[:,3:]]).astype(np.longdouble)
    assert np.array_equal(mesh['ids'],state['ids']) and np.array_equal(mesh['sites'],a[:,:2]) and np.array_equal(mesh['weights'],a[:,2])
    thermal=raw[:,17]-.5*(raw[:,13]+raw[:,16]);area=mesh['moments'][:,0]
    assert np.isfinite(raw).all() and np.isfinite(state['Psi']).all() and (thermal>0).all() and (area>0).all()
    pe=(gamma-1)*thermal/(P0*area)-1;mass=state['mass'].astype(np.longdouble);pressure={}
    for label,mask in dict(all=np.ones(n,bool),inner=r0<.4).items():
        pressure[label]=dict(RMS=float(np.sqrt(np.average(pe[mask]**2,weights=mass[mask]))),maximum=float(abs(pe[mask]).max()))
    change=raw-raw0;vector=np.r_[change[:,7:9].sum(0),change[:,17].sum(),(change[:,10]-change[:,11]).sum()]
    ledger=state['state'][-4:]-initial['state'][-4:];residual=vector-ledger
    # Same operational roundoff envelope as the running Ada continuous driver.
    step=max(1,int(round(float(state['time'])/.0025)))
    eps=np.finfo(float).eps;ops=128+4*n+4*step;g=ops*eps/(1-ops*eps)
    scale=np.r_[abs(raw[:,7:9]).sum(0)+abs(raw0[:,7:9]).sum(0),abs(raw[:,17]).sum()+abs(raw0[:,17]).sum(),abs(raw[:,10:12]).sum()+abs(raw0[:,10:12]).sum()]+abs(ledger)
    bound=g*scale;assert np.all(abs(residual)<=bound)
    combined=change[:,17].sum()+(state['Psi'].astype(np.longdouble)-initial['Psi'].astype(np.longdouble)).sum()
    energy_scale=abs(raw[:,17]).sum()+abs(raw0[:,17]).sum()+abs(state['Psi']).sum()+abs(initial['Psi']).sum()
    energy_bound=128*eps*(2*step+1)*energy_scale;assert abs(combined)<=energy_bound
    report=dict(status='passed_native_tentative_admissibility',time=float(state['time']),epoch=int(request.epoch),
        pressure=pressure,minimum_thermal=float(thermal.min()),P_E_L_ledger_residual=np.asarray(residual,float).tolist(),
        P_E_L_ledger_bound=np.asarray(bound,float).tolist(),combined_energy_change=float(combined),combined_energy_roundoff_bound=float(energy_bound),mass_bitwise_initial=True,accepted=False)
    folder=Path(bytes(request.output_dir).decode());folder.mkdir(exist_ok=False)
    write(folder/'admissibility.json',report)
    if request.outer:
        assert ctx['q48'],'Checkpoint requires original independent q48 density audit'
        fits=json.loads(Path(bytes(request.fits_path).decode()).read_text())
        assert [r['id'] for r in fits]==state['ids'].tolist() and all(r['status']=='passed' and r['Psi_error'] is not None for r in fits)
        np.savez(folder/'candidate.npz',**state)
        density=folder/'density';density.mkdir()
        with (density/'cells.jsonl').open('x') as f:
            for row in fits:f.write(json.dumps(row)+'\n')
        import hashlib
        write(density/'report.json',dict(status='passed_full_stage_Psi_density_preparation',state_sha256=hashlib.sha256((folder/'candidate.npz').read_bytes()).hexdigest(),signed_q48_Psi_error=sum(r['Psi_error'] for r in fits),absolute_q48_Psi_error=sum(abs(r['Psi_error']) for r in fits),maximum_q48_Psi_error=max(abs(r['Psi_error']) for r in fits),cells=n,time=float(state['time']),source='Native-owned fits with unchanged independent q48'))
        command(old/'project_steady_reference_on_full_curved_state_workers.py',['--state',folder/'candidate.npz','--geometry',ctx['geometry'],'--workers',8,'--output',folder/'reference'],folder/'reference.log')
        command(fast/'compare_full_curved_checkpoint_physics.py',['--state',folder/'candidate.npz','--geometry',ctx['geometry'],'--density',density,'--reference',folder/'reference','--output',folder/'physics'],folder/'physics.log')
        physical=json.loads((folder/'physics/report.json').read_text());gates={}
        for label,row in physical['regions'].items():
            if not row['cells']:continue
            p=row['relative_pressure_error'];v=row['fields']['velocity']['relative_RMS'];phase=row.get('velocity_direction_error_degrees')
            gates[label]=dict(pressure_RMS=p['RMS'],pressure_maximum=p['maximum'],velocity_relative_RMS=v,phase_RMS_degrees=phase['RMS'] if phase else None,passed=p['RMS']<.01 and p['maximum']<.1 and v is not None and v<.01 and (not phase or phase['RMS']<1.))
            if request.reserved[8]:
                density_error=row['relative_mass_reference_error']
                gates[label].update(density_RMS=density_error['RMS'],density_maximum=density_error['maximum'])
                gates[label]['passed']=gates[label]['passed'] and density_error['RMS']<.01 and density_error['maximum']<.1
        passed=all(row['passed'] for row in gates.values())
        write(folder/'physical_gate.json',dict(status='passed_checkpoint_physical_monitor' if passed else 'failed_checkpoint_physical_monitor',time=float(state['time']),regions=gates,criteria_source='Unchanged Astra 2026-10-19 monitoring; not a final long-time success tolerance'))
        assert passed,'Native candidate failed checkpoint physical monitor; no accepted commit'
    return folder/'admissibility.json'
