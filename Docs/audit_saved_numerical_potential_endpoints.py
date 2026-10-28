"""Saved numerical density Psi and endpoint prior/trace diagnostics.

No RHS, weight solve, update, analytical density or heat adjustment. Initial
theta is restored; each remaining saved state is fitted exactly once. q24/48
uses identical density objects, normalization and prior polygon parts.
"""
import argparse
import hashlib
import importlib.util
import json
import sys
import time
from pathlib import Path
import numpy as np

base=Path('/gpfs/kjhan/LagForce')
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args()
sys.path[:0]=[str(base/'Codex_compiled_midpoint_modules_20261018'),
    str(base/'Codex_face_state_batch_api_modules_20261018'),
    str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
from contact_prior_density import ContactPriorDensity,local_contact_priors,polygon_area
from density_quadrature import rule
from saved_contact_density_restore import restore_density
from embedded_laguerre_geometry import EmbeddedLaguerreGeometry
from run_canonical_av_branches_20261018 import read

def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def crosses(poly,radius):
    delta=np.roll(poly,-1,axis=0)-poly;length2=np.sum(delta*delta,axis=1)
    t=np.zeros(len(delta));good=length2>0
    t[good]=np.clip(-np.sum(poly[good]*delta[good],axis=1)/length2[good],0,1)
    minimum=np.linalg.norm(poly+t[:,None]*delta,axis=1).min()
    side=delta[:,0]*(-poly[:,1])-delta[:,1]*(-poly[:,0])
    if np.all(side>=0):minimum=0.
    return bool(minimum<=radius<=np.linalg.norm(poly,axis=1).max())

def stats(values,mask,weights=None):
    v=np.asarray(values)[mask]
    if not len(v):return dict(count=0)
    result=dict(count=len(v),minimum=float(v.min()),median=float(np.median(v)),
        p90=float(np.quantile(v,.9)),p99=float(np.quantile(v,.99)),maximum=float(v.max()))
    result['RMS']=float(np.sqrt(np.average(v*v,weights=None if weights is None else weights[mask])))
    return result

args.output.mkdir(exist_ok=False)
root=base/'Codex_compiled_midpoint_once_20261018'
checkpoints=dict(initial=root/'initial.npz',accepted=root/'accepted.npz',
    coarse=base/'Codex_contact_prior_evolution_dt1e4_20261018/step0010.npz',
    fine=base/'Codex_contact_prior_evolution_dt5e5_20261018/step0020.npz')
initial=np.load(checkpoints['initial']);ids=initial['ids'];mass=initial['mass']
a0=initial['state'][:-4].reshape(-1,20)
frozen=np.linalg.norm(a0[:,:2],axis=1)>=.60;active=~frozen
clip_path=base/'Codex_weighted_compression_geometry_20261017/audit.py'
potential_path=base/'Codex_weight_remap_moments_centered_20261017/audit.py'
clip=load('saved_clip',clip_path);potential=load('saved_potential',potential_path).potential
ic,t0=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');assert t0==0.
builder=EmbeddedLaguerreGeometry(np.column_stack([ic['x']-3,ic['y']-3]),ids,clip.clip)
manifest=dict(scope=__doc__,checkpoints={k:dict(path=str(p),sha256=sha(p)) for k,p in checkpoints.items()},
    sources={str(p):sha(p) for p in (Path(__file__),clip_path,potential_path)},
    restored_initial_archive={name:sha(base/'Codex_adjacency_prior_impact_20261018/stage0'/name)
        for name in ('common_inputs.npz','priors.json','fit_residuals.json')},
    canonical_initial_sha256=sha(base/'Codex_canonical_initial_20261017/run/khout.000000.dat'),
    environment_sha256=sha(base/'Codex_canonical_initial_20261017/run/environment.json'),
    runtime_modules={name:dict(path=str(m.__file__),sha256=sha(m.__file__))
        for name,m in list(sys.modules.items()) if getattr(m,'__file__',None)
        and Path(m.__file__).suffix=='.py' and Path(m.__file__).is_relative_to(base)},
    provider_calls=0,evolution_steps=0,quadrature_orders=[24,48],fit_order=12)
(args.output/'manifest.json').write_text(json.dumps(manifest,indent=2))
reports={}
for label,path in checkpoints.items():
    started=time.perf_counter();saved=np.load(path)
    assert np.array_equal(ids,saved['ids']) and np.array_equal(mass,saved['mass'])
    if 'frozen' in saved:assert np.array_equal(frozen,saved['frozen'].astype(bool))
    assert float(saved['time'])==dict(initial=0.,accepted=1e-4,coarse=.001,fine=.001)[label]
    a=saved['state'][:-4].reshape(-1,20);raw=np.column_stack([mass,a[:,3:]])
    polys,edges=builder(a[:,:2],a[:,2],float(saved['time']))
    priors=local_contact_priors(polys,raw,edges)
    if label=='initial':
        archive=base/'Codex_adjacency_prior_impact_20261018/stage0'
        common=np.load(archive/'common_inputs.npz')
        assert np.array_equal(common['state'],saved['state'])
        old_priors=json.loads((archive/'priors.json').read_text())['old']
        # Adjacency audit explicitly supplied the same working error arrays;
        # original provider uses default arrays. Their provenance label differs.
        # Require EVERY numeric/branch/other metadata field to remain exact.
        def without_error_label(meta):
            return {k:v for k,v in meta.items() if k!='coordinate_error_model'}
        assert [without_error_label(m) for m in priors]==[without_error_label(m) for m in old_priors]
        offsets=common['polygon_offsets'];vertices=common['polygon_vertices']
        old_polys=[vertices[lo:hi] for lo,hi in zip(offsets[:-1],offsets[1:])]
        assert all(np.array_equal(p,q) for p,q in zip(polys,old_polys))
        fits=json.loads((archive/'fit_residuals.json').read_text())
        densities=[restore_density(p,r,m,f['old']) for p,r,m,f in zip(polys,raw,priors,fits)]
        newfits=0
    else:
        densities=[ContactPriorDensity(p,r[0],r[1:3],r[3:7].reshape(2,2),m,12)
                   for p,r,m in zip(polys,raw,priors)]
        newfits=len(densities)
    # Archive immediately: expensive one-time fits can be reused after an audit failure.
    offsets=np.r_[0,np.cumsum([len(p) for p in polys])]
    np.savez(args.output/(label+'_reconstruction.npz'),ids=ids,state=saved['state'],mass=mass,
        polygon_offsets=offsets,polygon_vertices=np.concatenate(polys),
        theta=np.array([d.theta for d in densities]),logZ=np.array([d.logZ for d in densities]),
        origin=np.array([d.origin for d in densities]),scale=np.array([d.scale for d in densities]))
    (args.output/(label+'_priors.json')).write_text(json.dumps(priors,indent=2))
    (args.output/(label+'_fits.json')).write_text(json.dumps([dict(theta=d.theta.tolist(),
        iterations=d.iterations,feature_residual=d.residual) for d in densities],indent=2))
    integrals=[]
    for order in (24,48):
        values=[]
        for k,d in enumerate(densities):
            meta=d.cut_metadata
            # Original constructor drops zero-area parts. Preserve corresponding priors.
            if meta['flat']:part_priors=[meta['value']]
            else:
                from contact_prior_density import clip_halfplane
                n=np.asarray(meta['normal']);o=meta['offset']
                candidates=[(clip_halfplane(polys[k],-n,-o),meta['low']),
                            (clip_halfplane(polys[k],n,o),meta['high'])]
                part_priors=[rho for p,rho in candidates if polygon_area(p)>0]
            assert len(part_priors)==len(d.parts)
            total_mass=0.;Psi=0.
            for part,rho in zip(d.parts,part_priors):
                points,w=rule(part,order)
                dm=w*d._value_with_prior(points,np.full(len(w),rho))
                total_mass+=float(dm.sum());Psi+=float(dm@potential(points))
            values.append((total_mass,Psi))
            if (k+1)%400==0:print(json.dumps(dict(label=label,order=order,cells=k+1)),flush=True)
        integrals.append(np.array(values))
    q24,q48=integrals
    assert np.isfinite(q24).all() and np.isfinite(q48).all()
    assert np.all(q24[:,0]>0) and np.all(q48[:,0]>0)
    # Prior classification is independent of Kepler's analytical density values.
    boundary=np.zeros(len(ids),dtype=bool)
    for i,j,*_ in edges:
        if j>=0 and frozen[i]!=frozen[j]:boundary[i]=boundary[j]=True
    env=json.loads((base/'Codex_canonical_initial_20261017/run/environment.json').read_text())
    radii=[float(env['EUNHA_KEPLER_'+name]) for name in ('RIN','ROUT')]
    contact=np.array([any(crosses(p,r) for r in radii) for p in polys])
    regions=dict(all=np.ones(len(ids),dtype=bool),active=active,frozen=frozen,
        inner=np.linalg.norm(a0[:,:2],axis=1)<.4,contact=contact,active_boundary=boundary&active,
        frozen_boundary=boundary&frozen)
    contrast=np.array([m['high']/m['low'] for m in priors])
    theta=np.array([np.linalg.norm(d.theta) for d in densities])
    # Diagnostic quadrature only: split every face at both prior cuts, then
    # sample its interior with Gauss8. No face velocity, flux or RHS evaluation.
    gx,gw=np.polynomial.legendre.leggauss(8);gx=(gx+1)/2;gw/=2
    trace_rows=[]
    for edge_id,(i,j,p,q,n,length) in enumerate(edges):
        if j<0:continue
        cuts=[0.,1.];delta=q-p
        for cell in (i,j):
            meta=priors[cell]
            if meta['flat']:continue
            normal=np.asarray(meta['normal']);den=float(delta@normal)
            if den!=0:
                t=(meta['offset']-float(p@normal))/den
                if 0<t<1:cuts.append(float(t))
        cuts=sorted(set(cuts))
        for lo,hi in zip(cuts[:-1],cuts[1:]):
            points=p+(lo+(hi-lo)*gx)[:,None]*delta
            left=np.array([densities[i].value_trace(point,-n)[0] for point in points])
            right=np.array([densities[j].value_trace(point,n)[0] for point in points])
            jump=right-left;symmetric=2*abs(jump)/(left+right)
            trace_rows.extend(zip(np.full(8,edge_id),np.full(8,i),np.full(8,j),
                length*(hi-lo)*gw,left,right,jump,symmetric))
    traces=np.asarray(trace_rows);ti=traces[:,1].astype(int);tj=traces[:,2].astype(int)
    trace_regions={name:mask[ti]|mask[tj] for name,mask in regions.items()}
    np.savez(args.output/(label+'_diagnostics.npz'),ids=ids,active=active,q24=q24,q48=q48,
        contrast=contrast,theta_norm=theta,traces=traces,**{'region_'+k:v for k,v in regions.items()})
    result=dict(time=float(saved['time']),cells=len(ids),newton_fits=newfits,
        fit_residual_max=max(d.residual for d in densities),seconds=time.perf_counter()-started,
        numerical_Psi={name:dict(q24=float(q24[mask,1].sum()),q48=float(q48[mask,1].sum()),
            signed_q48_minus_q24=float((q48[mask,1]-q24[mask,1]).sum()),
            absolute_q48_minus_q24_sum=float(abs(q48[mask,1]-q24[mask,1]).sum())) for name,mask in regions.items()},
        numerical_mass_quadrature={name:dict(q48_minus_stored_max=float(abs(q48[mask,0]-mass[mask]).max(initial=0)),
            q48_minus_q24_absolute_sum=float(abs(q48[mask,0]-q24[mask,0]).sum())) for name,mask in regions.items()},
        cell_distributions={name:dict(prior_contrast=stats(contrast,mask),theta_norm=stats(theta,mask),
            nonflat=int(sum(not priors[k]['flat'] for k in np.flatnonzero(mask)))) for name,mask in regions.items()},
        face_distributions={name:dict(absolute_jump=stats(abs(traces[:,6]),mask,traces[:,3]),
            symmetric_relative_jump=stats(traces[:,7],mask,traces[:,3])) for name,mask in trace_regions.items()},
        trace_columns=['edge_id','i','j','length_weight','rho_left','rho_right','signed_jump','2absjump_over_sum'])
    reports[label]=result
    (args.output/(label+'_report.json')).write_text(json.dumps(result,indent=2))
    print(json.dumps(dict(label=label,completed=True,seconds=result['seconds'])),flush=True)

ledger_path=base/'Codex_saved_midpoint_Phi_flux_20261018/report.json'
ledger=json.loads(ledger_path.read_text());full=ledger['full_midpoint']
dPsi=reports['accepted']['numerical_Psi']['active']['q48']-reports['initial']['numerical_Psi']['active']['q48']
uncertainty=sum(reports[k]['numerical_Psi']['active']['absolute_q48_minus_q24_sum'] for k in ('initial','accepted'))
budget=dict(reconstructed_Psi_increment_q48=dPsi,
    transported_Psi_increment=full['transported_Psi_increment'],
    reconstructed_minus_transport=dPsi-full['transported_Psi_increment'],
    delta_hydro_plus_reconstructed_Psi_minus_boundary=full['actual_hydro_energy_change']+dPsi-full['hydro_plus_Phi_boundary_increment'],
    q24_q48_empirical_absolute_difference_sum=uncertainty)
report=dict(status='completed',manifest=manifest,states=reports,full_midpoint=budget,
    ledger_sha256=sha(ledger_path),limits='Diagnostic reconstructed Psi, not an evolved conservation law. q24/48 is empirical integration uncertainty, not a rigorous error bound. No long-time validation. Face distributions are Gauss8 length-weighted diagnostics, not original mass flux nodes.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(budget),flush=True)
