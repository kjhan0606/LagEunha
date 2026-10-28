"""Adjacency-only attribution: fixed old polygons/raw and common uncertainty.

No geometry provider, RHS, weight solve or evolution is changed.
"""
import argparse,hashlib,importlib.util,json,sys,time
from pathlib import Path
from fractions import Fraction as F
import numpy as np
base=Path('/gpfs/kjhan/LagForce');sys.path[:0]=[str(base/'Codex_compiled_midpoint_modules_20261018'),str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
from embedded_laguerre_geometry import EmbeddedLaguerreGeometry
from contact_prior_density import local_contact_priors,ContactPriorDensity,polygon_area
from geometry_uncertainty import working_vertex_error
spec=importlib.util.spec_from_file_location('clip',base/'Codex_weighted_compression_geometry_20261017/audit.py');clip=importlib.util.module_from_spec(spec);spec.loader.exec_module(clip)
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False);started=time.monotonic();reports={}
def sha(a):return hashlib.sha256(np.ascontiguousarray(a).tobytes()).hexdigest()
def emit(a):
    print(json.dumps(a),flush=True)
    with (args.output/'events.jsonl').open('a') as stream:stream.write(json.dumps(a)+'\n')
def signature(prior):
    if prior['flat']:return ('constant',F(prior['value']))
    low,high=F(prior['low']),F(prior['high'])
    if low==high:return ('constant',low)
    nx,ny=map(F,prior['normal']);offset=F(prior['offset']);lead=nx if nx else ny
    if not lead:raise ValueError('zero prior normal')
    if lead<0:nx=-nx;ny=-ny;offset=-offset;low,high=high,low;lead=-lead
    return ('cut',nx/lead,ny/lead,offset/lead,low,high)
def graph(edges,N):
    result=[set() for _ in range(N)]
    for i,j,*_ in edges:
        if j>=0:result[i].add(j);result[j].add(i)
    return result
def stencils(adj):
    out=[]
    for i,nb in enumerate(adj):
        s={i}|nb
        for j in nb:s.update(adj[j])
        out.append(s)
    return out
def fit_residual(fit,raw):
    p=fit.points;dm=fit.weighted_mass;measured=np.r_[dm.sum(),dm@p,(p.T@(dm[:,None]*p)).ravel()];target=raw[:7];delta=measured-target;scale=np.r_[raw[0],np.full(2,raw[0]*fit.scale),np.full(4,raw[0]*fit.area)]
    return dict(absolute=delta.tolist(),normalized=(delta/scale).tolist(),feature_residual=fit.residual,theta=fit.theta.tolist(),iterations=fit.iterations)

def evaluate_trace(folder,name,points,normals,weights,owners,global_pairs,fitA,fitB,face_ids=None):
    count=len(points);rhoAL=np.empty(count);rhoAR=np.empty(count);rhoBL=np.empty(count);rhoBR=np.empty(count);reused=0
    for i in range(len(fitA)):
        left=np.flatnonzero(owners[:,0]==i);right=np.flatnonzero(owners[:,1]==i)
        for mask,direction,a,b in [(left,-normals[left],rhoAL,rhoBL),(right,normals[right],rhoAR,rhoBR)]:
            if not len(mask):continue
            a[mask]=fitA[i].value_trace(points[mask],direction)
            if fitA[i] is fitB[i]:b[mask]=a[mask];reused+=len(mask)
            else:b[mask]=fitB[i].value_trace(points[mask],direction)
    boundary=owners[:,1]<0;rhoAR[boundary]=rhoAL[boundary];rhoBR[boundary]=rhoBL[boundary]
    assert all(np.isfinite(x).all() and np.all(x>0) for x in (rhoAL,rhoAR,rhoBL,rhoBR)) and np.all(weights>0)
    jumpA=abs(rhoAL-rhoAR)/np.maximum(rhoAL,rhoAR);jumpB=abs(rhoBL-rhoBR)/np.maximum(rhoBL,rhoBR)
    deltaL=abs(rhoAL-rhoBL);deltaR=abs(rhoAR-rhoBR);relL=deltaL/np.maximum(rhoAL,rhoBL);relR=deltaR/np.maximum(rhoAR,rhoBR)
    np.savez(folder/(name+'_trace.npz'),points=points,normals=normals,length_weights=weights,owners_local=owners,global_pairs=global_pairs,rho_old_L=rhoAL,rho_old_R=rhoAR,rho_new_L=rhoBL,rho_new_R=rhoBR,face_ids=face_ids if face_ids is not None else np.arange(count))
    return dict(nodes=count,total_length_weight=float(weights.sum()),old_relative_jump_max=float(jumpA.max()),new_relative_jump_max=float(jumpB.max()),old_relative_jump_length_sum=float(weights@jumpA),new_relative_jump_length_sum=float(weights@jumpB),maximum_relative_density_change=float(max(relL.max(),relR.max())),density_change_length_sum=float(weights@(deltaL+deltaR)),relative_density_change_length_sum=float(weights@(relL+relR)),unchanged_fit_sample_reuses=reused,boundary_mirror_nodes=int(boundary.sum()),scope='A/B at identical points and weights; density mirror R=L on the existing physical exterior BC. Geometry inventories reported separately.')

try:
    initial=np.load(base/'Codex_production_moment_geometry_20261018/stage0_inputs.npz');registry_initial=initial['xy'];full=base/'Codex_full_native_geometry_v2_20261018'
    for label in ('stage0','midpoint'):
        folder=args.output/label;folder.mkdir();saved=np.load(base/'Codex_compiled_midpoint_once_20261018'/label/'inputs.npz');state=saved['state'];mass=saved['mass'];ids=saved['ids'];a=state[:-4].reshape(-1,20);raw=np.column_stack([mass,a[:,3:]]);N=len(ids);lookup={int(g):i for i,g in enumerate(ids)}
        polygons,oldedges=EmbeddedLaguerreGeometry(registry_initial,ids,clip.clip)(a[:,:2],a[:,2],float(saved['time']));native=json.loads((full/label/'cpu1/exact.json').read_text());oldgraph=graph(oldedges,N);exactedges=[];external=[]
        for cell in native['cells']:
            i=lookup[cell['owner']]
            for witness in cell['witnesses']:
                if witness['added_id'] not in lookup:raise RuntimeError(('completion outside density-data patch',cell['owner'],witness))
            for edge in cell['edges']:
                assert F(edge['length_squared'])>0;gid=edge['neighbor']
                if gid in lookup:
                    j=lookup[gid]
                    if i<j:exactedges.append((i,j))
                else:
                    assert saved['frozen'][i];external.append((cell['owner'],gid))
        newgraph=graph(exactedges,N);sA=stencils(oldgraph);sB=stencils(newgraph)
        oldpairs={tuple(sorted((i,j))) for i,nb in enumerate(oldgraph) for j in nb};newpairs={tuple(sorted((i,j))) for i,nb in enumerate(newgraph) for j in nb}
        convert=lambda pair:[int(ids[i]) for i in pair]
        edge_difference=dict(old_only=[convert(p) for p in sorted(oldpairs-newpairs)],exact_only=[convert(p) for p in sorted(newpairs-oldpairs)])
        stencil_changes=[dict(global_id=int(ids[i]),old=sorted(map(int,ids[list(sA[i])])),new=sorted(map(int,ids[list(sB[i])])),added=sorted(map(int,ids[list(sB[i]-sA[i])])),removed=sorted(map(int,ids[list(sA[i]-sB[i])]))) for i in range(N) if sA[i]!=sB[i]]
        errors=[working_vertex_error(p,polygon_area(p)) for p in polygons];mass_error=np.zeros(N)
        priorA=local_contact_priors(polygons,raw,oldedges,vertex_errors=errors,mass_errors=mass_error);priorB=local_contact_priors(polygons,raw,exactedges,vertex_errors=errors,mass_errors=mass_error)
        for p,q in zip(priorA,priorB):
            for key in ('area_bounds','rho_bounds','maximum_vertex_error','mass_error'):assert p[key]==q[key]
        changed=[i for i in range(N) if signature(priorA[i])!=signature(priorB[i])];meta_changed=[i for i in range(N) if priorA[i]!=priorB[i]]
        (folder/'priors.json').write_text(json.dumps(dict(old=priorA,new=priorB,physical_change_global_ids=ids[changed].tolist(),metadata_change_global_ids=ids[meta_changed].tolist(),stencil_changes=stencil_changes,edge_difference=edge_difference),indent=2))
        np.savez(folder/'common_inputs.npz',state=state,raw=raw,ids=ids,polygon_offsets=np.r_[0,np.cumsum([len(p) for p in polygons])],polygon_vertices=np.concatenate(polygons),working_vertex_errors=np.concatenate(errors),mass_errors=mass_error)
        emit(dict(stage=label,operation='prior_comparison',two_hop_changes=len(stencil_changes),metadata_changes=len(meta_changed),physical_prior_changes=len(changed),edge_difference=edge_difference))
        fitA=[];fitB=[];residuals=[]
        for i,(poly,r,p,q) in enumerate(zip(polygons,raw,priorA,priorB)):
            try:
                fa=ContactPriorDensity(poly,r[0],r[1:3],r[3:7].reshape(2,2),p,12);fb=fa if i not in changed else ContactPriorDensity(poly,r[0],r[1:3],r[3:7].reshape(2,2),q,12)
            except BaseException:
                (folder/'fit_failure.json').write_text(json.dumps(dict(local_id=i,global_id=int(ids[i]),old=p,new=q),indent=2));np.savez(folder/'fit_failure.npz',polygon=poly,raw=r);raise
            fitA.append(fa);fitB.append(fb)
            ra=fit_residual(fa,r);rb=ra if fb is fa else fit_residual(fb,r);residuals.append(dict(global_id=int(ids[i]),old=ra,new=rb,reused=fb is fa))
        (folder/'fit_residuals.json').write_text(json.dumps(residuals,indent=2));emit(dict(stage=label,operation='fits_completed',baseline_fits=N,additional_new_fits=len(changed),seconds=time.monotonic()-started))
        # Inventory1: actual stored split face quadrature from the existing RHS.
        nodes=saved['node_inputs'];owners=saved['owners'];gp=np.column_stack([ids[owners[:,0]],np.where(owners[:,1]>=0,ids[np.maximum(owners[:,1],0)],-1)])
        oldtrace=evaluate_trace(folder,'old_saved_face_quadrature',nodes[:,:2],nodes[:,2:4],nodes[:,4],owners,gp,fitA,fitB)
        # Inventory2: native exact owner midpoint+separately rounded edge vector.
        canon=np.load(full/label/'canonical_owner_faces.npz');faces=canon['data'];z,w=np.polynomial.legendre.leggauss(8);points=[];normal=[];weight=[];owners2=[];globalpairs=[];faceid=[];duplicate=[]
        for k,row in enumerate(faces):
            gi,gj=int(row[0]),int(row[1]);i=lookup[gi];j=lookup.get(gj,-1);p=row[2:4]+.5*z[:,None]*row[4:6];points.extend(p);normal.extend(np.tile(row[7:9],(8,1)));weight.extend(.5*row[6]*w);owners2.extend([(i,j)]*8);globalpairs.extend([(gi,gj)]*8);faceid.extend([k]*8)
            unique=len(np.unique(p,axis=0))
            if unique<8:duplicate.append(dict(face_key=canon['keys'][k].tolist(),unique_points=unique,total_points=8,length=float(row[6])))
        exacttrace=evaluate_trace(folder,'exact_midpoint_vector_samples',np.array(points),np.array(normal),np.array(weight),np.array(owners2),np.array(globalpairs),fitA,fitB,np.array(faceid))
        report=dict(time=float(saved['time']),cells=N,edge_difference=edge_difference,two_hop_changed_cells=len(stencil_changes),metadata_changed_cells=len(meta_changed),physical_prior_changed_cells=len(changed),physical_prior_changed_global_ids=ids[changed].tolist(),old_nonflat=sum(not p['flat'] for p in priorA),new_nonflat=sum(not p['flat'] for p in priorB),baseline_fits=N,new_additional_fits=len(changed),unchanged_new_fit_object_aliases=sum(fa is fb for fa,fb in zip(fitA,fitB)),maximum_old_feature_residual=max(r['old']['feature_residual'] for r in residuals),maximum_new_feature_residual=max(r['new']['feature_residual'] for r in residuals),trace=dict(old_geometry=oldtrace,exact_geometry=exacttrace,rounded_duplicate_exact_faces=duplicate,exact_trace_scope='Point/one-sided trace samples on rounded midpoint+separate edge vectors. Not certified exact face integrals; repeated coordinates explicitly reported.'),external_exact_faces=len(external),external_prior_policy='Existing j<0 physical mirror boundary is not a density stencil neighbor; frozen reservoir cells within the patch ARE included on both graphs. No completion IDs outside supplied fit patch were discarded.',input_sha256=dict(raw=sha(raw),state=sha(state),ids=sha(ids),old_polygons=sha(np.concatenate(polygons)),common_working_vertex_errors=sha(np.concatenate(errors)),mass_error=sha(mass_error)),mass_error_scope='Stored raw treated as exact common input for this attribution only; not zero physical mass error.',coordinate_error_scope='Same working affine-roundoff arrays in A/B, not a certified bound on all upstream clipping errors.',scope=__doc__)
        reports[label]=report;(folder/'report.json').write_text(json.dumps(report,indent=2));emit(dict(stage=label,status='completed',report=report))
    summary=dict(status='completed',seconds=time.monotonic()-started,scope=__doc__,stages=reports,limits='Closure attribution at fixed old geometry only; no native geometry provider adopted, no complete RHS improvement or long Kepler stability claim.')
    (args.output/'summary.json').write_text(json.dumps(summary,indent=2));emit(dict(status='completed',seconds=summary['seconds']))
except BaseException:
    (args.output/'partial_report.json').write_text(json.dumps(reports,indent=2));raise
