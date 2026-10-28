"""Bounded .015->.02 map residual from preserved native base/half/next.

Exact rational interior polygon areas and face coefficients. Generator inverse
intervals include binary64 multiply/add rounding; arc cells remain uncertified.
No RHS, fits, root solve or evolution.
"""
import argparse,hashlib,json,math,struct
from fractions import Fraction as F
from pathlib import Path
import numpy as np


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):
    with np.load(p) as f:return {k:f[k] for k in f.files}
def rn_interval(v):
    assert math.isfinite(v)
    prev=float(np.nextafter(v,-np.inf));nxt=float(np.nextafter(v,np.inf))
    assert math.isfinite(prev) and math.isfinite(nxt)
    return (F(prev)+F(v))/2,(F(v)+F(nxt))/2

def floor_float(q):
    v=float(q);assert math.isfinite(v)
    if F(v)>q:v=float(np.nextafter(v,-np.inf))
    assert F(v)<=q
    return v

def ceil_float(q):
    v=float(q);assert math.isfinite(v)
    if F(v)<q:v=float(np.nextafter(v,np.inf))
    assert F(v)>=q
    return v

def inverse(a,b,h):
    lo,hi=rn_interval(b);pl=ceil_float(lo-F(a));pu=floor_float(hi-F(a));assert pl<=pu
    low=rn_interval(pl)[0]/h;high=rn_interval(pu)[1]/h
    # Also enclose FMA and possible zero/subnormal flush. This widens rather
    # than asserting an unobserved process MXCSR setting.
    low=min(low,(lo-F(a))/h);high=max(high,(hi-F(a))/h)
    if pl<=0<=pu:
        guard=F(float(np.finfo(float).tiny))/h
        low=min(low,-guard);high=max(high,guard)
    center=(F(b)-F(a))/h;assert low<=center<=high
    return low,high,center


def parse(cell):
    edges=cell['edges'];A=[tuple(F(v) for v in e['A']) for e in edges]
    B=[tuple(F(v) for v in e['B']) for e in edges]
    assert all(B[k]==A[(k+1)%len(A)] for k in range(len(A)))
    inside=all(sum(v*v for v in p)<=F(25,4) for p in A)
    return dict(owner=int(cell['owner']),edges=edges,A=A,B=B,inside=inside)

def area(cell):
    value=sum(a[0]*b[1]-a[1]*b[0] for a,b in zip(cell['A'],cell['B']))/2
    assert value>0
    return value


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--geometry-root',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(exist_ok=False)
    source=json.loads((args.geometry_root/'report.json').read_text());assert source['status']=='completed_three_saved_native_geometries'
    h=F(float(source['h']));assert float(h)==.005
    states={k:load(args.geometry_root/k/'state.npz') for k in ('base','half','next')}
    meshes={k:load(args.geometry_root/k/'geometry/mesh.npz') for k in states}
    raw={k:json.loads((args.geometry_root/k/'geometry/full_exact_supports.json').read_text()) for k in states}
    assert all(d['shared_exact_reversed'] for d in raw.values())
    polys={k:[parse(c) for c in d['cells']] for k,d in raw.items()}
    ids=states['base']['ids'];mass=states['base']['mass'];n=len(ids);lookup={int(v):i for i,v in enumerate(ids)}
    coords={k:states[k]['state'][:-4].reshape(n,20) for k in states}
    assert all(np.array_equal(ids,s['ids']) and np.array_equal(mass,s['mass']) for s in states.values())
    assert all([c['owner'] for c in d]==ids.tolist() for d in polys.values())
    base,half,nxt=coords['base'],coords['half'],coords['next'];z=states['next']['weight_rate_seed']
    assert np.array_equal((base[:,2]+float(h)*z).view(np.uint64),nxt[:,2].view(np.uint64))
    strict=np.array([all(polys[k][i]['inside'] and meshes[k]['arc_counts'][i]==0 for k in polys) for i in range(n)])
    g=[]
    with (args.output/'generator_intervals.jsonl').open('x') as f:
        for i in range(n):
            row=[inverse(float(base[i,j]),float(nxt[i,j]),h) for j in (0,1)];g.append(row)
            f.write(json.dumps(dict(id=int(ids[i]),lower=[str(v[0]) for v in row],upper=[str(v[1]) for v in row],center=[str(v[2]) for v in row]))+'\n')
    positions=[tuple(F(float(v)) for v in row[:2]) for row in half];weights=[F(float(v)) for v in half[:,2]];zf=[F(float(v)) for v in z]
    coefficients=[{} for _ in ids];constant=[F(0) for _ in ids];face_count=0
    def add(cell,site,vector,sign):
        if not strict[cell]:return
        old=coefficients[cell].setdefault(site,[F(0),F(0)])
        for axis in (0,1):old[axis]+=sign*vector[axis]
    for i,cell in enumerate(polys['half']):
        for edge,a,b in zip(cell['edges'],cell['A'],cell['B']):
            neighbor=int(edge['neighbor'])
            if neighbor<0 or int(ids[i])>neighbor:continue
            j=lookup[neighbor]
            if not strict[i] and not strict[j]:continue
            d=tuple(positions[j][k]-positions[i][k] for k in (0,1));e=tuple(b[k]-a[k] for k in (0,1));mid=tuple((a[k]+b[k])/2 for k in (0,1))
            norm2=sum(v*v for v in d);lam=(e[1]*d[0]-e[0]*d[1])/norm2
            assert lam>0 and (e[1],-e[0])==tuple(lam*v for v in d)
            plane=(sum(v*v for v in positions[j])-sum(v*v for v in positions[i])+weights[i]-weights[j])/2
            assert sum(d[k]*a[k] for k in (0,1))==plane and sum(d[k]*b[k] for k in (0,1))==plane
            ci=tuple(lam*(mid[k]-positions[i][k]) for k in (0,1));cj=tuple(lam*(positions[j][k]-mid[k]) for k in (0,1));c=lam*(zf[i]-zf[j])/2
            add(i,i,ci,1);add(i,j,cj,1);add(j,i,ci,-1);add(j,j,cj,-1)
            if strict[i]:constant[i]+=c
            if strict[j]:constant[j]-=c
            face_count+=1
    fields={k:np.full(n,np.nan) for k in ('delta_estimate','delta_lower','delta_upper','relative_delta','relative_lower','relative_upper','DeltaV','hVdot','bound_width','density_error_scale')}
    excludes=np.zeros(n,bool);exact_area=np.full(n,np.nan)
    with (args.output/'rational_bounds.jsonl').open('x') as f:
        for i in np.flatnonzero(strict):
            lo=hi=center=constant[i]
            for site,cs in coefficients[i].items():
                for axis,c in enumerate(cs):
                    a,b,m=g[site][axis];center+=c*m
                    lo+=c*(a if c>=0 else b);hi+=c*(b if c>=0 else a)
            A,B=area(polys['base'][i]),area(polys['next'][i]);dv=B-A
            dlo,dhi=dv-h*hi,dv-h*lo;estimate=dv-h*center;assert dlo<=estimate<=dhi
            excludes[i]=dlo>0 or dhi<0;exact_area[i]=float(A)
            vals=dict(delta_estimate=estimate,delta_lower=dlo,delta_upper=dhi,relative_delta=estimate/A,relative_lower=dlo/A,relative_upper=dhi/A,DeltaV=dv,hVdot=h*center,bound_width=dhi-dlo,density_error_scale=-estimate/A)
            for key,value in vals.items():fields[key][i]=floor_float(value) if key.endswith('lower') else ceil_float(value) if key.endswith('upper') or key=='bound_width' else float(value)
            f.write(json.dumps(dict(id=int(ids[i]),delta_lower=str(dlo),delta_upper=str(dhi),area_base=str(A),area_next=str(B),Vdot_center=str(center)))+'\n')
    # Empirical whole-disk check; no certification is claimed for arc cells.
    mesh=meshes['half'];gen=np.array([[float(v[2]) for v in row] for row in g],dtype=np.longdouble)
    ij=np.array([(lookup[int(a)],lookup[int(b)]) for a,b in mesh['face_keys']]);i,j=ij.T
    face=mesh['faces'].astype(np.longdouble);x=half[:,:2].astype(np.longdouble);d=x[j]-x[i];dg=gen[j]-gen[i]
    distance=np.sqrt((d*d).sum(1));wn=(.5*((x[j]+x[i]-2*face[:,:2])*dg).sum(1)+.5*(d*(gen[j]+gen[i])).sum(1)-.5*(z[j]-z[i]))/distance
    Vdot=np.zeros(n,np.longdouble);np.add.at(Vdot,i,face[:,4]*wn);np.add.at(Vdot,j,-face[:,4]*wn)
    empirical=meshes['next']['moments'][:,0].astype(np.longdouble)-meshes['base']['moments'][:,0].astype(np.longdouble)-np.longdouble(float(h))*Vdot
    # Same radial-envelope definition as the existing physical monitor.
    lo_radius=[];hi_radius=[]
    for cell in polys['half']:
        minimum=[];maximum=[];inside=True
        for a,b in zip(cell['A'],cell['B']):
            e=tuple(b[k]-a[k] for k in (0,1));norm=sum(v*v for v in e);normal=(e[1],-e[0])
            if sum(-normal[k]*a[k] for k in (0,1))>0:inside=False
            t=max(F(0),min(F(1),-sum(a[k]*e[k] for k in (0,1))/norm));point=tuple(a[k]+t*e[k] for k in (0,1))
            minimum.append(sum(v*v for v in point));maximum.extend(sum(v*v for v in p) for p in (a,b))
        if inside:minimum.append(F(0))
        lo_radius.append(math.sqrt(float(min(minimum))));hi_radius.append(min(2.5,math.sqrt(float(max(maximum)))))
    low=np.array(lo_radius);high=np.array(hi_radius);wall=mesh['arc_counts']>0;pad=32*np.finfo(float).eps*np.maximum(1,high)
    cin=(low-pad<=.5)&(high+pad>=.5);cout=(low-pad<=2)&(high+pad>=2)
    ic=load(Path('/gpfs/kjhan/LagForce/Codex_full_curved_kepler_IC_20261019/initial.npz'));r0=np.linalg.norm(ic['state'][:-4].reshape(n,20)[:,:2],axis=1)
    regions=dict(all=np.ones(n,bool),inner=r0<.4,cavity=(high<.5)&~cin,disk_core=(low>.5)&(high<2)&~cin&~cout,inner_contact=cin,outer_contact=cout,dynamic_buffer=(low>2)&~wall&~cout,wall=wall)
    summary={}
    for name,mask in regions.items():
        valid=mask&strict;row=dict(cells=int(mask.sum()),certified_cells=int(valid.sum()),geometry_bound_unavailable=int((mask&~strict).sum()),zero_excluded_cells=int((valid&excludes).sum()))
        if valid.any():
            w=mass[valid].astype(np.longdouble);v=fields['relative_delta'][valid];k=int(np.argmax(abs(v)))
            row.update(relative_delta_RMS=float(np.sqrt(np.sum(w*v*v)/w.sum())),relative_delta_signed_mean=float(np.sum(w*v)/w.sum()),maximum_relative_delta=float(abs(v[k])),maximum_id=int(ids[valid][k]),
                maximum_relative_bound_width=float(np.max(fields['bound_width'][valid]/exact_area[valid])),signed_delta_sum=float(fields['delta_estimate'][valid].astype(np.longdouble).sum()))
        summary[name]=row
    topology=np.array([set(e['neighbor'] for e in polys['base'][i]['edges'])!=set(e['neighbor'] for e in polys['next'][i]['edges']) for i in range(n)])
    np.savez(args.output/'map_residuals.npz',ids=ids,certified=strict,geometry_bound_unavailable=~strict,zero_excluded=excludes,topology_changed=topology,empirical_all_cell_delta=empirical,**fields,**regions)
    report=dict(status='completed_bounded_reconstructed_native_volume_map_diagnostic',start=.015,half=.0175,end=.02,h=float(h),source_geometry_report_sha256=sha(args.geometry_root/'report.json'),regions=summary,
        exact_canonical_faces_used=face_count,generator_contract='Exact rational RN preimages for multiply/add; also encloses FMA and zero-product FTZ/DAZ possibilities; no single inverse value used as exact stage generator',
        exact_geometry_contract='Interior original rational power polygons only; rational area and length/distance-cancelled affine face integral; global coefficients assembled before interval evaluation',
        arc_cells_certified=False,total_topology_changed=int(topology.sum()),
        empirical_global_Vdot_sum=float(Vdot.sum()),empirical_base_disk_area=float(meshes['base']['moments'][:,0].astype(np.longdouble).sum()),empirical_next_disk_area=float(meshes['next']['moments'][:,0].astype(np.longdouble).sum()),
        physical_providers=0,density_fits=0,root_solves=0,accepted_time_steps=0,
        limits='One native interval only, with reconstructed bounded generator, not original stored face-wn. Arc/unsupported geometry excluded from strict statistics. Region envelopes on preclip polygons are working masks; boundary cells are all uncertified. Map residual is not proof of a bug, dominant physical-error cause, time order or long-time stability.')
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(dict(status=report['status'],regions=summary)),flush=True)


if __name__=='__main__':main()
