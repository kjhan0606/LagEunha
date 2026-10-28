"""Potential transport identity versus finite-difference reconstructed potential.
Uses saved rates only; does not advance or overwrite a simulation.
"""
import argparse,hashlib,importlib.util,json,os,sys
from pathlib import Path
import numpy as np
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--polynomial-control',action='store_true');args=parser.parse_args()
base=Path('/gpfs/kjhan/LagForce');sys.path.insert(0,str(base));sys.path.insert(0,str(base/'Codex_coupled_relaxation_evolution_20261018'));sys.path.insert(0,str(base/'Codex_split_face_modules_20261018'))
import moment_ale_evolve as model
from affine_upwind_face_nodes import AffineUpwindNodes
from kinetic_factor_rhs import ALPHA
from run_canonical_av_branches_20261018 import read
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
geom=load('geom',base/'Codex_weight_remap_moments_centered_20261017/audit.py')
clip=load('clip',base/'Codex_weighted_compression_geometry_20261017/audit.py')
path=base/'Codex_continuous_initial_rhs_split_20261018/rhs_r0.55_f8.npz';saved=np.load(path)
ids=saved['ids'];mass=saved['mass'];y=saved['state'][:-4].reshape(-1,20);rates=saved['rates'][:-4].reshape(-1,20)
raw=np.column_stack([mass,y[:,3:]]);active=~saved['frozen'];factors=saved['covariance_factors']
ic,_=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');allx=np.column_stack([ic['x']-3,ic['y']-3]);x=allx[ids]
faces=np.loadtxt(base/'Codex_pair_hg_modes_t0_20261018/raw_faces.dat.rank0000');order=np.argsort(faces[:,0],kind='stable');owners=faces[order,0];lookup={int(k):i for i,k in enumerate(ids)}
os.environ['LAGFORCE_DENSITY_QUADRATURE_ORDER']='12';model.GAMMA=float(np.float32(1.666667))
polys=[];edges=[]
for local,i in enumerate(ids):
    ff=faces[order[np.searchsorted(owners,i):np.searchsorted(owners,i,side='right')]]
    vv=ff[:,6:8];polys.append(vv[np.argsort(np.arctan2(vv[:,1],vv[:,0]))]+x[local])
    for row in ff:
        j=lookup.get(int(row[1]),-1)
        if j>=0 and local>=j:continue
        a=x[local]+row[6:8];b=x[local]+row[8:10];delta=allx[int(row[1])]-x[local]
        edges.append((local,j,a,b,delta/np.linalg.norm(delta),np.linalg.norm(b-a)))
states=[model.decode(p,r) for p,r in zip(polys,raw)]
gravity=np.load(base/'Codex_continuous_initial_gravity_v2_20261018/gravity_q24.npz')['source']
root=base/('Codex_potential_transport_polynomial_20261018' if args.polynomial_control else 'Codex_potential_transport_closure_20261018');root.mkdir(exist_ok=False)
fluxes=[]
controls=[]
for nq in (8,16):
    provider=AffineUpwindNodes(x,mass,factors,rates[:,2],nq,lambda p:np.zeros(2));flow=np.zeros(len(ids));mass_rate=np.zeros(len(ids));polynomial_flux=np.zeros((len(ids),4))
    for edge_id,edge in enumerate(edges):
        i,j,a,b,n,length=edge
        if not active[i] and (j<0 or not active[j]):continue
        assert j>=0,'active patch touches external edge'
        nodes=provider(edge_id,edge,states,rates[:,:2]);ts=np.array([g for g,w in nodes]);dw=length*np.array([w for g,w in nodes]);pts=a+ts[:,None]*(b-a)
        d=np.linalg.norm(x[j]-x[i]);wn=((x[j]-x[i])@rates[j,:2]-(pts-x[i])@(rates[j,:2]-rates[i,:2]))/d+(rates[i,2]-rates[j,2])/(2*d)
        samples=[]
        for k in (i,j):
            s=states[k];vel=s['v']+(pts-s['c'])@s['B'].T
            # Broadcasting requires one point axis and one sample axis.
            normal=(vel@n)[:,None]+((ALPHA@factors[k].T/np.sqrt(mass[k]))@n)[None,:]-wn[:,None]
            samples.append(s['density'].value(pts)*.25*(np.maximum(normal,0) if k==i else np.minimum(normal,0)).sum(axis=1))
        fm=samples[0]+samples[1];amount=float(dw@fm);potential=float(dw@(fm*geom.potential(pts)))
        mass_rate[i]-=amount;mass_rate[j]+=amount;flow[i]-=potential;flow[j]+=potential
        features=np.column_stack([np.ones(len(pts)),pts,np.sum(pts*pts,axis=1)])
        flux=(dw*fm)@features;polynomial_flux[i]-=flux;polynomial_flux[j]+=flux
    volume=np.column_stack([np.zeros(len(ids)),raw[:,7:9],2*(raw[:,9]+raw[:,12])])
    expected=np.column_stack([np.zeros(len(ids)),rates[:,3:5],rates[:,5]+rates[:,8]])
    control=np.max(abs((expected-polynomial_flux-volume)[active]),axis=0)
    controls.append(dict(face_order=nq,polynomials=['1','x','y','r_squared'],maximum_derivative_defects=control.tolist()))
    fluxes.append(flow);print(json.dumps(dict(face_order=nq,active_mass_residual=float(np.max(abs(mass_rate[active]))))),flush=True)
if args.polynomial_control:
    report=dict(scope='Polynomial potential transport control at saved identical rates',results=controls,
        limits='Order8 uses the same quadrature as saved RHS; order16 is an independent quadrature check. Does not establish exact closure for spline potential.')
    (root/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report));sys.exit(0)
physical=fluxes[-1]-gravity[:,17]
def reconstructed_potential(eps):
    trial=y+eps*rates;positions=allx.copy();positions[ids]=trial[:,:2];weights=np.zeros(len(allx));weights[ids]=trial[:,2]
    cells=geom.polygons(positions,weights,ids,clip.clip);values=[]
    for i in np.flatnonzero(active):
        p=cells[i];rr=np.r_[mass[i],trial[i,3:]]
        density=model.Density(p,rr[0],rr[1:3],rr[3:7].reshape(2,2))
        pts,dv=geom.quadrature(p,24);values.append(float((dv*density.value(pts))@geom.potential(pts)))
    return np.array(values)
rows=[];previous=None
for eps in (1e-4,5e-5,2.5e-5):
    minus=reconstructed_potential(-eps);plus=reconstructed_potential(eps);derivative=(plus-minus)/(2*eps)
    defect=derivative-physical[active]
    row=dict(epsilon=eps,maximum_closure_defect=float(np.max(abs(defect))),closure_defect_l2=float(np.linalg.norm(defect)),
        summed_closure_defect=float(defect.sum()),physical_derivative_l2=float(np.linalg.norm(physical[active])))
    if previous is not None:row['successive_derivative_max_abs']=float(np.max(abs(derivative-previous)))
    previous=derivative;rows.append(row);np.savez(root/f'potential_eps{eps:g}.npz',ids=ids[active],minus=minus,plus=plus,derivative=derivative,physical=physical[active]);print(json.dumps(row),flush=True)
report=dict(scope=__doc__,source=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),results=rows,
    face8_to16_max_potential_flux_difference=float(np.max(abs((fluxes[1]-fluxes[0])[active]))),
    limits='Directional derivative of the fitted-density representation on perturbed Laguerre geometry. Not evolution. A nonzero defect indicates the potential functional is not exactly transported by the truncated moments; do not cancel by silently heating cells.')
(root/'report.json').write_text(json.dumps(report,indent=2))
