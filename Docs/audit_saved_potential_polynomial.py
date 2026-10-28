"""Saved degree-two potential identities and direct remainder decomposition.

No fitted state is recomputed; archived theta is restored. No RHS or evolution.
Local Taylor coefficients are frozen at each actual midpoint cell centroid.
"""
import argparse,hashlib,importlib.util,json,sys
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce')
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args()
sys.path[:0]=[str(base/'Codex_compiled_midpoint_modules_20261018'),
    str(base/'Codex_moment_flux_batch_modules_20261018'),str(base/'Codex_face_state_batch_api_modules_20261018'),
    str(base/'Codex_coupled_relaxation_evolution_20261018'),str(base)]
from saved_potential_taylor import taylor_coefficients,p2,raw_p2,raw_volume,integrate_remainder
from saved_contact_density_restore import restore_density
from audit_canonical_pair_history import gravity_coefficient
def load(name,path):
 spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
potential_path=base/'Codex_weight_remap_moments_centered_20261017/audit.py';potential=load('saved_potential',potential_path).potential
root=base/'Codex_compiled_midpoint_once_20261018';endpoints=base/'Codex_numerical_Psi_endpoints_v2_20261018'
mid=np.load(root/'midpoint/inputs.npz');N=len(mid['mass']);mass=mid['mass'];active=~mid['frozen'].astype(bool)
raw_mid=mid['cells'][:,:18];coeff=taylor_coefficients(raw_mid[:,1:3]/mass[:,None],potential,gravity_coefficient)
args.output.mkdir(exist_ok=False)
np.savez(args.output/'coefficients.npz',**coeff,ids=mid['ids'])
manifest=dict(scope=__doc__,provider_calls=0,newton_fits=0,source_sha256={str(p):sha(p) for p in (Path(__file__),potential_path,base/'Codex_moment_flux_batch_modules_20261018/saved_potential_taylor.py')})
manifest['runtime_modules']={name:dict(path=str(m.__file__),sha256=sha(m.__file__)) for name,m in list(sys.modules.items()) if getattr(m,'__file__',None) and Path(m.__file__).suffix=='.py' and Path(m.__file__).is_relative_to(base)}
manifest['fit_archives']={str(p):sha(p) for p in [endpoints/(label+suffix) for label in ('initial','accepted') for suffix in ('_reconstruction.npz','_priors.json','_fits.json')]+[base/'Codex_adjacency_prior_impact_20261018/midpoint'/name for name in ('common_inputs.npz','priors.json','fit_residuals.json')]}
(args.output/'manifest.json').write_text(json.dumps(manifest,indent=2))
results={};passed=True
def summary(v):
 v=np.asarray(v)[active]
 return dict(signed_sum=float(v.sum()),absolute_sum=float(abs(v).sum()),maximum=float(abs(v).max(initial=0)))
for label in ('stage0','midpoint'):
 data=np.load(root/label/'inputs.npz');reference=np.load(root/label/'reference.npz')
 nodes=data['node_inputs'];inc=data['incidence'];off=data['offsets'];gn=data['gravity_nodes'];goff=data['gravity_offsets']
 raw=data['cells'][:,:18];rate=reference['original_hydro'];wl=nodes[:,21:25];wr=nodes[:,25:29]
 # Preserve scalar grouping: sum(wl)+sum(wr), then multiply by scale.
 fm=nodes[:,4]*((((wl[:,0]+wl[:,1])+wl[:,2])+wl[:,3])+(((wr[:,0]+wr[:,1])+wr[:,2])+wr[:,3]))
 inventories={}
 for name in ('constant','x','y','xx_half','xy','yy_half','local_Taylor'):
  if name=='local_Taylor':c=coeff
  else:
   c=dict(center=np.zeros((N,2)),Phi=np.zeros(N),gradient=np.zeros((N,2)),H=np.zeros((N,2,2)))
   if name=='constant':c['Phi'][:]=1.
   elif name in ('x','y'):c['gradient'][:,0 if name=='x' else 1]=1.
   elif name=='xx_half':c['H'][:,0,0]=1.
   elif name=='yy_half':c['H'][:,1,1]=1.
   else:c['H'][:,0,1]=1.;c['H'][:,1,0]=1.
  face=[]
  for k in range(N):
   codes=inc[off[k]:off[k+1]];idx=abs(codes)-1;sign=np.where(codes<0,-1.,1.)
   face.append(np.sum(sign.astype(np.longdouble)*fm[idx]*p2(nodes[idx,:2],c,k)))
  face=np.array(face,dtype=np.longdouble);volume=raw_volume(raw,c);target=raw_p2(rate,c)
  defect=target-face-volume;bound=1e-14+1e-10*abs(target)
  okay=bool(np.all(abs(defect[active])<=bound[active]));passed &= okay
  inventories[name]=dict(passed=okay,identity_defect=summary(defect),mass_rate_residual=summary(rate[:,0]))
  if name=='local_Taylor':
   rem_face=[];rem_volume=[];vol_quad=[]
   for k in range(N):
    codes=inc[off[k]:off[k+1]];idx=abs(codes)-1;sign=np.where(codes<0,-1.,1.)
    point=nodes[idx,:2].astype(np.longdouble)
    rem_face.append(np.sum(sign*fm[idx]*(potential(point)-p2(point,c,k))))
    cells=gn[goff[k]:goff[k+1]].astype(np.longdouble);x=cells[:,:2];dm=cells[:,2];v=cells[:,3:5];gradP=c['gradient'][k]+(x-c['center'][k])@c['H'][k].T
    vol_quad.append(np.sum(dm*np.sum(v*gradP,axis=1)))
    rem_volume.append(np.sum(dm*np.sum(v*(-cells[:,5:7]-gradP),axis=1)))
   rem_face=np.array(rem_face);rem_volume=np.array(rem_volume);vol_quad=np.array(vol_quad)
   np.savez(args.output/(label+'_transport.npz'),face_P2=face,volume_P2_raw=volume,
    raw_P2_rate=target,identity_defect=defect,face_remainder=rem_face,volume_remainder=rem_volume,
    volume_P2_quadrature=vol_quad,gravity_energy_rate=reference['gravity'][:,17],mass_rate=rate[:,0])
   inventories[name].update(volume_quad_minus_raw=summary(vol_quad-volume),
    remainder_transport_rate=summary(rem_face+rem_volume))
 results[label]=inventories
 print(json.dumps(dict(stage=label,identities=inventories)),flush=True)

for label in ('initial','midpoint','accepted'):
 if label=='midpoint':
  folder=base/'Codex_adjacency_prior_impact_20261018/midpoint';common=np.load(folder/'common_inputs.npz')
  raw=common['raw'];offset=common['polygon_offsets'];vertices=common['polygon_vertices']
  priors=json.loads((folder/'priors.json').read_text())['old'];fits=[f['old'] for f in json.loads((folder/'fit_residuals.json').read_text())]
 else:
  saved=np.load(endpoints/(label+'_reconstruction.npz'));raw=np.column_stack([mass,saved['state'][:-4].reshape(-1,20)[:,3:]])
  offset=saved['polygon_offsets'];vertices=saved['polygon_vertices']
  priors=json.loads((endpoints/(label+'_priors.json')).read_text());fits=json.loads((endpoints/(label+'_fits.json')).read_text())
 polys=[vertices[a:b] for a,b in zip(offset[:-1],offset[1:])]
 densities=[restore_density(p,r,m,f) for p,r,m,f in zip(polys,raw,priors,fits)]
 integrals=integrate_remainder(densities,polys,coeff,potential)
 polynomial=raw_p2(raw,coeff)
 np.savez(args.output/(label+'_integrals.npz'),raw_P2=polynomial,integrals=integrals,
    fit_residual=np.array([d.residual for d in densities]))
 results[label+'_reconstruction']=dict(remainder_q48=summary(integrals[1,:,0]),
    remainder_q48_minus_q24=summary(integrals[1,:,0]-integrals[0,:,0]),
    P2_q48_minus_raw=summary(integrals[1,:,1]-polynomial),fit_residual_max=max(d.residual for d in densities))
 print(json.dumps(dict(endpoint=label,completed=True)),flush=True)
h=float(json.loads((root/'manifest.json').read_text())['dt'])
first=np.load(args.output/'initial_integrals.npz');last=np.load(args.output/'accepted_integrals.npz');transport=np.load(args.output/'midpoint_transport.npz')
polynomial_increment=last['raw_P2']-first['raw_P2'];remainder_increment=last['integrals'][1,:,0]-first['integrals'][1,:,0]
P2transport=h*(transport['face_P2']+transport['volume_P2_raw'])
Rtransport=h*(transport['face_remainder']+transport['volume_remainder'])
decomposition=dict(P2_increment_minus_transport=summary(polynomial_increment-P2transport),
    remainder_increment_minus_transport=summary(remainder_increment-Rtransport),
    volume_quadrature_minus_raw_increment=summary(h*(transport['volume_P2_quadrature']-transport['volume_P2_raw'])),
    combined_defect=summary(polynomial_increment+remainder_increment-P2transport-Rtransport),
    note='Each cell has its own fixed Taylor coefficients. Only actual Phi, not P2 alone, admits global shared-face cancellation. Fixed mass endpoints vs stored nonzero dm recorded separately.')
report=dict(status='passed_identities' if passed else 'identity_gate_failed',manifest=manifest,results=results,full_step=decomposition,
    limits='Longdouble remainder arithmetic reduces cancellation; density and its quadrature weights remain original binary64. No evolved Psi claim or epsilon sweep.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(decomposition),flush=True)
assert passed,'Polynomial moment identity failure preserved'
