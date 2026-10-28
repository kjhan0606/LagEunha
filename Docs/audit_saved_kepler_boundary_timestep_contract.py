"""Source/IC boundary contract and saved-stage timestep design diagnostics.

No timestep increase, domain extension, geometry solve or evolution.
"""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
base=Path('/gpfs/kjhan/LagForce');sys.path[:0]=[str(base)]
from run_canonical_av_branches_20261018 import read
from audit_canonical_pair_history import gravity_coefficient
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();args.output.mkdir(exist_ok=False)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def stats(value,mask):
 a=np.asarray(value)[mask];return dict(cells=len(a),minimum=float(a.min()),median=float(np.median(a)),p90=float(np.quantile(a,.9)),maximum=float(a.max())) if len(a) else dict(cells=0)
icpath=base/'Codex_canonical_initial_20261017/run/khout.000000.dat';ic,t=read(icpath);assert t==0.
envpath=icpath.parent/'environment.json';env=json.loads(envpath.read_text())
x=np.column_stack([ic['x']-float(env['EUNHA_KEPLER_CX']),ic['y']-float(env['EUNHA_KEPLER_CY'])]);radius=np.linalg.norm(x,axis=1)
velocity=np.column_stack([ic['vx'],ic['vy']]);expected=np.sqrt(gravity_coefficient(radius))[:,None]*np.column_stack([-x[:,1],x[:,0]])
rin=float(env['EUNHA_KEPLER_RIN']);rout=float(env['EUNHA_KEPLER_ROUT']);floor=float(env['EUNHA_KEPLER_FLOOR'])
rho_expected=np.where((radius>=rin)&(radius<=rout),1.,floor)
full=dict(cells=len(ic),bounds_minimum=x.min(0).tolist(),bounds_maximum=x.max(0).tolist(),
 radial_maximum=float(radius.max()),density_min=float(ic['den'].min()),density_max=float(ic['den'].max()),
 density_profile_error_max=float(abs(ic['den']-rho_expected).max()),pressure_min=float(ic['pressure'].min()),pressure_max=float(ic['pressure'].max()),
 circular_velocity_error_max=float(np.linalg.norm(velocity-expected,axis=1).max()),
 weights_zero=bool(np.all(ic['w2']==0)),outside_original_evolve_radius_cells=int((radius>=float(env['EUNHA_KEPLER_EVOLVE_R'])).sum()),
 within_disk_cells=int(((radius>=rin)&(radius<=rout)).sum()))
# Evaluate exact steady point-field limits at opposite box boundaries. This
# checks flux compatibility; it does not synthesize new numerical BC states.
g,w=np.polynomial.legendre.leggauss(48);coordinate=3*g
left=np.column_stack([np.full(len(g),-3.),coordinate]);right=left.copy();right[:,0]=3.
def point_state(points):
 r=np.linalg.norm(points,axis=1);rho=np.where((r>=rin)&(r<=rout),1.,floor)
 v=np.sqrt(gravity_coefficient(r))[:,None]*np.column_stack([-points[:,1],points[:,0]])
 P=float(ic['pressure'][0]);flux=np.column_stack([rho*v[:,0],rho*v[:,0]**2+P,rho*v[:,0]*v[:,1]])
 return rho,v,flux
rl,vl,fl=point_state(left);rr,vr,fr=point_state(right)
boundary=dict(original_source_BC='KH non-Sod periodic padding/post-stage; outer r>=2.5 has RK increments frozen in original Kepler source.',
 moment_patch_BC='Embedded nonperiodic power cells; j<0 uses zero-normal-speed mirror states; frozen r>=.60 raw storage masked. Not a physical long-time outer BC.',
 point_field_opposite_x_velocity_jump_max=float(np.linalg.norm(vr-vl,axis=1).max()),
 point_field_opposite_x_mass_flux_jump_max=float(abs(fr[:,0]-fl[:,0]).max()),
 point_field_opposite_x_tangential_momentum_flux_jump_max=float(abs(fr[:,2]-fl[:,2]).max()),
 point_field_box_normal_velocity_max=float(abs(vl[:,0]).max()),
 conclusion='The globally rotating floor has nonzero normal box velocity and an opposite-edge tangential momentum flux jump. Neither naive periodic identification nor stationary reflecting box walls reproduce this whole steady Kepler field. Dynamic buffer and an explicit compatible outer physical BC remain required.',
 full_builder_limit='Current EmbeddedLaguerreGeometry rejects unmatched bounding-box faces; it cannot simply promote all registry IDs to active boundary cells.')
root=base/'Codex_compiled_midpoint_once_20261018';h=float(json.loads((root/'manifest.json').read_text())['dt']);reports={}
for label,update in (('stage0','half_update'),('midpoint','full_update')):
 data=np.load(root/label/'inputs.npz');reference=np.load(root/label/'reference.npz');nodes=data['node_inputs'];owners=data['owners'];mass=data['mass'];N=len(mass)
 frozen=data['frozen'].astype(bool);active=~frozen;a=data['state'][:-4].reshape(-1,20);positions=a[:,:2];gen=data['generator'];z=data['weight_rates'];i=owners[:,0];j=owners[:,1];internal=j>=0
 normal=nodes[:,2:4];point=nodes[:,:2];scale=nodes[:,4];wn=np.zeros(len(nodes));ii=i[internal];jj=j[internal]
 delta=positions[jj]-positions[ii];distance=np.linalg.norm(delta,axis=1)
 wn[internal]=((delta[:,None,:]@gen[jj,:,None])[:,0,0]-((gen[jj]-gen[ii])[:,None,:]@(point[internal]-positions[ii])[:,:,None])[:,0,0])/distance+(z[ii]-z[jj])/(2*distance)
 vL=nodes[:,5:13].reshape(-1,4,2);vR=nodes[:,13:21].reshape(-1,4,2)
 relL=(vL@normal[:,:,None])[:,:,0]-wn[:,None];relR=(vR@normal[:,:,None])[:,:,0]-wn[:,None]
 gamma=float(np.float32(1.666667));cL=np.sqrt(gamma*(gamma-1)*nodes[:,29]);cR=np.sqrt(gamma*(gamma-1)*nodes[:,30])
 outgoing=np.bincount(i,weights=scale*nodes[:,21:25].sum(1),minlength=N)+np.bincount(j[internal],weights=-scale[internal]*nodes[internal,25:29].sum(1),minlength=N)
 relmax=np.zeros(N);signal=np.zeros(N);sound=np.zeros(N)
 np.maximum.at(relmax,i,abs(relL).max(1));np.maximum.at(relmax,j[internal],abs(relR[internal]).max(1))
 np.maximum.at(signal,i,abs(relL).max(1)+cL);np.maximum.at(signal,j[internal],abs(relR[internal]).max(1)+cR[internal]);np.maximum.at(sound,i,cL);np.maximum.at(sound,j[internal],cR[internal])
 archive=base/'Codex_adjacency_prior_impact_20261018'/label;common=np.load(archive/'common_inputs.npz');off=common['polygon_offsets'];verts=common['polygon_vertices'];polys=[verts[a:b] for a,b in zip(off[:-1],off[1:])]
 areas=[];widths=[];edge_min=[];hydraulic=[]
 for p in polys:
  d=np.roll(p,-1,axis=0)-p;lengths=np.linalg.norm(d,axis=1);good=lengths>0;normals=np.column_stack([d[good,1],-d[good,0]])/lengths[good,None]
  local=p-p[0];area=.5*np.sum(local[:,0]*np.roll(local[:,1],-1)-local[:,1]*np.roll(local[:,0],-1))
  widths.append(np.ptp(p@normals.T,axis=0).min());areas.append(area);edge_min.append(lengths[good].min());hydraulic.append(2*area/lengths.sum())
 widths=np.array(widths);hydraulic=np.array(hydraulic)
 B=reference['B'];strain=.5*(B+B.transpose(0,2,1));center=data['cells'][:,1:3]/mass[:,None]
 factor=np.load(root/update/'cuda.npy');Ac=factor[:,26:30].reshape(-1,2,2)
 quantities=dict(outgoing_mass_fraction=h*outgoing/mass,relative_sample_speed=relmax,sound_speed=sound,
  signal_Courant_hydraulic=h*signal/hydraulic,signal_Courant_minimum_width=h*signal/widths,
  omega_dt=h*np.sqrt(gravity_coefficient(np.linalg.norm(center,axis=1))),strain_norm_dt=h*np.linalg.norm(strain,ord=2,axis=(1,2)),
  covariance_generator_norm_dt=h*np.linalg.norm(reference['cov_A'],ord=2,axis=(1,2)),effective_factor_generator_norm_dt=h*np.linalg.norm(Ac,ord=2,axis=(1,2)),
  minimum_width=widths,hydraulic_length=hydraulic,minimum_encoded_edge=np.array(edge_min),area=np.array(areas))
 np.savez(args.output/(label+'.npz'),ids=data['ids'],active=active,**quantities)
 reports[label]={name:dict(active=stats(v,active),frozen=stats(v,frozen)) for name,v in quantities.items()}
factor_full=len(ic)/len(mass);stage_seconds=109.03531125560403;patch_days=10/h*2*stage_seconds/86400
cost=dict(patch_cells=len(mass),full_registry_cells=len(ic),cell_ratio=factor_full,patch_stage_seconds_observed=stage_seconds,
 fixed_dt=h,patch_t10_days_two_stages=patch_days,
 full_registry_t10_days_naive_cell_linear=patch_days*factor_full,
 full_registry_stage_seconds_naive_cell_linear=stage_seconds*factor_full,
 full_registry_peak_RSS_MiB_naive_linear=1721.34765625*factor_full,
 node_scaling_caveat='Only current patch node inventory known. Full-domain face count, adaptive face splits, priors and outer BC not built; cell-linear values are scenarios, not measurements or rigorous lower/upper bounds. Geometry global registry already contributes fixed cost. Native face-record and production diagnostic policies may change cost.')
paths=[Path(__file__),icpath,envpath,root/'manifest.json',base/'Codex_canonical_initial_20261017/kh.c',base/'Codex_canonical_initial_20261017/run/params.dat',base/'Codex_accepted_ledger_build_20261017/src/Exam/KH/util.c',base/'Codex_accepted_ledger_build_20261017/src/Exam/exam.c']
report=dict(status='completed_saved_contract_analysis',scope=__doc__,full_initial_registry=full,boundary_contract=boundary,timestep=reports,cost=cost,source_sha256={str(p):sha(p) for p in paths},limits='Dimensionless saved-state diagnostics do not guarantee positivity/realizability and do not authorize a larger timestep. No new full-domain topology, boundary extension or trajectory.')
(args.output/'report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report,indent=2))
