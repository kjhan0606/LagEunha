"""Actual t2 conditional moment patch: instantaneous hydro and boundary sensitivity.
No gravity or time evolution. Outside source sites remain fixed geometrically.
"""
import argparse,importlib.util,json,os,sys,time
from pathlib import Path
import numpy as np
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--failure-detail',action='store_true');parser.add_argument('--exact-pressure',action='store_true');parser.add_argument('--cutoff',type=float);args=parser.parse_args()
from scipy.spatial import cKDTree
base=Path('/gpfs/kjhan/LagForce');sys.path.insert(0,str(base/'Codex_coupled_relaxation_evolution_20261018'));sys.path.insert(0,str(base))
import moment_ale_evolve as model
from kinetic_factor_rhs import rhs
from run_canonical_av_branches_20261018 import read
def load(name,path):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m
geom=load('patch_geom',base/'Codex_weight_remap_moments_centered_20261017/audit.py')
clip=load('patch_clip',base/'Codex_weighted_compression_geometry_20261017/audit.py')
ic,t=read(base/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat')
allx=np.column_stack([ic['x']-3,ic['y']-3]);allw=.25*np.load(base/'Codex_evolving_weight_transport_20261017/weights_0.npz')['weights']
data=np.load(base/'Codex_checkpoint_covariance_factors_20261018/factors_w0.25.npz');ids=data['ids'];q=data['original_state'];mass=q[:,0];factors=list(data['factors'])
polys=geom.polygons(allx,allw,ids,clip.clip);lookup={int(v):i for i,v in enumerate(ids)};tree=cKDTree(allx);edges=[];max_match=0.;outside_faces=0
for i,(index,poly) in enumerate(zip(ids,polys)):
    _,neighbors=tree.query(allx[index],k=48);neighbors=neighbors[neighbors!=index]
    for a,b in zip(poly,np.roll(poly,-1,axis=0)):
        length=np.linalg.norm(b-a)
        if length<1e-14:continue
        mid=(a+b)/2;delta=allx[neighbors]-allx[index]
        err=abs(2*delta@(mid-allx[index])-np.sum(delta*delta,axis=1)-allw[index]+allw[neighbors])
        k=int(np.argmin(err));assert err[k]<1e-9
        max_match=max(max_match,float(err[k]));j=lookup.get(int(neighbors[k]),-1)
        if j<0:outside_faces+=1
        if j<0 or i<j:edges.append((i,j,a,b,np.array([b[1]-a[1],a[0]-b[0]])/length,length))
raw=np.column_stack([q[:,:3],q[:,3],q[:,4],q[:,4],q[:,5],q[:,6:12],q[:,12],q[:,13],q[:,13],q[:,14],q[:,15]])
y=np.r_[np.column_stack([allx[ids],allw[ids],raw[:,1:]]).ravel(),np.zeros(4)]
# Deliberately static geometry only for this instantaneous audit.
def saved_geometry(x,w,time):
    assert np.array_equal(x,allx[ids]) and np.array_equal(w,allw[ids])
    return polys,edges
model.geometry=saved_geometry;model.boundary_velocity=lambda p:np.zeros(2)
os.environ.update(LAGFORCE_QUADRATIC_BUDGET='projected',LAGFORCE_VELOCITY_FIT='radial_fixed_mean',LAGFORCE_VELOCITY_STENCIL='area',LAGFORCE_QUADRATIC_PRESSURE='0',LAGFORCE_KINETIC_SMOOTHING='0',LAGFORCE_PRESSURE_RIEMANN='hllc_davis',LAGFORCE_VECTOR_KINETIC_ROOT='0')
if args.exact_pressure:os.environ['LAGFORCE_PRESSURE_RIEMANN']='exact_euler_star'
name='Codex_checkpoint_patch_rhs'+('_exact' if args.exact_pressure else '')+('_failure' if args.failure_detail else '')+(f'_r{args.cutoff:g}' if args.cutoff is not None else '')+'_20261018'
root=base/name;root.mkdir(exist_ok=False);rows=[]
radius=np.linalg.norm(allx[ids],axis=1)
for quadratic in ((False,) if args.failure_detail else (False,True)):
    os.environ['LAGFORCE_QUADRATIC_VELOCITY']=str(int(quadratic))
    for cutoff in ((args.cutoff,) if args.cutoff is not None else ((.45,) if args.failure_detail else (.45,.5,.55))):
        frozen=radius>=cutoff;start=time.monotonic();row=dict(quadratic=quadratic,active_radius=cutoff,active_cells=int((~frozen).sum()),pressure_solver=os.environ['LAGFORCE_PRESSURE_RIEMANN'])
        try:
            d,diag=rhs(t,y,mass,covariance_factors=factors,frozen_mask=frozen)
            rates=d[:-4].reshape(-1,20);assert np.all(rates[frozen]==0)
            row.update(status='completed',mass_flux=diag['mass_flux'],minimum_pressure=diag['min_pressure'],generator_velocity_defect=diag['generator_velocity_defect'])
            np.savez(root/f'rhs_q{int(quadratic)}_r{cutoff:g}.npz',ids=ids,radius=radius,mass=mass,state=y,rates=d,pressure_rates=diag['pressure_rates'],frozen=frozen)
        except Exception as error:
            row.update(status='failed',error=repr(error))
            if args.failure_detail:
                frames=[];tb=error.__traceback__
                while tb:
                    frame=tb.tb_frame
                    if frame.f_code.co_name in ('rhs','fan','star'):
                        selected={}
                        for key in ('i','j','point','n','vl','vr','rl','rr','pl','pr','sl','sr','sm','ps','ql','qr','cl','cr','ul','ur','gamma'):
                            if key in frame.f_locals:
                                value=frame.f_locals[key]
                                selected[key]=value.tolist() if hasattr(value,'tolist') else value
                        frames.append(dict(function=frame.f_code.co_name,values=selected))
                    tb=tb.tb_next
                row['failure_frames']=frames
                np.savez(root/'source.npz',ids=ids,state=y,mass=mass,factors=factors,frozen=frozen)
        row['seconds']=time.monotonic()-start;rows.append(row);print(json.dumps(row),flush=True)
(root/'report.json').write_text(json.dumps(dict(scope=__doc__,source_time=t,cells=len(ids),shared_edges=len(edges),outside_faces=outside_faces,maximum_face_match_residual=max_match,results=rows,
    limits='Frozen exterior cells and fixed outside generators are artificial patch conditions. Compare common inner region before using this as a candidate evolution. No central gravity or physical stability conclusion.'),indent=2))
