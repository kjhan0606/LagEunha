"""Instantaneous hydro moment RHS on new continuous t0 projection.
Frozen outer patch, exact Euler pressure, no gravity or time evolution.
"""
import argparse,hashlib,json,os,sys,time
from pathlib import Path
import numpy as np
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--refine-faces',action='store_true')
parser.add_argument('--split-faces',action='store_true')
parser.add_argument('--quadratic-split',action='store_true')
parser.add_argument('--audit-covariance',action='store_true')
args=parser.parse_args()
if args.audit_covariance and not args.quadratic_split:parser.error('covariance audit requires quadratic split')
if args.quadratic_split:args.split_faces=True
if args.refine_faces and args.split_faces:parser.error('select one experiment')
base=Path('/gpfs/kjhan/LagForce');sys.path.insert(0,str(base));sys.path.insert(0,str(base/'Codex_coupled_relaxation_evolution_20261018'))
if args.split_faces:sys.path.insert(0,str(base/'Codex_split_face_modules_20261018'))
if args.quadratic_split:sys.path.insert(0,str(base/'Codex_quadratic_split_modules_20261018'))
import moment_ale_evolve as model
import fixedmass_hllc
from kinetic_factor_rhs import rhs
from run_canonical_av_branches_20261018 import read
source=base/'Codex_continuous_kepler_interface_split_20261018/moments_split_q16.npz'
data=np.load(source);ids=data['ids'];q=data['moments'];mass=q[:,0]
ic,t=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');assert t==0
allx=np.column_stack([ic['x']-3,ic['y']-3]);x=allx[ids];lookup={int(k):i for i,k in enumerate(ids)}
gamma=float(np.float32(1.666667))
# Explicit per-process EOS contract. Exact-star receives model.GAMMA directly;
# all kinetic energy and pressure-rate expressions use the same value.
model.GAMMA=gamma;fixedmass_hllc.GAMMA=gamma
os.environ.update(LAGFORCE_DENSITY_QUADRATURE_ORDER='12',LAGFORCE_QUADRATIC_VELOCITY='0',
    LAGFORCE_QUADRATIC_PRESSURE='0',LAGFORCE_KINETIC_SMOOTHING='0',
    LAGFORCE_PRESSURE_RIEMANN='exact_euler_star',LAGFORCE_VECTOR_KINETIC_ROOT='0')
if args.quadratic_split:
    os.environ.update(LAGFORCE_QUADRATIC_VELOCITY='1',LAGFORCE_QUADRATIC_BUDGET='projected',LAGFORCE_VELOCITY_FIT='radial_fixed_mean',LAGFORCE_VELOCITY_STENCIL='area')
faces=np.loadtxt(base/'Codex_pair_hg_modes_t0_20261018/raw_faces.dat.rank0000')
order=np.argsort(faces[:,0],kind='stable');owners=faces[order,0];polys=[];edges=[]
for local,i in enumerate(ids):
    ff=faces[order[np.searchsorted(owners,i):np.searchsorted(owners,i,side='right')]]
    v=ff[:,6:8];polys.append(v[np.argsort(np.arctan2(v[:,1],v[:,0]))]+x[local])
    for row in ff:
        other=int(row[1]);j=lookup.get(other,-1)
        if j>=0 and local>=j:continue
        a=x[local]+row[6:8];b=x[local]+row[8:10]
        delta=allx[other]-x[local];normal=delta/np.linalg.norm(delta)
        edges.append((local,j,a,b,normal,float(np.linalg.norm(b-a))))
raw=np.column_stack([q[:,:3],q[:,3],q[:,4],q[:,4],q[:,5],q[:,6:12],q[:,12],q[:,13],q[:,13],q[:,14],q[:,15]])
flat=np.r_[np.column_stack([x,np.zeros(len(ids)),raw[:,1:]]).ravel(),np.zeros(4)]
states=[model.decode(p,r) for p,r in zip(polys,raw)]
factors=[np.linalg.cholesky(s['R']) for s in states]
factor_error=max(float(np.max(abs(L@L.T-s['R']))) for L,s in zip(factors,states))
pressure_error=max(abs(s['P']/float(ic['pressure'][0])-1) for s in states)
def geometry(xx,ww,tt):
    assert np.array_equal(xx,x) and np.all(ww==0) and tt==0
    return polys,edges
model.geometry=geometry;model.boundary_velocity=lambda point:np.zeros(2)
root=base/('Codex_quadratic_covariance_tangent_audit_20261018' if args.audit_covariance else ('Codex_continuous_initial_rhs_quadratic_split_20261018' if args.quadratic_split else ('Codex_continuous_initial_rhs_split_20261018' if args.split_faces else ('Codex_continuous_initial_rhs_faces_20261018' if args.refine_faces else 'Codex_continuous_initial_rhs_20261018'))));root.mkdir(exist_ok=False)
rows=[];radius=np.linalg.norm(x,axis=1)
for cutoff in ((.55,) if args.refine_faces or args.split_faces else (.45,.5,.55)):
    for face_order in ((8,) if args.audit_covariance else ((16,32) if args.refine_faces else (4,8))):
        start=time.monotonic();frozen=radius>=cutoff
        row=dict(cutoff=cutoff,face_order=face_order,active_cells=int((~frozen).sum()))
        try:
            if args.split_faces:
                from affine_upwind_face_nodes import AffineUpwindNodes
                if args.quadratic_split:
                    from quadratic_upwind_face_nodes import QuadraticUpwindNodes as AffineUpwindNodes
                seed=np.load(base/'Codex_continuous_initial_rhs_faces_20261018/rhs_r0.55_f32.npz')
                assert np.array_equal(seed['state'],flat)
                z=seed['rates'][:-4].reshape(-1,20)[:,2].copy();history=[]
                for iteration in range(10):
                    provider=AffineUpwindNodes(x,mass,factors,z,face_order,model.boundary_velocity)
                    d,diag=rhs(0,flat,mass,covariance_factors=factors,frozen_mask=frozen,face_order=face_order,face_node_provider=provider)
                    newz=d[:-4].reshape(-1,20)[:,2];change=float(np.max(abs(newz-z)))
                    history.append(dict(iteration=iteration,weight_rate_change=change,split_faces=provider.split_faces,nodes=provider.nodes,mass_flux=diag['mass_flux']))
                    print(json.dumps(dict(face_order=face_order,split_iteration=history[-1])),flush=True)
                    if change<1e-11:break
                    z=newz.copy()
                else:raise RuntimeError('split quadrature/weight fixed point did not converge')
                row['split_history']=history
            else:
                d,diag=rhs(0,flat,mass,covariance_factors=factors,frozen_mask=frozen,face_order=face_order)
            rates=d[:-4].reshape(-1,20);assert np.all(rates[frozen]==0)
            if args.audit_covariance:
                direct=[];factor_rate=[]
                for i,(r,L) in enumerate(zip(raw,factors)):
                    m=mass[i];q0=r[1:3];p0=r[7:9]
                    v=np.r_[0.,rates[i,3:]];dq=v[1:3];dp=v[7:9]
                    S=r[3:7].reshape(2,2)-np.outer(q0,q0)/m
                    C=r[9:13].reshape(2,2)-np.outer(q0,p0)/m
                    B=np.linalg.solve(S,C).T
                    dS=v[3:7].reshape(2,2)-(np.outer(dq,q0)+np.outer(q0,dq))/m
                    dC=v[9:13].reshape(2,2)-(np.outer(dq,p0)+np.outer(q0,dp))/m
                    dT=v[13:17].reshape(2,2)-(np.outer(dp,p0)+np.outer(p0,dp))/m
                    direct.append(dT-B@dC-dC.T@B.T+B@dS@B.T)
                    A=diag['cov_A'][i];Qf=diag['cov_Qfactor'][i];R=L@L.T
                    factor_rate.append(A@R+R@A.T-L@diag['cov_K'][i]@L.T+Qf@Qf.T)
                direct=np.asarray(direct);factor_rate=np.asarray(factor_rate);err=factor_rate-direct
                active=~frozen;norm=float(np.linalg.norm(direct[active]))
                row['covariance_tangent_audit']=dict(maximum_absolute_difference=float(abs(err[active]).max()),
                    relative_L2_difference=float(np.linalg.norm(err[active])/norm),
                    direct_L2=norm,factor_L2=float(np.linalg.norm(factor_rate[active])),
                    maximum_implied_thermal_rate_change=float(.5*abs(np.trace(err[active],axis1=1,axis2=2)).max()))
                np.savez(root/'covariance_tangents.npz',ids=ids,frozen=frozen,direct=direct,factor=factor_rate)
            # Flat layout removes raw m; p columns9:11, D columns11:15,
            # Ehyd column19. Boundary ledger stores p2,E,L.
            closure=np.r_[rates[:,9:11].sum(axis=0),rates[:,19].sum(),(rates[:,12]-rates[:,13]).sum()]-d[-4:]
            row.update(status='completed',mass_flux=diag['mass_flux'],minimum_pressure=diag['min_pressure'],
                generator_velocity_defect=diag['generator_velocity_defect'],conservation_ledger_residual=closure.tolist(),
                maximum_active_pressure_rate=float(np.max(abs(diag['pressure_rates'][~frozen]))),
                weight_rate_range=float(np.ptp(rates[~frozen,2])))
            np.savez(root/f'rhs_r{cutoff:g}_f{face_order}.npz',ids=ids,state=flat,mass=mass,rates=d,
                pressure_rates=diag['pressure_rates'],frozen=frozen,covariance_factors=factors)
        except Exception as exc:row.update(status='failed',error=repr(exc))
        row['seconds']=time.monotonic()-start;rows.append(row);print(json.dumps(row),flush=True)
report=dict(scope=__doc__,source=str(source),sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
    gamma=gamma,maximum_relative_pressure_decode_error=pressure_error,maximum_factor_recomposition_error=factor_error,
    quadratic_velocity=args.quadratic_split,
    results=rows,limits='Instantaneous fixed geometry only. Frozen reservoir affects adjacent active cells. No gravity, time evolution, production restart, GPU or MPI test. Cholesky uses positive raw covariance without clipping.')
(root/'report.json').write_text(json.dumps(report,indent=2))
