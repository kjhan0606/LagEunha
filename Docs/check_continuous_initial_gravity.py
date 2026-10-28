"""Gravity moment-source quadrature on the new IC, plus saved hydro rates.
No time integration; potential energy transport is not implemented here.
"""
import argparse,hashlib,importlib.util,json,os,sys
from pathlib import Path
import numpy as np
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--quadratic',action='store_true');parser.add_argument('--export-backend',action='store_true');args=parser.parse_args()
base=Path('/gpfs/kjhan/LagForce');sys.path.insert(0,str(base));sys.path.insert(0,str(base/'Codex_coupled_relaxation_evolution_20261018'));sys.path.insert(0,str(base/'Codex_split_face_modules_20261018'))
import moment_ale_evolve as model
from moment_gravity_source import source
from audit_canonical_pair_history import gravity_coefficient
from run_canonical_av_branches_20261018 import read
spec=importlib.util.spec_from_file_location('geom',base/'Codex_weight_remap_moments_centered_20261017/audit.py');geom=importlib.util.module_from_spec(spec);spec.loader.exec_module(geom)
path=base/('Codex_continuous_initial_rhs_quadratic_split_20261018' if args.quadratic else 'Codex_continuous_initial_rhs_split_20261018')/'rhs_r0.55_f8.npz';saved=np.load(path)
ids=saved['ids'];mass=saved['mass'];flat=saved['state'];y=flat[:-4].reshape(-1,20);raw=np.column_stack([mass,y[:,3:]])
ic,_=read(base/'Codex_canonical_initial_20261017/run/khout.000000.dat');x=np.column_stack([ic['x']-3,ic['y']-3])
faces=np.loadtxt(base/'Codex_pair_hg_modes_t0_20261018/raw_faces.dat.rank0000');order=np.argsort(faces[:,0],kind='stable');owners=faces[order,0]
os.environ['LAGFORCE_DENSITY_QUADRATURE_ORDER']='12';model.GAMMA=float(np.float32(1.666667))
polys=[];states=[]
for i,r in zip(ids,raw):
    ff=faces[order[np.searchsorted(owners,i):np.searchsorted(owners,i,side='right')]]
    vv=ff[:,6:8];poly=vv[np.argsort(np.arctan2(vv[:,1],vv[:,0]))]+x[i]
    polys.append(poly);states.append(model.decode(poly,r))
if args.quadratic:
    import quadratic_velocity_reconstruction as quad
    os.environ.update(LAGFORCE_QUADRATIC_BUDGET='projected',LAGFORCE_VELOCITY_FIT='radial_fixed_mean',LAGFORCE_VELOCITY_STENCIL='area')
    reconstruction=quad.prepare(states,polys,[],mass,saved['covariance_factors'])
root=base/('Codex_gravity_quadrature_backends_20261018' if args.export_backend else ('Codex_continuous_quadratic_gravity_20261018' if args.quadratic else 'Codex_continuous_initial_gravity_v2_20261018'));root.mkdir(exist_ok=False)
rows=[];previous=None
for nq in ((12,) if args.export_backend else (12,24)):
    dg=[];centroid=[];Rrates=[];thermal=[];exactworks=[];Ridentity=[]
    samples=[]
    for s,p,r in zip(states,polys,raw):
        pts,dv=geom.quadrature(p,nq);dm=dv*s['density'].value(pts)
        rr=np.linalg.norm(pts,axis=1);g=gravity_coefficient(rr);acc=-g[:,None]*pts
        vel=s['v']+(pts-s['c'])@s['B'].T
        correction=quad.correction(s,pts) if args.quadratic else np.zeros_like(vel)
        vel+=correction
        if args.export_backend:samples.append(np.column_stack([pts,dm,vel,acc]))
        d=source(pts,dm,vel,acc);dg.append(d)
        thermal.append(d[17]-.5*(d[13]+d[16]))
        vc=np.sqrt(g)[:,None]*np.column_stack([-pts[:,1],pts[:,0]])
        exactworks.append(float(np.sum(dm*np.sum(vc*acc,axis=1))))
        c=s['c'];cg=-gravity_coefficient(np.array([np.linalg.norm(c)]))[0]*c
        centroid.append(r[0]*cg)
        # Quadratic mean can exchange energy with the stored residual budget.
        # Compare the differentiated moment identity with a direct integral.
        Cdot=d[9:13].reshape(2,2)-np.outer(c,d[7:9])
        Tdot=d[13:17].reshape(2,2)-np.outer(s['v'],d[7:9])-np.outer(d[7:9],s['v'])
        Rdot=Tdot-s['B']@Cdot-Cdot.T@s['B'].T;Rrates.append(Rdot)
        direct=correction.T@(dm[:,None]*acc)
        Ridentity.append(Rdot-direct-direct.T)
    dg=np.array(dg);active=~saved['frozen'];torque=dg[:,10]-dg[:,11]
    combined=saved['rates'].copy();rates=combined[:-4].reshape(-1,20);rates[active,3:]+=dg[active,1:]
    assert np.all(rates[~active]==0)
    np.savez(root/f'gravity_q{nq}.npz',ids=ids,source=dg,combined_rates=combined,centroid_force=centroid,residual_covariance_source=Rrates)
    norm=np.linalg.norm(dg[active,7:9]);diff=(np.array(centroid)-dg[:,7:9])[active]
    row=dict(order=nq,maximum_central_torque=float(np.max(abs(torque))),
        maximum_thermal_source=float(np.max(np.abs(thermal))),maximum_residual_covariance_source=float(np.max(abs(np.array(Rrates)))),
        maximum_analytic_circular_work=float(np.max(np.abs(exactworks))),maximum_reconstructed_gravity_work=float(np.max(abs(dg[:,17]))),
        maximum_residual_source_identity_error=float(np.max(abs(np.array(Ridentity)))),
        centroid_force_relative_l2=float(np.linalg.norm(diff)/norm),centroid_force_max_abs=float(np.max(abs(diff))))
    if previous is not None:row['successive_source_max_abs']=float(np.max(abs(dg-previous)));row['successive_active_force_relative_l2']=float(np.linalg.norm((dg-previous)[active,7:9])/norm)
    previous=dg;rows.append(row)
    if args.export_backend:
        # Explicit invalid fixtures must preserve the output sentinel.
        bad=samples[0][:1].copy();bad[0,2]=-1;samples.append(bad)
        bad=samples[0][:1].copy();bad[0,5]=np.nan;samples.append(bad)
        offsets=np.r_[0,np.cumsum([len(a) for a in samples])].astype('<u8')
        with (root/'samples.bin').open('xb') as f:np.concatenate(samples).astype('<f8').tofile(f)
        with (root/'offsets.bin').open('xb') as f:offsets.tofile(f)
out=dict(scope=__doc__,quadratic_velocity=args.quadratic,hydro_source=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),results=rows,
    limits='Gravity integrated over fitted smooth density and the selected conditional mean. Does not integrate original discontinuous density exactly. Combined rates exclude potential energy transport and moving-boundary validation. No full total-energy conservation or long-term stability claim.')
(root/'report.json').write_text(json.dumps(out,indent=2));print(json.dumps(out,indent=2))
