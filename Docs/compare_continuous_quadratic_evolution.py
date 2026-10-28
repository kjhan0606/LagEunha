"""Compare completed equal-time moving moment runs without inferring time order."""
import argparse,json,sys
from pathlib import Path
import numpy as np

base=Path('/gpfs/kjhan/LagForce')
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--three-level',action='store_true')
parser.add_argument('--active-radius',type=float,choices=(.55,.60,.65),default=.55)
parser.add_argument('--extended',action='store_true',help='Compare the two common-checkpoint continuations to t=1e-3')
parser.add_argument('--at-time',type=float,help='Compare accepted intermediate checkpoints, without requiring run completion')
args=parser.parse_args()
if args.at_time is not None and (not args.extended or not np.isfinite(args.at_time)):
    parser.error('--at-time requires --extended and a finite time')
if args.extended and (args.three_level or args.active_radius!=.60):
    parser.error('--extended requires --active-radius 0.60 and exactly two levels')
roots=[base/'Codex_continuous_quadratic_evolution_dt1e4_20261018',
       base/'Codex_continuous_quadratic_evolution_dt5e5_20261018']
if args.three_level:roots.append(base/'Codex_continuous_quadratic_evolution_dt25e6_20261018')
if args.active_radius!=.55:
    label=f'r{round(args.active_radius*100):03d}'
    roots=[base/f'Codex_quadratic_{label}_{tag}_20261018' for tag in
           (('dt1e4','dt5e5','dt25e6') if args.three_level else ('dt1e4','dt5e5'))]
if args.extended:
    roots=[base/f'Codex_quadratic_extended_{tag}_20261018' for tag in ('dt1e4','dt5e5')]
manifests=[json.loads((root/'manifest.json').read_text()) for root in roots]
assert all(m.get('active_radius',.55)==args.active_radius for m in manifests)
if args.at_time is not None:
    for m in manifests:
        span=args.at_time-m.get('start_time',0.)
        if span<=0 or args.at_time>m['end'] or abs(round(span/m['dt'])*m['dt']-span)>1e-14:
            parser.error('requested time must be a later common step boundary within both runs')
status=[]
for root in roots:
    if args.at_time is not None:
        events=[]
        with (root/'events.jsonl').open() as stream:
            for line in stream:
                try: row=json.loads(line)
                except json.JSONDecodeError: break  # Possibly an incomplete live write.
                if 'initial' in row or ('accepted' in row and row['accepted']['time']<=args.at_time+1e-14):events.append(row)
        accepted=[r['accepted'] for r in events if 'accepted' in r]
        ready=bool(accepted) and abs(accepted[-1]['time']-args.at_time)<1e-14
        status.append(dict(status='completed' if ready else 'pending',events=events,
                           completed_time=accepted[-1]['time'] if accepted else None))
        continue
    p=root/'summary.json'
    status.append(json.loads(p.read_text()) if p.exists() else {'status':'pending'})
if any(s['status']!='completed' for s in status):
    print(json.dumps(dict(scope=__doc__,ready=False,runs=[dict(root=str(r),status=s['status'],completed_time=s.get('completed_time')) for r,s in zip(roots,status)]),indent=2))
    sys.exit(2)
assert all(m['sha256']==manifests[0]['sha256'] and m['end']==manifests[0]['end'] for m in manifests)
assert all(m.get('start_time',0.)==manifests[0].get('start_time',0.) for m in manifests)
assert all(abs(manifests[i]['dt']/manifests[i+1]['dt']-2)<1e-14 for i in range(len(manifests)-1))
initial=[np.load(root/'initial.npz') for root in roots]
assert all(np.array_equal(initial[0][k],other[k]) for other in initial[1:] for k in ('state','mass','ids','factors'))
mass=initial[0]['mass'];rows=[];fields=[]
source=np.load(manifests[0]['input'])
assert np.array_equal(source['ids'],initial[0]['ids']) and np.array_equal(source['mass'],mass)
# The runner freezes by the original IC radius, not by the moved restart state.
initial_positions=source['state'][:-4].reshape(-1,20)[:,:2]
active=np.linalg.norm(initial_positions,axis=1)<args.active_radius
masks={'active':active,'common_inner':active&(np.linalg.norm(initial_positions,axis=1)<.4),
       'frozen_reservoir':~active}
if args.active_radius!=.55:masks['original_active']=~source['frozen']
for root,manifest,summary,ic in zip(roots,manifests,status,initial):
    target=manifest['end'] if args.at_time is None else args.at_time
    count=round((target-manifest.get('start_time',0.))/manifest['dt'])
    saved=np.load(root/f'step{count:04d}.npz')
    assert float(saved['time'])==summary['completed_time']
    assert abs(float(saved['time'])-target)<1e-14
    assert np.array_equal(saved['mass'],mass) and np.array_equal(saved['ids'],ic['ids'])
    a=saved['state'][:-4].reshape(-1,20);L=saved['factors']
    assert np.isfinite(a).all() and np.isfinite(L).all()
    assert np.array_equal(a[~active],ic['state'][:-4].reshape(-1,20)[~active])
    m=mass[:,None];q=a[:,3:5];p=a[:,9:11]
    S=a[:,5:9].reshape(-1,2,2)-q[:,:,None]*q[:,None,:]/m[:,:,None]
    C=a[:,11:15].reshape(-1,2,2)-q[:,:,None]*p[:,None,:]/m[:,:,None]
    Pi=a[:,15:19].reshape(-1,2,2)
    B=np.linalg.solve(S,C).transpose(0,2,1)
    R=L@L.transpose(0,2,1)
    reconstructed=p[:,:,None]*p[:,None,:]/m[:,:,None]+B@S@B.transpose(0,2,1)+R
    thermal=a[:,19]-.5*np.trace(Pi,axis1=1,axis2=2)
    assert np.all(thermal>0) and np.all(np.linalg.eigvalsh(S)>0) and np.all(np.linalg.eigvalsh(R)>0)
    ell=C[:,0,1]-C[:,1,0]
    fields.append(dict(generator=a[:,:2],weight=a[:,2],centroid=q/m,velocity=p/m,
        angular_momentum=ell,specific_thermal=thermal/mass,residual_covariance=R))
    accepted=[r['accepted'] for r in summary['events'] if 'accepted' in r]
    assert len(accepted)==count
    initial_diag=summary['events'][0]['initial'];final=accepted[-1]
    rows.append(dict(dt=manifest['dt'],steps=count,final=final,
        maximum_Pi_factor_mismatch=float(abs(Pi-reconstructed).max()),
        minimum_thermal=float(thermal.min()),
        reconstructed_energy_change=(final['active_hydro_energy']+final['active_reconstructed_potential'])-
            (initial_diag['active_hydro_energy']+initial_diag['active_reconstructed_potential']),
        energy_note='Open active patch: this energy change includes reservoir exchange; it is not a closed-system conservation error.'))
comparison={}
for region,mask in masks.items():
    comparison[region]={'cells':int(mask.sum()),'fields':{}}
    assert np.any(mask)
    for key in fields[0]:
        difference=(fields[1][key]-fields[0][key])[mask]
        squared=difference.reshape(int(mask.sum()),-1)**2
        comparison[region]['fields'][key]=dict(maximum_absolute_difference=float(abs(difference).max()),
            mass_weighted_RMS=float(np.sqrt(np.average(squared.sum(axis=1),weights=mass[mask]))))
        if args.three_level:
            fine=(fields[2][key]-fields[1][key])[mask]
            fine_rms=float(np.sqrt(np.average((fine.reshape(int(mask.sum()),-1)**2).sum(axis=1),weights=mass[mask])))
            coarse_rms=comparison[region]['fields'][key]['mass_weighted_RMS']
            # Report empirical ratios only. Topology transitions, spatial fit
            # errors and roundoff may prevent an asymptotic temporal regime.
            order=float(np.log2(coarse_rms/fine_rms)) if coarse_rms>0 and fine_rms>0 else None
            comparison[region]['fields'][key].update(fine_pair_mass_weighted_RMS=fine_rms,
                fine_pair_maximum_absolute_difference=float(abs(fine).max()),empirical_order=order)
print(json.dumps(dict(scope=__doc__,ready=True,runs=rows,same_time_differences=comparison,
    requested_checkpoint_time=args.at_time,full_run_completion_required=args.at_time is None,
    limits='Three-level ratios, when available, are empirical and do not prove asymptotic convergence. Two levels do not establish order. Very short stage test, not historical Kepler validation. Potential transport closure, frozen-boundary effects and production CUDA/MPI remain unresolved.'),indent=2))
