"""Controlled tangent-response probes of the original canonical executable.

Run on lageunha. Uses new directories; never changes source checkpoints or
production runs. Velocity perturbations leave thermal state and mass untouched.
This finite-time directional response is not an eigenvalue/stability proof.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import struct
import subprocess
import time

import numpy as np

BASE = Path('/gpfs/kjhan/LagForce')
SOURCE = BASE/'Codex_HG_constrained_20261017/omp8_t20/khout.000184.dat'
TEMPLATE = BASE/'Codex_entropy_step_audit_20261017/t2_half'
BINARY = BASE/'Codex_integrated_cpu_cuda_20261017/eunha2.exe'
EXPECTED_BINARY = 'ebb4c7dcbeee1a3b1e3da759b5ec2151f7f342775f2c4d0d78f65edfe4dc7697'
LAYOUT_PATH = BASE/'Codex_DPQ_backend_parity_20261017/layout.json'
FIELDS = ['x','y','mass','vx','vy','den','pressure','ie','volume','w2','w2ceil','csound']


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def read(p, layout):
    raw = bytearray(p.read_bytes())
    magic,n,size,stress = struct.unpack_from('<iiii',raw)
    assert magic == -1 and size == layout['size'] and stress == 1
    assert len(raw) == 16+n*size+52
    dtype = np.dtype({'names':['id']+FIELDS, 'formats':['<u8']+['<f8']*len(FIELDS),
                     'offsets':[layout['u4if']]+[layout[k] for k in FIELDS], 'itemsize':size})
    a = np.frombuffer(raw,dtype=dtype,count=n,offset=16)
    ids = a['id'] & layout['id_mask']
    assert np.array_equal(np.sort(ids),np.arange(n))
    t = struct.unpack_from('<d',raw,16+n*size)[0]
    return raw,a,np.argsort(ids),t


def prepare(root):
    assert sha(BINARY) == EXPECTED_BINARY
    layout = json.loads(LAYOUT_PATH.read_text())
    raw,a,order,t0 = read(SOURCE,layout)
    x,y = a['x']-3,a['y']-3
    r = np.hypot(x,y)
    theta = np.arctan2(y,x)
    sel = (r>.30)&(r<.50)&(a['den']<.1)
    window = np.zeros(len(a))
    window[sel] = np.sin(np.pi*(r[sel]-.30)/.20)**2
    tangential = np.column_stack((-y/r,x/r))
    seed = (window*a['csound']*np.cos(30*theta))[:,None]*tangential
    # Remove total linear/angular impulse using smooth compact support.
    basis = np.zeros((len(a),2,3))
    basis[:,0,0] = window
    basis[:,1,1] = window
    basis[:,0,2] = -window*y
    basis[:,1,2] = window*x
    moment = np.zeros_like(basis)
    moment[:,0,0] = 1
    moment[:,1,1] = 1
    moment[:,0,2] = -y
    moment[:,1,2] = x
    mass = a['mass']
    matrix = np.einsum('n,nij,nik->jk',mass,moment,basis)
    rhs = np.einsum('n,nij,ni->j',mass,moment,seed)
    seed -= np.einsum('nij,j->ni',basis,np.linalg.solve(matrix,rhs))
    assert np.isfinite(seed).all() and np.all(seed[~sel] == 0)
    root.mkdir(exist_ok=False)
    (root/'layout.json').write_text(json.dumps(layout,indent=2)+'\n')
    np.savez(root/'seed.npz',velocity=seed[order],selected=sel[order],theta=theta[order],radius=r[order])
    cases = []
    for dt, amplitudes in ((.001,(0.,.0005,-.0005,.00025,-.00025)),
                           (.0005,(0.,.00025,-.00025))):
        for eps in amplitudes:
            name = f'dt{dt:g}_eps{eps:+g}'
            run = root/name
            run.mkdir()
            payload,b,_,_ = read(SOURCE,layout)
            b['vx'] += eps*seed[:,0]
            b['vy'] += eps*seed[:,1]
            allowed = np.zeros(len(payload),dtype=bool)
            for off in (layout['vx'],layout['vy']):
                for byte in range(8): allowed[16+off+byte:16+len(a)*layout['size']:layout['size']] = True
            changed = np.frombuffer(payload,dtype='u1') != np.frombuffer(raw,dtype='u1')
            assert not np.any(changed & ~allowed)
            (run/SOURCE.name).write_bytes(payload)
            (run/'params.dat').write_bytes((TEMPLATE/'params.dat').read_bytes())
            env = json.loads((TEMPLATE/'environment.json').read_text())
            env.pop('GFS_TCF_SACC',None)
            # Stop between penultimate and last step; assert actual horizon later.
            env.update(HYDRO_TSTOP=str(t0+.01-.25*dt),GFS_MAX_DT=str(dt),EUNHA_DUMP_DT='1000000',
                       OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
            (run/'environment.json').write_text(json.dumps(env,indent=2)+'\n')
            v0 = np.column_stack((a['vx'],a['vy']))
            dv = np.column_stack((b['vx']-a['vx'],b['vy']-a['vy']))
            cases.append({'name':name,'dt':dt,'eps':eps,'steps':int(round(.01/dt)),
                          'start_sha256':sha(run/SOURCE.name),
                          'initial_impulse':np.einsum('n,nij,ni->j',mass,moment,dv).tolist(),
                          'initial_kinetic_change':float(np.sum(mass*(np.sum(v0*dv,axis=1)+.5*np.sum(dv*dv,axis=1)))),
                          'modified_only_velocity':True})
    manifest = {'source':str(SOURCE),'source_sha256':sha(SOURCE),'binary':str(BINARY),
                'binary_sha256':sha(BINARY),'t0':t0,'horizon':.01,'restart_step':184,
                'seed_mode':30,'seed_region':'.30 < current r < .50, current rho < .1',
                'perturbation':'compact tangential velocity mode scaled by local sound speed; projected to zero total Px, Py and Lz impulse; internal energy unchanged',
                'cases':cases,'production_runs_modified':False,
                'scope':'Finite-time directional derivative around already evolved state. Does not isolate a complete spectrum or establish long-time stability.'}
    (root/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps({'prepared':str(root),'cases':len(cases),'t0':t0}),flush=True)


def execute(root):
    manifest = json.loads((root/'manifest.json').read_text())
    assert sha(BINARY) == manifest['binary_sha256'] == EXPECTED_BINARY
    assert sha(SOURCE) == manifest['source_sha256']
    for case in manifest['cases']:
        run = root/case['name']
        env = {k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))}
        env.update(json.loads((run/'environment.json').read_text()))
        start = time.monotonic()
        with (run/'log').open('x') as log:
            proc = subprocess.Popen(['mpirun','-np','1',str(BINARY),'params.dat','184'],cwd=run,env=env,stdout=log,stderr=subprocess.STDOUT)
            (run/'runner_pid').write_text(str(proc.pid)+'\n')
            code = proc.wait()
        (run/'exit.code').write_text(str(code)+'\n')
        print(json.dumps({'case':case['name'],'exit':code,'wall_seconds':time.monotonic()-start}),flush=True)
        if code: raise RuntimeError(f'{case["name"]} exited {code}')


def analyze(root, partial=False):
    manifest = json.loads((root/'manifest.json').read_text())
    layout = json.loads((root/'layout.json').read_text())
    seed = np.load(root/'seed.npz')['velocity']
    _,base,order,_ = read(SOURCE,layout)
    base = base[order]
    mass = base['mass']
    rows = {}
    completed = []
    for case in manifest['cases']:
        run = root/case['name']
        if partial and not (run/'exit.code').exists():
            continue
        assert (run/'exit.code').read_text().strip() == '0'
        end = run/f'khout.{184+case["steps"]:06d}.dat'
        _,a,order,t = read(end,layout)
        a = a[order]
        assert abs(t-manifest['t0']-.01)<1e-12,(case['name'],t)
        assert np.array_equal(a['mass'],base['mass'])
        assert all(np.isfinite(a[f]).all() for f in FIELDS)
        assert np.all(a['pressure']>0)&np.all(a['volume']>0)&np.all(a['ie']>0)
        lines = (run/'log').read_text().splitlines()
        accepted = [dict(re.findall(r'(\w+)=([^\s]+)',line)) for line in lines if line.startswith('[ACCEPTED_STATE]')]
        assert len(accepted) == case['steps']
        assert all(q['finite']=='1' and q['retries']=='0' and abs(float(q['dt'])-case['dt'])<1e-14 for q in accepted)
        rk = [dict(re.findall(r'(\w+)=([^\s]+)',line)) for line in lines if line.startswith('[ACCEPTED_RK_BUDGET]')]
        # MPI stdout/stderr may interleave a [CPU ...] message inside one
        # RK4E record. Require a complete numeric token, allowing the next
        # bracketed logger tag as a separator; keep continuation lines.
        numeric = r'(\w+)=([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)(?=\s|\[|$)'
        floors = [dict(re.findall(numeric,block)) for block in '\n'.join(lines).split('[RK4E]')[1:]]
        assert len(rk)==len(floors)==case['steps']
        assert all(q['stage_mask_min']==q['stage_mask_max']=='15' for q in rk)
        for q in floors:
            for key in ('n_ie_le0','floor_inj','floor_cum','n_de','de_inj','de_cum','n_sfl','sfl_inj','sfl_cum'):
                assert float(q[key])==0,(case['name'],key,q[key])
        rows[case['dt'],case['eps']] = a
        case.update(end_sha256=sha(end),end_time=t,actual_steps=len(accepted),retries=0,
                    mass_bitwise_unchanged=True,positive_state=True,logged_energy_floor_injections_zero=True,
                    all_RK_stages_recorded=True,
                    maximum_accepted_energy_ledger_residual=float(max(abs(float(q['residual_E'])) for q in rk)),
                    final_total_energy_change=float(accepted[-1]['dE0']))
        completed.append(case['name'])
    r = np.hypot(base['x']-3,base['y']-3)
    masks = {'inner_background':(r<.5)&(base['den']<.1),
             'seed_region':np.load(root/'seed.npz')['selected'],
             'disk_inner_edge':(r>=.5)&(r<.6)&(base['den']>.1),
             'all_active':r<2.5}
    response = {}
    summaries = []
    for dt,eps in ((.001,.0005),(.001,.00025),(.0005,.00025)):
        if any(key not in rows for key in ((dt,eps),(dt,-eps),(dt,0.))):
            assert partial
            continue
        plus,minus,zero = rows[dt,eps],rows[dt,-eps],rows[dt,0.]
        derivative = {f:(plus[f]-minus[f])/(2*eps) for f in ('x','y','vx','vy','pressure','ie','volume')}
        response[dt,eps] = derivative
        dvel = np.column_stack((derivative['vx'],derivative['vy']))
        even = np.column_stack(((plus['vx']+minus['vx'])/2-zero['vx'],(plus['vy']+minus['vy'])/2-zero['vy']))
        record = {'dt':dt,'eps':eps,'regions':{}}
        # Fixed material angular labels, with local basis rotating with the
        # unperturbed trajectory. Do not select the largest mode after evolving.
        theta0 = np.arctan2(base['y']-3,base['x']-3)
        thetat = np.arctan2(zero['y']-3,zero['x']-3)
        tangent0 = np.column_stack((-np.sin(theta0),np.cos(theta0)))
        tangentt = np.column_stack((-np.sin(thetat),np.cos(thetat)))
        mask = masks['seed_region']
        phase = np.exp(-30j*theta0[mask])
        amp0 = np.sum(mass[mask]*np.sum(seed[mask]*tangent0[mask],axis=1)*phase)
        amp1 = np.sum(mass[mask]*np.sum(dvel[mask]*tangentt[mask],axis=1)*phase)
        gain = amp1/amp0
        record['fixed_material_mode30_corotating_basis'] = {'amplitude_gain':float(abs(gain)),
            'phase_change_rad':float(np.angle(gain)),
            'definition':'Tangential projection in the unperturbed endpoint local basis, Fourier phase from fixed t2 particle angular labels.'}
        for name,s in masks.items():
            w = mass[s]
            norm = lambda v:float(np.sqrt(np.sum(w[:,None]*v[s]**2)))
            norm0 = norm(seed)
            E0 = .5*np.sum(w*np.sum(seed[s]**2,axis=1))
            # A local acoustic norm, not a conserved perturbation energy for the full shearing flow.
            E1 = .5*np.sum(w*(np.sum(dvel[s]**2,axis=1)+derivative['pressure'][s]**2/(base['den'][s]*base['csound'][s])**2))
            record['regions'][name] = {'cells':int(s.sum()),'initial_velocity_norm':norm0,
                'velocity_response_norm':norm(dvel),'velocity_norm_gain':norm(dvel)/norm0 if norm0>0 else None,
                'velocity_projection_gain':float(np.sum(w[:,None]*seed[s]*dvel[s])/(2*E0)) if E0>0 else None,
                'acoustic_norm_energy_gain':float(E1/E0) if E0>0 else None,
                'pressure_response_rms_over_initial_pressure':float(np.sqrt(np.average((derivative['pressure'][s]/base['pressure'][s])**2,weights=w))),
                'unperturbed_pressure_change_rms_over_start_pressure':float(np.sqrt(np.average(((zero['pressure'][s]-base['pressure'][s])/base['pressure'][s])**2,weights=w))),
                'even_velocity_residual_norm':norm(even),
                'even_over_odd_velocity':norm(even)/(eps*norm(dvel)) if norm(dvel)>0 else None}
        summaries.append(record)
    comparisons = []
    for label,left,right in [('amplitude',(.001,.0005),(.001,.00025)),('time_step',(.001,.00025),(.0005,.00025))]:
        if left not in response or right not in response:
            assert partial
            continue
        metrics = {}
        for f in response[left]:
            a,b = response[left][f],response[right][f]
            s = masks['all_active']
            den = np.sum(mass[s]*b[s]**2)
            metrics[f] = float(np.sqrt(np.sum(mass[s]*(a[s]-b[s])**2)/den)) if den else None
        comparisons.append({'comparison':label,'relative_response_L2_difference':metrics})
    result = {'manifest':manifest,'completed_cases':completed,'all_cases_complete':len(completed)==len(manifest['cases']),
              'analysis_source_sha256':sha(Path(__file__)),
              'responses':summaries,'convergence_comparisons':comparisons,
              'scope':'One prescribed velocity direction over dt_total=.01. Acoustic norm uses start density/sound speed and is not an invariant. Positive gain does not prove an unstable eigenmode; negative gain does not exclude other growing modes. Topology/dissipation switches can make response nonsmooth.',
              'long_time_stability_verified':False}
    (root/('partial_result.json' if partial else 'result.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'responses':summaries,'comparisons':comparisons},indent=2),flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action',choices=['prepare','execute','analyze'])
    parser.add_argument('root',type=Path)
    parser.add_argument('--partial',action='store_true',help='Analyze only completed cases; never call missing data a pass.')
    args = parser.parse_args()
    if args.action == 'analyze': analyze(args.root,partial=args.partial)
    else:
        assert not args.partial
        globals()[args.action](args.root)
