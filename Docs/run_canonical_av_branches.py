"""Controlled t6 restart: original QBP compression blend versus chi disabled.

Runs on lageunha, one MPI rank and one OpenMP thread per independent case.
Original RK4/CUDA executable and all other physics settings are retained.
This is a local response diagnostic, not a long-time stability gate.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import re
import struct
import subprocess
import time

import numpy as np

BASE=Path('/gpfs/kjhan/LagForce')
SOURCE=BASE/'Codex_HG_constrained_20261017/omp8_t20/khout.000656.dat'
TEMPLATE=BASE/'Codex_early_growth_budget_20261017/t6_full'
BINARY=BASE/'Codex_integrated_cpu_cuda_20261017/eunha2.exe'
EXPECTED='ebb4c7dcbeee1a3b1e3da759b5ec2151f7f342775f2c4d0d78f65edfe4dc7697'
LAYOUT=BASE/'Codex_DPQ_backend_parity_20261017/layout.json'
FIELDS=['x','y','vx','vy','mass','volume','den','pressure','ie','w2']


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def read(p):
    layout=json.loads(LAYOUT.read_text());raw=p.read_bytes()
    magic,n,size,stress=struct.unpack_from('<iiii',raw)
    assert magic==-1 and stress==1 and size==layout['size'] and len(raw)==16+n*size+52
    dtype=np.dtype({'names':['id']+FIELDS,'formats':['<u8']+['<f8']*len(FIELDS),
                    'offsets':[layout['u4if']]+[layout[k] for k in FIELDS], 'itemsize':size})
    a=np.frombuffer(raw,dtype=dtype,count=n,offset=16)
    ids=a['id']&layout['id_mask'];order=np.argsort(ids)
    assert np.array_equal(ids[order],np.arange(n))
    return a[order],struct.unpack_from('<d',raw,16+n*size)[0]


def prepare(root):
    assert sha(BINARY)==EXPECTED
    a,t=read(SOURCE);assert abs(t-6.002110493120864)<1e-14
    root.mkdir(exist_ok=False)
    manifest={'source':str(SOURCE),'source_sha256':sha(SOURCE),'binary':str(BINARY),'binary_sha256':EXPECTED,
              'start_time':t,'horizon':.2,'dt_cap':.001,'dump_dt':.01,'cases':[],
              'scope':'Same evolved t6 state; local AV causal-response diagnostic, not an onset or long-time pass.',
              'storage_estimate_bytes':44*SOURCE.stat().st_size,
              'production_runs_modified':False,'fixed_mass':True,'weights_zero':True,'generator_steering':False}
    for chi in [1,0]:
        d=root/f'chi{chi}';d.mkdir()
        (d/SOURCE.name).write_bytes(SOURCE.read_bytes())
        (d/'params.dat').write_bytes((TEMPLATE/'params.dat').read_bytes())
        env=json.loads((TEMPLATE/'environment.json').read_text())
        env.update(GFS_QBP_CHI=str(chi),GFS_MAX_DT='.001',HYDRO_TSTOP=str(t+.2-.00025),
                   EUNHA_DUMP_DT='.01',EUNHA_DUMP_FINAL='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1')
        for key in ['GFS_FIXED_DT','GFS_TCF_SACC','GFS_TCF_DUMP']:
            assert key not in env
        (d/'environment.json').write_text(json.dumps(env,indent=2)+'\n')
        (d/'restart.step').write_text('656\n')
        manifest['cases'].append({'name':d.name,'chi':chi})
    (root/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(manifest),flush=True)


def execute(root):
    m=json.loads((root/'manifest.json').read_text())
    assert sha(BINARY)==m['binary_sha256']==EXPECTED and sha(SOURCE)==m['source_sha256']
    def run(case):
        d=root/case['name'];assert sha(d/SOURCE.name)==m['source_sha256']
        env={k:v for k,v in os.environ.items() if not k.startswith(('GFS_','EUNHA_','SEDOV_','LAGMFM_','HYDRO_','OMP_'))}
        env.update(json.loads((d/'environment.json').read_text()))
        begin=time.monotonic()
        with (d/'log').open('x') as log:
            p=subprocess.Popen(['mpirun','-np','1',str(BINARY),'params.dat','656'],cwd=d,env=env,stdout=log,stderr=subprocess.STDOUT)
            (d/'runner_pid').write_text(str(p.pid)+'\n')
            print(json.dumps({'case':case['name'],'runner_pid':p.pid,'started':True}),flush=True)
            code=p.wait()
        (d/'exit.code').write_text(str(code)+'\n')
        print(json.dumps({'case':case['name'],'exit':code,'wall_seconds':time.monotonic()-begin}),flush=True)
        return code
    with ThreadPoolExecutor(max_workers=2) as pool:
        codes=list(pool.map(run,m['cases']))
    assert codes==[0,0],codes


def analyze(root):
    from audit_canonical_pair_history import gravity_coefficient
    m=json.loads((root/'manifest.json').read_text());initial,t0=read(SOURCE)
    ic,_=read(BASE/'Codex_canonical_initial_20261017/run/khout.000000.dat')
    r0=np.hypot(ic['x']-3,ic['y']-3)
    masks={'inner':(r0<.5)&(ic['den']<.5),'inner_edge':(r0<.6)&(ic['den']>.5),
           'disk_core':(r0>=.6)&(r0<1.9)&(ic['den']>.5)}
    pairs=[(803,804),(804,933),(804,934),(804,805)]
    gamma=float(np.float32(float(re.search(r'define Gamma_hydro\s*=\s*([\d.eE+-]+)',(TEMPLATE/'params.dat').read_text())[1])))
    out={'scope':m['scope'],'manifest':m,'cases':{}}
    for case in m['cases']:
        d=root/case['name'];assert (d/'exit.code').read_text().strip()=='0'
        rows=[]
        for path in sorted(d.glob('khout.*.dat')):
            a,t=read(path);assert np.array_equal(a['mass'],initial['mass'])
            assert np.all(a['w2']==0)
            assert all(np.isfinite(a[k]).all() for k in FIELDS)
            assert all(np.all(a[k]>0) for k in ['volume','pressure','ie'])
            x=np.column_stack((a['x'],a['y']));v=np.column_stack((a['vx'],a['vy']))
            xc=x-3;radius=np.linalg.norm(xc,axis=1)
            vc=np.sqrt(gravity_coefficient(radius))[:,None]*np.column_stack((-xc[:,1],xc[:,0]))
            dv2=np.sum((v-vc)**2,axis=1)
            row={'time':t,'snapshot':str(path),'sha256':sha(path),'pairs':[],'regions':{}}
            for i,j in pairs:
                dx=x[j]-x[i];dv=v[j]-v[i];gap=float(np.linalg.norm(dx))
                h=float(np.sqrt(.5*(ic['volume'][i]+ic['volume'][j])))
                row['pairs'].append({'ids':[i,j],'gap_over_initial_h':gap/h,
                                     'separation_rate':float(dx@dv/gap),
                                     'pressure_over_t0':(a['pressure'][[i,j]]/ic['pressure'][[i,j]]).tolist(),
                                     'specific_angular_momentum':((x[[i,j],0]-3)*v[[i,j],1]-(x[[i,j],1]-3)*v[[i,j],0]).tolist()})
            ds=np.log(a['ie']/initial['ie'])+(gamma-1)*np.log(a['volume']/initial['volume'])
            for name,s in masks.items():
                w=a['mass'][s];err=a['pressure'][s]/ic['pressure'][s]-1
                row['regions'][name]={'pressure_error_mass_RMS':float(np.sqrt(np.average(err*err,weights=w))),
                                      'circular_velocity_error_mass_RMS':float(np.sqrt(np.average(dv2[s],weights=w))),
                                      'circular_velocity_perturbation_energy':float(.5*np.sum(w*dv2[s])),
                                      'radial_displacement_from_t0_mass_RMS':float(np.sqrt(np.average((radius[s]-r0[s])**2,weights=w))),
                                      'entropy_change_mass_mean':float(np.average(ds[s],weights=w)),
                                      'thermal_energy':float(a['ie'][s].sum())}
            rows.append(row)
        assert len(rows)>=20
        assert abs(rows[-1]['time']-t0-m['horizon'])<m['dt_cap'],rows[-1]['time']
        lines=(d/'log').read_text().splitlines()
        accepted=[dict(re.findall(r'(\w+)=([^\s]+)',s)) for s in lines if s.startswith('[ACCEPTED_STATE]')]
        assert accepted and all(s['finite']=='1' and s['retries']=='0' for s in accepted)
        rk=[dict(re.findall(r'(\w+)=([^\s]+)',s)) for s in lines if s.startswith('[ACCEPTED_RK_BUDGET]')]
        assert len(rk)==len(accepted) and all(s['stage_mask_min']==s['stage_mask_max']=='15' for s in rk)
        floors=[dict(re.findall(r'(\w+)=([^\s]+)',s)) for s in lines if s.startswith('[RK4E]')]
        assert floors
        # MPI/runtime output can split the long RK4E line. Inspect each named
        # injection field across the full log, requiring one per accepted step.
        for key in ['floor_inj','de_inj','sfl_inj']:
            values=re.findall(r'\b'+key+r'=([^\s]+)', '\n'.join(lines))
            assert len(values)==len(accepted),(key,len(values),len(accepted))
            assert all(float(value)==0 for value in values),key
        out['cases'][case['name']]={'accepted_steps':len(accepted),'RK4_four_stages_verified':True,
                                  'logged_energy_floor_injections_zero':True,
                                  'mass_unchanged':True,'max_abs_step_energy_budget_residual':max(abs(float(s['residual_E'])) for s in rk),
                                  'history':rows}
    with (root/'analysis.json').open('x') as f:json.dump(out,f,indent=2);f.write('\n')
    print(json.dumps({name:{'steps':v['accepted_steps'],'final':v['history'][-1]} for name,v in out['cases'].items()},indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['prepare','execute','analyze']);ap.add_argument('root',type=Path)
    args=ap.parse_args();globals()[args.mode](args.root)
