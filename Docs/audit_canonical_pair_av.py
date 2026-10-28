"""Inspect AV around the t6 closest pair using saved faces and accepted steps.

HG acceleration is a dt->0 estimate from completed steps, not a direct dump.
No simulation is launched; all pressure, mass, velocity states remain intact.
"""
import argparse
import hashlib
import json
import struct
from pathlib import Path

import numpy as np
from audit_canonical_pair_history import gravity_coefficient


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();base=Path('/gpfs/kjhan/LagForce')
    layout=json.loads((base/'Codex_DPQ_backend_parity_20261017/layout.json').read_text())
    fields=['x','y','vx','vy','mass','volume','pressure','ie','w2']
    dtype=np.dtype({'names':['id']+fields,'formats':['<u8']+['<f8']*len(fields),
                    'offsets':[layout['u4if']]+[layout[k] for k in fields], 'itemsize':layout['size']})
    def read(path):
        raw=path.read_bytes();magic,n,size,stress=struct.unpack_from('<iiii',raw)
        assert magic==-1 and stress==1 and size==dtype.itemsize and len(raw)==16+n*size+52
        a=np.frombuffer(raw,dtype=dtype,count=n,offset=16)
        ids=a['id']&layout['id_mask'];order=np.argsort(ids)
        assert np.array_equal(ids[order],np.arange(n))
        return a[order],struct.unpack_from('<d',raw,16+n*size)[0],hashlib.sha256(raw).hexdigest()
    root=base/'Codex_TCF_rhs_identity_20261017/t6'
    a,t,sha=read(root/'khout.000656.dat');assert np.all(a['w2']==0)
    rhs=np.loadtxt(root/'rhs.dat');rhs=rhs[np.argsort(rhs[:,0])]
    faces=np.loadtxt(root/'rhs.dat.faces')
    assert np.array_equal(rhs[:,0],np.arange(len(a))) and np.isfinite(faces).all()
    for k,col in [('x',1),('y',2),('vx',3),('vy',4),('mass',7)]:
        assert np.max(abs(a[k]-rhs[:,col]))<1e-12
    x=rhs[:,1:3]-3;v=rhs[:,3:5];mass=rhs[:,7]
    grav=-gravity_coefficient(np.linalg.norm(x,axis=1))[:,None]*x
    ahydro=rhs[:,8:10]/mass[:,None]
    ii=faces[:,0].astype(int);jj=faces[:,1].astype(int)
    fq=-(faces[:,4]*faces[:,14])[:,None]*faces[:,6:8]
    aq=np.zeros_like(x);np.add.at(aq,ii,fq);aq/=mass[:,None]
    apress=ahydro-aq
    pi,pj=803,804;n=(x[pj]-x[pi])/np.linalg.norm(x[pj]-x[pi])
    gamma=float(np.float32(1.666667))
    # Use each run's actual gamma setting rather than an assumed ideal value.
    import re
    gamma=float(np.float32(float(re.search(r'define Gamma_hydro\s*=\s*([\d.eE+-]+)',(root/'params.dat').read_text())[1])))
    cs=np.sqrt(gamma*rhs[:,5]*rhs[:,6]/mass)
    face_rows=[]
    for k in np.flatnonzero((ii==pi)|(ii==pj)):
        i,j=ii[k],jj[k];q=faces[k,14];pb=.5*(rhs[i,5]+rhs[j,5])
        du=faces[k,9];qin=faces[k,18]
        mm=0. if qin*du<=0 else np.copysign(min(abs(qin),abs(du)),qin)
        assert abs(q-faces[k,15]*mm)<3e-6*max(abs(q),abs(faces[k,15]*mm),1e-20)
        heat=.5*faces[k,4]*q*du
        face_rows.append({'i':int(i),'j':int(j),'distance':float(faces[k,5]),
                          'cell_velocity_jump_over_max_sound_speed':float(du/max(cs[i],cs[j])),
                          'QBP_input_over_max_sound_speed':float(qin/max(cs[i],cs[j])),
                          'Q_over_mean_pressure':float(q/pb),'Q':float(q),
                          'compression_blend_chi':float(faces[k,20]),
                          'face_density_factors':faces[k,16:18].tolist(),
                          'AV_acceleration_along_803_to_804':float(fq[k]@n/mass[i]),
                          'AV_acceleration_magnitude':float(np.linalg.norm(fq[k])/mass[i]),
                          'heat_rate':float(heat),'heat_rate_over_cell_U':float(heat/a['ie'][i])})
    env=json.loads((root/'environment.json').read_text())
    assert env.get('GFS_TCF_QLAW','1')=='1' and env['GFS_TCF_FACE_Z']=='1'
    snapshots=[];acc=[];dts=[];envs=[];binary=[]
    for case,expected in [('t6_full',.001),('t6_half',.0005),('t6_quarter',.00025)]:
        d=base/'Codex_early_growth_budget_20261017'/case
        assert (d/'exit.code').read_text().strip()=='0'
        start,ts,hs=read(d/'khout.000656.dat');end,te,he=read(d/'khout.000657.dat')
        assert hs==sha and np.array_equal(start['mass'],end['mass'])
        dt=te-ts;assert abs(dt-expected)<3e-15
        log=(d/'log').read_text();accepted=[s for s in log.splitlines() if s.startswith('[ACCEPTED_STATE]')]
        assert len(accepted)==1
        record=dict(re.findall(r'(\w+)=([^\s]+)',accepted[0]))
        assert record['finite']=='1' and record['retries']=='0'
        acc.append((np.column_stack((end['vx'],end['vy']))-v)/dt);dts.append(dt)
        binary.append((d/'binary.sha256').read_text().split()[0]);envs.append(json.loads((d/'environment.json').read_text()))
        snapshots.append({'run':str(d),'dt':dt,'start_sha256':hs,'end_sha256':he,'binary_sha256':binary[-1]})
    assert len(set(binary))==1
    for key in envs[0]:
        if key not in ['GFS_MAX_DT','HYDRO_TSTOP']:
            assert all(e[key]==envs[0][key] for e in envs)
    # First-order extrapolation of accepted mean acceleration. Remaining
    # error includes finite dt, operator splitting and all untracked sources.
    extrap=[2*acc[1]-acc[0],2*acc[2]-acc[1]]
    inferred=[z-grav-ahydro for z in extrap]
    pair_delta=lambda z:float(n@(z[pj]-z[pi]))
    cells=[]
    for i in [pi,pj]:
        cells.append({'id':i,'pressure':float(rhs[i,5]),'density':float(mass[i]/rhs[i,6]),
                      'sound_speed':float(cs[i]),'TCF_heat_rate_over_U':float(rhs[i,11]/a['ie'][i]),
                      'pressure_acceleration':apress[i].tolist(),'TCF_acceleration':aq[i].tolist(),
                      'gravity_acceleration':grav[i].tolist(),
                      'HG_plus_unaccounted_acceleration_extrapolated_coarse':inferred[0][i].tolist(),
                      'HG_plus_unaccounted_acceleration_extrapolated_fine':inferred[1][i].tolist(),
                      'extrapolation_difference_norm':float(np.linalg.norm(inferred[1][i]-inferred[0][i]))})
    out={'scope':'t6 pair 803/804, existing diagnostics and accepted steps; no AV settings changed.',
         'source_time':t,'source_sha256':sha,'weights_all_zero':True,
         'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'rhs_sha256':hashlib.sha256((root/'rhs.dat').read_bytes()).hexdigest(),
         'faces_sha256':hashlib.sha256((root/'rhs.dat.faces').read_bytes()).hexdigest(),
         'TCF_law':'Q=Z_face*minmod(QBP_input,cell_velocity_jump); default QLAW=1; positive A*Q*du heat.',
         'face_compression_ratio_definition':'Normal pair velocity jump / max(cell sound speeds); not a shock Mach number or shock detection.',
         'cells':cells,'faces':face_rows,'accepted_steps':snapshots,
         'pair_radial_relative_acceleration':{'pressure':pair_delta(apress),'TCF':pair_delta(aq),'gravity':pair_delta(grav),
             'HG_plus_unaccounted_extrapolated_coarse':pair_delta(inferred[0]),
             'HG_plus_unaccounted_extrapolated_fine':pair_delta(inferred[1]),
             'accepted_mean_by_dt':[pair_delta(z) for z in acc]},
         'limits':['HG estimate is total accepted acceleration minus saved pre-HG hydro and initial gravity; not independently instrumented HG.',
                   'No causal inference about earlier shocks from a single-time face audit.',
                   'A local TCF contribution can accelerate a selected pair while the global operator dissipates kinetic energy.',
                   'TCF faces printed at about seven significant digits; no bitwise-force assertion.']}
    with args.output.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__':main()
