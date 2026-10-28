"""Same-geometry actual/circular AV and frozen-state chi=0 counterfactual.

Reconstructs face impedance for the changed minmod argument; merely rescaling
the original Q would not implement the code's nonlinear dissipation law.
Only existing diagnostic outputs are read. No production setting is changed.
"""
import argparse
import hashlib
import json
import re
import struct
from pathlib import Path

import numpy as np
from audit_canonical_pair_history import gravity_coefficient


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args();base=Path('/gpfs/kjhan/LagForce')
    layout=json.loads((base/'Codex_DPQ_backend_parity_20261017/layout.json').read_text())
    fields=['x','y','den','mass']
    dtype=np.dtype({'names':['id']+fields,'formats':['<u8']+['<f8']*len(fields),
                    'offsets':[layout['u4if']]+[layout[k] for k in fields], 'itemsize':layout['size']})
    initial_path=base/'Codex_canonical_initial_20261017/run/khout.000000.dat'
    raw=initial_path.read_bytes();_,n,size,_=struct.unpack_from('<iiii',raw)
    assert size==dtype.itemsize and len(raw)==16+n*size+52
    initial=np.frombuffer(raw,dtype=dtype,count=n,offset=16)
    ids=initial['id']&layout['id_mask'];order=np.argsort(ids)
    assert np.array_equal(ids[order],np.arange(n));initial=initial[order]
    r0=np.hypot(initial['x']-3,initial['y']-3)
    masks={'initial_inner':(r0<.5)&(initial['den']<.5),
           'initial_inner_edge':(r0<.6)&(initial['den']>.5),
           'initial_disk_core':(r0>=.6)&(r0<1.9)&(initial['den']>.5)}
    cases={
        't2':(base/'Codex_t2_onset_rhs_20261017/t2',base/'Codex_t2_onset_rhs_20261017/circular'),
        't4':(base/'Codex_TCF_rhs_identity_20261017/t4',base/'Codex_TCF_circular_rhs_20261017/t4'),
        't6':(base/'Codex_TCF_rhs_identity_20261017/t6',base/'Codex_TCF_circular_rhs_20261017/t6')}
    out={'scope':'Saved single-RHS diagnostics at t2/t4/t6. Frozen-state AV activation audit; no evolved counterfactual.',
         'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
         'initial_sha256':hashlib.sha256(raw).hexdigest(),
         'limits':['Actual/circular heat difference is nonlinear and not an additive causal budget.',
                   'chi=0 is reconstructed offline with nonlinear impedance recomputed; no production adoption.',
                   'Face inputs have ~7 significant digits; reconstruction tolerances are reported.',
                   'Circular velocity at current perturbed geometry is a diagnostic, not the original evolved solution.',
                   'TCF only; HG and long-time response excluded.'], 'cases':{}}
    for case,roots in cases.items():
        tables=[];f_tables=[];row={};gamma=None
        for label,root in zip(['actual','circular'],roots):
            assert (root/'exit.code').read_text().strip()=='0'
            env=json.loads((root/'environment.json').read_text())
            assert env['GFS_TCF_FACE_Z']=='1' and env['GFS_TCF_QBP']=='1' and env['GFS_QBP_CHI']=='1'
            assert env.get('GFS_TCF_QLAW','1')=='1'
            gamma=float(np.float32(float(re.search(r'define Gamma_hydro\s*=\s*([\d.eE+-]+)',(root/'params.dat').read_text())[1])))
            a=np.loadtxt(root/'rhs.dat');a=a[np.argsort(a[:,0])]
            f=np.loadtxt(root/'rhs.dat.faces');f=f[np.lexsort((f[:,7],f[:,6],f[:,1],f[:,0]))]
            assert a.shape==(n,25) and f.shape[1]==21
            assert np.array_equal(a[:,0],np.arange(n)) and np.isfinite(a).all() and np.isfinite(f).all()
            tables.append(a);f_tables.append(f)
            i,j=f[:,0].astype(int),f[:,1].astype(int)
            rho=a[:,7]/a[:,6];U=a[:,5]*a[:,6]/(gamma-1)
            rhoL,rhoR=rho[i]*f[:,16],rho[j]*f[:,17]
            pL,pR=a[i,5]*f[:,16]**gamma,a[j,5]*f[:,17]**gamma
            cM=np.maximum(np.sqrt(gamma*pL/rhoL),np.sqrt(gamma*pR/rhoR))
            mu=rhoL*rhoR/(rhoL+rhoR)
            def reconstruct(qin):
                du=f[:,9]
                mm=np.where(qin*du>0,np.copysign(np.minimum(abs(qin),abs(du)),qin),0.)
                ZD=mu*(cM+np.maximum(mm,0))
                ah=(gamma+1)*abs(mm)/8
                ZH=mu*(ah+np.sqrt(ah*ah+cM*cM))
                return np.maximum(ZD,ZH)*mm
            q=reconstruct(f[:,18]);q0=reconstruct(f[:,19])
            relative_q_error=float(abs(q-f[:,14]).sum()/max(abs(f[:,14]).sum(),1e-300))
            assert relative_q_error<3e-6,(case,label,relative_q_error)
            blend=(1-f[:,20])*f[:,19]+f[:,20]*f[:,9]
            blend_error=float(abs(blend-f[:,18]).sum()/max(abs(f[:,18]).sum(),1e-300))
            assert blend_error<3e-6
            sources={'rhs':str(root/'rhs.dat'),'rhs_sha256':hashlib.sha256((root/'rhs.dat').read_bytes()).hexdigest(),
                     'faces_sha256':hashlib.sha256((root/'rhs.dat.faces').read_bytes()).hexdigest(),
                     'binary':(root/'binary.sha256').read_text().strip()}
            models={};heat_arrays={};force_arrays={}
            for mode,Q in [('saved_chi',f[:,14]),('chi_zero',q0)]:
                heat=np.zeros(n);force=np.zeros((n,2))
                np.add.at(heat,i,.5*f[:,4]*Q*f[:,9])
                np.add.at(force,i,-(f[:,4]*Q)[:,None]*f[:,6:8])
                assert np.min(heat)>-1e-20
                models[mode]={'regions':{},'tracked_cells':{}}
                for name,m in masks.items():
                    fm=m[i];w=initial['mass'][m]
                    models[mode]['regions'][name]={
                        'heat_sum':float(heat[m].sum()),
                        'heat_over_U_mass_mean':float(np.average(heat[m]/U[m],weights=w)),
                        'AV_acceleration_mass_RMS':float(np.sqrt(np.average(np.sum(force[m]**2,axis=1)/a[m,7]**2,weights=w))),
                        'chi_quantiles':np.quantile(f[fm,20],[.5,.9,.99,1]).tolist()}
                for k in [803,804,805,933,934]:
                    models[mode]['tracked_cells'][str(k)]={'heat_over_U':float(heat[k]/U[k]),
                                                           'acceleration':(force[k]/a[k,7]).tolist()}
                heat_arrays[mode]=heat;force_arrays[mode]=force
            residual=float(abs(heat_arrays['saved_chi']-a[:,11]).sum()/max(abs(a[:,11]).sum(),1e-300))
            assert residual<3e-6
            powers={}
            if label=='actual':
                x=a[:,1:3]-3;r=np.linalg.norm(x,axis=1);u=r/.5
                omega=np.sqrt(gravity_coefficient(r))
                jx=np.column_stack((-x[:,1],x[:,0]));vc=omega[:,None]*jx
                delta=a[:,3:5]-vc
                # l=(1/r) d(k)/dr for g=-k(r)x. Then
                # delta.Dv_c[delta] = l/(2 omega) (x.delta)(Jx.delta).
                ell=-3/np.maximum(r,1e-30)**5
                lo=u<=.5;mid=(u>.5)&(u<1)
                ell[lo]=(-384/5+96*u[lo])/.5**5
                qmid=u[mid]
                ell[mid]=(1/(5*qmid**5)-48/qmid+384/5-32*qmid)/.5**5
                shear=-a[:,7]*ell/(2*omega)*np.sum(x*delta,axis=1)*np.sum(jx*delta,axis=1)
                pressure_force=a[:,8:10]-force_arrays['saved_chi']
                pressure_work=np.sum(delta*pressure_force,axis=1)
                def circular(xx):
                    om=np.sqrt(gravity_coefficient(np.linalg.norm(xx,axis=1)))
                    return om[:,None]*np.column_stack((-xx[:,1],xx[:,0]))
                fd_errors=[]
                active=r<2.4
                for epsilon in [1e-5,1e-6]:
                    fd=(circular(x+epsilon*delta)-circular(x-epsilon*delta))/(2*epsilon)
                    fd_shear=-a[:,7]*np.sum(delta*fd,axis=1)
                    err=float(abs(fd_shear[active]-shear[active]).sum()/max(abs(shear[active]).sum(),1e-300))
                    assert err<2e-5,(case,epsilon,err)
                    fd_errors.append({'epsilon':epsilon,'relative_L1_error':err})
                for name,mask in masks.items():
                    pw=float(pressure_work[mask].sum());sw=float(shear[mask].sum())
                    qw={mode:float(np.sum(delta[mask]*ff[mask])) for mode,ff in force_arrays.items()}
                    powers[name]={'kinetic_energy_relative_to_circular':float(.5*np.sum(a[mask,7,None]*delta[mask]**2)),
                                  'pressure_work':pw,'differential_rotation_work':sw,'TCF_work':qw,
                                  'total_without_HG':{mode:pw+sw+work for mode,work in qw.items()}}
                powers['derivative_check']=fd_errors
            target_faces=[]
            for k in np.flatnonzero((i==804)&np.isin(j,[803,805,933,934])):
                pb=.5*(a[i[k],5]+a[j[k],5])
                target_faces.append({'i':804,'j':int(j[k]),'chi':float(f[k,20]),
                                     'raw_QBP_delta':float(f[k,19]),'blended_QBP_input':float(f[k,18]),
                                     'cell_velocity_jump':float(f[k,9]),'Q_saved_over_Pmean':float(f[k,14]/pb),
                                     'Q_chi_zero_over_Pmean':float(q0[k]/pb)})
            row[label]={'sources':sources,'Q_reconstruction_relative_L1_error':relative_q_error,
                        'blend_reconstruction_relative_L1_error':blend_error,'heat_reconstruction_relative_L1_error':residual,
                        'models':models,'tracked_faces':target_faces,
                        'actual_velocity_perturbation_power':powers}
        a,b=tables;f,g=f_tables
        assert np.array_equal(a[:,[0,1,2,5,6,7]],b[:,[0,1,2,5,6,7]])
        assert np.array_equal(f[:,:8],g[:,:8]) and np.array_equal(f[:,16:18],g[:,16:18])
        x=b[:,1:3]-3;om=np.sqrt(gravity_coefficient(np.linalg.norm(x,axis=1)))
        vc=om[:,None]*np.column_stack((-x[:,1],x[:,0]))
        active=np.linalg.norm(x,axis=1)<2.4
        error=float(abs(vc[active]-b[active,3:5]).max());assert error<1e-12
        row['same_geometry_and_thermodynamic_state']=True
        row['circular_velocity_max_abs_error']=error
        row['ratios']={}
        for name in masks:
            actual=row['actual']['models']['saved_chi']['regions'][name]['heat_sum']
            actual0=row['actual']['models']['chi_zero']['regions'][name]['heat_sum']
            circ=row['circular']['models']['saved_chi']['regions'][name]['heat_sum']
            circ0=row['circular']['models']['chi_zero']['regions'][name]['heat_sum']
            row['ratios'][name]={'circular_over_actual_heat':circ/actual,
                                 'actual_chi_zero_over_saved_heat':actual0/actual,
                                 'circular_chi_zero_over_saved_heat':circ0/circ}
        out['cases'][case]=row
    with args.output.open('x') as f:json.dump(out,f,indent=2);f.write('\n')
    print(json.dumps({k:{'ratios':v['ratios'],'actual_target':v['actual']['tracked_faces'],
                          'circular_target':v['circular']['tracked_faces']} for k,v in out['cases'].items()},indent=2))


if __name__=='__main__':main()
