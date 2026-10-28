"""Resolve saved pre-HG hydrodynamic RHS along actual close-pair directions.

Uses existing dumps only. TCF face records have about seven significant digits;
pressure is the total hydro force minus reconstructed Q. HG is excluded.
"""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from audit_canonical_pair_history import gravity_coefficient


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--history',type=Path,required=True)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    base=Path('/gpfs/kjhan/LagForce')
    history=json.loads((args.history/'result.json').read_text())
    data=np.load(args.history/'pair_history.npz')
    pairs=data['pairs']; h=data['h0']; times=data['times']
    cases={'t2':base/'Codex_t2_onset_rhs_20261017/t2',
           't4':base/'Codex_TCF_rhs_identity_20261017/t4',
           't6':base/'Codex_TCF_rhs_identity_20261017/t6',
           't10':base/'Codex_TCF_rhs_t10_20261017/t10'}
    layout=json.loads((base/'Codex_DPQ_backend_parity_20261017/layout.json').read_text())
    fields=['x','y','vx','vy','mass','pressure','volume']
    dtype=np.dtype({'names':['id']+fields,'formats':['<u8']+['<f8']*len(fields),
                    'offsets':[layout['u4if']]+[layout[k] for k in fields], 'itemsize':layout['size']})
    result={'scope':'Existing pre-HG saved RHS, one rank; no new evolution.',
            'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            'formula':'d_double_dot=n.(a_j-a_i)+|v_relative_transverse|^2/d',
            'limits':['HG excluded; no total applied-step acceleration claim.',
                      'TCF reconstruction limited by printed face precision (~7 significant digits).',
                      'Pressure force inferred as saved hydro force minus reconstructed TCF Q.',
                      'Instantaneous acceleration does not establish a mode stiffness or continuous binding.'],
            'cases':{}}
    for name,root in cases.items():
        assert (root/'exit.code').read_text().strip()=='0'
        rhs=np.loadtxt(root/'rhs.dat');rhs=rhs[np.argsort(rhs[:,0])]
        faces=np.loadtxt(root/'rhs.dat.faces')
        assert rhs.shape==(36864,25) and faces.shape[1]==21
        assert np.array_equal(rhs[:,0],np.arange(len(rhs)))
        assert np.isfinite(rhs).all() and np.isfinite(faces).all()
        index=int(np.argmin(abs(times-float(name[1:]))))
        path=Path(history['sources'][index]['path'])
        raw=path.read_bytes()
        assert hashlib.sha256(raw).hexdigest()==history['sources'][index]['sha256']
        a=np.frombuffer(raw,dtype=dtype,count=len(rhs),offset=16)
        a=a[np.argsort(a['id']&layout['id_mask'])]
        match={k:float(np.max(abs(a[k]-rhs[:,col]))) for k,col in
               [('x',1),('y',2),('vx',3),('vy',4),('pressure',5),('mass',7)]}
        # If a split kick has changed velocity, do not silently identify the
        # force dump as the same snapshot state.
        assert max(match.values())<1e-12,(name,match)
        ii=faces[:,0].astype(int);jj=faces[:,1].astype(int)
        assert np.all((ii>=0)&(ii<len(rhs))&(jj>=0)&(jj<len(rhs)))
        force_q=np.zeros((len(rhs),2))
        np.add.at(force_q,ii,-(faces[:,4]*faces[:,14])[:,None]*faces[:,6:8])
        heat=np.zeros(len(rhs))
        np.add.at(heat,ii,.5*faces[:,4]*faces[:,14]*faces[:,9])
        active=np.hypot(rhs[:,1]-3,rhs[:,2]-3)<2.4
        heat_res=float(abs(heat[active]-rhs[active,11]).sum()/max(abs(rhs[active,11]).sum(),1e-300))
        assert heat_res<3e-6,(name,heat_res)
        mass=rhs[:,7];x=rhs[:,1:3]-3;v=rhs[:,3:5]
        grav=-gravity_coefficient(np.linalg.norm(x,axis=1))[:,None]*x
        aq=force_q/mass[:,None];apress=(rhs[:,8:10]-force_q)/mass[:,None]
        selected=np.argsort(data['separation'][index]/h)[:12]
        records=[]
        for pair_index in selected:
            i,j=pairs[pair_index]
            dx=x[j]-x[i];dv=v[j]-v[i];d=float(np.linalg.norm(dx));n=dx/d
            vd=float(dv@n);vt2=max(0.,float(dv@dv)-vd*vd)
            acc={'gravity':grav[j]-grav[i],'pressure':apress[j]-apress[i],'TCF_Q':aq[j]-aq[i]}
            radial={key:float(n@val) for key,val in acc.items()}
            kin=vt2/d
            mu=mass[i]*mass[j]/(mass[i]+mass[j])
            work={key:float(mu*dv@val) for key,val in acc.items()}
            direct=(ii==i)&(jj==j)
            reverse=(ii==j)&(jj==i)
            pb=.5*(rhs[i,5]+rhs[j,5])
            # For this unweighted Voronoi run, the pressure geometry term
            # tangent to the shared face has zero projection along x_j-x_i.
            # Evaluate both directed records to retain printed asymmetry.
            fp_i=np.sum(-(faces[direct,4]*pb)[:,None]*faces[direct,6:8],axis=0)
            fp_j=np.sum(-(faces[reverse,4]*pb)[:,None]*faces[reverse,6:8],axis=0)
            fq_i=np.sum(-(faces[direct,4]*faces[direct,14])[:,None]*faces[direct,6:8],axis=0)
            fq_j=np.sum(-(faces[reverse,4]*faces[reverse,14])[:,None]*faces[reverse,6:8],axis=0)
            direct_pressure=float(n@(fp_j/mass[j]-fp_i/mass[i]))
            direct_q=float(n@(fq_j/mass[j]-fq_i/mass[i]))
            records.append({'ids':[int(i),int(j)],'gap':d,
                            'gap_over_initial_h':float(d/h[pair_index]),
                            'separation_rate':vd,'relative_speed':float(np.linalg.norm(dv)),
                            'radial_relative_acceleration':radial,
                            'transverse_kinematic_term':kin,
                            'separation_second_derivative_without_HG':float(sum(radial.values())+kin),
                            'relative_kinetic_power':work,
                            'pressure_force_interpretation':'Total pressure acceleration from all faces, not direct pair attraction.',
                            'direct_shared_face_pressure_radial_relative_acceleration':direct_pressure,
                            'other_faces_pressure_radial_relative_acceleration':radial['pressure']-direct_pressure,
                            'direct_shared_face_TCF_radial_relative_acceleration':direct_q,
                            'other_faces_TCF_radial_relative_acceleration':radial['TCF_Q']-direct_q,
                            'cell_pressure':rhs[[i,j],5].tolist(),
                            'cell_density':(mass[[i,j]]/rhs[[i,j],6]).tolist(),
                            'cell_pressure_acceleration_magnitude':np.linalg.norm(apress[[i,j]],axis=1).tolist(),
                            'cell_gravity_acceleration_magnitude':np.linalg.norm(grav[[i,j]],axis=1).tolist(),
                            'cell_pressure_specific_torque':(x[[i,j],0]*apress[[i,j],1]-x[[i,j],1]*apress[[i,j],0]).tolist(),
                            'share_recorded_face':bool(direct.any()),
                            'direct_faces_Q':faces[direct,14].tolist(),
                            'direct_faces_applied_impedance':faces[direct,15].tolist(),
                            'direct_faces_compressive_velocity_jump':faces[direct,9].tolist(),
                            'direct_faces_reconstructed_velocity_jump':faces[direct,8].tolist()})
        result['cases'][name]={'time':float(times[index]),'source':str(root),
                              'rhs_sha256':hashlib.sha256((root/'rhs.dat').read_bytes()).hexdigest(),
                              'face_sha256':hashlib.sha256((root/'rhs.dat.faces').read_bytes()).hexdigest(),
                              'snapshot_match_max_abs':match,'TCF_heat_reconstruction_relative_L1':heat_res,
                              'closest_pairs':records}
    with args.output.open('x') as f:json.dump(result,f,indent=2);f.write('\n')
    print(json.dumps({k:{'time':v['time'],'closest':v['closest_pairs'][0]} for k,v in result['cases'].items()},indent=2))


if __name__=='__main__':main()
