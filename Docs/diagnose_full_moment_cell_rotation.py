"""Saved-cell rotation/objectivity diagnostics with no added energy or RHS.

ell/trace(S) is a specified mass-metric rigid-rotation projection. It need
not equal half the affine vorticity in an anisotropic deforming cell. This
projection does not uniquely identify all physical spin versus shear.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from compare_full_moment_checkpoint_errors import load,REGIONS


def rotation(raw):
    a=np.asarray(raw,dtype=np.longdouble);m=a[:,0];q=a[:,1:3];p=a[:,7:9]
    S=a[:,3:7].reshape(-1,2,2)-q[:,:,None]*q[:,None,:]/m[:,None,None]
    C=a[:,9:13].reshape(-1,2,2)-q[:,:,None]*p[:,None,:]/m[:,None,None]
    determinant=S[:,0,0]*S[:,1,1]-S[:,0,1]*S[:,1,0]
    assert (m>0).all() and (S[:,0,0]>0).all() and (determinant>0).all()
    inverse=np.empty_like(S);inverse[:,0,0]=S[:,1,1];inverse[:,1,1]=S[:,0,0]
    inverse[:,0,1]=-S[:,0,1];inverse[:,1,0]=-S[:,1,0];inverse/=determinant[:,None,None]
    B=(inverse@C).transpose(0,2,1)
    inertia=S[:,0,0]+S[:,1,1];ell=C[:,0,1]-C[:,1,0];omega=ell/inertia
    rigid=np.zeros_like(B);rigid[:,0,1]=-omega;rigid[:,1,0]=omega
    remainder=B-rigid
    trace=lambda x:np.trace(x,axis1=1,axis2=2)
    affine=.5*trace(B@S@B.transpose(0,2,1));rigid_energy=.5*omega*omega*inertia
    remainder_energy=.5*trace(remainder@S@remainder.transpose(0,2,1))
    return dict(ell=ell,inertia=inertia,angular_momentum_rotation_rate=omega,
        half_affine_vorticity=.5*(B[:,1,0]-B[:,0,1]),rigid_projection_energy=rigid_energy,
        affine_energy=affine,orthogonal_affine_remainder_energy=remainder_energy,
        energy_projection_identity_residual=affine-rigid_energy-remainder_energy,
        angular_projection_identity_residual=ell-omega*inertia)


def fixture_raw(S,B):
    raw=np.zeros((1,18));raw[:,0]=1.;raw[:,3:7]=np.asarray(S).reshape(1,4)
    raw[:,9:13]=(S@B.T).reshape(1,4)
    Pi=B@S@B.T+np.eye(2)*.01;raw[:,13:17]=Pi.reshape(1,4);raw[:,17]=.5*np.trace(Pi)+.1
    return raw


def fixtures():
    S=np.array([[2.,.3],[.3,1.]])
    B=np.array([[0.,-2.1],[2.1,0.]])
    r=rotation(fixture_raw(S,B));assert abs(r['angular_momentum_rotation_rate'][0]-2.1)<1e-13
    assert abs(r['half_affine_vorticity'][0]-2.1)<1e-13
    assert r['orthogonal_affine_remainder_energy'][0]<1e-26
    theta=.31;Q=np.array([[np.cos(theta),-np.sin(theta)],[np.sin(theta),np.cos(theta)]])
    transformed=rotation(fixture_raw(Q@S@Q.T,Q@B@Q.T))
    for key in ('ell','inertia','angular_momentum_rotation_rate','half_affine_vorticity','rigid_projection_energy','affine_energy'):
        assert abs(r[key][0]-transformed[key][0])<1e-13
    strain=rotation(fixture_raw(S,np.array([[1.,.5],[.5,-1.]])))
    assert abs(strain['half_affine_vorticity'][0])<1e-15
    assert abs(strain['angular_momentum_rotation_rate'][0])>.01
    assert abs(strain['energy_projection_identity_residual'][0])<1e-13
    return dict(rigid_rotation=True,rotated_frame_invariance=True,
                anisotropic_strain_not_misidentified_as_vorticity=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--checkpoint',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();assert not args.output.exists();checks=fixtures()
    data,physical,source=load(args.checkpoint)
    actual=rotation(data['actual_raw']);reference=rotation(data['reference_raw'])
    mass=np.asarray(data['actual_mass'],dtype=np.longdouble);regions={}
    for name in REGIONS:
        mask=data[name]
        if not mask.any():continue
        w=mass[mask];rms=lambda x:float(np.sqrt(np.sum(w*np.asarray(x,dtype=np.longdouble)**2)/w.sum()))
        region=dict(cells=int(mask.sum()))
        for key in ('ell','angular_momentum_rotation_rate','half_affine_vorticity'):
            region[key]=dict(actual_RMS=rms(actual[key][mask]),reference_RMS=rms(reference[key][mask]),
                error_RMS=rms(actual[key][mask]-reference[key][mask]))
        region['rigid_projection_energy_error_over_reference_thermal_RMS']=rms(
            (actual['rigid_projection_energy'][mask]-reference['rigid_projection_energy'][mask])/data['reference_thermal'][mask])
        region['maximum_energy_projection_identity_residual']=float(abs(actual['energy_projection_identity_residual'][mask]).max())
        region['maximum_angular_projection_identity_residual']=float(abs(actual['angular_projection_identity_residual'][mask]).max())
        region['projection_rate_minus_half_vorticity_RMS']=rms((actual['angular_momentum_rotation_rate']-actual['half_affine_vorticity'])[mask])
        regions[name]=region
    args.output.mkdir()
    np.savez(args.output/'rotation.npz',ids=data['ids'],**{'actual_'+k:v for k,v in actual.items()},**{'reference_'+k:v for k,v in reference.items()})
    report=dict(status='completed_saved_cell_rotation_diagnostics',time=physical['time'],source=source,
        fixtures=checks,regions=regions,physical_providers=0,accepted_time_steps=0,
        bookkeeping='Rigid projection energy is PART of affine kinetic energy; never add it to stored E. ell/trace(S) and half-vorticity are separately reported.',
        limits='Saved-state diagnostics only. No independent spin DOF, imposed viscosity or steering added. Torque evolution, objective transport and Kepler shear-error suppression still require actual operator/evolution validation.')
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(status=report['status'],time=report['time'],all=regions['all'])),flush=True)


if __name__=='__main__':main()
