"""Saved-checkpoint kinetic energy bookkeeping and subcell error scales.

No new model, state update, RHS or acceptance threshold. The minimum energy
compatible with intrinsic angular momentum is a lower bound, not a unique
rotation/shear decomposition. R energy is unresolved kinetic covariance;
it is not an additional energy added to the already stored total E.
"""
import argparse
import json
from pathlib import Path
import numpy as np
from compare_full_moment_checkpoint_errors import load,REGIONS


def energies(data,prefix):
    raw=data[prefix+'raw'];mass=raw[:,0]
    assert raw.shape==(len(mass),18) and (mass>0).all()
    q=raw[:,1:3];p=raw[:,7:9]
    Q=raw[:,3:7].reshape(-1,2,2);D=raw[:,9:13].reshape(-1,2,2)
    S=Q-q[:,:,None]*q[:,None,:]/mass[:,None,None]
    C=D-q[:,:,None]*p[:,None,:]/mass[:,None,None]
    assert (np.linalg.eigvalsh(S)>0).all()
    B=np.linalg.solve(S,C).transpose(0,2,1)
    affine=.5*np.trace(B@S@B.transpose(0,2,1),axis1=1,axis2=2)
    unresolved=.5*np.trace(data[prefix+'R'],axis1=1,axis2=2)
    center=.5*np.sum(p*p,axis=1)/mass
    ell=C[:,0,1]-C[:,1,0]
    minimum_rotation=ell*ell/(2*np.trace(S,axis1=1,axis2=2))
    thermal=data[prefix+'thermal']
    assert (thermal>0).all() and (unresolved>0).all()
    balance=raw[:,17].astype(np.longdouble)-center.astype(np.longdouble)-affine.astype(np.longdouble)-unresolved.astype(np.longdouble)-thermal.astype(np.longdouble)
    return dict(center_kinetic=center,affine_kinetic=affine,unresolved_kinetic=unresolved,
        intrinsic_angular_momentum=ell,minimum_rotation_kinetic=minimum_rotation,
        thermal=thermal,stored_total_energy=raw[:,17],energy_decomposition_residual=balance,
        minimum_rotation_excess_over_affine=np.maximum(0,minimum_rotation-affine))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--checkpoint',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    data,report,provenance=load(args.checkpoint)
    actual=energies(data,'actual_');reference=energies(data,'reference_')
    weights=data['actual_mass'].astype(np.longdouble);regions={}
    for name in REGIONS:
        mask=data[name]
        if not mask.any():continue
        w=weights[mask]
        def rms(x):
            x=np.asarray(x,dtype=np.longdouble)
            return float(np.sqrt(np.sum(w*x*x)/w.sum()))
        terms={}
        for key in ('center_kinetic','affine_kinetic','unresolved_kinetic','minimum_rotation_kinetic'):
            error=actual[key][mask].astype(np.longdouble)-reference[key][mask].astype(np.longdouble)
            denominator=reference['thermal'][mask].astype(np.longdouble)
            normalized=error/denominator
            terms[key]=dict(error_RMS=rms(error),error_divided_by_reference_thermal_RMS=rms(normalized),
                error_divided_by_reference_thermal_maximum=float(abs(normalized).max()),
                actual_sum=float(actual[key][mask].astype(np.longdouble).sum()),
                reference_sum=float(reference[key][mask].astype(np.longdouble).sum()))
        regions[name]=dict(cells=int(mask.sum()),terms=terms,
            maximum_absolute_energy_decomposition_residual=float(abs(actual['energy_decomposition_residual'][mask]).max()),
            maximum_minimum_rotation_excess_over_affine=float(actual['minimum_rotation_excess_over_affine'][mask].max()),
            reference_thermal_sum=float(reference['thermal'][mask].astype(np.longdouble).sum()))
    result=dict(status='completed_saved_subcell_energy_diagnostics',time=report['time'],source=provenance,
        regions=regions,physical_providers=0,accepted_time_steps=0,new_acceptance_thresholds=False,
        bookkeeping='E = |p|^2/(2m) + trace(B S B^T)/2 + trace(R)/2 + thermal; the ell^2/(2 trace(S)) lower bound belongs inside affine kinetic energy, never added again',
        limits='Subcell differences normalized by thermal energy for scale only; no unique spin/shear split, certified rounding bound, convergence order or long-time accuracy is claimed.')
    args.output.mkdir(exist_ok=False)
    np.savez(args.output/'energies.npz',ids=data['ids'],**{'actual_'+k:v for k,v in actual.items()},**{'reference_'+k:v for k,v in reference.items()})
    (args.output/'report.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(dict(status=result['status'],time=result['time'],all=regions['all'])),flush=True)


if __name__=='__main__':main()
