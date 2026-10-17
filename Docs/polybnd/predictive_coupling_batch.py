"""Autonomous basic gates: saved-state full HLLC and positive kinetic coupling."""
import json,time
from pathlib import Path
import numpy as np
from adaptive_coupled_mass import atomic
from adaptive_moment_batch import setup
import moment_ale_evolve as model
import kinetic_factor_rhs as kinetic
from predictive_boundary_density import BoundaryDensity


def affine_rhs_errors(dy,mass):
    a=dy[:-4].reshape(-1,20);dq=a[:,3:5];dQ=a[:,5:9].reshape(-1,2,2)
    u=np.array([.2,0.]);B=np.array([[0.,1.1],[0.,0.]])
    expectedp=dq@B.T
    expectedD=np.einsum('ni,j->nij',dq,u)+dQ@B.T
    bdq=dq@B.T
    expectedPi=np.einsum('ni,j->nij',bdq,u)+np.einsum('i,nj->nij',u,bdq)+B@dQ@B.T
    return dict(p=float(np.max(abs(a[:,9:11]-expectedp)/mass[:,None])),
        D=float(np.max(abs(a[:,11:15].reshape(-1,2,2)-expectedD)/mass[:,None,None])),
        Pi=float(np.max(abs(a[:,15:19].reshape(-1,2,2)-expectedPi)/mass[:,None,None])))


def metrics(dy,diag,mass):
    a=dy[:-4].reshape(-1,20)
    ledger=np.r_[np.sum(a[:,9:11],axis=0),np.sum(a[:,19]),np.sum(a[:,12]-a[:,13])]-dy[-4:]
    return dict(mass_flux=float(diag['mass_flux']),closure_moment_error=float(diag['closure_moment_error']),
        density_min=float(diag['closure_min_density']),density_max=float(diag['closure_max_density']),
        closure_volume_cap_excess=float(diag['closure_volume_cap_excess']),pressure_rate=float(max(abs(diag['pressure_rates']))),
        affine_rhs_errors=affine_rhs_errors(dy,mass),ledger_rate_error=float(max(abs(ledger))),
        generator_velocity_defect=float(diag['generator_velocity_defect']),
        max_weight_rate=float(max(abs(a[:,2]))))


def saved_case(name):
    cp=json.loads((Path('inputs')/name/'checkpoint.json').read_text());cfg=cp['config']
    if cfg['kind']!='shear_contact' or cp['t']!=.01:raise RuntimeError('unexpected saved input')
    rebuilt,_=setup('shear_contact',cfg['n'],cfg['amplitude']);mass=np.asarray(cp['mass'])
    np.testing.assert_allclose(mass,rebuilt,rtol=1e-13,atol=1e-15)
    return cp,mass,np.asarray(cp['state']),[np.asarray(f) for f in cp['factors']]


def static(name,mode):
    path=Path('static_'+name+'_'+mode+'.json')
    if path.exists():
        old=json.loads(path.read_text())
        if old['passed']:return old
        raise RuntimeError('failed static result requires a diagnosed revision')
    started=time.monotonic();cp,mass,y,L=saved_case(name);provider=BoundaryDensity(mode)
    args=dict(face_density_provider=provider,face_segments=64,face_order=4)
    full,df=model.rhs(cp['t'],y,mass,**args)
    primary,dk=kinetic.rhs(cp['t'],y,mass,covariance_factors=L,exact_bin_transport=True,**args)
    low,dl=kinetic.rhs(cp['t'],y,mass,covariance_factors=L,exact_bin_transport=True,**dict(args,face_order=2))
    difference=float(np.max(abs(primary-low)/np.maximum(1.,abs(primary))))
    vals=dict(full_hllc=metrics(full,df,mass),positive_kinetic=metrics(primary,dk,mass))
    passed=all(v['mass_flux']<1e-10 and v['closure_moment_error']<1e-8 and v['density_min']>0
        and v['closure_volume_cap_excess']<1e-8 and v['pressure_rate']<1e-8 and max(v['affine_rhs_errors'].values())<1e-8
        and v['ledger_rate_error']<1e-10 and v['generator_velocity_defect']<1e-12 for v in vals.values()) and difference<1e-8
    result=dict(case=name,mode=mode,passed=passed,metrics=vals,split_order2_vs4_difference=difference,
        seconds=time.monotonic()-started,note='No exact density/face data supplied. Full HLLC and current positive kinetic closure both re-solve weights. Segmented face profile; kinetic positive-part roots split exactly.')
    atomic(path,result);return result


def main():
    results={}
    for name in ('n3_r10','n3_r100','n4_r100'):
        for mode in ('neighbor','moment'):
            key=name+'_'+mode
            try:results[key]=static(name,mode)
            except Exception as exc:results[key]=dict(passed=False,error=repr(exc))
            atomic('static_progress.json',results)
            print('static',key,'passed',results[key]['passed'],flush=True)
    report=dict(completed=True,passed=all(x['passed'] for x in results.values()),cases=results)
    atomic('static_report.json',report)
    if not report['passed']:raise SystemExit(2)
    import sys
    if '--static-only' in sys.argv:return
    from predictive_evolution_batch import run_basic_sequence
    run_basic_sequence()


if __name__=='__main__':main()
