"""Physical map-residual scales, potential remainder and gauge fixtures."""
import argparse
import ctypes as C
import hashlib
import json
from pathlib import Path

import numpy as np
from check_full_moment_coupled_state import fixture


class Potential(C.Structure):
    _fields_=[('origin',C.c_longdouble*2),('scale',C.c_longdouble),
              ('reference',C.c_longdouble*6),('projection',C.c_longdouble*6),
              ('sigma',C.c_longdouble),('rank_error',C.c_longdouble),('dependent',C.c_int)]


class Residual(C.Structure):
    _fields_=[('maximum',C.c_longdouble*10),('cell',C.c_longlong*10),('dependent_cells',C.c_int)]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--library',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(exist_ok=False)
    lib=C.CDLL(str(args.library))
    lib.full_moment_coupled_residual.argtypes=[C.c_int,C.c_double]+[C.c_void_p]*10
    lib.full_moment_coupled_residual.restype=C.c_int
    lib.full_moment_coupled_converged.argtypes=[C.POINTER(Residual)];lib.full_moment_coupled_converged.restype=C.c_int
    s=fixture([.7,-.4],[.3,.8],np.array([[.3,.04],[.04,.2]]),np.array([[0.,-1.2],[1.2,0.]]),np.array([[.02,.004],[.004,.03]]))
    base=np.r_[s,s,[0.,0.,0.,0.]];mass=np.ones(2);area=np.ones(2);psi=np.array([-1.2,-1.2],dtype=np.longdouble)
    parameters=(Potential*2)()
    for row in parameters:
        row.origin[:]=(.7,-.4);row.scale=1.;row.reference[:]=(-2.,.3,.2,.01,.02,.03)
        row.projection[:]=row.reference[:];row.sigma=1e-6;row.rank_error=1e-10;row.dependent=0
    def run(mapped,potential=parameters,mapped_psi=None):
        mapped_psi=psi.copy() if mapped_psi is None else mapped_psi
        result=Residual();failed=C.c_int(-1)
        ptr=lambda x:x.ctypes.data_as(C.c_void_p)
        code=lib.full_moment_coupled_residual(2,float(np.float32(1.666667)),ptr(base),ptr(mass),ptr(area),ptr(base),ptr(mapped),ptr(psi),ptr(mapped_psi),C.cast(potential,C.c_void_p),C.byref(result),C.byref(failed))
        return code,bool(lib.full_moment_coupled_converged(C.byref(result))) if code==0 else False,result
    cases={}
    code,passed,result=run(base.copy());assert code==0 and passed and max(result.maximum)==0
    cases['identity']=True
    mapped=base.copy();mapped[2]+=.125;mapped[22]+=.125
    code,passed,result=run(mapped);assert code==0 and passed and result.maximum[1]==0
    cases['common_weight_gauge_removed']=True
    mapped=base.copy();mapped[2]+=1e-5
    code,passed,result=run(mapped);assert code==0 and not passed and result.maximum[1]>1e-6
    cases['relative_weight_error_detected']=True
    mapped=base.copy();mapped[0]+=1e-7
    code,passed,result=run(mapped);assert code==0 and not passed and result.maximum[0]>1e-8
    cases['generator_error_detected']=True
    mapped=base.copy();mapped[19]+=1e-8
    code,passed,result=run(mapped);assert code==0 and not passed and result.maximum[7]>1e-8
    cases['cold_thermal_scale_retained']=True
    mapped_psi=psi.copy();mapped_psi[0]+=np.longdouble('2e-14')
    code,passed,result=run(base.copy(),mapped_psi=mapped_psi)
    assert code==0 and not passed and result.maximum[8]>1e-8
    cases['independent_Psi_error_detected']=True
    mapped=base.copy();mapped[3]+=1e-12
    dq=np.longdouble(mapped[3])-np.longdouble(base[3]);ox=np.longdouble(.7);oy=np.longdouble(-.4)
    correction=np.longdouble(.3)*dq+np.longdouble(.01)*(-2*ox*dq)+np.longdouble(.02)*(-oy*dq)
    mapped_psi=psi.copy();mapped_psi[0]+=correction
    code,passed,result=run(mapped,mapped_psi=mapped_psi)
    assert code==0 and passed and result.maximum[8]<1e-8
    assert abs(mapped_psi[0]-psi[0])/np.longdouble('1e-6')>1e-8
    cases['projected_raw_feature_change_removed_from_Psi_residual']=True
    for row in parameters:
        row.dependent=1;row.reference[:]=(-1.2,0.,0.,0.,0.,0.);row.projection[:]=row.reference[:]
        row.sigma=0.;row.rank_error=0.
    code,passed,result=run(base.copy());assert code==0 and passed and result.dependent_cells==2
    cases['declared_dependent_potential_direct_constraint']=True
    mapped_psi=psi.copy();mapped_psi[0]+=np.longdouble('1e-12')
    code,passed,result=run(base.copy(),mapped_psi=mapped_psi)
    assert code==0 and not passed and result.maximum[9]>1
    cases['inconsistent_dependent_potential_rejected']=True
    for row in parameters:row.dependent=0
    code,passed,result=run(base.copy());assert code!=0 and not passed
    cases['independent_rank_failure_not_hidden_by_floor']=True
    report=dict(status='passed_native_coupled_residual_fixtures',cases=cases,physical_providers=0,accepted_time_steps=0,
        library_sha256=hashlib.sha256(args.library.read_bytes()).hexdigest(),
        normal_moment_tolerance=1e-10,thermal_and_independent_Psi_remainder_tolerance=1e-8,
        limits='Residual evaluator only; no nonlinear solve, conditioning certification or physical evolution.')
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)


if __name__=='__main__':main()
