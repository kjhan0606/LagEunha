"""Native arithmetic-midpoint/centered-state fixtures, zero physical providers."""
import argparse
import ctypes as C
from decimal import Decimal,localcontext
import hashlib
import json
from pathlib import Path

import numpy as np


class Centered(C.Structure):
    _fields_=[('center',C.c_longdouble*2),('velocity',C.c_longdouble*2),
              ('S',C.c_longdouble*4),('C',C.c_longdouble*4),('R',C.c_longdouble*4),
              ('thermal',C.c_longdouble)]


def pointer(a):return a.ctypes.data_as(C.c_void_p)


def decimal_R(a,m):
    with localcontext() as context:
        context.prec=80;d=lambda x:Decimal.from_float(float(x));m=d(m)
        q=list(map(d,a[3:5]));p=list(map(d,a[9:11]))
        S=[[d(a[5+2*i+j])-q[i]*q[j]/m for j in range(2)] for i in range(2)]
        B=[[d(a[11+2*i+j])-q[i]*p[j]/m for j in range(2)] for i in range(2)]
        determinant=S[0][0]*S[1][1]-S[0][1]*S[1][0]
        inv=[[S[1][1]/determinant,-S[0][1]/determinant],[-S[1][0]/determinant,S[0][0]/determinant]]
        return [np.longdouble(str(d(a[15+2*i+j])-p[i]*p[j]/m-
            sum(B[k][i]*inv[k][l]*B[l][j] for k in range(2) for l in range(2))))
            for i in range(2) for j in range(2)]


def fixture(center,velocity,S,B,R,m=1.):
    q=m*np.array(center);p=m*np.array(velocity);Q=S+np.outer(q,q)/m
    mixed=S@B.T;D=mixed+np.outer(q,p)/m
    Pi=np.outer(p,p)/m+B@S@B.T+R;Pi[1,0]=Pi[0,1];Q[1,0]=Q[0,1]
    a=np.r_[center,0.,q,Q.ravel(),p,D.ravel(),Pi.ravel(),.5*np.trace(Pi)+.02]
    assert a.shape==(20,);return a


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--library',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(exist_ok=False)
    lib=C.CDLL(str(args.library));lib.full_moment_coupled_midpoint.argtypes=[C.c_int]+[C.c_void_p]*9
    lib.full_moment_coupled_midpoint.restype=C.c_int
    lib.full_moment_centered_decode.argtypes=[C.c_void_p,C.c_double,C.POINTER(Centered)]
    lib.full_moment_centered_decode.restype=C.c_int
    def midpoint(a,b,m,psi_a,psi_b):
        n=len(m);state=np.full(20*n+4,-123.);psi=np.full(n,np.longdouble(-123));L=np.full((n,2,2),-123.)
        failed=np.array([-1],dtype=np.int32)
        code=lib.full_moment_coupled_midpoint(n,pointer(a),pointer(b),pointer(m),pointer(psi_a),pointer(psi_b),pointer(state),pointer(psi),pointer(L),pointer(failed))
        return code,int(failed[0]),state,psi,L
    S=np.array([[.3,.04],[.04,.2]]);R=np.array([[.02,.004],[.004,.03]])
    spin=np.array([[0.,-1.2],[1.2,0.]])
    a=fixture([.7,-.4],[.3,.8],S,spin,R)
    b=fixture([.75,-.35],[.25,.85],S,spin,R)
    ba=np.r_[a,[0.,0.,0.,0.]];bb=np.r_[b,[0.,0.,0.,0.]]
    m=np.array([1.]);pa=np.array([-1.2],dtype=np.longdouble);pb=np.array([-1.1],dtype=np.longdouble)
    code,failed,state,psi,L=midpoint(ba,bb,m,pa,pb)
    assert code==0 and failed==-1
    assert np.array_equal(state,.5*(ba+bb)) and np.array_equal(psi,.5*(pa+pb))
    exact=np.array(decimal_R(state[:20],1.)).reshape(2,2)
    diff=np.max(abs(L[0].astype(np.longdouble)@L[0].astype(np.longdouble).T-exact))
    assert diff<np.longdouble('1e-16')
    failures={}
    for label in ('nonfinite','asymmetric_Q','asymmetric_Pi','nonpositive_S','nonpositive_R','nonpositive_thermal','invalid_mass'):
        bad=ba.copy();bad_m=m.copy()
        if label=='nonfinite':bad[0]=np.nan
        if label=='asymmetric_Q':bad[6]=np.nextafter(bad[6],np.inf)
        if label=='asymmetric_Pi':bad[16]=np.nextafter(bad[16],np.inf)
        if label=='nonpositive_S':bad[5]=-1.
        if label=='nonpositive_R':bad[15]=-1.
        if label=='nonpositive_thermal':bad[19]=0.
        if label=='invalid_mass':bad_m[0]=0.
        code,cell,out,out_psi,out_L=midpoint(bad,bad,bad_m,pa,pa)
        assert code!=0 and np.all(out==-123.) and np.all(out_psi==-123.) and np.all(out_L==-123.)
        failures[label]=dict(status=code,failed_cell=cell,outputs_unchanged=True)
    source=Path('/gpfs/kjhan/LagForce/Codex_full_disk_continuous_two_dt_to_t1_20261019/fine/checkpoints/fine_0010/accepted.npz')
    source_sha=hashlib.sha256(source.read_bytes()).hexdigest()
    with np.load(source) as f:
        saved=f['state'];mass=f['mass'];potential=f['Psi'];old_L=f['factors']
    code,cell,mid,mid_psi,mid_L=midpoint(saved,saved,mass,potential,potential)
    assert code==0,(code,cell)
    assert np.array_equal(mid,saved) and np.array_equal(mid_psi,potential)
    selected=np.unique(np.linspace(0,len(mass)-1,32,dtype=int));maximum=np.longdouble(0)
    for i in selected:
        expected=np.array(decimal_R(mid[20*i:20*i+20],mass[i])).reshape(2,2)
        actual=mid_L[i].astype(np.longdouble)@mid_L[i].astype(np.longdouble).T
        maximum=max(maximum,np.max(abs(actual-expected))/np.longdouble(mass[i]))
    assert maximum<np.longdouble('1e-15')
    assert hashlib.sha256(source.read_bytes()).hexdigest()==source_sha
    report=dict(status='passed_native_coupled_midpoint_state_fixtures',source=str(source),source_sha256=source_sha,
        library_sha256=hashlib.sha256(args.library.read_bytes()).hexdigest(),saved_cells=len(mass),
        identical_endpoint_raw_and_Psi_unchanged=True,offcenter_rotating_nondiagonal_fixture_passed=True,
        decimal_80digit_schur_complement_cells=len(selected),maximum_decimal_R_error_per_mass=float(maximum),
        negative_fixtures=failures,physical_providers=0,accepted_time_steps=0,
        limits='State construction only. No nonlinear midpoint convergence, covariance-map time-order proof, provider or trajectory test. Restored L need not equal stored L bitwise after raw-moment rounding.')
    (args.output/'report.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report),flush=True)


if __name__=='__main__':main()
