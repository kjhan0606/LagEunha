"""Instantaneous perturbation-energy work from canonical C force dumps.

For vc=r*Omega*e_phi and central gravity balancing vc:
 d(1/2 m |v-vc|^2)/dt = (v-vc).Fhydro - m*r*Omega'*dvr*dvt.
HG splitting, frozen cells and finite-step changes are outside this RHS.
"""
import argparse
import json
from pathlib import Path
import numpy as np

ap=argparse.ArgumentParser()
ap.add_argument('root')
ap.add_argument('--output',required=True)
ap.add_argument('--circular-root')
ap.add_argument('--cases', nargs='+', default=['t4', 't6'])
args=ap.parse_args()
out={}
for case in args.cases:
    root=Path(args.root)/case
    c=np.loadtxt(root/'rhs.dat'); c=c[np.argsort(c[:,0])]
    f=np.loadtxt(root/'rhs.dat.faces')
    ids=c[:,0].astype(np.int64)
    owner=np.searchsorted(ids,f[:,0].astype(np.int64))
    assert np.all(owner<len(ids)) and np.array_equal(ids[owner],f[:,0])
    fq=np.zeros((len(c),2))
    np.add.at(fq[:,0],owner,-f[:,4]*f[:,14]*f[:,6])
    np.add.at(fq[:,1],owner,-f[:,4]*f[:,14]*f[:,7])
    x,y=c[:,1]-3,c[:,2]-3
    r=np.sqrt(x*x+y*y); u=r/.5
    assert np.all(r>0)
    mass=np.ones_like(r); derivative=np.zeros_like(r)
    lo=u<=.5; mid=(u>.5)&(u<1)
    q=u[lo]; mass[lo]=q**3*(32/3-192/5*q*q+32*q**3)
    derivative[lo]=(32*q*q-192*q**4+192*q**5)/.5
    q=u[mid]; mass[mid]=-1/15+q**3*(64/3-48*q+192/5*q*q-32/3*q**3)
    derivative[mid]=(64*q*q-192*q**3+192*q**4-64*q**5)/.5
    omega=np.sqrt(mass/r**3)
    romegaprime=.5*omega*(r*derivative/mass-3)
    dv=np.column_stack((c[:,3]+y*omega,c[:,4]-x*omega))
    dvr=(x*dv[:,0]+y*dv[:,1])/r
    dvt=(-y*dv[:,0]+x*dv[:,1])/r
    pressure_force=c[:,8:10]-fq
    uniform_force=-2.7777796321447566e-5*c[:,15:17]
    wq=np.sum(dv*fq,axis=1)
    wp=np.sum(dv*pressure_force,axis=1)
    ws=-c[:,7]*romegaprime*dvr*dvt
    ep=.5*c[:,7]*np.sum(dv*dv,axis=1)
    # Independent directional derivative of perturbation kinetic energy.
    # Move x along v and v along the instantaneous hydro+gravity RHS.
    def circular(xx,yy):
        rr=np.sqrt(xx*xx+yy*yy); uu=rr/.5
        mm=np.ones_like(rr)
        low=uu<=.5; middle=(uu>.5)&(uu<1)
        q=uu[low]; mm[low]=q**3*(32/3-192/5*q*q+32*q**3)
        q=uu[middle]; mm[middle]=-1/15+q**3*(64/3-48*q+192/5*q*q-32/3*q**3)
        oo=np.sqrt(mm/rr**3)
        return np.column_stack((-yy*oo,xx*oo))
    vel=c[:,3:5]
    acc=c[:,8:10]/c[:,7,None]-omega[:,None]**2*np.column_stack((x,y))
    checks=[]
    active=r<2.5
    predicted=wp+wq+ws
    for eps in [1e-5,1e-6]:
        plus=vel+eps*acc-circular(x+eps*vel[:,0],y+eps*vel[:,1])
        minus=vel-eps*acc-circular(x-eps*vel[:,0],y-eps*vel[:,1])
        numeric=.25*c[:,7]*(np.sum(plus*plus,axis=1)-np.sum(minus*minus,axis=1))/eps
        relative=float(np.abs(numeric[active]-predicted[active]).sum()/np.abs(predicted[active]).sum())
        checks.append({'epsilon':eps,'relative_L1_error':relative})
    assert checks[-1]['relative_L1_error']<1e-5
    out[case]={}
    out[case]['directional_derivative_checks']=checks
    for name,mask in [('inner',r<.5),('disk',(r>=.5)&(r<=2)),('outer_active',(r>2)&(r<2.5))]:
        out[case][name]={'cells':int(mask.sum()),'perturbation_kinetic_energy':float(ep[mask].sum()),
                        'pressure_work':float(wp[mask].sum()),'Q_work':float(wq[mask].sum()),
                        'differential_rotation_work':float(ws[mask].sum()),
                        'total_RHS':float((wp+wq+ws)[mask].sum()),
                        'Q_heat':float(c[mask,11].sum()),
                        'uniform_pressure_work':float(np.sum(dv[mask]*uniform_force[mask])),
                        'uniform_pressure_force_L1':float(np.linalg.norm(uniform_force[mask],axis=1).sum()),
                        'actual_pressure_force_L1':float(np.linalg.norm(pressure_force[mask],axis=1).sum())}
    if args.circular_root:
        circ=np.loadtxt(Path(args.circular_root)/case/'rhs.dat')
        circ=circ[np.argsort(circ[:,0])]
        assert np.array_equal(c[:,[0,1,2,5,6,7]],circ[:,[0,1,2,5,6,7]])
        # dU=-P*dV; kinetic pressure work has the opposite sign: +P*dV.
        volume_work=c[:,5]*(c[:,12]-circ[:,12])
        out[case]['global_virtual_work']={
            'pressure_force_work':float(wp.sum()),
            'positive_P_delta_volume_rate':float(volume_work.sum()),
            'difference':float(wp.sum()-volume_work.sum()),
            'relative_difference':float(abs(wp.sum()-volume_work.sum())/max(abs(wp.sum()),abs(volume_work.sum()),1e-300)),
            'scope':'Whole dump including frozen cells; delta velocity actual minus circular. Regional equality is not required because faces transfer energy across region boundaries. Face-output rounding limits force reconstruction.'}
result={'cases':out,'scope':'Instantaneous material sum for cells currently in each region, not derivative of a fixed Eulerian-region integral. Excludes HG split kicks. Q force reconstructed from 7-digit directed-face output; pressure force is total minus Q. Canonical unweighted TCF, no other viscosity.'}
Path(args.output).write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
