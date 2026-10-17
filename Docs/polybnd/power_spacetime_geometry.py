"""Static 2D geometry identity gate: fixed sites, affine weight path, fixed box.
Enumerate topology events, integrate shared swept-face geometric moments.
No fluid time evolution, pressure update, or moving-generator claim.
"""
import itertools,json
import numpy as np
from weighted_mass_patch import power_cells
from positive_transport_patch import rule

DOMAIN=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])


def constraints(x,w0,dw,i):
    normals=[];b0=[];b1=[]
    for j in range(len(x)):
        if j==i:continue
        normals.append(2*(x[j]-x[i]));b0.append(x[j]@x[j]-x[i]@x[i]+w0[i]-w0[j]);b1.append(dw[i]-dw[j])
    normals.extend([[-1.,0.],[1.,0.],[0.,-1.],[0.,1.]])
    b0.extend([0.,1.,0.,1.]);b1.extend([0.,0.,0.,0.])
    return np.array(normals),np.array(b0),np.array(b1)


def event_times(x,w0,w1):
    candidates=[0.,1.];dw=w1-w0
    for i in range(len(x)):
        A,b,c=constraints(x,w0,dw,i)
        for j,k in itertools.combinations(range(len(A)),2):
            M=A[[j,k]]
            if abs(np.linalg.det(M))<1e-13:continue
            v0=np.linalg.solve(M,b[[j,k]]);v1=np.linalg.solve(M,c[[j,k]])
            for l in range(len(A)):
                if l in (j,k):continue
                slope=A[l]@v1-c[l]
                if abs(slope)<1e-14:continue
                t=(b[l]-A[l]@v0)/slope
                if 1e-12<t<1-1e-12 and np.max(A@(v0+t*v1)-b-t*c)<1e-11:
                    candidates.append(float(t))
    ordered=sorted(candidates);out=[ordered[0]]
    for t in ordered[1:]:
        if t-out[-1]>1e-10:out.append(t)
    return np.array(out)


def basis(p):
    p=np.asarray(p)
    return np.column_stack([np.ones(len(p)),p[:,0],p[:,1],p[:,0]**2,p[:,0]*p[:,1],p[:,1]**2])


def moments(cells):
    result=[]
    for cell in cells:
        pts,q=rule(cell);result.append(q@basis(pts))
    return np.array(result)


def rate(x,w,dw):
    cells=power_cells(x,w,DOMAIN);r=np.zeros((len(x),6));topology=[]
    for i,cell in enumerate(cells):
        if cell.area<=0:raise RuntimeError('empty cell')
        p=np.asarray(cell.exterior.coords[:-1]);delta=x-x[i]
        bound=np.sum(delta*delta,axis=1)+w[i]-w
        for a,b in zip(p,np.roll(p,-1,axis=0)):
            length=np.linalg.norm(b-a)
            if length<1e-13:continue
            residual=abs(2*delta@((a+b)/2-x[i])-bound);residual[i]=np.inf;j=int(np.argmin(residual))
            if residual[j]>1e-10 or i>=j:continue
            topology.append((i,j))
            speed=(dw[i]-dw[j])/(2*np.linalg.norm(x[j]-x[i]))
            pts=np.array([a+u*(b-a) for u in (.5-.5/np.sqrt(3),.5+.5/np.sqrt(3))])
            swept=length*speed*np.mean(basis(pts),axis=0)
            r[i]+=swept;r[j]-=swept
    return r,topology


def integrate(x,w0,w1,split=True):
    times=event_times(x,w0,w1) if split else np.array([0.,1.]);dw=w1-w0
    total=np.zeros((len(x),6));topologies=[]
    for left,right in zip(times,times[1:]):
        mid=(left+right)/2;half=(right-left)/2
        topologies.append(rate(x,w0+mid*dw,dw)[1])
        for q in (-1/np.sqrt(3),1/np.sqrt(3)):
            total+=half*rate(x,w0+(mid+half*q)*dw,dw)[0]
    actual=moments(power_cells(x,w1,DOMAIN))-moments(power_cells(x,w0,DOMAIN))
    return dict(events=times.tolist(),topologies=topologies,
                area_error=float(np.max(abs(total[:,0]-actual[:,0]))),
                first_moment_error=float(np.max(abs(total[:,1:3]-actual[:,1:3]))),
                second_moment_error=float(np.max(abs(total[:,3:]-actual[:,3:]))),
                global_exchange_error=float(np.max(abs(total.sum(axis=0)))))


def probe():
    x=np.array([[.25,.25],[.75,.25],[.25,.75],[.75,.75]])
    results=[]
    for crossing in (.2,.5,.8):
        w0=np.array([.04*crossing,0.,0.,0.]);w1=np.array([-.04*(1-crossing),0.,0.,0.])
        results.append(dict(expected_event=crossing,split=integrate(x,w0,w1),unsplit=integrate(x,w0,w1,False)))
    return results

if __name__=='__main__':print(json.dumps(probe(),indent=2))
