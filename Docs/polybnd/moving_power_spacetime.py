"""Geometry-only gate: linear physical generator trajectories and weight paths.
Fixed unit box. Polynomial topology events; shared face moment integration.
This is not a fluid evolution or an assertion of fixed cell mass.
"""
import itertools,json
import numpy as np
from numpy.polynomial import Polynomial as Poly
from weighted_mass_patch import power_cells
from power_spacetime_geometry import DOMAIN,basis,moments


def constraints(x,v,w,dw,i):
    rows=[]
    for j in range(len(x)):
        if j==i:continue
        A=[Poly([2*(x[j,k]-x[i,k]),2*(v[j,k]-v[i,k])]) for k in range(2)]
        b=Poly([x[j]@x[j]-x[i]@x[i]+w[i]-w[j],2*(x[j]@v[j]-x[i]@v[i])+dw[i]-dw[j],v[j]@v[j]-v[i]@v[i]])
        rows.append(A+[b])
    rows.extend([[Poly([a]),Poly([b]),Poly([c])] for a,b,c in [(-1,0,0),(1,0,1),(0,-1,0),(0,1,1)]])
    return rows


def roots(p):
    p=p.trim(tol=1e-13)
    if p.degree()==0:return []
    return [float(z.real) for z in p.roots() if abs(z.imag)<1e-9 and 1e-11<z.real<1-1e-11]


def events(x,v,w,dw,constraint_builder=None):
    candidates=[0.,1.];singular_splits=0
    for i in range(len(x)):
        rows=(constraints if constraint_builder is None else constraint_builder)(x,v,w,dw,i)
        for j,k in itertools.combinations(range(len(rows)),2):
            a,b=rows[j],rows[k]
            rr=roots(a[0]*b[1]-a[1]*b[0]);candidates.extend(rr);singular_splits+=len(rr)
        for ids in itertools.combinations(range(len(rows)),3):
            a,b,c=[rows[j] for j in ids]
            det=a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0])
            for t in roots(det):
                values=np.array([[p(t) for p in row] for row in rows]);A=values[:,:2];bounds=values[:,2]
                for j,k in itertools.combinations(ids,2):
                    if abs(np.linalg.det(A[[j,k]]))<1e-11:continue
                    point=np.linalg.solve(A[[j,k]],bounds[[j,k]])
                    if np.max(A@point-bounds)<1e-9:
                        candidates.append(t);break
    ordered=sorted(candidates);out=[ordered[0]]
    for t in ordered[1:]:
        if t-out[-1]>1e-9:out.append(t)
    return np.array(out),singular_splits


def rate(x,v,w,dw):
    cells=power_cells(x,w,DOMAIN);r=np.zeros((len(x),6));topology=[]
    for i,cell in enumerate(cells):
        if cell.area<=0:raise RuntimeError('empty cell')
        p=np.asarray(cell.exterior.coords[:-1]);delta=x-x[i];bound=np.sum(delta*delta,axis=1)+w[i]-w
        for a,b in zip(p,np.roll(p,-1,axis=0)):
            length=np.linalg.norm(b-a)
            if length<1e-13:continue
            residual=abs(2*delta@((a+b)/2-x[i])-bound);residual[i]=np.inf;j=int(np.argmin(residual))
            if residual[j]>1e-9 or i>=j:continue
            topology.append((i,j));distance=np.linalg.norm(x[j]-x[i])
            pts=np.array([a+s*(b-a) for s in (.5-.5/np.sqrt(3),.5+.5/np.sqrt(3))])
            speed=(x[j]@v[j]-x[i]@v[i]-pts@(v[j]-v[i])+.5*(dw[i]-dw[j]))/distance
            exchange=length/2*(speed@basis(pts));r[i]+=exchange;r[j]-=exchange
    return r,sorted(topology)


def integrate(x,v,w,dw,order=8,split=True):
    ts,singular=events(x,v,w,dw) if split else (np.array([0.,1.]),0)
    q,weights=np.polynomial.legendre.leggauss(order);total=np.zeros((len(x),6));tops=[]
    for left,right in zip(ts,ts[1:]):
        mid=(left+right)/2;half=(right-left)/2
        tops.append(rate(x+mid*v,v,w+mid*dw,dw)[1])
        for z,a in zip(q,weights):
            t=mid+half*z;total+=half*a*rate(x+t*v,v,w+t*dw,dw)[0]
    actual=moments(power_cells(x+v,w+dw,DOMAIN))-moments(power_cells(x,w,DOMAIN))
    error=np.max(abs(total-actual),axis=0)
    return dict(order=order,split=split,events=ts.tolist(),topologies=tops,singular_splits=singular,
                area_error=float(error[0]),first_moment_error=float(max(error[1:3])),
                second_moment_error=float(max(error[3:])),global_exchange_error=float(np.max(abs(total.sum(0)))))


def probe():
    x=np.array([[.25,.25],[.75,.25],[.25,.75],[.75,.75]])
    result=[]
    for duration in (.02,.1,.3):
        # Exact generator displacement in v_phys=(.2+1.1*y,0): no steering.
        v=duration*np.column_stack([.2+1.1*x[:,1],np.zeros(4)])
        w=np.array([.02,0.,0.,0.]);dw=np.array([-.04,0.,0.,0.])
        for order,split in [(2,True),(4,True),(8,True),(16,True),(16,False)]:
            result.append(dict(duration=duration,**integrate(x,v,w,dw,order,split)))
    return result

if __name__=='__main__':print(json.dumps(probe(),indent=2))
