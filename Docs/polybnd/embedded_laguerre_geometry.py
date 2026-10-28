"""Moving Laguerre patch embedded in fixed exterior generator positions.

Global regular triangulation supplies clipping neighbors; local face records
are rebuilt after every motion/weight change. Research 2-D geometry only.
"""
import numpy as np
from scipy.spatial import ConvexHull


class EmbeddedLaguerreGeometry:
    def __init__(self, positions, ids, clip, bounds=(-3.,3.,-3.,3.)):
        self.positions=np.asarray(positions,dtype=float).copy()
        self.ids=np.asarray(ids,dtype=int).copy()
        if self.positions.ndim!=2 or self.positions.shape[1]!=2 or not np.isfinite(self.positions).all():
            raise ValueError('finite 2-D positions required')
        if len(np.unique(self.ids))!=len(self.ids) or np.any(self.ids<0) or np.any(self.ids>=len(self.positions)):
            raise ValueError('unique valid patch ids required')
        self.lookup={int(k):i for i,k in enumerate(self.ids)}
        self.clip=clip;self.bounds=bounds

    def __call__(self,x,w,t):
        x=np.asarray(x,dtype=float);w=np.asarray(w,dtype=float)
        if x.shape!=(len(self.ids),2) or w.shape!=(len(self.ids),) or not np.isfinite(x).all() or not np.isfinite(w).all():
            raise ValueError('finite patch positions and weights required')
        positions=self.positions.copy();positions[self.ids]=x
        weights=np.zeros(len(positions));weights[self.ids]=w
        lifted=np.column_stack([positions,np.sum(positions**2,axis=1)-weights])
        hull=ConvexHull(lifted)
        tri=hull.simplices[hull.equations[:,2]<-1e-12]
        pairs=np.unique(np.sort(np.concatenate([tri[:,[0,1]],tri[:,[0,2]],tri[:,[1,2]]]),axis=1),axis=0)
        neighbors=[[] for _ in positions]
        for i,j in pairs:neighbors[i].append(j);neighbors[j].append(i)
        lo,hi,bottom,top=self.bounds
        box=np.array([[lo,bottom],[hi,bottom],[hi,top],[lo,top]])
        polys=[];edges=[]
        for local,global_id in enumerate(self.ids):
            nb=np.array(sorted(neighbors[global_id]),dtype=int)
            if not len(nb):raise RuntimeError('hidden/empty Laguerre cell')
            delta=positions[nb]-x[local];distance=np.linalg.norm(delta,axis=1)
            if np.any(distance==0):raise RuntimeError('coincident generators')
            rhs=np.sum(delta**2,axis=1)+w[local]-weights[nb]
            p=box-x[local]
            for d,b in zip(delta,rhs):
                p=self.clip(p,2*d,b)
                if len(p)<3:raise RuntimeError('empty Laguerre cell')
            signed=.5*np.sum(p[:,0]*np.roll(p[:,1],-1)-p[:,1]*np.roll(p[:,0],-1))
            if signed<=0:raise RuntimeError('nonpositive cell area')
            polys.append(p+x[local])
            for a,b in zip(p,np.roll(p,-1,axis=0)):
                length=float(np.linalg.norm(b-a))
                if length==0:continue
                residual=abs(2*(delta@((a+b)/2))-rhs)/(2*distance)
                k=int(np.argmin(residual))
                if residual[k]>1e-10:
                    raise RuntimeError('patch touches bounding box or unmatched face')
                j=self.lookup.get(int(nb[k]),-1)
                if j>=0 and local>=j:continue
                normal=delta[k]/distance[k]
                if max(abs(2*np.dot(delta[k],a)-rhs[k]),abs(2*np.dot(delta[k],b)-rhs[k]))/(2*distance[k])>1e-10:
                    raise RuntimeError('face endpoints do not share power plane')
                edges.append((local,j,a+x[local],b+x[local],normal,length))
        return polys,edges
