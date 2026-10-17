"""Exact positive-part integrals for affine normal transport in constant bins."""
import numpy as np


def positive_integral(a,b):
    """Mean max(linear(a,b),0), and positive fraction/subgradient."""
    a=np.asarray(a);b=np.asarray(b);value=np.zeros_like(a);fraction=np.zeros_like(a)
    both=(a>0)&(b>0);value[both]=(a[both]+b[both])/2;fraction[both]=1.
    start=(a>0)&(b<=0);end=(b>0)&(a<=0)
    value[start]=a[start]**2/(2*(a[start]-b[start]));fraction[start]=a[start]/(a[start]-b[start])
    value[end]=b[end]**2/(2*(b[end]-a[end]));fraction[end]=b[end]/(b[end]-a[end])
    fraction[(a==0)&(b==0)]=.5
    return value,fraction


def arrays(bins):
    return {k:np.asarray([r[k] for r in bins]) for k in ('i','j','scale','distance','sl0','sl1','sr0','sr1','rl','rr')}


def evaluate(data,count,z):
    i=data['i'];j=data['j'];internal=j>=0;shift=np.zeros(len(i))
    shift[internal]=(z[i[internal]]-z[j[internal]])/(2*data['distance'][internal])
    l0=data['sl0']-shift[:,None];l1=data['sl1']-shift[:,None]
    r0=data['sr0']-shift[:,None];r1=data['sr1']-shift[:,None]
    lp,lf=positive_integral(l0,l1);rp,rf=positive_integral(r0,r1)
    negative=(r0+r1)/2-rp
    flux=data['scale']*(data['rl']*.25*lp.sum(1)+data['rr']*.25*negative.sum(1))
    rho=data['rl']*.25*lf.sum(1)+data['rr']*.25*(1-rf).sum(1)
    residual=np.zeros(count);np.add.at(residual,i,flux);np.add.at(residual,j[internal],-flux[internal])
    lap=np.zeros((count,count));ii=i[internal];jj=j[internal]
    c=data['scale'][internal]*rho[internal]/(2*data['distance'][internal])
    np.add.at(lap,(ii,ii),c);np.add.at(lap,(jj,jj),c);np.add.at(lap,(ii,jj),-c);np.add.at(lap,(jj,ii),-c)
    return residual,lap


def solve(bins,count):
    data=arrays(bins);z=np.zeros(count)
    for iteration in range(50):
        residual,lap=evaluate(data,count,z);norm=float(max(abs(residual)))
        if norm<2e-13:return z,norm,iteration+1
        step=np.r_[np.linalg.lstsq(lap[:-1,:-1],residual[:-1],rcond=1e-13)[0],0.]
        for bt in range(30):
            trial=z+step*2.**(-bt)
            if max(abs(evaluate(data,count,trial)[0]))<norm:z=trial;break
        else:raise RuntimeError('exact-bin kinetic mass line search failed')
    raise RuntimeError('exact-bin kinetic mass solve failed')


def cuts(record,z):
    shift=(z[record['i']]-z[record['j']])/(2*record['distance']) if record['j']>=0 else 0.
    values=[0.,1.]
    for a,b in zip(np.r_[record['sl0'],record['sr0']]-shift,np.r_[record['sl1'],record['sr1']]-shift):
        if a*b<0:values.append(float(a/(a-b)))
    return np.unique(values)


def make_bins(edges,x,states,factors,mass,gen,trace,segments,boundary_velocity,velocity,alpha):
    bins=[]
    for edge_id,(i,j,a,b,n,length) in enumerate(edges):
        d=np.linalg.norm(x[j]-x[i]) if j>=0 else 0.
        def endpoint(g):
            point=a+g*(b-a);vl=velocity(states[i],point);vls=vl+alpha@factors[i].T/np.sqrt(mass[i])
            if j>=0:
                vr=velocity(states[j],point);vrs=vr+alpha@factors[j].T/np.sqrt(mass[j])
                wn=((x[j]-x[i])@gen[j]-(gen[j]-gen[i])@(point-x[i]))/d
            else:
                wn=boundary_velocity(point)@n;vrs=vls-2*((vls@n-wn)[:,None]*n)
            return vls@n-wn,vrs@n-wn
        for k in range(segments):
            sl0,sr0=endpoint(k/segments);sl1,sr1=endpoint((k+1)/segments)
            mid=(k+.5)/segments;rl=trace(i,edge_id,mid);rr=trace(j,edge_id,mid) if j>=0 else rl
            bins.append(dict(i=i,j=j,scale=length/segments,distance=d,sl0=sl0,sl1=sl1,sr0=sr0,sr1=sr1,rl=rl,rr=rr))
    return bins


def make_bin_arrays(edges,x,states,factors,mass,gen,trace,segments,boundary_velocity,alpha):
    """Vectorized equivalent for affine velocity/boundary fields."""
    packed={k:[] for k in ('i','j','scale','distance','sl0','sl1','sr0','sr1','rl','rr')}
    parameters=np.arange(segments+1)/segments
    for edge_id,(i,j,a,b,n,length) in enumerate(edges):
        points=a+parameters[:,None]*(b-a)
        sl=states[i];vl=sl['v']+(points-sl['c'])@sl['B'].T
        vls=vl[:,None,:]+(alpha@factors[i].T/np.sqrt(mass[i]))[None,:,:]
        if j>=0:
            d=np.linalg.norm(x[j]-x[i]);sr=states[j];vr=sr['v']+(points-sr['c'])@sr['B'].T
            vrs=vr[:,None,:]+(alpha@factors[j].T/np.sqrt(mass[j]))[None,:,:]
            wn=((x[j]-x[i])@gen[j]-(points-x[i])@(gen[j]-gen[i]))/d
        else:
            d=0.;wa=boundary_velocity(a)@n;wb=boundary_velocity(b)@n;wn=wa+parameters*(wb-wa)
            vrs=vls-2*((vls@n-wn[:,None])[:,:,None]*n)
        left=vls@n-wn[:,None];right=vrs@n-wn[:,None]
        values=dict(i=np.full(segments,i),j=np.full(segments,j),scale=np.full(segments,length/segments),distance=np.full(segments,d),
            sl0=left[:-1],sl1=left[1:],sr0=right[:-1],sr1=right[1:],
            rl=np.array([trace(i,edge_id,(k+.5)/segments) for k in range(segments)]),
            rr=np.array([trace(j if j>=0 else i,edge_id,(k+.5)/segments) for k in range(segments)]))
        for key,value in values.items():packed[key].append(value)
    return {key:np.concatenate(value) for key,value in packed.items()}
