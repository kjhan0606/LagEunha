"""Independent Python reference and normal-range bounds for live-stage audit."""
import numpy as np
from scipy.linalg import expm
from moment_factor_error_bounds import check,make_bounds
eps=np.finfo(float).eps

def factor_reference(base,linear,A,K,Q,offsets,h,raw,frozen):
    n=len(raw);out=np.zeros((n,63))
    for i in range(n):
        out[i,:4]=base[i].ravel();out[i,4:22]=raw[i]
        if frozen[i]:continue
        U,s,Vh=np.linalg.svd(linear[i]);scaled=s/s.max()
        ratio=scaled[:,None]*scaled[None,:]/(scaled[:,None]**2+scaled[None,:]**2)
        loss=U@((Vh@K[i]@Vh.T)*ratio)@U.T;Ac=A[i]-loss
        eh=expm(.5*h*Ac);ef=expm(h*Ac)
        columns=np.column_stack([ef@base[i],np.sqrt(h)*(eh@Q[offsets[i]:offsets[i+1]].T)])
        _,r=np.linalg.qr(columns.T,mode='reduced');newL=(np.where(np.diag(r)<0,-1.,1.)[:,None]*r).T
        R=newL@newL.T;m=raw[i,0];q=raw[i,1:3];p=raw[i,7:9]
        S=raw[i,3:7].reshape(2,2)-np.outer(q,q)/m;C=raw[i,9:13].reshape(2,2)-np.outer(q,p)/m
        B=np.linalg.solve(S,C).T;pi=np.outer(p,p)/m+B@S@B.T+R;pi[1,0]=pi[0,1]
        out[i,:4]=newL.ravel();out[i,17:21]=pi.ravel()
        d=out[i,22:62];d[:4]=loss.ravel();d[4:8]=Ac.ravel();d[8:12]=eh.ravel();d[12:16]=ef.ravel();d[16:20]=R.ravel();d[20:24]=B.ravel();d[24:28]=(pi-raw[i,13:17].reshape(2,2)).ravel()
        d[28]=raw[i,17]-.5*(raw[i,13]+raw[i,16]);d[29]=raw[i,17]-.5*np.trace(pi);d[30]=d[29]-d[28];d[39]=np.linalg.norm(columns@columns.T,np.inf)
    return out

def audit_factor(inputs,reference,outputs):
    base,linear,A,K,Q,offset,h,raw,frozen=inputs;n=len(raw)
    deck=np.column_stack([base.reshape(n,4),linear.reshape(n,4),A.reshape(n,4),K.reshape(n,4),raw,np.full(n,h),frozen])
    tol,metrics=make_bounds(deck,reference,Q,offset,n);result={};frozen=frozen.astype(bool)
    for name,out in outputs.items():
        assert np.isfinite(out).all() and not np.any(out[:,62])
        row={}
        for label,sl in [('L',slice(0,4)),('raw',slice(4,22)),('loss',slice(22,26)),('Ac',slice(26,30)),('Ehalf',slice(30,34)),('Efull',slice(34,38)),('R',slice(38,42)),('B',slice(42,46)),('deltaPi',slice(46,50)),('thermal',slice(50,53))]:
            row[label]=check(label,out[:,sl],reference[:,sl],tol[:,sl])
        assert np.array_equal(out[frozen,:4].view(np.uint64),base[frozen].reshape(-1,4).view(np.uint64))
        assert np.array_equal(out[frozen,4:22].view(np.uint64),raw[frozen].view(np.uint64))
        for k in list(range(13))+[17]:assert np.array_equal(out[:,4+k].view(np.uint64),raw[:,k].view(np.uint64))
        assert np.array_equal(out[:,18].view(np.uint64),out[:,19].view(np.uint64))
        assert np.all(out[~frozen,0]>0) and np.all(out[~frozen,3]>0) and not np.any(out[~frozen,1])
        assert np.all(out[~frozen,51]>0) and np.max(out[:,53])<128*eps
        for i in np.flatnonzero(~frozen):
            metric=metrics[i];L=out[i,:4].reshape(2,2);U=metric['U'];thin=np.linalg.solve(L,U).T
            assert np.linalg.norm(thin.T@thin-np.eye(2),np.inf)<64*metric['qr_gamma']*metric['condition_U']**2+64*metric['exp_bound']
            assert out[i,58]<64*metric['qr_gamma']*metric['condition_U']**2
            check('UUT',L@L.T,U@U.T,2*tol[i,38:42].reshape(2,2))
        row['maximum_QR_orthogonality']=float(out[:,58].max());row['minimum_thermal']=float(out[~frozen,51].min());result[name]=row
    return result

def stage_reference(data,captured,gravity,covariance_source,covariance_generator):
    raw=data['cells'][:,:18];pv=data['cells'][:,18];n=len(raw);c=raw[:,1:3]/raw[:,0,None]
    F=np.zeros((n,2));H=np.zeros((n,2,2));dr=np.zeros((n,18));work=np.zeros(n);K=np.zeros((n,2,2));Q=np.zeros((len(data['incidence']),16))
    a=data['node_inputs'];extra=data['extra'];off=data['offsets'];inc=data['incidence'];frozen=data['frozen'].astype(bool)
    for i in range(n):
        for t in range(off[i],off[i+1]):
            code=int(inc[t]);node=abs(code)-1;side=0 if code<0 else 1;sign=-1 if side==0 else 1
            v=a[node];dr[i]+=sign*captured['expected'][node];traction=v[4]*v[31]*v[2:4]
            F[i]+=sign*traction;H[i]+=sign*np.outer(v[:2]-c[i],traction);work[i]+=sign*v[4]*v[31]*v[32]
            xi=extra[node,8*side:8*side+8].reshape(4,2);loss=v[21:25] if side==0 else -v[25:29]
            K[i]+=v[4]/raw[i,0]*(xi.T@(loss[:,None]*xi))
            incoming=-v[4]*v[25:29] if side==0 else v[4]*v[21:25]
            velocities=v[13:21].reshape(4,2) if side==0 else v[5:13].reshape(4,2)
            qi=np.sqrt(incoming)[:,None]*(velocities-extra[node,16+2*side:18+2*side]);Q[t,:8]=qi.ravel()
            if not frozen[i]:Q[t,8:]=qi.ravel()
    src=np.zeros_like(dr);B=[];G=[];Ag=[];pressure=np.zeros_like(dr)
    for i,r in enumerate(raw):
        m=r[0];p=r[7:9];D=r[9:13].reshape(2,2);S=r[3:7].reshape(2,2)-np.outer(r[1:3],r[1:3])/m
        hh=H[i]+pv[i]*np.eye(2);A=np.linalg.solve(S,hh).T;a0=F[i]/m-A@c[i]
        src[i,1:3]=p;src[i,3:7]=(D+D.T).ravel();src[i,9:13]=(r[13:17].reshape(2,2)+pv[i]*np.eye(2)).ravel();src[i,13:17]=(np.outer(a0,p)+np.outer(p,a0)+A@D+D.T@A.T).ravel()
        B.append(np.linalg.solve(S,D-np.outer(r[1:3],p)/m).T)
        g=covariance_source(r,gravity[i]);G.append(g);Ag.append(covariance_generator(data['factors'][i],g))
        pressure[i,7:9]=F[i];pressure[i,9:13]=(hh+np.outer(c[i],F[i])).ravel();pressure[i,13:17]=src[i,13:17];pressure[i,17]=work[i]
    unmask=dr+src;masked=unmask.copy();masked[frozen]=0.;pressure[frozen]=0.;B=np.array(B);G=np.array(G);Ag=np.array(Ag)
    cov=-B+Ag;cov[frozen]=0.;km=K.copy();km[frozen]=0.
    total=masked.copy();total[~frozen]+=gravity[~frozen]
    return dict(raw=total,unmasked_hydro=unmask,hydro=masked,pressure=pressure,source=src,F=F,H=H,K=km,K_unmasked=K,Q=Q,gravity=gravity,B=B,G=G,Ag=Ag,cov_A=cov)

def audit_stage(data,reference,outputs):
    r=data['cells'][:,:18];n=len(r);a=data['node_inputs'];off=data['offsets'];inc=data['incidence'];fr=data['frozen'].astype(bool);gn=data['gravity_nodes'];goff=data['gravity_offsets']
    # Absolute primitive node bounds before cancellation.
    vl=a[:,5:13].reshape(-1,4,2);vr=a[:,13:21].reshape(-1,4,2);wl=abs(a[:,21:25]);wr=abs(a[:,25:29]);x=abs(a[:,:2]);fm=wl.sum(1)+wr.sum(1)
    fp=np.einsum('nk,nkj->nj',wl,abs(vl))+np.einsum('nk,nkj->nj',wr,abs(vr))+abs(a[:,31,None])*abs(a[:,2:4])
    pi=np.einsum('nk,nki,nkj->nij',wl,abs(vl),abs(vl))+np.einsum('nk,nki,nkj->nij',wr,abs(vr),abs(vr))
    energy=(wl*(.5*(vl*vl).sum(2)+abs(a[:,29,None]))).sum(1)+(wr*(.5*(vr*vr).sum(2)+abs(a[:,30,None]))).sum(1)+abs(a[:,31]*a[:,32])
    nb=np.column_stack([fm,x*fm[:,None],(x[:,:,None]*x[:,None,:]*fm[:,None,None]).reshape(-1,4),fp,(x[:,:,None]*fp[:,None,:]).reshape(-1,4),pi.reshape(-1,4),energy])*a[:,4,None]
    fb=np.zeros((n,18));fa=np.zeros(n);ha=np.zeros(n);wa=np.zeros(n);qbound=np.zeros_like(reference['Q']);ka=np.zeros((n,2,2));incoming=np.zeros((n,2,2));c=r[:,1:3]/r[:,0,None]
    for i in range(n):
        for t in range(off[i],off[i+1]):
            code=int(inc[t]);k=abs(code)-1;side=0 if code<0 else 1;v=a[k];fb[i]+=nb[k]
            traction=abs(v[4]*v[31]*v[2:4]);fa[i]+=np.linalg.norm(traction);ha[i]+=np.linalg.norm(v[:2]-c[i])*np.linalg.norm(traction)
            wa[i]+=abs(v[4]*v[31]*v[32])
            xi=data['extra'][k,8*side:8*side+8].reshape(4,2);loss=v[21:25] if side==0 else -v[25:29];ka[i]+=v[4]/r[i,0]*(abs(xi).T@(loss[:,None]*abs(xi)))
            incoming_rate=-v[4]*v[25:29] if side==0 else v[4]*v[21:25]
            vv=v[13:21].reshape(4,2) if side==0 else v[5:13].reshape(4,2)
            dq=32*eps*np.sqrt(incoming_rate)[:,None]*(abs(vv)+abs(data['extra'][k,16+2*side:18+2*side]))
            qbound[t,:8]=dq.ravel()
            if not fr[i]:qbound[t,8:]=dq.ravel()
            qq=abs(reference['Q'][t,:8].reshape(4,2));incoming[i]+=dq.T@qq+qq.T@dq+dq.T@dq+(64*eps+8*(off[i+1]-off[i])*eps)*(qq.T@qq)
    gc=np.diff(goff);primitive=np.zeros((len(gn),18));dm=gn[:,2];xx=abs(gn[:,:2]);v=abs(gn[:,3:5]);acc=abs(gn[:,5:7])
    primitive[:,7:9]=dm[:,None]*acc;primitive[:,9:13]=(dm[:,None,None]*xx[:,:,None]*acc[:,None,:]).reshape(-1,4)
    primitive[:,13:17]=(dm[:,None,None]*(v[:,:,None]*acc[:,None,:]+acc[:,:,None]*v[:,None,:])).reshape(-1,4);primitive[:,17]=dm*(v*acc).sum(1)
    ga=np.add.reduceat(primitive,goff[:-1],axis=0);gt=(64*eps+2*gc[:,None]*eps)*ga
    rt=(64*eps+2*np.diff(off)[:,None]*eps)*fb;kt=(64*eps+8*np.diff(off)[:,None,None]*eps)*ka
    covtol=np.zeros((n,2,2));btol=np.zeros((n,2,2));Gtol=np.zeros((n,2,2));Agtol=np.zeros((n,2,2));pt=np.zeros((n,18))
    for i,rr in enumerate(r):
        m=rr[0];q=rr[1:3];p=rr[7:9];D=rr[9:13].reshape(2,2);S=rr[3:7].reshape(2,2)-np.outer(q,q)/m;C=D-np.outer(q,p)/m
        sn=np.linalg.norm(S,2);inv=1/np.linalg.eigvalsh(S).min();B=reference['B'][i];bn=np.linalg.norm(B,2)
        ds=16*eps*np.linalg.norm(abs(rr[3:7].reshape(2,2))+abs(np.outer(q,q)/m));dc=16*eps*np.linalg.norm(abs(D)+abs(np.outer(q,p)/m));eta=inv*ds;assert eta<1
        db=inv*(dc+ds*bn+128*eps*(sn*bn+np.linalg.norm(C)))/(1-eta);btol[i]=db
        count=off[i+1]-off[i];dF=(32+2*count)*eps*fa[i]
        dH=(32+2*count)*eps*ha[i]+16*eps*(np.linalg.norm(c[i])+np.max(np.linalg.norm(a[np.abs(inc[off[i]:off[i+1]])-1,:2],axis=1)))*fa[i]
        H=reference['H'][i]+data['cells'][i,18]*np.eye(2);A=np.linalg.solve(S,H).T;an=np.linalg.norm(A);a0=reference['F'][i]/m-A@c[i]
        dA=inv*(dH+ds*an+128*eps*(sn*an+np.linalg.norm(H)))/(1-eta);da0=dF/m+dA*np.linalg.norm(c[i])+32*eps*(np.linalg.norm(reference['F'][i])/m+an*np.linalg.norm(c[i]))
        dpip=2*(da0*np.linalg.norm(p)+dA*np.linalg.norm(D))+128*eps*(np.linalg.norm(a0)*np.linalg.norm(p)+an*np.linalg.norm(D))
        pt[i,7:9]=dF;pt[i,9:13]=dH+np.linalg.norm(c[i])*dF+64*eps*(np.linalg.norm(H)+np.linalg.norm(c[i])*np.linalg.norm(reference['F'][i])+abs(data['cells'][i,18]));pt[i,13:17]=dpip;pt[i,17]=(32+2*count)*eps*wa[i]
        rt[i,1:3]+=16*eps*abs(p);rt[i,3:7]+=16*eps*(abs(D)+abs(D.T)).ravel();rt[i,9:13]+=16*eps*(abs(rr[13:17])+abs(data['cells'][i,18])*np.eye(2).ravel());rt[i,13:17]+=dpip
        g=reference['gravity'][i];gp=g[7:9];gC=g[9:13].reshape(2,2)-np.outer(q,gp)/m;gT=g[13:17].reshape(2,2)-(np.outer(gp,p)+np.outer(p,gp))/m
        dgC=np.linalg.norm(gt[i,9:13])+np.linalg.norm(q)*np.linalg.norm(gt[i,7:9])/m+32*eps*(np.linalg.norm(g[9:13])+np.linalg.norm(q)*np.linalg.norm(gp)/m)
        dgT=np.linalg.norm(gt[i,13:17])+2*np.linalg.norm(p)*np.linalg.norm(gt[i,7:9])/m+32*eps*(np.linalg.norm(g[13:17])+2*np.linalg.norm(p)*np.linalg.norm(gp)/m)
        dG=dgT+2*(db*np.linalg.norm(gC)+bn*dgC+db*dgC)+128*eps*(np.linalg.norm(gT)+2*bn*np.linalg.norm(gC));Gtol[i]=dG
        L=data['factors'][i];R=L@L.T;lmin=np.linalg.eigvalsh(R).min();assert lmin>0;Ag=reference['Ag'][i]
        dR=32*eps*np.linalg.norm(abs(L)@abs(L.T));assert dR<lmin
        dag=(dG+2*dR*np.linalg.norm(Ag)+128*eps*(2*np.linalg.norm(R)*np.linalg.norm(Ag)+np.linalg.norm(reference['G'][i])))/(2*(lmin-dR))
        Agtol[i]=dag;covtol[i]=db+dag+16*eps*(bn+np.linalg.norm(Ag))
    total_tol=rt+gt+16*eps*(abs(reference['hydro'])+abs(reference['gravity']));total_tol[fr]=0.;kt[fr]=0.;covtol[fr]=0.;pt[fr]=0.
    ht=rt.copy();ht[fr]=0.
    if 'original_hydro' in reference:check('original Python RHS replay',reference['hydro'],reference['original_hydro'],ht)
    result={}
    for name,out in outputs.items():
        st=out['stage'];qq=out['Q'];row={}
        for label,actual,target,tol in [('raw',st[:,:18],reference['raw'],total_tol),('gravity',st[:,30:48],reference['gravity'],gt),('B',st[:,18:22],reference['B'].reshape(-1,4),btol.reshape(-1,4)),('K',st[:,22:26],reference['K'].reshape(-1,4),kt.reshape(-1,4)),('Q',qq,reference['Q'],qbound),('G',st[:,48:52],reference['G'].reshape(-1,4),Gtol.reshape(-1,4)),('Ag',st[:,52:56],reference['Ag'].reshape(-1,4),Agtol.reshape(-1,4)),('cov_A',st[:,56:60],reference['cov_A'].reshape(-1,4),covtol.reshape(-1,4))]:row[label]=check(label,actual,target,tol)
        assert np.isfinite(st).all() and not np.any(st[:,93]);assert np.max(st[:,61])<128*eps
        assert not np.any(st[fr,:18]) and not np.any(st[fr,22:30]) and not np.any(st[fr,56:60])
        row['pressure']=check('pressure',out['hydro'][:,36:54],reference['pressure'],pt)
        for i in range(n):
            q=qq[off[i]:off[i+1],:8].reshape(-1,2);check('incoming',q.T@q,st[i,67:71].reshape(2,2),incoming[i]+64*eps*(abs(q).T@abs(q)))
            refq=reference['Q'][off[i]:off[i+1],8:].reshape(-1,2)
            tolerance=incoming[i]+64*eps*(abs(refq).T@abs(refq)) if not fr[i] else np.zeros((2,2))
            check('masked incoming',st[i,26:30].reshape(2,2),refq.T@refq,tolerance)
        total=np.r_[st[:,7:9].sum(0),st[:,17].sum(),(st[:,10]-st[:,11]).sum()]
        primitive_total=fb.sum(0)+ga.sum(0);bound=(128+2*n)*eps*np.r_[primitive_total[7:9],primitive_total[17],primitive_total[10]+primitive_total[11]]
        bound+=np.r_[total_tol[:,7:9].sum(0),total_tol[:,17].sum(),(total_tol[:,10]+total_tol[:,11]).sum()]
        row['global_ledger']=check('stage ledger',total,out['ledgers'].reshape(3,4).sum(0),bound);row['global_ledger_bound']=bound.tolist()
        row['maximum_sylvester_backward']=float(st[:,61].max());result[name]=row
    return result
