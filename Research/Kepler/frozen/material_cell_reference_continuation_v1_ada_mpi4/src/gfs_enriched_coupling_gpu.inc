/* Actual native CUDA polygons + canonical shared-face Euler moment RHS.
 * The only nonlinear solve here is the stage weight RATE. No weights,
 * particles, evolved moments or checkpoints are advanced by this audit.
 */
#include <map>
#include <set>
#include <vector>
#include <climits>
#include "gfs_ale_fan.h"

struct GfsCoupledPair {
    int i,j,edge,gi,gj,count,offset,status;
    double x1,y1,x2,y2,nx,ny,inv2d,t0,first[2];int ns[2];
};
struct GfsCoupledJac {int i,j;double a;};
struct GfsCoupledVelocity {double x,y;};

__global__ static void gfs_ec_velocity(const GfsEnrichAudit *a,int n,
    const GfsEnrichConfig *cfg,GfsCoupledVelocity *v,unsigned *bad){
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;
    double u[4];
    if(!gfs_enrich_evaluate(&a[i].basis,cfg,&a[i].state,0,0,u) ||
       !(u[0]>0) || !(gfs_dg_pressure(u,cfg->field.ic.gamma)>0)){atomicAdd(bad,1u);return;}
    v[i].x=u[1]/u[0];v[i].y=u[2]/u[0];
}

__global__ static void gfs_ec_pair_geometry(ParticleSoA p,const PolygonOut *poly,
    const GfsEnrichAudit *a,const GfsEnrichConfig *cfg,GfsCoupledPair *pairs,int m){
    int k=blockIdx.x*blockDim.x+threadIdx.x;if(k>=m)return;
    GfsCoupledPair q=pairs[k];q.status=-1;q.count=0;
    double x[GPU_MAX_CORNERS],y[GPU_MAX_CORNERS];
    if(!gfs_moment_support_polygon(p,poly+q.i,q.i,x,y)){pairs[k]=q;return;}
    int e=q.edge,j=(e+1)%poly[q.i].n_corners;
    q.x1=x[e];q.y1=y[e];q.x2=x[j];q.y2=y[j];
    double ex=q.x2-q.x1,ey=q.y2-q.y1,len=hypot(ex,ey);
    if(len==0){q.status=1;pairs[k]=q;return;}
    double dx=p.x[q.j]-p.x[q.i],dy=p.y[q.j]-p.y[q.i],d=hypot(dx,dy);
    if(!(d>0) || !(ey*dx-ex*dy>0)){q.status=-2;pairs[k]=q;return;}
    q.nx=dx/d;q.ny=dy/d;q.inv2d=.5/d;
    double ox=cfg->field.ic.cx-a[q.gi].basis.generator[0]-q.x1;
    double oy=cfg->field.ic.cy-a[q.gi].basis.generator[1]-q.y1;
    double t=fmax(0.,fmin(1.,(ox*ex+oy*ey)/(len*len)));
    double distance=hypot(t*ex-ox,t*ey-oy);
    if(distance>=len){q.t0=0;q.ns[0]=0;q.ns[1]=1;q.first[0]=0;q.first[1]=1;}
    else{
        q.t0=t;
        for(int side=0;side<2;side++){
            double delta=fabs((side?1.:0.)-t);q.ns[side]=0;q.first[side]=0;if(delta==0)continue;
            q.first[side]=fmin(1.,(distance+1.e-6)/(len*delta));q.ns[side]=1;
            while(ldexp(q.first[side],2*(q.ns[side]-1))<1 && q.ns[side]<GFS_ENRICH_SEGMENTS)q.ns[side]++;
            if(ldexp(q.first[side],2*(q.ns[side]-1))<1){q.status=-3;pairs[k]=q;return;}
        }
    }
    q.count=cfg->na*(q.ns[0]+q.ns[1]);q.status=1;pairs[k]=q;
}

__device__ static void gfs_ec_line(const GfsCoupledPair &p,const GfsEnrichConfig *cfg,
    int q,double *x,double *y,double *weight){
    int node=q%cfg->na,seg=q/cfg->na,side=seg<p.ns[0]?0:1;if(side)seg-=p.ns[0];
    double lo=seg?ldexp(p.first[side],2*(seg-1)):0,hi=fmin(1.,ldexp(p.first[side],2*seg));
    double delta=(side?1.:0.)-p.t0,t=p.t0+delta*(lo+(hi-lo)*cfg->an[node]);
    *x=p.x1+t*(p.x2-p.x1);*y=p.y1+t*(p.y2-p.y1);
    *weight=hypot(p.x2-p.x1,p.y2-p.y1)*fabs(delta)*(hi-lo)*cfg->aw[node];
}

__global__ static void gfs_ec_points(ParticleSoA p,const GfsEnrichAudit *a,
    const GfsEnrichConfig *cfg,const GfsCoupledVelocity *velocity,
    const GfsCoupledPair *pairs,int m,GfsWeightPoint *points,unsigned *bad,unsigned *hll){
    int k=blockIdx.x*blockDim.x+threadIdx.x;if(k>=m)return;
    GfsCoupledPair edge=pairs[k];
    double dx=p.x[edge.j]-p.x[edge.i],dy=p.y[edge.j]-p.y[edge.i];
    GfsCoupledVelocity vi=velocity[edge.gi],vj=velocity[edge.gj];
    for(int q=0;q<edge.count;q++){
        GfsWeightPoint z={};z.i=edge.gi;z.j=edge.gj;z.inv2d=edge.inv2d;
        gfs_ec_line(edge,cfg,q,&z.xi,&z.yi,&z.weight);z.xj=z.xi-dx;z.yj=z.yi-dy;
        double u[4],v[4];
        if(!gfs_enrich_evaluate(&a[z.i].basis,cfg,&a[z.i].state,z.xi,z.yi,u) ||
           !gfs_enrich_evaluate(&a[z.j].basis,cfg,&a[z.j].state,z.xj,z.yj,v) || !(u[0]>0) || !(v[0]>0)){
            atomicAdd(bad,1u);points[edge.offset+q]=z;continue;
        }
        GfsAlePrimitive l={u[0],gfs_dg_pressure(u,cfg->field.ic.gamma),u[1]/u[0],u[2]/u[0]};
        GfsAlePrimitive r={v[0],gfs_dg_pressure(v,cfg->field.ic.gamma),v[1]/v[0],v[2]/v[0]};
        z.fan=gfs_ale_fan_build(l,r,edge.nx,edge.ny,cfg->field.ic.gamma);
        z.w0=vj.x*edge.nx+vj.y*edge.ny-2*edge.inv2d*((vj.x-vi.x)*z.xi+(vj.y-vi.y)*z.yi);
        if(!z.fan.valid || !(z.weight>0))atomicAdd(bad,1u);if(z.fan.hll)atomicAdd(hll,1u);
        points[edge.offset+q]=z;
    }
}

__global__ static void gfs_ec_residual(const GfsWeightPoint *points,int m,const double *z,
    double *R,double *diag,GfsCoupledJac *jac){
    int k=blockIdx.x*blockDim.x+threadIdx.x;if(k>=m)return;const GfsWeightPoint &p=points[k];
    GfsAleFanSample f=gfs_weight_flux(&p,z[p.i],z[p.j]);double flux=p.weight*f.flux[0];
    double a=p.weight*p.inv2d*f.density;jac[k]={p.i,p.j,a};
    atomicAdd(R+p.i,-flux);atomicAdd(R+p.j,flux);atomicAdd(diag+p.i,a);atomicAdd(diag+p.j,a);
}

__global__ static void gfs_ec_jvp(const GfsCoupledJac *a,int m,const double *v,double *out){
    int k=blockIdx.x*blockDim.x+threadIdx.x;if(k>=m)return;
    GfsCoupledJac q=a[k];double t=q.a*(v[q.i]-v[q.j]);atomicAdd(out+q.i,t);atomicAdd(out+q.j,-t);
}

__global__ static void gfs_ec_volume(ParticleSoA p,const PolygonOut *polygon,
    const int *global_index,int np,const GfsEnrichAudit *state,const GfsEnrichConfig *cfg,
    double *rhs,unsigned *bad){
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=np)return;int gi=global_index[i];
    const GfsEnrichAudit &a=state[gi];double x[GPU_MAX_CORNERS],y[GPU_MAX_CORNERS];
    if(!gfs_moment_support_polygon(p,polygon+i,i,x,y)){atomicAdd(bad,1u);return;}
    GfsMomentPolygon poly={x,y,polygon[i].n_corners,0};GfsEnrichQuadrature rule;
    if(!gfs_enrich_quadrature_init(poly,cfg,&a.basis,&rule)){atomicAdd(bad,1u);return;}
    double rate[4][GFS_ENRICH_N]={};
    for(int q=0;q<gfs_enrich_nquad(poly,cfg,&rule);q++){
        double xx,yy,w;if(!gfs_enrich_quad(poly,cfg,&rule,q,&xx,&yy,&w)){atomicAdd(bad,1u);return;}
        if(w==0)continue;double gx=0,gy=0;
        if(cfg->field.ic.kind==2){
            double X=a.basis.generator[0]+xx-cfg->field.ic.cx,Y=a.basis.generator[1]+yy-cfg->field.ic.cy;
            double inverse=pow(X*X+Y*Y+cfg->field.ic.eps*cfg->field.ic.eps,-1.5);gx=-X*inverse;gy=-Y*inverse;
        }
        if(!gfs_enrich_volume_point(&a.basis,cfg,&a.state,xx,yy,gx,gy,w,rate)){atomicAdd(bad,1u);return;}
    }
    for(int c=0;c<4;c++)for(int k=0;k<GFS_ENRICH_N;k++)rhs[gi*40+c*10+k]=rate[c][k];
}

__device__ static void gfs_ec_state_from_q(const GfsEnrichAudit &a,const double *q,double u[4]){
    for(int c=0;c<4;c++){
        u[c]=a.state.c[c][0];for(int k=1;k<GFS_ENRICH_N;k++)u[c]=fma(a.state.c[c][k],q[k]-a.basis.mean[k],u[c]);
    }
}

__global__ static void gfs_ec_face_rhs(const GfsCoupledPair *pair,int m,const GfsWeightPoint *point,
    const double *z,const GfsEnrichAudit *a,const GfsEnrichConfig *cfg,
    double *rhs,double *geometry,double *pair_out,unsigned *bad){
    int k=blockIdx.x*blockDim.x+threadIdx.x;if(k>=m)return;const GfsCoupledPair &edge=pair[k];
    double ri[4][10]={},rj[4][10]={},gi[4][10]={},gj[4][10]={},ledger[10]={};
    for(int q=0;q<edge.count;q++){
        const GfsWeightPoint &p=point[edge.offset+q];GfsAleFanSample f=gfs_weight_flux(&p,z[p.i],z[p.j]);
        double qi[10],qj[10],ui[4],uj[4],wn=p.w0+p.inv2d*(z[p.i]-z[p.j]);
        if(!f.valid || !gfs_enrich_q(&a[p.i].basis,cfg,10,p.xi,p.yi,qi,0,0) ||
           !gfs_enrich_q(&a[p.j].basis,cfg,10,p.xj,p.yj,qj,0,0)){atomicAdd(bad,1u);return;}
        gfs_ec_state_from_q(a[p.i],qi,ui);gfs_ec_state_from_q(a[p.j],qj,uj);
        for(int c=0;c<4;c++){
            ledger[c]+=p.weight*f.flux[c];
            for(int j=0;j<10;j++){
                double ti=j?qi[j]-a[p.i].basis.mean[j]:1,tj=j?qj[j]-a[p.j].basis.mean[j]:1;
                ri[c][j]-=p.weight*f.flux[c]*ti;rj[c][j]+=p.weight*f.flux[c]*tj;
                gi[c][j]+=p.weight*wn*ui[c]*ti;gj[c][j]-=p.weight*wn*uj[c]*tj;
            }
        }
        ledger[4]+=p.weight*wn;ledger[5]+=p.weight*p.inv2d*f.density;
        ledger[6]+=p.weight*gfs_ale_fan_sample(&p.fan,p.w0).flux[0];ledger[7]+=p.weight;
        double Xi=a[p.i].basis.generator[0]+p.xi-cfg->field.ic.cx,Yi=a[p.i].basis.generator[1]+p.yi-cfg->field.ic.cy;
        double Xj=a[p.j].basis.generator[0]+p.xj-cfg->field.ic.cx,Yj=a[p.j].basis.generator[1]+p.yj-cfg->field.ic.cy;
        ledger[8]+=p.weight*(Xi*f.flux[2]-Yi*f.flux[1]);ledger[9]+=p.weight*(Xj*f.flux[2]-Yj*f.flux[1]);
    }
    for(int c=0;c<4;c++)for(int j=0;j<10;j++){
        atomicAdd(rhs+edge.gi*40+c*10+j,ri[c][j]);atomicAdd(rhs+edge.gj*40+c*10+j,rj[c][j]);
        atomicAdd(geometry+edge.gi*40+c*10+j,gi[c][j]);atomicAdd(geometry+edge.gj*40+c*10+j,gj[c][j]);
    }
    for(int j=0;j<10;j++)pair_out[k*10+j]=ledger[j];
}

template<class T> static T *gfs_ec_alloc(size_t n){T *p=NULL;if(n)CUDA_CHECK(cudaMalloc(&p,n*sizeof(T)));return p;}
static double gfs_ec_dot(const double *a,const double *b,int n){long double s=0;for(int i=0;i<n;i++)s+=(long double)a[i]*b[i];return (double)s;}

struct GfsCoupledWorkspace {
    int n,m;cudaStream_t stream;GfsEnrichSum sum;void *user;
    GfsWeightPoint *points;GfsCoupledJac *jac;double *z,*v,*R,*diag;
    std::vector<double> mu;
    void residual(const double *x,double *r,double *d){
        CUDA_CHECK(cudaMemcpyAsync(z,x,n*sizeof(double),cudaMemcpyHostToDevice,stream));
        CUDA_CHECK(cudaMemsetAsync(R,0,n*sizeof(double),stream));CUDA_CHECK(cudaMemsetAsync(diag,0,n*sizeof(double),stream));
        if(m)gfs_ec_residual<<<(m+255)/256,256,0,stream>>>(points,m,z,R,diag,jac);
        CUDA_CHECK(cudaGetLastError());CUDA_CHECK(cudaMemcpyAsync(r,R,n*sizeof(double),cudaMemcpyDeviceToHost,stream));
        CUDA_CHECK(cudaMemcpyAsync(d,diag,n*sizeof(double),cudaMemcpyDeviceToHost,stream));CUDA_CHECK(cudaStreamSynchronize(stream));
        sum(user,r,n);sum(user,d,n);for(int i=0;i<n;i++){r[i]+=mu[i]*x[i];d[i]+=mu[i];}
    }
    void product(const double *x,double *out){
        CUDA_CHECK(cudaMemcpyAsync(v,x,n*sizeof(double),cudaMemcpyHostToDevice,stream));CUDA_CHECK(cudaMemsetAsync(R,0,n*sizeof(double),stream));
        if(m)gfs_ec_jvp<<<(m+255)/256,256,0,stream>>>(jac,m,v,R);
        CUDA_CHECK(cudaGetLastError());CUDA_CHECK(cudaMemcpyAsync(out,R,n*sizeof(double),cudaMemcpyDeviceToHost,stream));CUDA_CHECK(cudaStreamSynchronize(stream));
        sum(user,out,n);for(int i=0;i<n;i++)out[i]+=mu[i]*x[i];
    }
    double residual_norm(const double *r){double q=0;for(int i=0;i<n;i++)q=fmax(q,fabs(r[i])/mu[i]);return q;}
    double merit(const double *r){long double q=0;for(int i=0;i<n;i++)q+=(long double)r[i]*r[i]/mu[i];return (double)q;}
    int cg(const double *b,const double *d,double *solution){
        std::vector<double> r(b,b+n),s(n),p(n),Ap(n);for(int i=0;i<n;i++){solution[i]=0;s[i]=p[i]=r[i]/d[i];}
        double rz=gfs_ec_dot(r.data(),s.data(),n),initial=gfs_ec_dot(r.data(),r.data(),n);
        if(initial==0)return 0;
        for(int iteration=1;iteration<=300;iteration++){
            product(p.data(),Ap.data());double den=gfs_ec_dot(p.data(),Ap.data(),n);
            if(!(den>0) || !isfinite(den))return -1;double alpha=rz/den;
            for(int i=0;i<n;i++){solution[i]+=alpha*p[i];r[i]-=alpha*Ap[i];}
            if(gfs_ec_dot(r.data(),r.data(),n)<=1.e-24*initial)return iteration;
            for(int i=0;i<n;i++)s[i]=r[i]/d[i];double next=gfs_ec_dot(r.data(),s.data(),n),beta=next/rz;
            for(int i=0;i<n;i++)p[i]=s[i]+beta*p[i];rz=next;
        }return -1;
    }
};

extern "C" int gpu_enriched_coupled_rhs(GPUContext *ctx,int np,int n,int rank,
    const unsigned long long *ids,const GfsEnrichAudit *initial,const GfsEnrichConfig *cfg,
    const double perturb[2],GfsEnrichSum sum,void *user,double *wdot,double *rhs,double *geometry,
    double report[EC_NREPORT]){
    if(n<=0 || np<0)return 1;int nt=ctx->d_parts.n_total;
    std::map<unsigned long long,int> index;for(int i=0;i<n;i++)if(!index.insert({ids[i],i}).second)return 1;
    std::vector<int> local_index(np);std::vector<GfsEnrichAudit> state(initial,initial+n);
    for(int i=0;i<n;i++)for(int k=0;k<10;k++){
        double rho=state[i].state.c[0][k],px=state[i].state.c[1][k],py=state[i].state.c[2][k];
        state[i].state.c[1][k]+=perturb[0]*rho;state[i].state.c[2][k]+=perturb[1]*rho;
        state[i].state.c[3][k]+=perturb[0]*px+perturb[1]*py+.5*(perturb[0]*perturb[0]+perturb[1]*perturb[1])*rho;
    }
    CUDA_CHECK(cudaMemcpy(ctx->h_parts.indx,ctx->d_parts.indx,nt*sizeof(long long),cudaMemcpyDeviceToHost));
    std::vector<PolygonOut> polygons(np);if(np)CUDA_CHECK(cudaMemcpy(polygons.data(),ctx->d_polygons,np*sizeof(PolygonOut),cudaMemcpyDeviceToHost));
    std::vector<GfsCoupledPair> pairs;std::set<std::pair<unsigned long long,unsigned long long> > seen;
    for(int i=0;i<np;i++){
        auto owner=index.find((unsigned long long)ctx->h_parts.indx[i]);if(owner==index.end())return 1;local_index[i]=owner->second;
        for(int edge=0;edge<polygons[i].n_corners;edge++){
            int j=polygons[i].edge_neigh[edge];if(j<0 || j>=nt)return 1;
            unsigned long long a=ctx->h_parts.indx[i],b=ctx->h_parts.indx[j];if(a>=b)continue;
            auto neighbor=index.find(b);if(neighbor==index.end() || !seen.insert({a,b}).second)return 1;
            GfsCoupledPair p={};p.i=i;p.j=j;p.edge=edge;p.gi=owner->second;p.gj=neighbor->second;pairs.push_back(p);
        }
    }
    cudaStream_t stream=(cudaStream_t)ctx->stream;int m=(int)pairs.size();
    GfsCoupledPair *dp=gfs_ec_alloc<GfsCoupledPair>(m);GfsEnrichAudit *da=gfs_ec_alloc<GfsEnrichAudit>(n);
    GfsEnrichConfig *dc=gfs_ec_alloc<GfsEnrichConfig>(1);GfsCoupledVelocity *dv=gfs_ec_alloc<GfsCoupledVelocity>(n);
    unsigned *bad=gfs_ec_alloc<unsigned>(1),*hll=gfs_ec_alloc<unsigned>(1);int *dlocal=gfs_ec_alloc<int>(np);
    CUDA_CHECK(cudaMemcpyAsync(da,state.data(),n*sizeof(*da),cudaMemcpyHostToDevice,stream));
    CUDA_CHECK(cudaMemcpyAsync(dc,cfg,sizeof(*dc),cudaMemcpyHostToDevice,stream));
    if(m)CUDA_CHECK(cudaMemcpyAsync(dp,pairs.data(),m*sizeof(*dp),cudaMemcpyHostToDevice,stream));
    if(np)CUDA_CHECK(cudaMemcpyAsync(dlocal,local_index.data(),np*sizeof(int),cudaMemcpyHostToDevice,stream));
    CUDA_CHECK(cudaMemsetAsync(bad,0,sizeof(unsigned),stream));CUDA_CHECK(cudaMemsetAsync(hll,0,sizeof(unsigned),stream));
    gfs_ec_velocity<<<(n+127)/128,128,0,stream>>>(da,n,dc,dv,bad);
    if(m)gfs_ec_pair_geometry<<<(m+127)/128,128,0,stream>>>(ctx->d_parts,ctx->d_polygons,da,dc,dp,m);
    CUDA_CHECK(cudaGetLastError());if(m)CUDA_CHECK(cudaMemcpyAsync(pairs.data(),dp,m*sizeof(*dp),cudaMemcpyDeviceToHost,stream));
    CUDA_CHECK(cudaStreamSynchronize(stream));size_t count=0;
    for(int i=0;i<m;i++){
        if(pairs[i].status!=1 || pairs[i].count<0){fprintf(stderr,"[EC] invalid pair geometry=%d\n",i);return 1;}
        pairs[i].offset=(int)count;count+=pairs[i].count;if(count>INT_MAX)return 1;
    }
    if(m)CUDA_CHECK(cudaMemcpyAsync(dp,pairs.data(),m*sizeof(*dp),cudaMemcpyHostToDevice,stream));
    GfsCoupledWorkspace w={};w.n=n;w.m=(int)count;w.stream=stream;w.sum=sum;w.user=user;
    w.points=gfs_ec_alloc<GfsWeightPoint>(count);w.jac=gfs_ec_alloc<GfsCoupledJac>(count);
    w.z=gfs_ec_alloc<double>(n);w.v=gfs_ec_alloc<double>(n);w.R=gfs_ec_alloc<double>(n);w.diag=gfs_ec_alloc<double>(n);w.mu.resize(n);
    for(int i=0;i<n;i++){w.mu[i]=initial[i].integral[0]/initial[i].basis.geometry.area;if(!(w.mu[i]>0))return 1;wdot[i]=0;}
    if(m)gfs_ec_points<<<(m+63)/64,64,0,stream>>>(ctx->d_parts,da,dc,dv,dp,m,w.points,bad,hll);
    unsigned failures=0,fallback=0;CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpyAsync(&failures,bad,sizeof(unsigned),cudaMemcpyDeviceToHost,stream));
    CUDA_CHECK(cudaMemcpyAsync(&fallback,hll,sizeof(unsigned),cudaMemcpyDeviceToHost,stream));CUDA_CHECK(cudaStreamSynchronize(stream));
    double totals[4]={(double)count,(double)m,(double)failures,(double)fallback};sum(user,totals,4);
    if(totals[2]){fprintf(stderr,"[EC] invalid physical states/fans=%g\n",totals[2]);return 1;}
    report[EC_POINTS]=totals[0];report[EC_PAIRS]=totals[1];report[EC_HLL_POINTS]=totals[3];
    std::vector<double> r(n),d(n),direction(n),minus(n),trial(n),tr(n),td(n);
    int newton=0,cg=0,backtracks=0;bool converged=false;
    for(;newton<40;newton++){
        w.residual(wdot,r.data(),d.data());double scale=1;for(int i=0;i<n;i++)scale=fmax(scale,1+fabs(wdot[i]));
        if(w.residual_norm(r.data())<1.e-12*scale){converged=true;break;}
        for(int i=0;i<n;i++)minus[i]=-r[i];int iterations=w.cg(minus.data(),d.data(),direction.data());
        if(iterations<0){fprintf(stderr,"[EC] CG failed\n");return 1;}cg+=iterations;
        double factor=1,merit=w.merit(r.data());bool accepted=false;
        for(int retry=0;retry<32;retry++){
            for(int i=0;i<n;i++)trial[i]=wdot[i]+factor*direction[i];w.residual(trial.data(),tr.data(),td.data());
            if(w.merit(tr.data())<merit*(1-1.e-4*factor)){
                for(int i=0;i<n;i++)wdot[i]=trial[i];accepted=true;break;
            }factor*=.5;backtracks++;
        }
        if(!accepted){fprintf(stderr,"[EC] nonlinear line search failed\n");return 1;}
    }
    if(!converged){fprintf(stderr,"[EC] Newton budget exhausted\n");return 1;}
    w.residual(wdot,r.data(),d.data());report[EC_RESIDUAL]=w.residual_norm(r.data());
    report[EC_NEWTON]=newton;report[EC_CG]=cg;report[EC_BACKTRACKS]=backtracks;
    std::vector<double> u(n),v(n),Ju(n),Jv(n);for(int i=0;i<n;i++){u[i]=sin(.37*(i+1));v[i]=cos(.29*(i+1));}
    w.product(u.data(),Ju.data());w.product(v.data(),Jv.data());
    double uv=gfs_ec_dot(u.data(),Jv.data(),n),vu=gfs_ec_dot(v.data(),Ju.data(),n);
    report[EC_SYMMETRY]=fabs(uv-vu)/(1+fabs(uv)+fabs(vu));
    long double bound=0;for(int i=0;i<n;i++)bound+=(long double)w.mu[i]*u[i]*u[i];
    report[EC_RAYLEIGH]=gfs_ec_dot(u.data(),Ju.data(),n)/(double)bound;
    double epsilon=1.e-5;for(int i=0;i<n;i++)trial[i]=wdot[i]+epsilon*u[i];w.residual(trial.data(),tr.data(),td.data());
    for(int i=0;i<n;i++)trial[i]=wdot[i]-epsilon*u[i];w.residual(trial.data(),minus.data(),td.data());
    for(int i=0;i<n;i++)tr[i]=(tr[i]-minus[i])/(2*epsilon);
    report[EC_SECANT]=gfs_ec_dot(u.data(),tr.data(),n)/(double)bound;
    w.residual(wdot,r.data(),d.data());
    double *drhs=gfs_ec_alloc<double>((size_t)n*40),*dgeo=gfs_ec_alloc<double>((size_t)n*40),*dledger=gfs_ec_alloc<double>((size_t)m*10);
    CUDA_CHECK(cudaMemsetAsync(drhs,0,(size_t)n*40*sizeof(double),stream));CUDA_CHECK(cudaMemsetAsync(dgeo,0,(size_t)n*40*sizeof(double),stream));
    if(np)gfs_ec_volume<<<(np+63)/64,64,0,stream>>>(ctx->d_parts,ctx->d_polygons,dlocal,np,da,dc,drhs,bad);
    if(m)gfs_ec_face_rhs<<<(m+63)/64,64,0,stream>>>(dp,m,w.points,w.z,da,dc,drhs,dgeo,dledger,bad);
    CUDA_CHECK(cudaGetLastError());std::vector<double> ledger((size_t)m*10);
    CUDA_CHECK(cudaMemcpyAsync(rhs,drhs,(size_t)n*40*sizeof(double),cudaMemcpyDeviceToHost,stream));
    CUDA_CHECK(cudaMemcpyAsync(geometry,dgeo,(size_t)n*40*sizeof(double),cudaMemcpyDeviceToHost,stream));
    if(m)CUDA_CHECK(cudaMemcpyAsync(ledger.data(),dledger,(size_t)m*10*sizeof(double),cudaMemcpyDeviceToHost,stream));
    CUDA_CHECK(cudaMemcpyAsync(&failures,bad,sizeof(unsigned),cudaMemcpyDeviceToHost,stream));CUDA_CHECK(cudaStreamSynchronize(stream));
    double fail=(double)failures;sum(user,&fail,1);if(fail){fprintf(stderr,"[EC] full Euler RHS invalid=%g\n",fail);return 1;}
    sum(user,rhs,n*40);sum(user,geometry,n*40);long double mass_sum=0,gauge=0;double max_evolution=0;
    for(int i=0;i<n;i++){
        mass_sum+=rhs[i*40];gauge+=(long double)w.mu[i]*wdot[i];
        const GfsEnrichBasis &b=initial[i].basis;double rate=sqrt(2*b.scale[3]/b.scale[0])/fmin(b.geometry.l00,b.geometry.l11);
        for(int c=0;c<4;c++)for(int j=0;j<10;j++)max_evolution=fmax(max_evolution,
            fabs(rhs[i*40+c*10+j]-geometry[i*40+c*10+j])/(b.geometry.area*b.scale[c]*rate));
    }
    report[EC_MASS_SUM]=(double)mass_sum;report[EC_GAUGE]=(double)gauge;report[EC_MAX_RHS_MINUS_GEOMETRY]=max_evolution;
    char path[160];snprintf(path,sizeof(path),"enriched_pairs_r%d.bin",rank);FILE *out=fopen(path,"wx");if(!out)return 1;
    unsigned long long header[4]={0x4746534543503031ULL,1,(unsigned long long)m,144};fwrite(header,sizeof(header),1,out);
    for(int k=0;k<m;k++){
        const GfsCoupledPair &p=pairs[k];unsigned long long key[2]={ids[p.gi],ids[p.gj]};fwrite(key,sizeof(key),1,out);
        fwrite(initial[p.gi].basis.generator,sizeof(double),2,out);double corners[4]={p.x1,p.y1,p.x2,p.y2};fwrite(corners,sizeof(corners),1,out);
        fwrite(ledger.data()+k*10,sizeof(double),10,out);
    }
    if(ferror(out) || fclose(out))return 1;
    CUDA_CHECK(cudaFree(dp));CUDA_CHECK(cudaFree(da));CUDA_CHECK(cudaFree(dc));CUDA_CHECK(cudaFree(dv));
    CUDA_CHECK(cudaFree(bad));CUDA_CHECK(cudaFree(hll));CUDA_CHECK(cudaFree(dlocal));
    CUDA_CHECK(cudaFree(w.points));CUDA_CHECK(cudaFree(w.jac));CUDA_CHECK(cudaFree(w.z));CUDA_CHECK(cudaFree(w.v));
    CUDA_CHECK(cudaFree(w.R));CUDA_CHECK(cudaFree(w.diag));CUDA_CHECK(cudaFree(drhs));CUDA_CHECK(cudaFree(dgeo));CUDA_CHECK(cudaFree(dledger));
    return 0;
}
