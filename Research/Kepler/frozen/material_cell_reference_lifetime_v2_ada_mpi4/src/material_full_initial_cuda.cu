#include <cuda_runtime.h>
#include <stdint.h>
#include "material_full_initial_batch.h"
using DD=MFDoubleDouble;
__constant__ GfsEnrichConfig mf_config;
__device__ static DD dd(double a){return {a,0};}
__device__ static DD add(DD a,DD b){double s=a.hi+b.hi,z=s-a.hi,e=(a.hi-(s-z))+(b.hi-z);e+=a.lo+b.lo;double h=s+e;return {h,e-(h-s)};}
__device__ static DD neg(DD a){return {-a.hi,-a.lo};}
__device__ static DD sub(DD a,DD b){return add(a,neg(b));}
__device__ static DD mul(DD a,DD b){double p=a.hi*b.hi,e=fma(a.hi,b.hi,-p);e+=a.hi*b.lo+a.lo*b.hi;e+=a.lo*b.lo;double h=p+e;return {h,e-(h-p)};}
__device__ static DD divide(DD a,DD b){double q=a.hi/b.hi;DD r=sub(a,mul(dd(q),b));return add(dd(q),dd((r.hi+r.lo)/b.hi));}
__device__ static double rounded(DD a){return a.hi+a.lo;}
__device__ static int active_edges(const STTriangle&t,const GfsEnrichQuadrature&r,int*edge){
    int n=0;for(int k=0;k<3;k++){int j=(k+1)%3;double bx=t.vertex[k][0]-r.anchor[0],by=t.vertex[k][1]-r.anchor[1],cx=t.vertex[j][0]-r.anchor[0],cy=t.vertex[j][1]-r.anchor[1],jac=bx*cy-by*cx;
        if(k==r.anchor_edge||(bx==0&&by==0)||(cx==0&&cy==0))jac=0;if(jac<0)return -1;if(jac>0)edge[n++]=k;}return n;
}
__device__ static int point(const STTriangle&t,const GfsEnrichQuadrature&r,size_t s,int slot,const int*edge,STPoint*q){
    *q={};q->triangle=s;q->boundary=-1;int a=slot%mf_config.na,radial=(slot/mf_config.na)%mf_config.nr,sec=(slot/(mf_config.na*mf_config.nr))%r.ns,e=edge[slot/(mf_config.na*mf_config.nr*r.ns)];
    double X[3],Y[3];for(int k=0;k<3;k++){X[k]=t.vertex[k][0];Y[k]=t.vertex[k][1];}GfsMomentPolygon p={X,Y,3,0};int original=((e*r.ns+sec)*mf_config.nr+radial)*mf_config.na+a;
    if(!gfs_enrich_quad(p,&mf_config,&r,original,&q->X,&q->Y,&q->weight)||!(q->weight>0))return 0;
    double u[4];if(!gfs_enrich_reference(&mf_config,q->X+12,q->Y+12,u,0,0))return 0;q->rho=u[0];
    DD ux=sub(dd(X[1]),dd(X[0])),uy=sub(dd(Y[1]),dd(Y[0])),vx=sub(dd(X[2]),dd(X[0])),vy=sub(dd(Y[2]),dd(Y[0])),dx=sub(dd(q->X),dd(X[0])),dy=sub(dd(q->Y),dd(Y[0])),det=sub(mul(ux,vy),mul(uy,vx));
    q->l[1]=rounded(divide(sub(mul(dx,vy),mul(dy,vx)),det));q->l[2]=rounded(divide(sub(mul(ux,dy),mul(uy,dx)),det));q->l[0]=1-q->l[1]-q->l[2];return q->l[0]>0&&q->l[1]>0&&q->l[2]>0;
}
// Optional positive original-point subset, used only by evolved operators.
struct MFRuleView {const uint64_t*offset=nullptr;const uint32_t*index=nullptr;const double*weight=nullptr;};
__device__ static int evolved_count(MFRuleView v,size_t s,int original){uint64_t n=v.offset?v.offset[s+1]-v.offset[s]:0;return n?int(n):original;}
__device__ static int evolved_point(const STTriangle&t,const GfsEnrichQuadrature&r,size_t s,int slot,const int*edge,MFRuleView v,STPoint*q){
    uint64_t n=v.offset?v.offset[s+1]-v.offset[s]:0;
    if(n&&(!v.index||!v.weight||uint64_t(slot)>=n))return 0;
    int original=n?int(v.index[v.offset[s]+slot]):slot;
    if(!point(t,r,s,original,edge,q))return 0;
    if(n){q->weight=v.weight[v.offset[s]+slot];if(!(q->weight>0)||!isfinite(q->weight))return 0;}
    return 1;
}
__global__ static void mass_and_initial_gravity(const STTriangle*t,const STCoefficient*c,
        const GfsEnrichQuadrature*rule,size_t start,size_t count,MFInitialBlock*out){
    size_t at=blockIdx.x;if(at>=count)return;size_t s=start+at;const auto&T=t[s];const auto&R=rule[s];int tid=threadIdx.x;
    __shared__ double Nvalue[64][ST_MAX],force[2][64][ST_MAX],weight[64],density[64],metric[64][MF_INITIAL_METRICS],Jvalue[64];
    __shared__ int edge[3],ne,nq,bad[64],bad_kind[64];
    if(tid==0){out[at]={};out[at].minJ=DBL_MAX;out[at].bad_point=-1;ne=active_edges(T,R,edge);nq=ne>0?ne*R.ns*mf_config.na*mf_config.nr:0;if(ne<=0){out[at].bad_point=0;out[at].bad_kind=1;}}
    __syncthreads();if(ne<=0)return;
    int row=0,col=0,pair=tid,nt=T.n*(T.n+1)/2;if(pair<nt){while(pair>=T.n-row){pair-=T.n-row;row++;}col=row+pair;}
    DD mass={0,0},fx={0,0},fy={0,0},sum={0,0};
    for(int base=0;base<nq;base+=64){
        if(tid<64){int qindex=base+tid;bad[tid]=0;bad_kind[tid]=0;weight[tid]=0;density[tid]=0;Jvalue[tid]=DBL_MAX;for(int j=0;j<ST_MAX;j++){Nvalue[tid][j]=force[0][tid][j]=force[1][tid][j]=0;}for(int j=0;j<MF_INITIAL_METRICS;j++)metric[tid][j]=0;
            if(qindex<nq){STPoint q;STJet N[ST_MAX];GfsMaterialPoint12 p;
                if(!point(T,R,s,qindex,edge,&q)){bad[tid]=1;bad_kind[tid]=2;}
                else if(!st_shape(T,q,N)||!st_point(T,c,N,&p)){bad[tid]=1;bad_kind[tid]=3;}
                else{weight[tid]=q.weight;density[tid]=q.rho;Jvalue[tid]=p.J;for(int j=0;j<T.n;j++)Nvalue[tid][j]=N[j].v;
                    double z=p.position[0]*p.position[0]+p.position[1]*p.position[1]+mf_config.field.ic.eps*mf_config.field.ic.eps,g=pow(z,-1.5),gx=-p.position[0]*g,gy=-p.position[1]*g,sigma[3]={0,0,0};
                    for(int begin=0;begin<T.n;begin+=12){int n=min(12,T.n-begin);double A[12],Ax[12],Ay[12],F[2][12]={},work=0,gwork=0;for(int j=0;j<n;j++){A[j]=N[begin+j].v;Ax[j]=N[begin+j].x;Ay[j]=N[begin+j].y;}if(!gfs_material12_stress_at_point(n,&p,A,Ax,Ay,q.rho,q.weight,sigma,gx,gy,F,&work,&gwork)){bad[tid]=1;bad_kind[tid]=4;}for(int j=0;j<n;j++){force[0][tid][begin+j]=F[0][j];force[1][tid][begin+j]=F[1][j];}}
                    double dm=q.weight*q.rho,x=p.position[0],y=p.position[1],vx=p.velocity[0],vy=p.velocity[1],omega=pow(q.X*q.X+q.Y*q.Y+mf_config.field.ic.eps*mf_config.field.ic.eps,-.75),v0x=-q.Y*omega,v0y=q.X*omega;
                    double value[MF_INITIAL_METRICS]={dm,q.weight,q.weight*p.J,q.weight*p.J*(p.grad_v[0]+p.grad_v[3]),q.weight*mf_config.field.ic.pressure/(mf_config.field.ic.gamma-1),.5*dm*(vx*vx+vy*vy),dm*(x*vy-y*vx),-dm/sqrt(z),dm*(vx*gx+vy*gy),dm*(x*gy-y*gx),dm*x,dm*y,dm*vx,dm*vy,dm*(x*x+y*y),.5*dm*((vx-v0x)*(vx+v0x)+(vy-v0y)*(vy+v0y))};for(int j=0;j<MF_INITIAL_METRICS;j++)metric[tid][j]=value[j];
                }
            }
        }
        __syncthreads();int end=min(64,nq-base);
        if(tid<nt)for(int k=0;k<end;k++)mass=add(mass,mul(mul(mul(dd(weight[k]),dd(density[k])),dd(Nvalue[k][row])),dd(Nvalue[k][col])));
        if(tid<T.n)for(int k=0;k<end;k++){fx=add(fx,dd(force[0][k][tid]));fy=add(fy,dd(force[1][k][tid]));}
        if(tid<MF_INITIAL_METRICS)for(int k=0;k<end;k++)sum=add(sum,dd(metric[k][tid]));
        if(tid==0)for(int k=0;k<end;k++){if(bad[k]&&out[at].bad_point<0){out[at].bad_point=base+k;out[at].bad_kind=bad_kind[k];}if(weight[k]>0)out[at].points++;out[at].minJ=fmin(out[at].minJ,Jvalue[k]);}
        __syncthreads();
    }
    if(tid<nt)out[at].mass[row*ST_MAX+col]=out[at].mass[col*ST_MAX+row]=mass;
    if(tid<T.n){out[at].force[0][tid]=fx;out[at].force[1][tid]=fy;}if(tid<MF_INITIAL_METRICS)out[at].metric[tid]=sum;
}
struct MFCachedPoint;
struct MFCubatureAux;
struct Context {STTriangle*t=nullptr;STCoefficient*c=nullptr;STCoefficient*acceleration=nullptr;GfsEnrichQuadrature*rule=nullptr;MFInitialBlock*out=nullptr;cudaEvent_t begin=nullptr,end=nullptr;size_t nt=0,nc=0;MFCurvatureBlock*curvature=nullptr;STCoefficient*generating=nullptr;MFEvolvedCell*cell=nullptr;MFEvolvedBlock*evolved=nullptr;size_t cells=0;MFCachedPoint*cache=nullptr;uint32_t*cache_prefix=nullptr;size_t cache_capacity=0,cache_start=0,cache_count=0;MFCubatureAux*cubature_aux=nullptr;DD*cubature_output=nullptr;size_t cubature_pairs=0;uint64_t*reduced_offset=nullptr;uint32_t*reduced_index=nullptr;double*reduced_weight=nullptr;bool reduced_active=false,reduced_reference_only=false;};
static MFRuleView evolved_rule_view(const Context*s){return s->reduced_active?MFRuleView{s->reduced_offset,s->reduced_index,s->reduced_weight}:MFRuleView{};}
extern "C" void mf_initial_destroy(void*ptr){auto*s=static_cast<Context*>(ptr);if(!s)return;cudaFree(s->t);cudaFree(s->c);cudaFree(s->acceleration);cudaFree(s->curvature);cudaFree(s->generating);cudaFree(s->cell);cudaFree(s->evolved);cudaFree(s->cache);cudaFree(s->cache_prefix);cudaFree(s->cubature_aux);cudaFree(s->cubature_output);cudaFree(s->reduced_offset);cudaFree(s->reduced_index);cudaFree(s->reduced_weight);cudaFree(s->rule);cudaFree(s->out);if(s->begin)cudaEventDestroy(s->begin);if(s->end)cudaEventDestroy(s->end);delete s;}
extern "C" int mf_initial_create(const STTriangle*t,size_t nt,const STCoefficient*c,size_t nc,const GfsEnrichQuadrature*rule,const GfsEnrichConfig*cfg,void**output){
    *output=nullptr;auto*s=new Context;s->nt=nt;s->nc=nc;cudaError_t error=cudaSuccess;
    #define TRY(call) do{if(error==cudaSuccess)error=(call);}while(0)
    TRY(cudaMalloc((void**)&s->t,nt*sizeof(*t)));TRY(cudaMalloc((void**)&s->c,nc*sizeof(*c)));TRY(cudaMalloc((void**)&s->rule,nt*sizeof(*rule)));TRY(cudaMalloc((void**)&s->out,MF_BATCH_TRIANGLES*sizeof(*s->out)));
    TRY(cudaMemcpy(s->t,t,nt*sizeof(*t),cudaMemcpyHostToDevice));TRY(cudaMemcpy(s->c,c,nc*sizeof(*c),cudaMemcpyHostToDevice));TRY(cudaMemcpy(s->rule,rule,nt*sizeof(*rule),cudaMemcpyHostToDevice));TRY(cudaMemcpyToSymbol(mf_config,cfg,sizeof(*cfg)));TRY(cudaEventCreate(&s->begin));TRY(cudaEventCreate(&s->end));
    #undef TRY
    if(error!=cudaSuccess){mf_initial_destroy(s);return error;}*output=s;return 0;
}
extern "C" int mf_initial_batch(void*ptr,size_t start,size_t count,MFInitialBlock*out,double*seconds){
    auto*s=static_cast<Context*>(ptr);if(count>MF_BATCH_TRIANGLES||start+count>s->nt)return cudaErrorInvalidValue;cudaError_t error;
    #define TRY(call) do{error=(call);if(error!=cudaSuccess)return error;}while(0)
    TRY(cudaEventRecord(s->begin));mass_and_initial_gravity<<<count,256>>>(s->t,s->c,s->rule,start,count,s->out);TRY(cudaGetLastError());TRY(cudaEventRecord(s->end));TRY(cudaEventSynchronize(s->end));float ms=0;TRY(cudaEventElapsedTime(&ms,s->begin,s->end));*seconds=ms/1000.;TRY(cudaMemcpy(out,s->out,count*sizeof(*out),cudaMemcpyDeviceToHost));
    #undef TRY
    return 0;
}
__global__ static void curvature_kernel(const STTriangle*t,const STCoefficient*c,const STCoefficient*acc,
        const GfsEnrichQuadrature*rule,size_t start,size_t count,MFCurvatureBlock*out){
    size_t at=blockIdx.x;if(at>=count)return;size_t s=start+at;const auto&T=t[s];const auto&R=rule[s];int tid=threadIdx.x;
    __shared__ double metric[64][MF_CURV_METRICS],chi[64],jac[64],velocity_error[64];
    __shared__ int edge[3],ne,nq,bad[64];
    if(tid==0){out[at]={};out[at].minJ=DBL_MAX;out[at].bad_point=out[at].max_chi_point=-1;ne=active_edges(T,R,edge);nq=ne>0?ne*R.ns*mf_config.na*mf_config.nr:0;if(ne<=0)out[at].bad_point=0;}
    __syncthreads();if(ne<=0)return;DD sum={0,0};
    for(int base=0;base<nq;base+=64){
        if(tid<64){for(int j=0;j<MF_CURV_METRICS;j++)metric[tid][j]=0;chi[tid]=velocity_error[tid]=0;jac[tid]=DBL_MAX;bad[tid]=0;int qindex=base+tid;
            if(qindex<nq){STPoint q;STJet N[ST_MAX];GfsMaterialPoint12 p;if(!point(T,R,s,qindex,edge,&q)||!st_shape(T,q,N)||!st_point(T,c,N,&p))bad[tid]=1;
                else{double acceleration[2]={},D[4]={};for(int j=0;j<T.n;j++)for(int d=0;d<2;d++){double a=acc[T.id[j]].v[d];acceleration[d]+=a*N[j].v;D[2*d]+=a*N[j].x;D[2*d+1]+=a*N[j].y;}
                    double x=q.X,y=q.Y,z=x*x+y*y+mf_config.field.ic.eps*mf_config.field.ic.eps,o=pow(z,-.75),o2=pow(z,-1.5),B[4]={1.5*x*y*o/z,-o+1.5*y*y*o/z,o-1.5*x*x*o/z,-1.5*x*y*o/z};
                    double exact=-2*o2+3*(x*x+y*y)*o2/z-(B[0]*B[0]+2*B[1]*B[2]+B[3]*B[3]);double current=D[0]+D[3]-(p.grad_v[0]*p.grad_v[0]+2*p.grad_v[1]*p.grad_v[2]+p.grad_v[3]*p.grad_v[3]);
                    chi[tid]=mf_config.field.ic.gamma*fabs(current-exact)/o2;jac[tid]=p.J;velocity_error[tid]=fmax(fabs(p.velocity[0]+y*o),fabs(p.velocity[1]-x*o));double dm=q.weight*q.rho,X=p.position[0],Y=p.position[1],vx=p.velocity[0],vy=p.velocity[1],r2=X*X+Y*Y+mf_config.field.ic.eps*mf_config.field.ic.eps,g=pow(r2,-1.5),gx=-X*g,gy=-Y*g;
                    double value[MF_CURV_METRICS]={q.weight,q.weight*chi[tid],q.weight*p.J,q.weight*p.J*(p.grad_v[0]+p.grad_v[3]),dm,dm*(vx*acceleration[0]+vy*acceleration[1]),dm*(X*acceleration[1]-Y*acceleration[0]),dm*(vx*gx+vy*gy),dm*(X*gy-Y*gx),.5*dm*(vx*vx+vy*vy),dm*X,dm*Y};for(int j=0;j<MF_CURV_METRICS;j++)metric[tid][j]=value[j];
                    if(!isfinite(chi[tid])||!isfinite(value[5])||!isfinite(value[6]))bad[tid]=1;
                }
            }
        }
        __syncthreads();int end=min(64,nq-base);if(tid<MF_CURV_METRICS)for(int k=0;k<end;k++)sum=add(sum,dd(metric[k][tid]));
        if(tid==0)for(int k=0;k<end;k++){if(bad[k]&&out[at].bad_point<0)out[at].bad_point=base+k;if(chi[k]>out[at].max_chi){out[at].max_chi=chi[k];out[at].max_chi_point=base+k;}out[at].minJ=fmin(out[at].minJ,jac[k]);out[at].max_velocity_error=fmax(out[at].max_velocity_error,velocity_error[k]);if(!bad[k])out[at].points++;}
        __syncthreads();
    }
    if(tid<MF_CURV_METRICS)out[at].metric[tid]=sum;
}
extern "C" int mf_curvature_upload(void*ptr,const STCoefficient*a,size_t nc){
    auto*s=static_cast<Context*>(ptr);if(nc!=s->nc||s->acceleration||s->curvature)return cudaErrorInvalidValue;cudaError_t e=cudaMalloc((void**)&s->acceleration,nc*sizeof(*a));if(e==cudaSuccess)e=cudaMalloc((void**)&s->curvature,MF_BATCH_TRIANGLES*sizeof(MFCurvatureBlock));if(e==cudaSuccess)e=cudaMemcpy(s->acceleration,a,nc*sizeof(*a),cudaMemcpyHostToDevice);return e;
}
extern "C" int mf_curvature_batch(void*ptr,size_t start,size_t count,MFCurvatureBlock*out,double*seconds){
    auto*s=static_cast<Context*>(ptr);if(count>MF_BATCH_TRIANGLES||start+count>s->nt||!s->acceleration)return cudaErrorInvalidValue;cudaError_t e;
    #define TRY(call) do{e=(call);if(e!=cudaSuccess)return e;}while(0)
    TRY(cudaEventRecord(s->begin));curvature_kernel<<<count,256>>>(s->t,s->c,s->acceleration,s->rule,start,count,s->curvature);TRY(cudaGetLastError());TRY(cudaEventRecord(s->end));TRY(cudaEventSynchronize(s->end));float ms=0;TRY(cudaEventElapsedTime(&ms,s->begin,s->end));*seconds=ms/1000.;TRY(cudaMemcpy(out,s->curvature,count*sizeof(*out),cudaMemcpyDeviceToHost));
    #undef TRY
    return 0;
}
#include "material_full_evolved_cuda.inc"
#include "material_full_gpu_mass.inc"
#include "material_full_mass_resident.inc"
#include "material_full_cached_cuda.inc"

#include "material_full_cubature_gpu.inc"

#include "material_full_reduced_rule.inc"
