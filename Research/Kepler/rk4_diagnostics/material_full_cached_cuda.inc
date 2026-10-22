// Bounded transient cache of the SAME original point/shape/field evaluations.
// Complete source cells share each tile, so canonical delta and mean pressure
// can be formed before reusing the tile. No quadrature order is changed.
struct MFCachedPoint {STJet N[ST_MAX];GfsMaterialPoint12 p;double weight,rho,j0;int bad;};
static_assert(sizeof(MFCachedPoint)==568,"declared bounded point cache ABI");
__global__ static void mf_cache_points(const STTriangle*t,const STCoefficient*c,
        const STCoefficient*generating,const GfsEnrichQuadrature*rule,
        const MFEvolvedCell*cell,size_t start,size_t count,const uint32_t*prefix,MFRuleView view,MFCachedPoint*cache){
    size_t at=blockIdx.x;if(at>=count)return;size_t s=start+at;const auto&T=t[s];const auto&R=rule[s];
    __shared__ int edge[3],ne,nq;int tid=threadIdx.x;
    if(tid==0){ne=active_edges(T,R,edge);nq=evolved_count(view,s,ne>0?ne*R.ns*mf_config.na*mf_config.nr:0);}__syncthreads();
    size_t n=prefix[at+1]-prefix[at];
    for(size_t qindex=tid;qindex<n;qindex+=blockDim.x){auto&out=cache[prefix[at]+qindex];out.bad=0;
        if(ne<=0||size_t(nq)!=n){out.bad=1;continue;}STPoint q;GfsMaterialPoint12 p0;
        if(!evolved_point(T,R,s,qindex,edge,view,&q)){out.bad=2;continue;}
        if(!st_shape(T,q,out.N)||!st_point(T,c,out.N,&out.p)){out.bad=3;continue;}
        if(!st_point(T,generating,out.N,&p0)){out.bad=4;continue;}
        out.weight=q.weight;out.rho=q.rho;out.j0=p0.J+cell[T.cell].generating_delta;
    }
}
__global__ static void evolved_cached_kernel(const STTriangle*t,const STCoefficient*c,
        const STCoefficient*generating,const GfsEnrichQuadrature*rule,
        const MFEvolvedCell*cell,const MFCachedPoint*cache,const uint32_t*prefix,int pass,size_t start,size_t count,MFRuleView view,MFEvolvedBlock*out){
    size_t at=blockIdx.x;if(at>=count)return;size_t s=start+at;const auto&T=t[s];const auto&R=rule[s];int tid=threadIdx.x;
    __shared__ double force[2][64][ST_MAX],metric[64][MF_EVOL_METRICS];
    __shared__ double jac[64],jhat[64],heat[64],grad[64],length[64],sound[64],omega[64];
    __shared__ int edge[3],ne,nq,bad[64];
    if(tid==0){out[at]={};out[at].minJ=out[at].min_j=out[at].min_U=out[at].min_length=DBL_MAX;out[at].bad_point=-1;ne=active_edges(T,R,edge);nq=evolved_count(view,s,ne>0?ne*R.ns*mf_config.na*mf_config.nr:0);if(ne<=0){out[at].bad_point=0;out[at].bad_kind=1;}}
    __syncthreads();if(ne<=0)return;DD sum={0,0},fx={0,0},fy={0,0};
    for(int base=0;base<nq;base+=64){
        if(tid<64){for(int k=0;k<MF_EVOL_METRICS;k++)metric[tid][k]=0;for(int d=0;d<2;d++)for(int k=0;k<ST_MAX;k++)force[d][tid][k]=0;
            jac[tid]=jhat[tid]=heat[tid]=length[tid]=DBL_MAX;grad[tid]=sound[tid]=omega[tid]=0;bad[tid]=0;int qindex=base+tid;
            if(qindex<nq){const auto&cached=cache[prefix[at]+qindex];STPoint q={};const STJet*N=cached.N;GfsMaterialPoint12 p={};
                if(cached.bad)bad[tid]=cached.bad;
                else{q.weight=cached.weight;q.rho=cached.rho;p=cached.p;const auto&C=cell[T.cell];double div=p.grad_v[0]+p.grad_v[3],jd=p.J*div,dm=q.weight*q.rho;
                    jac[tid]=p.J;metric[tid][0]=q.weight;metric[tid][1]=q.weight*p.J;metric[tid][2]=q.weight*jd;metric[tid][3]=dm;
                    double x=p.position[0],y=p.position[1],vx=p.velocity[0],vy=p.velocity[1],z=x*x+y*y+mf_config.field.ic.eps*mf_config.field.ic.eps;
                    double g=pow(z,-1.5),gx=-x*g,gy=-y*g;
                    grad[tid]=sqrt(p.grad_v[0]*p.grad_v[0]+p.grad_v[1]*p.grad_v[1]+p.grad_v[2]*p.grad_v[2]+p.grad_v[3]*p.grad_v[3]);omega[tid]=sqrt(g);
                    double refh=DBL_MAX;for(int k=0;k<3;k++)refh=fmin(refh,1/hypot(T.grad[k][0],T.grad[k][1]));
                    double tr=p.F[0]*p.F[0]+p.F[1]*p.F[1]+p.F[2]*p.F[2]+p.F[3]*p.F[3],maxsv=sqrt(.5*(tr+sqrt(fmax(0.,tr*tr-4*p.J*p.J))));length[tid]=refh*p.J/maxsv;
                    if(pass>0){
                        if(!(cached.j0>0))bad[tid]=4;
                        else{double j0=cached.j0,j=p.J+C.delta;
                            if(!(j0>0)||!(j>0))bad[tid]=5;
                            else{double gamma=mf_config.field.ic.gamma,P0=mf_config.field.ic.pressure,lr=log1p((j-j0)/j0),dp=P0*expm1(-gamma*lr),P=P0+dp,U0=P0*q.weight*j0/(gamma-1),U=U0*exp((1-gamma)*lr),Udot=-P*q.weight*(jd+C.delta_dot);
                                if(!(P>0)||!(U>0)||!isfinite(Udot)||!isfinite(lr))bad[tid]=6;
                                else{jhat[tid]=j;heat[tid]=U;sound[tid]=sqrt(gamma*P*j/q.rho);
                                    metric[tid][4]=q.weight*dp;metric[tid][5]=U;metric[tid][6]=Udot;metric[tid][7]=q.weight*j;
                                    metric[tid][8]=.5*dm*(vx*vx+vy*vy);metric[tid][9]=dm*(x*vy-y*vx);metric[tid][10]=-dm/sqrt(z);metric[tid][11]=dm*(vx*gx+vy*gy);metric[tid][12]=dm*(x*gy-y*gx);metric[tid][13]=dm*x;metric[tid][14]=dm*y;metric[tid][15]=dm*vx;metric[tid][16]=dm*vy;metric[tid][17]=dm*(x*x+y*y);metric[tid][18]=2*dm*(x*vx+y*vy);
                                    if(pass==2){double residual=dp-C.mean_deltaP,sigma[3]={-residual,0,-residual};for(int begin=0;begin<T.n;begin+=12){int n=min(12,T.n-begin);double A[12],Ax[12],Ay[12],F[2][12]={},work=0,gwork=0;for(int k=0;k<n;k++){A[k]=N[begin+k].v;Ax[k]=N[begin+k].x;Ay[k]=N[begin+k].y;}
                                            if(!gfs_material12_stress_at_point(n,&p,A,Ax,Ay,q.rho,q.weight,sigma,gx,gy,F,&work,&gwork))bad[tid]=7;
                                            for(int d=0;d<2;d++)for(int k=0;k<n;k++)force[d][tid][begin+k]=F[d][k];}}
                                }
                            }
                        }
                    }
                }
            }
        }
        __syncthreads();int n=min(64,nq-base);
        if(tid<MF_EVOL_METRICS)for(int k=0;k<n;k++)sum=add(sum,dd(metric[k][tid]));
        if(tid<T.n&&pass==2)for(int k=0;k<n;k++){fx=add(fx,dd(force[0][k][tid]));fy=add(fy,dd(force[1][k][tid]));}
        if(tid==0)for(int k=0;k<n;k++){if(bad[k]&&out[at].bad_point<0){out[at].bad_point=base+k;out[at].bad_kind=bad[k];}if(!bad[k])out[at].points++;out[at].minJ=fmin(out[at].minJ,jac[k]);out[at].max_J=fmax(out[at].max_J,jac[k]);out[at].min_j=fmin(out[at].min_j,jhat[k]);out[at].min_U=fmin(out[at].min_U,heat[k]);out[at].max_grad_v=fmax(out[at].max_grad_v,grad[k]);out[at].min_length=fmin(out[at].min_length,length[k]);out[at].max_sound=fmax(out[at].max_sound,sound[k]);out[at].max_omega=fmax(out[at].max_omega,omega[k]);}
        __syncthreads();
    }
    if(tid<MF_EVOL_METRICS)out[at].metric[tid]=sum;if(tid<T.n&&pass==2){out[at].force[0][tid]=fx;out[at].force[1][tid]=fy;}
}
extern "C" int mf_evolved_cache_reserve(void*ptr,size_t capacity){
    auto*s=static_cast<Context*>(ptr);if(!s||!s->generating||s->cache||!capacity)return cudaErrorInvalidValue;
    cudaError_t e=cudaMalloc((void**)&s->cache,capacity*sizeof(MFCachedPoint));if(e==cudaSuccess)e=cudaMalloc((void**)&s->cache_prefix,(MF_BATCH_TRIANGLES+1)*sizeof(uint32_t));if(e==cudaSuccess)s->cache_capacity=capacity;return e;
}
extern "C" int mf_evolved_cache_build(void*ptr,size_t start,size_t count,const uint32_t*prefix,MFEvolvedBlock*out,double*seconds){
    auto*s=static_cast<Context*>(ptr);if(!s||!s->cache||!count||count>MF_BATCH_TRIANGLES||start+count>s->nt||prefix[0]!=0||prefix[count]>s->cache_capacity)return cudaErrorInvalidValue;cudaError_t e;
    #define TRY(call) do{e=(call);if(e!=cudaSuccess)return e;}while(0)
    TRY(cudaMemcpy(s->cache_prefix,prefix,(count+1)*sizeof(uint32_t),cudaMemcpyHostToDevice));TRY(cudaEventRecord(s->begin));
    mf_cache_points<<<count,64>>>(s->t,s->c,s->generating,s->rule,s->cell,start,count,s->cache_prefix,evolved_rule_view(s),s->cache);TRY(cudaGetLastError());
    evolved_cached_kernel<<<count,128>>>(s->t,s->c,s->generating,s->rule,s->cell,s->cache,s->cache_prefix,0,start,count,evolved_rule_view(s),s->evolved);TRY(cudaGetLastError());TRY(cudaEventRecord(s->end));TRY(cudaEventSynchronize(s->end));float ms;TRY(cudaEventElapsedTime(&ms,s->begin,s->end));*seconds=ms/1000.;TRY(cudaMemcpy(out,s->evolved,count*sizeof(*out),cudaMemcpyDeviceToHost));s->cache_start=start;s->cache_count=count;
    #undef TRY
    return 0;
}
extern "C" int mf_evolved_cache_batch(void*ptr,int pass,size_t start,size_t count,MFEvolvedBlock*out,double*seconds){
    auto*s=static_cast<Context*>(ptr);if(!s||pass<1||pass>2||!s->cache||start!=s->cache_start||count!=s->cache_count)return cudaErrorInvalidValue;cudaError_t e;
    #define TRY(call) do{e=(call);if(e!=cudaSuccess)return e;}while(0)
    TRY(cudaEventRecord(s->begin));evolved_cached_kernel<<<count,128>>>(s->t,s->c,s->generating,s->rule,s->cell,s->cache,s->cache_prefix,pass,start,count,evolved_rule_view(s),s->evolved);TRY(cudaGetLastError());TRY(cudaEventRecord(s->end));TRY(cudaEventSynchronize(s->end));float ms;TRY(cudaEventElapsedTime(&ms,s->begin,s->end));*seconds=ms/1000.;TRY(cudaMemcpy(out,s->evolved,count*sizeof(*out),cudaMemcpyDeviceToHost));
    #undef TRY
    return 0;
}
