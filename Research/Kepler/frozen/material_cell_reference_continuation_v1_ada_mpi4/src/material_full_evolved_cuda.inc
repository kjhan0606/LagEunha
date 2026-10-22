// Included in the native CUDA translation unit: same fixed shape functions,
// positive quadrature and original symmetric stress operator as initialization.
__global__ static void evolved_kernel(const STTriangle*t,const STCoefficient*c,
        const STCoefficient*generating,const GfsEnrichQuadrature*rule,
        const MFEvolvedCell*cell,int pass,size_t start,size_t count,MFRuleView view,MFEvolvedBlock*out){
    size_t at=blockIdx.x;if(at>=count)return;size_t s=start+at;const auto&T=t[s];const auto&R=rule[s];int tid=threadIdx.x;
    __shared__ double force[2][64][ST_MAX],metric[64][MF_EVOL_METRICS];
    __shared__ double jac[64],jhat[64],heat[64],grad[64],length[64],sound[64],omega[64];
    __shared__ int edge[3],ne,nq,bad[64];
    if(tid==0){out[at]={};out[at].minJ=out[at].min_j=out[at].min_U=out[at].min_length=DBL_MAX;out[at].bad_point=-1;ne=active_edges(T,R,edge);nq=evolved_count(view,s,ne>0?ne*R.ns*mf_config.na*mf_config.nr:0);if(ne<=0){out[at].bad_point=0;out[at].bad_kind=1;}}
    __syncthreads();if(ne<=0)return;DD sum={0,0},fx={0,0},fy={0,0};
    for(int base=0;base<nq;base+=64){
        if(tid<64){for(int k=0;k<MF_EVOL_METRICS;k++)metric[tid][k]=0;for(int d=0;d<2;d++)for(int k=0;k<ST_MAX;k++)force[d][tid][k]=0;
            jac[tid]=jhat[tid]=heat[tid]=length[tid]=DBL_MAX;grad[tid]=sound[tid]=omega[tid]=0;bad[tid]=0;int qindex=base+tid;
            if(qindex<nq){STPoint q;STJet N[ST_MAX];GfsMaterialPoint12 p,p0;
                if(!evolved_point(T,R,s,qindex,edge,view,&q))bad[tid]=2;
                else if(!st_shape(T,q,N)||!st_point(T,c,N,&p))bad[tid]=3;
                else{const auto&C=cell[T.cell];double div=p.grad_v[0]+p.grad_v[3],jd=p.J*div,dm=q.weight*q.rho;
                    jac[tid]=p.J;metric[tid][0]=q.weight;metric[tid][1]=q.weight*p.J;metric[tid][2]=q.weight*jd;metric[tid][3]=dm;
                    double x=p.position[0],y=p.position[1],vx=p.velocity[0],vy=p.velocity[1],z=x*x+y*y+mf_config.field.ic.eps*mf_config.field.ic.eps;
                    double g=pow(z,-1.5),gx=-x*g,gy=-y*g;
                    grad[tid]=sqrt(p.grad_v[0]*p.grad_v[0]+p.grad_v[1]*p.grad_v[1]+p.grad_v[2]*p.grad_v[2]+p.grad_v[3]*p.grad_v[3]);omega[tid]=sqrt(g);
                    double refh=DBL_MAX;for(int k=0;k<3;k++)refh=fmin(refh,1/hypot(T.grad[k][0],T.grad[k][1]));
                    double tr=p.F[0]*p.F[0]+p.F[1]*p.F[1]+p.F[2]*p.F[2]+p.F[3]*p.F[3],maxsv=sqrt(.5*(tr+sqrt(fmax(0.,tr*tr-4*p.J*p.J))));length[tid]=refh*p.J/maxsv;
                    if(pass>0){
                        if(!st_point(T,generating,N,&p0))bad[tid]=4;
                        else{double j0=p0.J+C.generating_delta,j=p.J+C.delta;
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
extern "C" int mf_evolved_bind(void*ptr,const STCoefficient*g,size_t nc,size_t cells){
    auto*s=static_cast<Context*>(ptr);if(nc!=s->nc||s->generating||!cells)return cudaErrorInvalidValue;s->cells=cells;cudaError_t e;
    #define TRY(call) do{e=(call);if(e!=cudaSuccess)return e;}while(0)
    TRY(cudaMalloc((void**)&s->generating,nc*sizeof(*g)));TRY(cudaMalloc((void**)&s->cell,cells*sizeof(MFEvolvedCell)));TRY(cudaMalloc((void**)&s->evolved,MF_BATCH_TRIANGLES*sizeof(MFEvolvedBlock)));TRY(cudaMemcpy(s->generating,g,nc*sizeof(*g),cudaMemcpyHostToDevice));
    #undef TRY
    return 0;
}
extern "C" int mf_evolved_stage(void*ptr,const STCoefficient*c,size_t nc,const MFEvolvedCell*cell,size_t cells){
    auto*s=static_cast<Context*>(ptr);if(nc!=s->nc||cells!=s->cells||!s->generating)return cudaErrorInvalidValue;cudaError_t e=cudaMemcpy(s->c,c,nc*sizeof(*c),cudaMemcpyHostToDevice);if(e==cudaSuccess)e=cudaMemcpy(s->cell,cell,cells*sizeof(*cell),cudaMemcpyHostToDevice);return e;
}
extern "C" int mf_evolved_batch(void*ptr,int pass,size_t start,size_t count,MFEvolvedBlock*out,double*seconds){
    auto*s=static_cast<Context*>(ptr);if(pass<0||pass>2||!count||count>MF_BATCH_TRIANGLES||start+count>s->nt||!s->generating)return cudaErrorInvalidValue;cudaError_t e;
    #define TRY(call) do{e=(call);if(e!=cudaSuccess)return e;}while(0)
    TRY(cudaEventRecord(s->begin));evolved_kernel<<<count,128>>>(s->t,s->c,s->generating,s->rule,s->cell,pass,start,count,evolved_rule_view(s),s->evolved);TRY(cudaGetLastError());TRY(cudaEventRecord(s->end));TRY(cudaEventSynchronize(s->end));float ms;TRY(cudaEventElapsedTime(&ms,s->begin,s->end));*seconds=ms/1000.;TRY(cudaMemcpy(out,s->evolved,count*sizeof(*out),cudaMemcpyDeviceToHost));
    #undef TRY
    return 0;
}
extern "C" int mf_evolved_context_range(void*ptr,const MFEvolvedCell*cell,size_t first,size_t count){
    auto*s=static_cast<Context*>(ptr);if(!s||!s->generating||first>s->cells||count>s->cells-first)return cudaErrorInvalidValue;
    if(!count)return 0;
    return cudaMemcpy(s->cell+first,cell+first,count*sizeof(*cell),cudaMemcpyHostToDevice);
}
