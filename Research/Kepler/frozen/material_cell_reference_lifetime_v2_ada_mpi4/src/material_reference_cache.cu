// All original physical kernels and exact cell solve kernels are unchanged.
#include "material_cell_condensation_gpu.cu"
#include "material_reference_cache.h"
struct MFRPoint {GfsMaterialPoint12 p;int bad;};
static_assert(sizeof(MFRPoint)==112,"separate current point ABI");
struct MFRContext {
    Context*physical=nullptr;size_t point_capacity=0,double_capacity=0,total_doubles=0;
    bool resident=false;double*host=nullptr,*payload=nullptr;MFRPoint*point=nullptr;
    uint32_t*prefix=nullptr;uint64_t*double_prefix=nullptr;int*bad=nullptr;
    size_t start=0,count=0;uint64_t offset=0,bytes=0,generated_doubles=0;
};
__device__ static void mfr_shape(const STTriangle&t,const double*fixed,STJet*N){
    // Identical original P3 evaluation order. Enriched functions are restored
    // from their original double values/reference gradients, without rounding.
    STJet l[3];for(int j=0;j<3;j++)l[j]={fixed[j],t.grad[j][0],t.grad[j][1]};
    for(int a=0;a<10;a++){int index[3];gfs_material_p3_index(a,index);N[a]=STJet(1);for(int k=0;k<3;k++)for(int j=0;j<index[k];j++)N[a]=N[a]*(STJet(3)*l[k]-STJet(j))/STJet(index[k]-j);}
    for(int a=10;a<t.n;a++){size_t i=6+3*(a-10);N[a]={fixed[i],fixed[i+1],fixed[i+2]};}
}
__global__ static void mfr_generate_kernel(const STTriangle*t,const STCoefficient*generating,const GfsEnrichQuadrature*rule,
        const MFEvolvedCell*cell,size_t start,size_t count,const uint32_t*prefix,const uint64_t*double_prefix,MFRuleView view,double*payload,int*bad){
    size_t at=blockIdx.x;if(at>=count)return;size_t s=start+at;const auto&T=t[s];const auto&R=rule[s];
    __shared__ int edge[3],ne,nq;int tid=threadIdx.x;
    if(tid==0){bad[at]=0;ne=active_edges(T,R,edge);nq=evolved_count(view,s,ne>0?ne*R.ns*mf_config.na*mf_config.nr:0);}__syncthreads();
    size_t n=prefix[at+1]-prefix[at],stride=6+3*(T.n-10);
    if(ne<=0||size_t(nq)!=n||double_prefix[at+1]-double_prefix[at]!=n*stride){if(tid==0)bad[at]=1;return;}
    for(size_t qindex=tid;qindex<n;qindex+=blockDim.x){
        STPoint q;STJet N[ST_MAX];GfsMaterialPoint12 p0;int error=0;
        if(!evolved_point(T,R,s,qindex,edge,view,&q))error=2;
        else if(!st_shape(T,q,N)||!st_point(T,generating,N,&p0))error=3;
        else{
            double*out=payload+double_prefix[at]+qindex*stride,j0=p0.J+cell[T.cell].generating_delta;
            if(!(j0>0)||!isfinite(j0))error=4;
            else{for(int j=0;j<3;j++)out[j]=q.l[j];out[3]=q.weight;out[4]=q.rho;out[5]=j0;
                for(int a=10;a<T.n;a++){size_t i=6+3*(a-10);out[i]=N[a].v;out[i+1]=N[a].x;out[i+2]=N[a].y;}
                for(size_t i=0;i<stride;i++)if(!isfinite(out[i]))error=5;}
        }
        if(error)atomicExch(bad+at,error);
    }
}
__global__ static void mfr_current_points(const STTriangle*t,const STCoefficient*c,size_t start,size_t count,
        const uint32_t*prefix,const uint64_t*double_prefix,const double*payload,MFRPoint*cache){
    size_t at=blockIdx.x;if(at>=count)return;const auto&T=t[start+at];size_t stride=6+3*(T.n-10),n=prefix[at+1]-prefix[at];
    for(size_t qindex=threadIdx.x;qindex<n;qindex+=blockDim.x){auto&out=cache[prefix[at]+qindex];const double*fixed=payload+double_prefix[at]+qindex*stride;STJet N[ST_MAX];mfr_shape(T,fixed,N);out.bad=0;if(!st_point(T,c,N,&out.p))out.bad=3;}
}
// Generated from the verified original cached physical kernel: only its
// point loader and j0/weight/rho storage references differ.
#include "material_reference_cached_kernel.inc"
extern "C" void mfr_destroy(void*ptr){auto*s=static_cast<MFRContext*>(ptr);if(!s)return;cudaFree(s->payload);cudaFree(s->point);cudaFree(s->prefix);cudaFree(s->double_prefix);cudaFree(s->bad);delete s;}
extern "C" uint64_t mfr_device_bytes(void*ptr){auto*s=static_cast<MFRContext*>(ptr);return s?s->bytes:0;}
extern "C" int mfr_create(void*ctx,size_t points,size_t doubles,size_t total,int resident,double*host,void**out){
    *out=nullptr;auto*p=static_cast<Context*>(ctx);if(!p||!p->generating||!points||!doubles||!total||(!resident&&!host))return cudaErrorInvalidValue;
    auto*s=new MFRContext;s->physical=p;s->point_capacity=points;s->double_capacity=doubles;s->total_doubles=total;s->resident=resident;s->host=host;cudaError_t e=cudaSuccess;
    #define MFR_ALLOC(field,n) do{if(e==cudaSuccess){size_t bytes=(n)*sizeof(*s->field);e=cudaMalloc((void**)&s->field,bytes);if(e==cudaSuccess)s->bytes+=bytes;}}while(0)
    MFR_ALLOC(payload,resident?total:doubles);MFR_ALLOC(point,points);MFR_ALLOC(prefix,MF_BATCH_TRIANGLES+1);MFR_ALLOC(double_prefix,MF_BATCH_TRIANGLES+1);MFR_ALLOC(bad,MF_BATCH_TRIANGLES);
    #undef MFR_ALLOC
    if(e!=cudaSuccess){mfr_destroy(s);return e;}*out=s;return 0;
}
static int mfr_prefixes(MFRContext*s,Context*p,size_t start,size_t count,uint64_t offset,const uint32_t*prefix,const uint64_t*doubles){
    if(!s||s->physical!=p||!count||count>MF_BATCH_TRIANGLES||start+count>p->nt||prefix[0]||doubles[0]||prefix[count]>s->point_capacity||doubles[count]>s->double_capacity||offset>s->total_doubles||doubles[count]>s->total_doubles-offset)return cudaErrorInvalidValue;
    auto e=cudaMemcpy(s->prefix,prefix,(count+1)*sizeof(uint32_t),cudaMemcpyHostToDevice);if(e==cudaSuccess)e=cudaMemcpy(s->double_prefix,doubles,(count+1)*sizeof(uint64_t),cudaMemcpyHostToDevice);return e;
}
extern "C" int mfr_generate(void*ctx,void*ptr,size_t start,size_t count,uint64_t offset,const uint32_t*prefix,const uint64_t*doubles,double*seconds){
    auto*p=static_cast<Context*>(ctx);auto*s=static_cast<MFRContext*>(ptr);cudaError_t e;
    #define MFR_TRY(call) do{e=(call);if(e!=cudaSuccess)return e;}while(0)
    if(!p||!s||p->reduced_reference_only||offset!=s->generated_doubles)return cudaErrorInvalidValue;
    int error=mfr_prefixes(s,p,start,count,offset,prefix,doubles);if(error)return error;double*payload=s->payload+(s->resident?offset:0);
    MFR_TRY(cudaEventRecord(p->begin));mfr_generate_kernel<<<count,64>>>(p->t,p->generating,p->rule,p->cell,start,count,s->prefix,s->double_prefix,evolved_rule_view(p),payload,s->bad);MFR_TRY(cudaGetLastError());
    int bad[MF_BATCH_TRIANGLES];MFR_TRY(cudaMemcpy(bad,s->bad,count*sizeof(int),cudaMemcpyDeviceToHost));for(size_t i=0;i<count;i++)if(bad[i])return cudaErrorInvalidValue;
    if(!s->resident)MFR_TRY(cudaMemcpy(s->host+offset,s->payload,doubles[count]*sizeof(double),cudaMemcpyDeviceToHost));
    MFR_TRY(cudaEventRecord(p->end));MFR_TRY(cudaEventSynchronize(p->end));float ms;MFR_TRY(cudaEventElapsedTime(&ms,p->begin,p->end));*seconds=ms/1000.;s->generated_doubles+=doubles[count];return 0;
    #undef MFR_TRY
}
extern "C" int mfr_build(void*ctx,void*ptr,size_t start,size_t count,uint64_t offset,const uint32_t*prefix,const uint64_t*doubles,MFEvolvedBlock*out,double*seconds){
    auto*p=static_cast<Context*>(ctx);auto*s=static_cast<MFRContext*>(ptr);cudaError_t e;
    #define MFR_TRY(call) do{e=(call);if(e!=cudaSuccess)return e;}while(0)
    int error=mfr_prefixes(s,p,start,count,offset,prefix,doubles);if(error)return error;double*payload=s->payload+(s->resident?offset:0);
    MFR_TRY(cudaEventRecord(p->begin));if(!s->resident)MFR_TRY(cudaMemcpy(s->payload,s->host+offset,doubles[count]*sizeof(double),cudaMemcpyHostToDevice));
    mfr_current_points<<<count,64>>>(p->t,p->c,start,count,s->prefix,s->double_prefix,payload,s->point);MFR_TRY(cudaGetLastError());
    mfr_evolved_cached_kernel<<<count,128>>>(p->t,p->c,p->generating,p->rule,p->cell,s->point,s->prefix,s->double_prefix,payload,0,start,count,evolved_rule_view(p),p->evolved);MFR_TRY(cudaGetLastError());
    MFR_TRY(cudaEventRecord(p->end));MFR_TRY(cudaEventSynchronize(p->end));float ms;MFR_TRY(cudaEventElapsedTime(&ms,p->begin,p->end));*seconds=ms/1000.;MFR_TRY(cudaMemcpy(out,p->evolved,count*sizeof(*out),cudaMemcpyDeviceToHost));s->start=start;s->count=count;s->offset=offset;return 0;
    #undef MFR_TRY
}
extern "C" int mfr_batch(void*ctx,void*ptr,int pass,size_t start,size_t count,MFEvolvedBlock*out,double*seconds){
    auto*p=static_cast<Context*>(ctx);auto*s=static_cast<MFRContext*>(ptr);if(!s||s->physical!=p||pass<1||pass>2||s->start!=start||s->count!=count)return cudaErrorInvalidValue;cudaError_t e;
    #define MFR_TRY(call) do{e=(call);if(e!=cudaSuccess)return e;}while(0)
    double*payload=s->payload+(s->resident?s->offset:0);MFR_TRY(cudaEventRecord(p->begin));
    mfr_evolved_cached_kernel<<<count,128>>>(p->t,p->c,p->generating,p->rule,p->cell,s->point,s->prefix,s->double_prefix,payload,pass,start,count,evolved_rule_view(p),p->evolved);MFR_TRY(cudaGetLastError());
    MFR_TRY(cudaEventRecord(p->end));MFR_TRY(cudaEventSynchronize(p->end));float ms;MFR_TRY(cudaEventElapsedTime(&ms,p->begin,p->end));*seconds=ms/1000.;MFR_TRY(cudaMemcpy(out,p->evolved,count*sizeof(*out),cudaMemcpyDeviceToHost));return 0;
    #undef MFR_TRY
}

// Only the completed cache's owning context can enter reference-only mode.
// The caller additionally proves the same-stage direct witness on the host.
extern "C" int mfr_release_generation_inputs(void*ctx,void*ptr,uint64_t expected_nodes,uint64_t*freed){
    auto*p=static_cast<Context*>(ctx);auto*s=static_cast<MFRContext*>(ptr);
    if(!p||!s||s->physical!=p||!freed||!expected_nodes||p->reduced_reference_only||
       !p->reduced_active||!p->reduced_offset||!p->reduced_index||!p->reduced_weight||
       s->generated_doubles!=s->total_doubles)return cudaErrorInvalidValue;
    uint64_t nodes=0;auto e=cudaMemcpy(&nodes,p->reduced_offset+p->nt,sizeof(nodes),cudaMemcpyDeviceToHost);
    if(e!=cudaSuccess)return e;if(nodes!=expected_nodes)return cudaErrorInvalidValue;
    e=cudaDeviceSynchronize();if(e!=cudaSuccess)return e;
    // Block every original direct selected-point entry before releasing inputs.
    p->reduced_reference_only=true;
    e=cudaFree(p->reduced_index);if(e!=cudaSuccess)return e;p->reduced_index=nullptr;
    e=cudaFree(p->reduced_weight);if(e!=cudaSuccess)return e;p->reduced_weight=nullptr;
    *freed=nodes*(sizeof(uint32_t)+sizeof(double));return 0;
}
