// Reuse the original native kernel; static geometry/points stay on device.
#include "material_geometric_pressure_cuda.cu"
struct GPPersistent {
    STTriangle *t=nullptr; STCoefficient *c=nullptr; STPoint *q=nullptr;
    GPInput *input=nullptr; GPOutput *out=nullptr;
    size_t nc=0,nq=0; cudaEvent_t begin=nullptr,end=nullptr;
};
extern "C" void gp_persistent_destroy(void *ptr) {
    auto *s=static_cast<GPPersistent*>(ptr); if(!s)return;
    cudaFree(s->t);cudaFree(s->c);cudaFree(s->q);
    cudaFree(s->input);cudaFree(s->out);
    if(s->begin)cudaEventDestroy(s->begin);
    if(s->end)cudaEventDestroy(s->end);delete s;
}
extern "C" int gp_persistent_create(const STTriangle*t,size_t nt,
        size_t nc,const STPoint*q,size_t nq,void **context) {
    *context=nullptr;auto *s=new GPPersistent;s->nc=nc;s->nq=nq;
    cudaError_t e=cudaSetDevice(0);
    #define GP_TRY(call) do{if(e==cudaSuccess)e=(call);}while(0)
    GP_TRY(cudaMalloc((void**)&s->t,nt*sizeof(*t)));
    GP_TRY(cudaMalloc((void**)&s->c,nc*sizeof(*s->c)));
    GP_TRY(cudaMalloc((void**)&s->q,nq*sizeof(*q)));
    GP_TRY(cudaMalloc((void**)&s->input,nq*sizeof(*s->input)));
    GP_TRY(cudaMalloc((void**)&s->out,nq*sizeof(*s->out)));
    GP_TRY(cudaMemcpy(s->t,t,nt*sizeof(*t),cudaMemcpyHostToDevice));
    GP_TRY(cudaMemcpy(s->q,q,nq*sizeof(*q),cudaMemcpyHostToDevice));
    GP_TRY(cudaEventCreate(&s->begin));GP_TRY(cudaEventCreate(&s->end));
    #undef GP_TRY
    if(e!=cudaSuccess){gp_persistent_destroy(s);return e;}
    *context=s;return 0;
}
extern "C" int gp_persistent_evaluate(void*context,const STCoefficient*c,
        const STPhysics*p,const GPInput*input,GPOutput*out,double*seconds) {
    auto*s=static_cast<GPPersistent*>(context);cudaError_t e=cudaSuccess;
    #define GP_TRY(call) do{e=(call);if(e!=cudaSuccess)return e;}while(0)
    GP_TRY(cudaMemcpy(s->c,c,s->nc*sizeof(*c),cudaMemcpyHostToDevice));
    GP_TRY(cudaMemcpy(s->input,input,s->nq*sizeof(*input),cudaMemcpyHostToDevice));
    GP_TRY(cudaEventRecord(s->begin));
    evaluate<<<(s->nq+127)/128,128>>>(s->t,s->c,s->q,s->nq,*p,s->input,s->out);
    GP_TRY(cudaGetLastError());GP_TRY(cudaEventRecord(s->end));
    GP_TRY(cudaEventSynchronize(s->end));float ms=0;
    GP_TRY(cudaEventElapsedTime(&ms,s->begin,s->end));*seconds=ms/1000.;
    GP_TRY(cudaMemcpy(out,s->out,s->nq*sizeof(*out),cudaMemcpyDeviceToHost));
    #undef GP_TRY
    return 0;
}
