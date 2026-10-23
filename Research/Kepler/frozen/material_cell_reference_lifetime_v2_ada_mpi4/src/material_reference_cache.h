#ifndef MATERIAL_REFERENCE_CACHE_H
#define MATERIAL_REFERENCE_CACHE_H
#include "material_full_initial_batch.h"
extern "C" int mfr_create(void*,size_t,size_t,size_t,int,double*,void**);
extern "C" void mfr_destroy(void*);
extern "C" int mfr_release_generation_inputs(void*,void*,uint64_t,uint64_t*);
extern "C" uint64_t mfr_device_bytes(void*);
extern "C" int mfr_generate(void*,void*,size_t,size_t,uint64_t,const uint32_t*,const uint64_t*,double*);
extern "C" int mfr_build(void*,void*,size_t,size_t,uint64_t,const uint32_t*,const uint64_t*,MFEvolvedBlock*,double*);
extern "C" int mfr_batch(void*,void*,int,size_t,size_t,MFEvolvedBlock*,double*);
namespace reference_cache {
struct Binding {
    void*self=nullptr;
    int(*build)(void*,void*,size_t,size_t,const uint32_t*,MFEvolvedBlock*,double*)=nullptr;
    int(*batch)(void*,void*,int,size_t,size_t,MFEvolvedBlock*,double*)=nullptr;
};
inline Binding active;
inline int build(void*ctx,size_t first,size_t count,const uint32_t*prefix,MFEvolvedBlock*out,double*seconds){return active.self?active.build(active.self,ctx,first,count,prefix,out,seconds):mf_evolved_cache_build(ctx,first,count,prefix,out,seconds);}
inline int batch(void*ctx,int pass,size_t first,size_t count,MFEvolvedBlock*out,double*seconds){return active.self?active.batch(active.self,ctx,pass,first,count,out,seconds):mf_evolved_cache_batch(ctx,pass,first,count,out,seconds);}
}
#endif
