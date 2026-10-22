#ifndef MATERIAL_FULL_GPU_MASS_H
#define MATERIAL_FULL_GPU_MASS_H
#include "material_full_initial_batch.h"
extern "C" int mf_mass_create(const STTriangle*,size_t,size_t,const MFDoubleDouble*,
    const MFDoubleDouble*,const MFDoubleDouble*,const uint64_t*,const uint32_t*,
    const uint8_t*,size_t,void**);
extern "C" int mf_mass_apply(void*,const MFDoubleDouble*,MFDoubleDouble*,int,double*);
extern "C" void mf_mass_destroy(void*);
#endif
