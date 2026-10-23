#ifndef MATERIAL_COMPACT_GPU_H
#define MATERIAL_COMPACT_GPU_H
#include "material_full_mass_resident.h"
extern "C" int cm_create(size_t,const uint64_t*,const uint32_t*,const MFDoubleDouble*,const MFDoubleDouble*,
    const uint64_t*,size_t,const uint64_t*,size_t,const uint64_t*,size_t,void**);
extern "C" void cm_destroy(void*);
extern "C" int cm_begin(void*,const MFDoubleDouble*);
extern "C" int cm_product(void*,int);
extern "C" int cm_pack(void*,int,int,MFDoubleDouble*);
extern "C" int cm_scatter(void*,int,int,size_t,size_t,int,const MFDoubleDouble*);
extern "C" int cm_update(void*,int,MFDoubleDouble);
extern "C" int cm_dot(void*,int,int,MFDoubleDouble*);
extern "C" int cm_finish(void*,MFDoubleDouble*);
#endif
