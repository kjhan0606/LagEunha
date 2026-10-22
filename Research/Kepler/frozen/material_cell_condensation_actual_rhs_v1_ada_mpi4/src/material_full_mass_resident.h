#ifndef MATERIAL_FULL_MASS_RESIDENT_H
#define MATERIAL_FULL_MASS_RESIDENT_H
#include "material_full_initial_batch.h"
enum {MF_CG_X=0,MF_CG_R=1,MF_CG_Z=2,MF_CG_P=3,MF_CG_AP=4};
extern "C" int mf_resident_create(void*,const uint64_t*,size_t,const uint64_t*,size_t,const uint64_t*,size_t,void**);
extern "C" void mf_resident_destroy(void*);
extern "C" int mf_resident_begin(void*,const MFDoubleDouble*);
extern "C" int mf_resident_product(void*,int);
extern "C" int mf_resident_pack(void*,int,int,MFDoubleDouble*);
extern "C" int mf_resident_scatter(void*,int,int,size_t,size_t,int,const MFDoubleDouble*);
extern "C" int mf_resident_update(void*,int,MFDoubleDouble);
extern "C" int mf_resident_dot(void*,int,int,MFDoubleDouble*);
extern "C" int mf_resident_finish(void*,MFDoubleDouble*);
#endif
