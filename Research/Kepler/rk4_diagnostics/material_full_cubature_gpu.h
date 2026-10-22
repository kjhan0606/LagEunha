#ifndef MATERIAL_FULL_CUBATURE_GPU_H
#define MATERIAL_FULL_CUBATURE_GPU_H
#include "material_full_initial_batch.h"
extern "C" int mf_cubature_reserve(void*,size_t);
extern "C" int mf_cubature_tile(void*);
extern "C" int mf_cubature_rows(void*,size_t,size_t,size_t,MFDoubleDouble*,double*);
extern "C" int mf_cubature_weights(void*,size_t,double*,size_t);
#endif
