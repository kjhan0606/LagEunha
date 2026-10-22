#ifndef MATERIAL_CELL_CONDENSATION_GPU_H
#define MATERIAL_CELL_CONDENSATION_GPU_H
#include "material_full_mass_resident.h"
struct MCCell {uint64_t I,G,L,IG;uint32_t ni,ng;};
extern "C" int mcc_create(size_t,size_t,const MCCell*,size_t,const uint32_t*,size_t,const uint32_t*,size_t,
    const uint32_t*,const uint8_t*,const uint64_t*,const uint64_t*,const MFDoubleDouble*,size_t,const MFDoubleDouble*,size_t,void**);
extern "C" void mcc_destroy(void*);
extern "C" uint64_t mcc_device_bytes(void*);
extern "C" int mcc_rhs(void*,const MFDoubleDouble*,const MFDoubleDouble*,MFDoubleDouble*,MFDoubleDouble*);
extern "C" int mcc_recover(void*,const MFDoubleDouble*,const MFDoubleDouble*,MFDoubleDouble*,MFDoubleDouble*);
#endif
