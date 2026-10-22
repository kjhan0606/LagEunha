#ifndef MATERIAL_COMPACT_P3_CUDA_H
#define MATERIAL_COMPACT_P3_CUDA_H
#include "material_compact_p3.h"
#include <cuda_runtime.h>
namespace compact_p3 {
struct Sample {uint64_t triangle;double l[3];};
struct SampleValue {
    GfsMaterialPoint12 map;
    double rho_reference,energy_reference;
    int valid;
};
}
// All pointers are caller-owned DEVICE buffers, belonging to one immutable
// representation epoch. Asynchronous launch; caller checks stream completion.
// Values are reference densities, not a new entropy/EOS or a force operator.
extern "C" int compact_p3_sample_device(const compact_p3::Triangle*,size_t,
    const STCoefficient*,const compact_p3::ThermalDensity*,
    const compact_p3::Sample*,size_t,compact_p3::SampleValue*,cudaStream_t);
#endif
