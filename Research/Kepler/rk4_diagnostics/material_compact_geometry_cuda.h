#ifndef MATERIAL_COMPACT_GEOMETRY_CUDA_H
#define MATERIAL_COMPACT_GEOMETRY_CUDA_H
#include "material_compact_geometry.h"
#include <cuda_runtime.h>
namespace compact_geometry {
struct Sample {uint64_t patch;double l[3];};
struct SampleValue {GfsMaterialPoint12 map;double rho_reference=0,energy_reference=0;int valid=0;};
struct StressInput {Sample point;double weight,rho_reference,sigma[3],gravity[2];};
struct BoundaryInput {Sample point;double normal[2],weight,pressure;};
}
// Async caller-owned device buffers. component geometry DOFs share the SAME
// coefficient buffer as polynomial DOFs; they cannot be frozen by this API.
extern "C" int compact_geometry_sample_device(const compact_geometry::Patch*,size_t,
    const STCoefficient*,size_t,const compact_p3::ThermalDensity*,
    const compact_geometry::Sample*,size_t,compact_geometry::SampleValue*,cudaStream_t);
extern "C" int compact_geometry_stress_device(const compact_geometry::Patch*,size_t,
    const STCoefficient*,size_t,const compact_geometry::StressInput*,size_t,
    compact_geometry::StressValue*,cudaStream_t);
extern "C" int compact_geometry_exterior_device(const compact_geometry::Patch*,size_t,
    const STCoefficient*,size_t,const compact_geometry::BoundaryInput*,size_t,
    compact_geometry::BoundaryValue*,cudaStream_t);
#endif
