#ifndef MATERIAL_COMPACT_THERMO_CUDA_H
#define MATERIAL_COMPACT_THERMO_CUDA_H
#include "material_compact_thermo.h"
#include <cuda_runtime.h>
// Original three-pass decomposition, with geometry-completed P3 and positive
// reference rho/e fields. Output retains the ORIGINAL MFEvolvedBlock ABI.
// pass0: Q1,QJ,QJdot,Qrho. pass1: heat/pressure/current moments.
// pass2: same moments plus residual-pressure/gravity forces. Host cell
// context uses CURRENT boundary V/Vdot and the SAME positive functional;
// mean-pressure/exterior force must be added from compact_volume.
// Generating coefficients/delta and projected e0 are immutable in an epoch.
// All buffers are caller-owned DEVICE buffers; launch is asynchronous.
extern "C" int compact_material_batch_device(const compact_geometry::Patch*,size_t,
    const STCoefficient*,const STCoefficient*,size_t,const compact_p3::ThermalDensity*,
    const MFEvolvedCell*,size_t,const uint64_t*,
    const compact_thermo::QuadraturePoint*,size_t,size_t,size_t,int,
    compact_thermo::Config,MFEvolvedBlock*,cudaStream_t);
#endif
