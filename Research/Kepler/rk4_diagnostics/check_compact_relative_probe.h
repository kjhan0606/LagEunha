#ifndef CHECK_COMPACT_RELATIVE_PROBE_H
#define CHECK_COMPACT_RELATIVE_PROBE_H
#include "material_compact_geometry.h"
#include <cuda_runtime.h>
struct RelativeInput {double t,lo,hi,B[2],C[2],l[3];};
struct RelativeOutput {STJet relative[2],mode[2];int valid;};
extern "C" int compact_relative_probe(const compact_geometry::Patch*,const RelativeInput*,size_t,RelativeOutput*);
#endif
