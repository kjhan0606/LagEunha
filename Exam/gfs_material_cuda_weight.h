#ifndef GFS_MATERIAL_CUDA_WEIGHT_H
#define GFS_MATERIAL_CUDA_WEIGHT_H
#include "gfs_material_trace_csr2d.h"
#include "gfs_material_exterior.h"
#ifdef __cplusplus
extern "C" {
#endif
int gfs_material_cuda_exterior_pcg(int,int,int,const GfsMaterialExteriorEdge*,const double*,double*,int*,double*);
int gfs_material_cuda_weight_pcg(int,int,const GfsTraceDirected2D*,const unsigned char*,const double*,const double*,double*,int*,double*);
#ifdef __cplusplus
}
#endif
#endif
