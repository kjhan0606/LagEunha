#ifndef GFS_MATERIAL_CUDA_PROJECT_H
#define GFS_MATERIAL_CUDA_PROJECT_H
#include "gfs_material_cuda_rhs.h"
#ifdef __cplusplus
extern "C" {
#endif
/* One cell/scalar projection function shared by C and CUDA. No state cache. */
/* Only n_owned rows are computed; all n_sample physical states are uploaded.
 * Prescribed Riemann face slots are not gradient samples. Failure leaves no
 * usable output; callers must discard the entire batch. */
int gfs_cuda_project_gradients_samples(int n_owned,int n_sample,int directed_faces,
 const int*,const GfsTraceDirected2D*,const GfsTraceDirected2D*,
 const GfsFvPrimitive*,const unsigned char*,double*);
/* Closed/periodic API: n_sample == n_owned. */
int gfs_cuda_project_gradients(int,int,const int*,const GfsTraceDirected2D*,
 const GfsTraceDirected2D*,const GfsFvPrimitive*,const unsigned char*,double*);
#ifdef __cplusplus
}
#endif
#endif
