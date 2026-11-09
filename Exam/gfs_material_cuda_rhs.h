#ifndef GFS_MATERIAL_CUDA_RHS_H
#define GFS_MATERIAL_CUDA_RHS_H
#include "gfs_material_trace_csr2d.h"
#include "gfs_material_face_response.h"
#include "gfs_material_face_normal2d.h"
#ifdef __cplusplus
extern "C" {
#endif
/* Same status as CPU: 1 valid, 2 vacuum, 0 invalid/backend failure. */
int gfs_material_cuda_rhs2d(const GfsTraceCsr2D*,const GfsFvPrimitive*,double,const double*,double*);
/* Optional prescribed wall speeds; trace limiter must use this same array.
 * Output scalar records total external wall power, separate from H. */
int gfs_material_cuda_rhs_moving_walls2d(const GfsTraceCsr2D*,const GfsFvPrimitive*,double,const double*,const double*,double*,double*);
/* Optional same-call face response output; never reused across states. */
int gfs_material_cuda_rhs_cache2d(const GfsTraceCsr2D*,const GfsFvPrimitive*,double,const double*,double*,GfsMaterialFaceResponse*);
/* Direct-FV owner pressure shifts after reconstruction; NULL preserves old path. */
int gfs_material_cuda_owner_rhs_cache2d(const GfsTraceCsr2D*,const GfsFvPrimitive*,double,const double*,double*,GfsMaterialFaceResponse*,const double*);
/* Batched frozen linear actions; output rows are contiguous (Dg, C^Tg). */
int gfs_material_cuda_actions2d(const GfsTraceCsr2D*,const GfsMaterialFaceResponse*,const double*,int,const double*,double*,double*);
#ifdef __cplusplus
}
#endif
#endif
