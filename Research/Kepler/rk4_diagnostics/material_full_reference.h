#ifndef MATERIAL_FULL_REFERENCE_H
#define MATERIAL_FULL_REFERENCE_H
#include <stdint.h>
// Native geometry input to material state construction, not a particle restart.
typedef struct {
    uint64_t id,n;
    double x,y,mass,volume,vx,vy,U,P;
} MFCell;
typedef struct {
    uint64_t neighbor;
    double nx,ny,ax,ay,bx,by;
} MFFace;
#ifdef __cplusplus
extern "C" {
#endif
int gfs_material_prepare_reference(void *gpu_context,int np,int fortran_comm,
    double gamma,double pressure,double eps,double box);
int gfs_material_initialize_full(const char *retained_reference,int fortran_comm,
    double gamma,double pressure,double eps,double box);
#ifdef __cplusplus
}
#endif
#endif
