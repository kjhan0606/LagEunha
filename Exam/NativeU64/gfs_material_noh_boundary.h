#ifndef GFS_MATERIAL_NOH_BOUNDARY_H
#define GFS_MATERIAL_NOH_BOUNDARY_H
/* Include after eunha.h. Single-rank qualification of the existing square
 * Noh IC (rho0=1, inward speed=1, p0=1e-6). No periodic exterior aliases. */
#include "gfs_material_stage_time.h"
#include "../gfs_noh_exterior.h"
static inline int gfs_material_noh_native_mode(void){
 const char *gate=getenv("GFS_MATERIAL_NOH_NATIVE"),*ic=getenv("EUNHA_IC");
 return gate&&!strcmp(gate,"1")&&ic&&!strcmp(ic,"noh2d");
}
static inline int gfs_material_noh_configuration(SimParameters*s,GfsNohExterior*b){
 double t;int phase;double L=Xmax_HydroExam(s);
 if(!b||!gfs_material_noh_native_mode()||NID(s)!=1||NX(s)<16||NX(s)>2048||NX(s)%2||NY(s)!=NX(s)||
    VORO_NP(s)!=(ptrdiff_t)NX(s)*NY(s)||Xmin_HydroExam(s)!=0||Ymin_HydroExam(s)!=0||
    !(L>0)||!isfinite(L)||Ymax_HydroExam(s)!=L||!(GAS_GAMMA(s)>1)||
    !gfs_material_native_evaluation_time(&t,&phase))return 0;
 *b=(GfsNohExterior){2,{.5*L,.5*L,0},1,1,1e-6,GAS_GAMMA(s),t};return 1;
}
int gfs_material_native_boundary_power(double *power);
#endif
