#ifndef GFS_CYLINDER_MATERIAL_BOUNDARY_H
#define GFS_CYLINDER_MATERIAL_BOUNDARY_H
/* Include after eunha.h. Boundary state belongs to RunCylinder, not CUDA. */
static inline int gfs_material_cylinder_mode(void){
 const char *a=getenv("GFS_MATERIAL_CYLINDER_NATIVE"),*b=getenv("EUNHA_IC");
 return a&&!strcmp(a,"1")&&b&&!strcmp(b,"cylinder");
}
typedef struct {
 int n,nt,*alias,*image;
 long long *owned_ids;
 double *lowx,*lowy,center[2],radius,lx,ly;
} GfsCylinderMaterialSites;
const GfsCylinderMaterialSites *gfs_material_cylinder_sites(int,int,const long long*);
void gfs_material_cylinder_padding(SimParameters*,double);
int gfs_material_cylinder_search(SimParameters*,double,double*,double*,double*,double*);
#ifdef GFS_MATERIAL_SEARCH_EXTENT
int gfs_material_cylinder_search_extent(SimParameters*,double*,double*,double*,double*,const double*);
#endif
void gfs_material_cylinder_accept(double,double,double);
/* Owned file must already be complete and represent this accepted time. */
void gfs_material_cylinder_checkpoint_boundary(SimParameters*,const char*,double);
int gfs_material_cylinder_refresh(SimParameters*,double);
int gfs_material_native_wall_power(double*);
void cyl_postStage(SimParameters*);
#endif
