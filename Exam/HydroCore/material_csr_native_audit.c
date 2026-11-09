#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../exam_gpu.h"
#include "../gfs_fv_native.h"
#include "../gfs_material_csr2d.h"
#include "../gfs_material_csr_native.inc"
int main(void){
 for(int mode=0;mode<3;mode++){
  int n=8,nc=n*n;GfsFvMesh m={0};m.dim=2;m.n=nc;m.L[0]=m.L[1]=1;m.gamma=1.4;m.bc[0]=m.bc[1]=GFS_FV_REFLECT;m.cell=calloc(nc,sizeof(*m.cell));
  for(int j=0;j<n;j++)for(int i=0;i<n;i++){
   int c=j*n+i;m.cell[c].q[0]=(i+.5+(mode?.07*sin(3*c):0))/n;m.cell[c].q[1]=(j+.5+(mode?.07*cos(5*c):0))/n;m.cell[c].W=mode==2?.001*sin(c)/(n*n):0;
  }
  if(gfs_fv_geometry(&m))return 1;
  FaceCSR f={0};ParticleSoA p={0};f.n_particles=nc;p.n_total=nc;
  f.face_offset=calloc(nc+1,sizeof(int));p.x=malloc(nc*sizeof(double));p.y=malloc(nc*sizeof(double));
  for(int c=0;c<nc;c++){f.face_offset[c+1]=f.face_offset[c]+m.geometry[c].faces;p.x[c]=m.cell[c].q[0];p.y[c]=m.cell[c].q[1];}
  int nf=f.n_faces_total=f.face_offset[nc];
  f.c1x=malloc(nf*sizeof(double));f.c1y=malloc(nf*sizeof(double));f.c2x=malloc(nf*sizeof(double));f.c2y=malloc(nf*sizeof(double));f.neighbor_idx=malloc(nf*sizeof(int));f.is_ghost=malloc(nf*sizeof(int));
  for(int c=0;c<nc;c++)for(int j=0;j<m.geometry[c].faces;j++){
   GfsFvNativeFace *face=&m.geometry[c].face[j];int k=f.face_offset[c]+j;
   if(face->count!=3)return 2;
   GfsFvPoint *pt=&m.geometry[c].point[face->point];
   double dx=(pt[2].x[0]-pt[0].x[0])/sqrt(.6),dy=(pt[2].x[1]-pt[0].x[1])/sqrt(.6);
   f.c1x[k]=pt[1].x[0]-.5*dx;f.c2x[k]=pt[1].x[0]+.5*dx;f.c1y[k]=pt[1].x[1]-.5*dy;f.c2y[k]=pt[1].x[1]+.5*dy;
   f.neighbor_idx[k]=face->j;f.is_ghost[k]=face->mirror!=0;
  }
  GfsMaterialCsr2D imported={0};if(!gfs_material_import_facecsr(&f,&p,&imported)){fprintf(stderr,"import failed mode %d\n",mode);return 3;}
  double err=0,area=0;
  for(int c=0;c<nc;c++){double a=gfs_mm_area(imported.offset[c+1]-imported.offset[c],imported.vertex+imported.offset[c],imported.x);area+=a;err=fmax(err,fabs(a-m.geometry[c].V));}
  printf("mode %d: native cells %d shared nodes %d total area %.17g maximum cell area error %.3g merged coordinate distance %.3g\n",mode,nc,imported.nodes,area,err,imported.max_merge_distance);
  if(fabs(area-1)>1e-12||err>1e-12)return 4;
  gfs_material_csr_free(&imported);
  /* A falsely declared boundary on a shared interior edge must be rejected. */
  int badface=-1;for(int k=0;k<nf;k++)if(!f.is_ghost[k]){badface=k;break;}
  if(badface<0)return 5;
  f.is_ghost[badface]=1;
  if(gfs_material_import_facecsr(&f,&p,&imported))return 6;
  if(imported.x||imported.vertex||imported.offset)return 7;
  f.is_ghost[badface]=0;
  free(f.face_offset);free(f.c1x);free(f.c1y);free(f.c2x);free(f.c2y);free(f.neighbor_idx);free(f.is_ghost);free(p.x);free(p.y);gfs_fv_geometry_free(&m);free(m.cell);
 }
 return 0;
}
