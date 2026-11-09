/* Runs the original Cylinder padding. Finite differences of that code are
 * compared with the shared kinematic adapter, not a duplicate ghost generator.
 * No time integration or change to a physical state is performed here. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "eunha.h"
#include "voro.h"
#include "nnost.h"
#include "gnnost.h"
#include "exam.h"
#include "exam2d.h"
#include "Cylinder/cylinder.h"
#include "gfs_material_cylinder_boundary2d.h"
#include "gfs_material_power_face_motion.h"
#include "gfs_material_face_response.h"

void *my_malloc(size_t n){return malloc(n);}
void my_free(void*p){free(p);}
void paddingTreeVorork4Particles(SimParameters*s,postype width){
 (void)s;(void)width;fprintf(stderr,"Unexpected MPI padding in serial audit\n");abort();
}
#define CHECK(x) do{if(!(x)){fprintf(stderr,"FAIL %d: %s\n",__LINE__,#x);return 1;}}while(0)
static int sample(SimParameters*s,double h,CylinderPaddingOwner**map,
 treevorostressrk4particletype**ghost,int*count){
 treevorostressrk4particletype*p=(void*)VORORK4_TBP(s);int n=VORO_NP(s);
 for(int i=0;i<n;i++){p[i].x+=h*p[i].vx;p[i].y+=h*p[i].vy;p[i].weight.W+=h*(.03+.001*i);}
 paddingCylinderParticlesMapped(s,.2,map,count);
 *ghost=(void*)VORORK4_TBPP(s);
 for(int i=0;i<n;i++){p[i].x-=h*p[i].vx;p[i].y-=h*p[i].vy;p[i].weight.W-=h*(.03+.001*i);}
 VORORK4_TBPP(s)=NULL;VORO_NPAD(s)=0;return 1;
}
int main(void){
 SimParameters*s=calloc(1,sizeof(*s));CHECK(s);
 TVORORK4_DDINFO(s)=calloc(1,sizeof(*TVORORK4_DDINFO(s)));CHECK(TVORORK4_DDINFO(s));
 TVORORK4_DDINFO(s)[0].n_size=sizeof(treevorostressrk4particletype);
 NID(s)=1;NDDINFO(s)=0;MPI_COMM(s)=MPI_COMM_WORLD;
 CYL_CX(s)=CYL_CY(s)=4;CYL_R(s)=1;CYL_RHO(s)=1;CYL_P(s)=.08;CYL_UINF(s)=1;
 GAS_GAMMA(s)=1.4;SIMBOX(s).x.max=SIMBOX(s).y.max=8;
 const double xy[][2]={{5.1,4.1},{4.1,5.1},{2.9,3.9},{3.9,2.9},
  {.1,.1},{7.9,.1},{.1,7.9},{7.9,7.9},{4,.1},{4,7.9},{.1,4},{7.9,4},{6,6}};
 int n=sizeof(xy)/sizeof(xy[0]);VORO_NP(s)=VORO_TNP(s)=n;
 treevorostressrk4particletype*p=calloc(n,sizeof(*p));CHECK(p);VORORK4_TBP(s)=(void*)p;
 for(int i=0;i<n;i++){
  p[i].x=xy[i][0];p[i].y=xy[i][1];p[i].vx=.3+i*.02;p[i].vy=.4+i*.03;
  p[i].mass=.01;p[i].volume=.01;p[i].pressure=.08;p[i].den=1;
  p[i].u4if.indx=19+7*i;p[i].weight.W=.01+.001*i;
 }
 CylinderPaddingOwner *map;treevorostressrk4particletype*g;int ng;
 CHECK(sample(s,0,&map,&g,&ng));int counts[5]={0},corners=0;
 GfsMaterialCylinderImage2D*image=calloc(ng,sizeof(*image));CHECK(image);
 double center[2]={4,4},max_position=0,max_velocity=0,max_weight=0;
 for(int k=0;k<ng;k++){
  CylinderPaddingOwner*m=map+k;int j=m->source_index;CHECK(j>=0&&j<n);
  double x[2]={p[j].x,p[j].y},v[2]={p[j].vx,p[j].vy},gv[2]={g[k].vx,g[k].vy};
  CHECK(gfs_material_cylinder_image2d(m,n,(long long)PINDX(p+j),x,v,p[j].weight.W,.03+.001*j,gv,center,1,8,8,image+k));
  for(int a=0;a<2;a++)max_position=fmax(max_position,fabs(image[k].position[a]+image[k].low[a]-(a?g[k].y:g[k].x)));
  CHECK(image[k].weight==g[k].weight.W&&image[k].weight_owner==j);
  CHECK((image[k].wall_owner>=0)==(m->kind==CYL_IMAGE_RADIAL));
  CHECK((image[k].periodic_owner>=0)==(m->kind==CYL_IMAGE_PERIODIC));
  CHECK(image[k].fluid_velocity[0]==g[k].vx&&image[k].fluid_velocity[1]==g[k].vy);
  if(m->kind>=CYL_IMAGE_INFLOW){CHECK(image[k].mesh_velocity[0]==-p[j].vx);if(m->y_shift)corners++;}
  counts[m->kind]++;
  GfsMaterialCylinderImage2D unused;
  CHECK(!gfs_material_cylinder_image2d(m,n,(long long)PINDX(p+j)+1,x,v,p[j].weight.W,0,gv,center,1,8,8,&unused));
 }
 const double hs[]={1e-4,1e-5,1e-6};
 for(int hidx=0;hidx<3;hidx++){
  double h=hs[hidx];CylinderPaddingOwner *mp,*mm;treevorostressrk4particletype *gp,*gm;int np,nm;
  CHECK(sample(s,h,&mp,&gp,&np)&&sample(s,-h,&mm,&gm,&nm)&&np==ng&&nm==ng);
  for(int k=0;k<ng;k++){
   CHECK(mp[k].source_id==map[k].source_id&&mm[k].source_id==map[k].source_id&&mp[k].kind==map[k].kind&&mm[k].kind==map[k].kind);
   max_velocity=fmax(max_velocity,fabs((gp[k].x-gm[k].x)/(2*h)-image[k].mesh_velocity[0]));
   max_velocity=fmax(max_velocity,fabs((gp[k].y-gm[k].y)/(2*h)-image[k].mesh_velocity[1]));
   max_weight=fmax(max_weight,fabs((gp[k].weight.W-gm[k].weight.W)/(2*h)-image[k].weight_rate));
  }
  free(mp);free(mm);free(gp);free(gm);
 }
 CHECK(corners==4&&counts[1]>0&&counts[2]>0&&counts[3]>0&&counts[4]>0);
 CHECK(max_position<1e-13&&max_velocity<1e-8&&max_weight<1e-9);

 /* A stationary open mirror is not a moving material boundary. Test the
  * actual right ghost and the common geometric face-rate formula. */
 int k=-1;for(int a=0;a<ng;a++)if(map[a].kind==CYL_IMAGE_OUTFLOW&&map[a].y_shift==0){k=a;break;}CHECK(k>=0);
 int j=map[k].source_index;double v[2]={p[j].vx,p[j].vy};
 double d[2]={image[k].position[0]-p[j].x,0},c[2]={.5*d[0],0},face_rate;
 GfsMaterialPowerFaceMotion f;CHECK(gfs_material_power_face_motion(2,d,c,1,&f));
 CHECK(gfs_material_power_face_motion_rate(2,&f,v,image[k].mesh_velocity,.03+.001*j,image[k].weight_rate,&face_rate));
 CHECK(fabs(face_rate)<1e-14);
 double target[2];
 for(int different=0;different<2;different++){
  GfsFvPrimitive a={1,.08,{1,0,0}},b=a;GfsMaterialFaceResponse left,right;
  double nl[3]={-1,0,0},nr[3]={1,0,0};
  CHECK(gfs_material_face_response(2,a,b,nl,1.4,&left)==1);
  a.u[0]=b.u[0]=different?1.2:1;
  CHECK(gfs_material_face_response(2,a,b,nr,1.4,&right)==1);
  target[different]=left.ustar+right.ustar;
 }
 CHECK(fabs(target[0])<1e-14&&fabs(target[1]-.2)<1e-14);
 printf("{\"ghosts\":%d,\"kinds\":[%d,%d,%d,%d,%d],\"open_corners\":%d,\"max_position_error\":%.17g,\"max_coordinate_derivative_error\":%.17g,\"max_weight_derivative_error\":%.17g,\"stationary_open_face_rate\":%.17g,\"balanced_total_volume_target\":%.17g,\"unbalanced_total_volume_target\":%.17g,\"source_id_mismatch_rejected\":true,\"time_evolution\":false}\n",ng,counts[0],counts[1],counts[2],counts[3],counts[4],corners,max_position,max_velocity,max_weight,face_rate,target[0],target[1]);
 free(image);free(map);free(g);free(p);free(TVORORK4_DDINFO(s));free(s);return 0;
}
