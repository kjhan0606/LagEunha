#include <stdio.h>
#include <string.h>
#include <math.h>
#include "../gfs_material_mesh2d.h"
int main(void){
 const int cells[2][4]={{0,1,4,3},{1,2,5,4}};
 double x[6][2]={{0,0},{1,0},{2,0},{0,1},{1,1},{2,1}},v[6][2],trial[6][2];
 const double mass[2]={1,2};double dt=.001,maxerr=0;
 for(int step=0;step<1000;step++){
  /* Exact steady shear: affine map, fixed volume, shared material boundary.
     Prescribed-motion geometry control, NOT a hydrodynamic solution. */
  for(int i=0;i<6;i++){v[i][0]=x[i][1];v[i][1]=0;for(int a=0;a<2;a++)trial[i][a]=x[i][a]+dt*v[i][a];}
  for(int c=0;c<2;c++){
   if(!gfs_mm_convex_step(4,cells[c],x,v,dt))return 1;
   double sweep=0;for(int k=0;k<4;k++)sweep+=gfs_mm_edge_sweep(cells[c][k],cells[c][(k+1)%4],x,v,dt);
   double a=gfs_mm_area(4,cells[c],x),b=gfs_mm_area(4,cells[c],trial);
   if(!(b>0))return 2;
   maxerr=fmax(maxerr,fabs(b-a-sweep));
   if(fabs(mass[c]/b-mass[c])>1e-11)return 3;
  }
  if(gfs_mm_edge_sweep(1,4,x,v,dt)+gfs_mm_edge_sweep(4,1,x,v,dt)!=0)return 4;
  memcpy(x,trial,sizeof(x));
 }
 /* Nonzero compression/rotation: test quadratic area change too. */
 for(int i=0;i<6;i++){
  v[i][0]=-.2*x[i][0]+.3*x[i][1];
  v[i][1]=-.1*x[i][0]-.15*x[i][1];
  for(int a=0;a<2;a++)trial[i][a]=x[i][a]+.1*v[i][a];
 }
 for(int c=0;c<2;c++){
  if(!gfs_mm_convex_step(4,cells[c],x,v,.1))return 6;
  double sweep=0;for(int k=0;k<4;k++)sweep+=gfs_mm_edge_sweep(cells[c][k],cells[c][(k+1)%4],x,v,.1);
  maxerr=fmax(maxerr,fabs(gfs_mm_area(4,cells[c],trial)-gfs_mm_area(4,cells[c],x)-sweep));
 }
 /* Deliberately invert a cell: reject before modifying accepted geometry. */
 memset(v,0,sizeof(v));v[4][1]=-3;
 if(gfs_mm_convex_step(4,cells[0],x,v,1))return 5;
 printf("PASS 1000 shared-vertex shear geometry steps to t=1; area-sweep error %.17g; inversion rejected; masses unchanged\n",maxerr);
 return maxerr>1e-12;
}
