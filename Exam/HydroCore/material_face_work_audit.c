#include <stdio.h>
#include <math.h>
#include "../gfs_face_work.h"
static double area(double x[4][2]) {
 double a=0; for(int i=0;i<4;i++){int j=(i+1)%4;a+=x[i][0]*x[j][1]-x[j][0]*x[i][1];}return .5*a;
}
int main(void){
 /* Independent shoelace derivative versus integrated edge velocity for a
    convex, nonrectangular cell under affine motion. */
 double x[4][2]={{0,0},{1.2,.1},{1,1.1},{-.2,.8}};
 double v[4][2],xp[4][2],xm[4][2];
 double h=1e-5,flux=0,power=0,P=2.3;
 for(int i=0;i<4;i++){
  v[i][0]=.7+.2*x[i][0]+.6*x[i][1];
  v[i][1]=-.4-.3*x[i][0]-.1*x[i][1];
  for(int k=0;k<2;k++){xp[i][k]=x[i][k]+h*v[i][k];xm[i][k]=x[i][k]-h*v[i][k];}
 }
 for(int i=0;i<4;i++){
  int j=(i+1)%4;double sx=x[j][1]-x[i][1],sy=x[i][0]-x[j][0];
  double wx=.5*(v[i][0]+v[j][0]),wy=.5*(v[i][1]+v[j][1]);
  flux+=gfs_face_sweep2d(sx,sy,wx,wy);
  power+=gfs_face_work2d(P,sx,sy,0,0,wx,wy);
  double a=gfs_face_work2d(P,sx,sy,.1,-.2,wx,wy);
  double b=gfs_face_work2d(P,-sx,-sy,-.1,.2,wx,wy);
  if(fabs(a+b)>1e-14)return 1;
 }
 double fd=(area(xp)-area(xm))/(2*h);
 if(fabs(fd-flux)>1e-9||fabs(power+P*fd)>1e-9)return 2;
 /* Equal and opposite sweep errors cancel volume, not pressure-weighted work. */
 double defect=gfs_face_work_defect2d(2,1,0,0,0,.1,0,0,0)
              +gfs_face_work_defect2d(1,-1,0,0,0,.1,0,0,0);
 if(fabs(defect+.1)>1e-14)return 3;
 printf("PASS shared-face cancellation; area derivative error %.3g; pressure work error %.3g; nonuniform-pressure defect %.17g\n",fabs(fd-flux),fabs(power+P*fd),defect);
 return 0;
}
