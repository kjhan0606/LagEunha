/* Deterministic geometric gates; no hydrodynamic verification is implied. */
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "laguerre.h"
static int check(const char *name,double value,double expected) {
  double error=fabs(value-expected)/expected;
  printf("%s volume=%.17g expected=%.17g relative=%.17g\n",name,value,expected,error);fflush(stdout);
  return !isfinite(value)||error>2e-12;
}
int main(int argc,char **argv) {
  Vertex3D v[1024];int failures=0;
  if(argc>1 && strcmp(argv[1],"contracts")==0) {
    Point3D center={0},input[6]={0},neighbors[6]={0};int work[1024];
    for(int i=0;i<6;i++) { input[i].index=i+1;double sign=i%2?1:-1;
      if(i/2==0) input[i].x=sign;if(i/2==1) input[i].y=sign;if(i/2==2) input[i].z=sign;
    }
    input[1].weightSquared=1.5;
    int n=constructLaguerreCell3D(&center,input,neighbors,6,v,1024,4,0,work);
    if(n<4) return 1;
    failures+=check("generator outside cell",computePolyhedronVolume3D(v,n),.25);
    input[1].weightSquared=3.;
    n=constructLaguerreCell3D(&center,input,neighbors,6,v,1024,4,0,work);
    printf("empty cell return=%d expected=0\n",n);failures+=(n!=0);
    Vertex3D small[8];Point3D one={0};one.index=1;one.x=one.y=one.z=1.;
    n=constructLaguerreCell3D(&center,&one,neighbors,1,small,8,4,0,work);
    printf("insufficient capacity return=%d expected=-1\n",n);failures+=(n!=-1);
    one.x=one.y=one.z=0.;
    n=constructLaguerreCell3D(&center,&one,neighbors,1,v,1024,4,0,work);
    printf("coincident site return=%d expected=-1\n",n);failures+=(n!=-1);
    return failures?1:0;
  }
  if(argc>1 && strcmp(argv[1],"clipped")==0) {
    Point3D center={0},input[6]={0},neighbors[6]={0};int work[1024];
    for(int i=0;i<6;i++) { input[i].index=i+1;double sign=i%2?1:-1;
      if(i/2==0) input[i].x=sign;if(i/2==1) input[i].y=sign;if(i/2==2) input[i].z=sign;
    }
    int n=constructLaguerreCell3D(&center,input,neighbors,6,v,1024,4,0,work);
    printf("clipped vertices=%d\n",n);fflush(stdout);
    return check("six-neighbor cell",computePolyhedronVolume3D(v,n),1.);
  }
  for(int k=0;k<3;k++) {
    double L=k==0?1.:(k==1?2.:.125);int n=initBoundingBox3D(v,1024,L);
    failures+=check("cube",computePolyhedronVolume3D(v,n),L*L*L);
    failures+=check("origin face pyramid",computePyramidVolume3D(v,v[0].adjacent[0]),L*L*L/6);
    for(int i=0;i<n;i++) { double x=v[i].x,y=v[i].y,z=v[i].z;
      v[i].x=3+2*x+.3*y-.2*z;v[i].y=-7+3*y+.4*z;v[i].z=11+.5*z;
    }
    failures+=check("translated affine cell",computePolyhedronVolume3D(v,n),3*L*L*L);
  }
  memset(v,0,sizeof(v));double xyz[4][3]={{0,0,0},{1,0,0},{0,2,0},{0,0,3}};
  for(int i=0;i<4;i++) {v[i].status=STATUS_ACTIVE;v[i].x=xyz[i][0];v[i].y=xyz[i][1];v[i].z=xyz[i][2];int k=0;
    for(int j=0;j<4;j++) if(j!=i) v[i].neighborSites[k++]=j;
  }
  failures+=check("tetrahedron",computePolyhedronVolume3D(v,4),1.);
  return failures?1:0;
}
