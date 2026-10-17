#include <stdlib.h>
#include "laguerre.h"
int cell(int n,double w0,double *data,double *output,double *volume){
 Point3D c={0};c.weightSquared=w0;Point3D *input=calloc(n,sizeof(*input)),*neighbors=calloc(n,sizeof(*neighbors));
 Vertex3D *v=calloc(512,sizeof(*v));int *work=calloc(512,sizeof(int));
 for(int i=0;i<n;i++){input[i].x=data[4*i];input[i].y=data[4*i+1];input[i].z=data[4*i+2];input[i].weightSquared=data[4*i+3];input[i].index=i+1;}
 int count=constructLaguerreCell3D(&c,input,neighbors,n,v,512,4.,0,work);
 if(count>0){*volume=computePolyhedronVolume3D(v,count);for(int i=0;i<count;i++){output[3*i]=v[i].x;output[3*i+1]=v[i].y;output[3*i+2]=v[i].z;}}
 free(input);free(neighbors);free(v);free(work);return count;
}
