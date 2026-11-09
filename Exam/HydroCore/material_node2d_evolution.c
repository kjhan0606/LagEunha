/* Coupled time-integration control for the shared native C/CUDA operators.
 * Structured quadrilateral connectivity is a test fixture, NOT a replacement
 * for production Voronoi/Laguerre geometry or its MPI ownership. First order
 * cell states; SSPRK2. No claim of completed KH/RT or PLM integration.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../gfs_material_mesh2d.h"
#include "../gfs_material_node2d.h"
#include "../gfs_material_reconstruction2d.h"
#include "../exam_gpu.h"
#include "../gfs_fv_native.h"
#include "../gfs_material_csr2d.h"
#include "../gfs_material_csr_native.inc"
#define NC 1024
#define NV 1089
#define GAMMA 1.4
static int n,nc,nv,ids[NC][64],corners[NC],free_x[NV],free_y[NV],cnt[NV],owner[NV][8];
static double mass[NC];
static int order=1,neighbors[NC][64],neighbor_count[NC];
static double center[NC][2],primitive[NC][4],gradient[NC][4][2];
typedef struct {double x[NV][2],p[NC][2],E[NC];} State;
static State current,stage,next,k1,k2;
static GfsMaterialHalfEdge edges[NV][8];
static double vtx[NV][2],zero[NV][2];
static double max_residual=0,min_heat=0;
static int primitives(const State *s,int c,double *rho,double *p,double u[2]){
 double V=gfs_mm_area(corners[c],ids[c],s->x);if(!(V>0))return 0;
 u[0]=s->p[c][0]/mass[c];u[1]=s->p[c][1]/mass[c];*rho=mass[c]/V;
 *p=(GAMMA-1)*(s->E[c]-.5*mass[c]*(u[0]*u[0]+u[1]*u[1]))/V;
 return *p>0&&isfinite(*p);
}
static int valid(const State *s){
 double area=0;
 for(int c=0;c<nc;c++){
  double rho,p,u[2];if(!primitives(s,c,&rho,&p,u)||!gfs_mm_convex_step(corners[c],ids[c],s->x,zero,0))return 0;
  area+=gfs_mm_area(corners[c],ids[c],s->x);
 }
 return fabs(area-1)<1e-10;
}
static int reconstruction(const State *s){
 for(int c=0;c<nc;c++){
  double rho,p,u[2];if(!primitives(s,c,&rho,&p,u))return 0;
  primitive[c][0]=rho;primitive[c][1]=p;primitive[c][2]=u[0];primitive[c][3]=u[1];
  int o=ids[c][0];double area2=0,cx=0,cy=0;
  for(int k=1;k<corners[c]-1;k++){
   int i=ids[c][k],j=ids[c][k+1];double ax=s->x[i][0]-s->x[o][0],ay=s->x[i][1]-s->x[o][1],bx=s->x[j][0]-s->x[o][0],by=s->x[j][1]-s->x[o][1];
   double cross=ax*by-ay*bx;area2+=cross;cx+=cross*(ax+bx);cy+=cross*(ay+by);
  }
  if(!(area2>0))return 0;
  center[c][0]=s->x[o][0]+cx/(3*area2);center[c][1]=s->x[o][1]+cy/(3*area2);
 }
 memset(gradient,0,sizeof(gradient));if(order==1)return 1;
 for(int c=0;c<nc;c++){
  double dx[64],dy[64],dq[64],dvx[64],dvy[64],vx[64],vy[64],condition=0;
  int count=neighbor_count[c];
  for(int k=0;k<count;k++){int j=neighbors[c][k];dx[k]=center[j][0]-center[c][0];dy[k]=center[j][1]-center[c][1];dvx[k]=primitive[j][2]-primitive[c][2];dvy[k]=primitive[j][3]-primitive[c][3];}
  for(int k=0;k<corners[c];k++){vx[k]=s->x[ids[c][k]][0]-center[c][0];vy[k]=s->x[ids[c][k]][1]-center[c][1];}
  for(int f=0;f<4;f++){
   double low=primitive[c][f],high=low;
   for(int k=0;k<count;k++){int j=neighbors[c][k];dq[k]=primitive[j][f]-primitive[c][f];low=fmin(low,primitive[j][f]);high=fmax(high,primitive[j][f]);}
   if(!gfs_material_fit2d(count,dx,dy,dq,gradient[c][f],&condition))return 0;
   if(f<2){double alpha=gfs_material_scalar_alpha2d(corners[c],vx,vy,gradient[c][f][0],gradient[c][f][1],primitive[c][f],low,high);gradient[c][f][0]*=alpha;gradient[c][f][1]*=alpha;}
  }
  double alpha=gfs_velocity_midpoint_alpha(count,dx,dy,dvx,dvy,gradient[c][2][0],gradient[c][2][1],gradient[c][3][0],gradient[c][3][1],condition,hypot(primitive[c][2],primitive[c][3]));
  if(alpha<0)return 0;
  for(int f=2;f<4;f++)for(int a=0;a<2;a++)gradient[c][f][a]*=alpha;
 }
 return 1;
}
static int rhs(const State *s,State *rate){
 memset(rate,0,sizeof(*rate));memset(cnt,0,sizeof(cnt));
 if(!reconstruction(s))return 0;
 for(int c=0;c<nc;c++){
  double rho,p,u[2];if(!primitives(s,c,&rho,&p,u))return 0;
  for(int k=0;k<corners[c];k++){
   int a=ids[c][k],b=ids[c][(k+1)%corners[c]];
   double sx=s->x[b][1]-s->x[a][1],sy=s->x[a][0]-s->x[b][0],l=hypot(sx,sy);
   if(!(l>0))return 0;
   for(int t=0;t<2;t++){
    int node=t?b:a,idx=cnt[node]++;if(idx>=8)return 0;
    GfsMaterialHalfEdge *e=&edges[node][idx];owner[node][idx]=c;
    e->n[0]=sx/l;e->n[1]=sy/l;e->length=.5*l;
    double dx=s->x[node][0]-center[c][0],dy=s->x[node][1]-center[c][1],q[4];
    for(int f=0;f<4;f++)q[f]=primitive[c][f]+gradient[c][f][0]*dx+gradient[c][f][1]*dy;
    e->rho=q[0];e->p=q[1];e->u[0]=q[2];e->u[1]=q[3];
   }
  }
 }
 for(int node=0;node<nv;node++){
  int fx=free_x[node],fy=free_y[node],it;
  double v[2]={0,0},res;
  for(int k=0;k<cnt[node];k++)for(int a=0;a<2;a++)v[a]+=edges[node][k].u[a]/cnt[node];
  if(!fx)v[0]=0;
  if(!fy)v[1]=0;
  if(!gfs_material_node_solve(cnt[node],edges[node],GAMMA,fx,fy,v,&it,&res))return 0;
  max_residual=fmax(max_residual,res);
  for(int a=0;a<2;a++)rate->x[node][a]=v[a];
  for(int k=0;k<cnt[node];k++){
   double f[2],de,heat;
   if(!gfs_material_halfedge_rate(&edges[node][k],GAMMA,v,f,&de,&heat))return 0;
   int c=owner[node][k];rate->p[c][0]+=f[0];rate->p[c][1]+=f[1];rate->E[c]+=de;
   min_heat=fmin(min_heat,heat);
  }
 }
 return 1;
}
static void forward(State *out,const State *a,const State *rate,double dt){
 *out=*a;
 for(int i=0;i<nv;i++)for(int k=0;k<2;k++)out->x[i][k]+=dt*rate->x[i][k];
 for(int i=0;i<nc;i++){out->E[i]+=dt*rate->E[i];for(int k=0;k<2;k++)out->p[i][k]+=dt*rate->p[i][k];}
}
static int path_valid(const State *a,const State *rate,double dt){
 for(int i=0;i<nv;i++)for(int k=0;k<2;k++)vtx[i][k]=rate->x[i][k];
 for(int c=0;c<nc;c++)if(!gfs_mm_convex_step(corners[c],ids[c],a->x,vtx,dt))return 0;
 return 1;
}
static int initialize(const char *kind,int native){
 double rho,p;
 if(native){
  GfsFvMesh m={0};m.dim=2;m.n=nc;m.L[0]=m.L[1]=1;m.gamma=GAMMA;m.bc[0]=m.bc[1]=GFS_FV_REFLECT;m.cell=calloc(nc,sizeof(*m.cell));
  for(int j=0;j<n;j++)for(int i=0;i<n;i++){
   int c=j*n+i;m.cell[c].q[0]=(i+.5+.07*sin(3*c))/n;m.cell[c].q[1]=(j+.5+.07*cos(5*c))/n;
   m.cell[c].W=native==2?.001*sin(c)/(n*n):0;
  }
  if(gfs_fv_geometry(&m))return 0;
  FaceCSR f={0};ParticleSoA part={0};f.n_particles=nc;part.n_total=nc;
  f.face_offset=calloc(nc+1,sizeof(int));part.x=malloc(nc*sizeof(double));part.y=malloc(nc*sizeof(double));
  for(int c=0;c<nc;c++){f.face_offset[c+1]=f.face_offset[c]+m.geometry[c].faces;part.x[c]=m.cell[c].q[0];part.y[c]=m.cell[c].q[1];}
  int nf=f.n_faces_total=f.face_offset[nc];
  f.c1x=malloc(nf*sizeof(double));f.c1y=malloc(nf*sizeof(double));f.c2x=malloc(nf*sizeof(double));f.c2y=malloc(nf*sizeof(double));f.neighbor_idx=malloc(nf*sizeof(int));f.is_ghost=malloc(nf*sizeof(int));
  for(int c=0;c<nc;c++)for(int j=0;j<m.geometry[c].faces;j++){
   GfsFvNativeFace *face=&m.geometry[c].face[j];int k=f.face_offset[c]+j;if(face->count!=3)return 0;
   GfsFvPoint *pt=&m.geometry[c].point[face->point];
   double dx=(pt[2].x[0]-pt[0].x[0])/sqrt(.6),dy=(pt[2].x[1]-pt[0].x[1])/sqrt(.6);
   f.c1x[k]=pt[1].x[0]-.5*dx;f.c2x[k]=pt[1].x[0]+.5*dx;f.c1y[k]=pt[1].x[1]-.5*dy;f.c2y[k]=pt[1].x[1]+.5*dy;
   f.neighbor_idx[k]=face->j;f.is_ghost[k]=face->mirror!=0;
  }
  GfsMaterialCsr2D mesh={0};if(!gfs_material_import_facecsr(&f,&part,&mesh))return 0;
  nv=mesh.nodes;if(nv>NV)return 0;
  for(int i=0;i<nv;i++)for(int a=0;a<2;a++)current.x[i][a]=mesh.x[i][a];
  for(int c=0;c<nc;c++){
   corners[c]=mesh.offset[c+1]-mesh.offset[c];if(corners[c]>64)return 0;
   for(int k=0;k<corners[c];k++)ids[c][k]=mesh.vertex[mesh.offset[c]+k];
   double x=m.cell[c].q[0]+m.geometry[c].centroid[0]-.5,y=m.cell[c].q[1]+m.geometry[c].centroid[1]-.5;
   rho=strcmp(kind,"contact")==0?(x<0?1:.125):1;p=strcmp(kind,"pulse")==0?1+.2*exp(-50*(x*x+y*y)):1;
   mass[c]=rho*m.geometry[c].V;current.E[c]=p*m.geometry[c].V/(GAMMA-1);
  }
  gfs_material_csr_free(&mesh);free(f.face_offset);free(f.c1x);free(f.c1y);free(f.c2x);free(f.c2y);free(f.neighbor_idx);free(f.is_ghost);free(part.x);free(part.y);gfs_fv_geometry_free(&m);free(m.cell);
 }else{
  nv=(n+1)*(n+1);
  for(int j=0;j<=n;j++)for(int i=0;i<=n;i++){int v=j*(n+1)+i;current.x[v][0]=(double)i/n;current.x[v][1]=(double)j/n;}
  for(int j=0;j<n;j++)for(int i=0;i<n;i++){
   int c=j*n+i,a=j*(n+1)+i;corners[c]=4;ids[c][0]=a;ids[c][1]=a+1;ids[c][2]=a+n+2;ids[c][3]=a+n+1;
   double x=(i+.5)/n-.5,y=(j+.5)/n-.5,volume=1./(n*n);
   rho=strcmp(kind,"contact")==0?(x<0?1:.125):1;p=strcmp(kind,"pulse")==0?1+.2*exp(-50*(x*x+y*y)):1;
   mass[c]=rho*volume;current.E[c]=p*volume/(GAMMA-1);
  }
 }
 for(int i=0;i<nv;i++){
  free_x[i]=current.x[i][0]>1e-12&&current.x[i][0]<1-1e-12;
  free_y[i]=current.x[i][1]>1e-12&&current.x[i][1]<1-1e-12;
 }
 int incident_count[NV]={0},incident[NV][8];
 for(int c=0;c<nc;c++)for(int k=0;k<corners[c];k++){
  int v=ids[c][k];if(incident_count[v]>=8)return 0;incident[v][incident_count[v]++]=c;
 }
 for(int c=0;c<nc;c++)for(int k=0;k<corners[c];k++){
  int v=ids[c][k];
  for(int j=0;j<incident_count[v];j++){
   int d=incident[v][j];if(d==c)continue;
   int pos=0;while(pos<neighbor_count[c]&&neighbors[c][pos]<d)pos++;
   if(pos<neighbor_count[c]&&neighbors[c][pos]==d)continue;
   if(neighbor_count[c]>=64)return 0;
   for(int a=neighbor_count[c];a>pos;a--)neighbors[c][a]=neighbors[c][a-1];
   neighbors[c][pos]=d;neighbor_count[c]++;
  }
 }
 return 1;
}
int main(int argc,char **argv){
 n=argc>1?atoi(argv[1]):8;const char *kind=argc>2?argv[2]:"pulse";int native=argc>3?atoi(argv[3]):0;
 order=argc>6?atoi(argv[6]):1;
 if((order!=1&&order!=2)||n<2||n>16)return 1;
 nc=n*n;double target=.1;
 if(!initialize(kind,native))return 9;
 double E0=0;for(int c=0;c<nc;c++)E0+=current.E[c];
 if(!valid(&current))return 2;
 double factor=argc>5?atof(argv[5]):1;
 if(!(factor>0&&factor<=1))return 10;
 double t=0,dt=factor*.05/n/sqrt(1.4*1.2/.125);int steps=0,rejected=0;
 while(t<target){
  double h=fmin(dt,target-t);if(!rhs(&current,&k1))return 3;int accepted=0,failure=0;
  for(int tries=0;tries<40;tries++){
   failure=1;
   if(path_valid(&current,&k1,h)){
    forward(&stage,&current,&k1,h);
    failure=2;
    if(valid(&stage)){
     failure=3;
     if(rhs(&stage,&k2)){
      failure=4;
      if(path_valid(&stage,&k2,h)){
     forward(&next,&stage,&k2,h);
     for(int i=0;i<nv;i++)for(int k=0;k<2;k++)next.x[i][k]=.5*(current.x[i][k]+next.x[i][k]);
     for(int i=0;i<nc;i++){next.E[i]=.5*(current.E[i]+next.E[i]);for(int k=0;k<2;k++)next.p[i][k]=.5*(current.p[i][k]+next.p[i][k]);}
     failure=5;
     if(valid(&next)){accepted=1;break;}
      }
     }
    }
   }
   rejected++;h*=.5;
  }
  if(!accepted||h<1e-14){
   fprintf(stderr,"STOP time=%.17g step=%d trial_dt=%.17g failure_stage=%d rejected=%d\n",t,steps,h,failure,rejected);
   if(argc>4){FILE *f=fopen(argv[4],"w");if(f){fprintf(f,"# time %.17g cells %d nodes %d\n",t,nc,nv);for(int i=0;i<nv;i++)fprintf(f,"v %d %.17g %.17g %.17g %.17g\n",i,current.x[i][0],current.x[i][1],k1.x[i][0],k1.x[i][1]);for(int c=0;c<nc;c++){fprintf(f,"c %d %.17g %.17g %.17g %.17g %d",c,mass[c],current.p[c][0],current.p[c][1],current.E[c],corners[c]);for(int k=0;k<corners[c];k++)fprintf(f," %d",ids[c][k]);fputc('\n',f);}fclose(f);}}
   return 4;
  }
  current=next;t+=h;steps++;if(steps>100000)return 5;
 }
 double E=0,minP=1e100,maxP=0,maxSpeed=0,area=0;
 for(int c=0;c<nc;c++){double rho,p,u[2];if(!primitives(&current,c,&rho,&p,u))return 6;E+=current.E[c];area+=gfs_mm_area(corners[c],ids[c],current.x);minP=fmin(minP,p);maxP=fmax(maxP,p);maxSpeed=fmax(maxSpeed,hypot(u[0],u[1]));}
 printf("{\"case\":\"%s\",\"N\":%d,\"time\":%.17g,\"steps\":%d,\"rejected\":%d,\"energy_relative_error\":%.17g,\"area_error\":%.17g,\"min_pressure\":%.17g,\"max_pressure\":%.17g,\"max_speed\":%.17g,\"node_residual\":%.17g,\"min_dissipative_heat\":%.17g}\n",kind,n,t,steps,rejected,(E-E0)/E0,area-1,minP,maxP,maxSpeed,max_residual,min_heat);
 if(fabs((E-E0)/E0)>1e-10||min_heat < -1e-12)return 7;
 if(strcmp(kind,"contact")==0&&(maxSpeed>1e-12||fabs(maxP-1)>1e-12||fabs(minP-1)>1e-12))return 8;
 return 0;
}
