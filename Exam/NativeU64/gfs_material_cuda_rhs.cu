#include <cuda_runtime.h>
#include <stdio.h>
#include "../gfs_material_cuda_rhs.h"
#ifdef GFS_MATERIAL_DIRECT_FV
#include "../gfs_material_direct_flux.h"
#endif
#ifdef GFS_MATERIAL_STAR_ENTROPY
#ifndef GFS_MATERIAL_DIRECT_FV
#error Material star entropy candidate requires direct force/volume/work
#endif
#include "../gfs_material_star_entropy.h"
#endif
/* Reuses the qualified split face-response/scatter helpers. Reflection and
 * gamma match native CPU, rather than historical fixture-only constants. */
__global__ static void material_faces(const GfsFvPrimitive*c,const GfsTraceFace2D*f,const GfsTraceEntry*entry,const double*alpha,double gamma,int nf,double*out,int*status,GfsMaterialFaceResponse*cache,const double*owner_shift,double*face_flux,const double*wall_normal,double*wall_power){
 int k=blockIdx.x*blockDim.x+threadIdx.x;if(k>=nf)return;
 if(wall_power)wall_power[k]=0;
 const GfsTraceFace2D face=f[k];GfsFvPrimitive q[2]={};double rev[2]={0};int start[2]={face.sl,face.sr},len[2]={face.nl,face.nr};
 for(int side=0;side<2;side++)for(int j=start[side];j<start[side]+len[side];j++){
  const GfsTraceEntry e=entry[j];GfsFvPrimitive a=c[e.cell];int owner=side?face.j:face.i;
  GfsMaterialTraceWeights w=gfs_material_entry_weights(&e,owner,alpha?alpha+3*owner:NULL);
  q[side].rho+=w.rho*a.rho;q[side].p+=w.pressure*a.p;
  for(int axis=0;axis<2;axis++)q[side].u[axis]+=w.velocity*gfs_trace_velocity_sign(&e,axis)*a.u[axis];
  rev[side]+=w.reversible*(face.nx*a.u[0]+face.ny*gfs_trace_velocity_sign(&e,1)*a.u[1]);
 }
#ifdef GFS_MATERIAL_GRADIENT_PROJECT
 for(int side=0;side<2;side++){q[side].rho=face.thermo[2*side];q[side].p=face.thermo[2*side+1];}
#endif
 if(owner_shift){q[0].p-=owner_shift[face.i];q[1].p-=owner_shift[face.j];}
 double normal[3]={face.nx,face.ny,0};GfsMaterialFaceResponse r;GfsMaterialFaceExchange x;
#ifdef GFS_MATERIAL_DIRECT_WALLS
 if(face.i==face.j){
  GfsMaterialDirectFlux wall;
  if(owner_shift||cache||!gfs_material_direct_moving_wall_flux(2,q[0],normal,face.area,gamma,wall_normal?wall_normal[k]:0,&wall)){status[k]=0;return;}
  if(face_flux){face_flux[4*k]=wall.force[0];face_flux[4*k+1]=wall.force[1];face_flux[4*k+2]=wall.volume;face_flux[4*k+3]=wall.energy;}
  else {atomicAdd(out+4*face.i,wall.force[0]);atomicAdd(out+4*face.i+1,wall.force[1]);
   if(wall_normal){atomicAdd(out+4*face.i+2,wall.volume);atomicAdd(out+4*face.i+3,wall.energy);}}
  if(wall_power)wall_power[k]=wall.energy;status[k]=1;return;
 }
#endif
 int st=gfs_material_face_response(2,q[0],q[1],normal,gamma,&r);status[k]=st;if(st!=1)return;
 #ifdef GFS_MATERIAL_STAR_ENTROPY
 double entropy_theta;
 if(!gfs_material_star_entropy_limit(2,c[face.i],c[face.j],normal,gamma,
     &r.pstar,&r.ustar,&entropy_theta)){status[k]=0;return;}
 #endif
 if(cache)cache[k]=r;
#ifdef GFS_MATERIAL_DIRECT_FV
 GfsMaterialDirectFlux flux;if(!gfs_material_direct_flux(2,normal,face.area,&r,&flux)){status[k]=0;return;}
 if(face_flux){face_flux[4*k]=flux.force[0];face_flux[4*k+1]=flux.force[1];face_flux[4*k+2]=flux.volume;face_flux[4*k+3]=flux.energy;}
 else {for(int side=0;side<2;side++){int i=side?face.j:face.i;double sign=side?-1:1;atomicAdd(out+4*i,sign*flux.force[0]);atomicAdd(out+4*i+1,sign*flux.force[1]);atomicAdd(out+4*i+2,sign*flux.volume);atomicAdd(out+4*i+3,sign*flux.energy);}}
 return;
#endif
#ifdef GFS_MATERIAL_SPLIT_REVERSIBLE
 r.mean_speed=r.alpha_left*rev[0]+r.alpha_right*rev[1];
 if(cache)cache[k]=r;
#endif
 if(!gfs_material_face_exchange(2,normal,face.area,c[face.i].p,c[face.j].p,&r,&x)){status[k]=0;return;}
 for(int side=0;side<2;side++)for(int j=start[side];j<start[side]+len[side];j++){
  const GfsTraceEntry e=entry[j];int owner=side?face.j:face.i;
  GfsMaterialTraceWeights w=gfs_material_entry_weights(&e,owner,alpha?alpha+3*owner:NULL);
  for(int axis=0;axis<2;axis++)atomicAdd(out+4*e.cell+axis,gfs_trace_velocity_sign(&e,axis)*gfs_material_scatter_force(w,&x,side,axis));
  atomicAdd(out+4*e.cell+2,w.pressure*x.side_volume[side]);
 }
 atomicAdd(out+4*face.i+2,x.owner_volume[0]);atomicAdd(out+4*face.j+2,x.owner_volume[1]);
 atomicAdd(out+4*face.i+3,x.owner_heat[0]);atomicAdd(out+4*face.j+3,x.owner_heat[1]);
}
#ifdef GFS_MATERIAL_DIRECT_FV
/* Each cell traverses incident faces in canonical face-index order, the
 * same order as the host direct loop. No floating-point atomic reduction. */
static __global__ void material_direct_gather(int n,const int*off,const int*occ,const double*flux,double*out){
 int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;double sum[4]={0};
 for(int k=off[i];k<off[i+1];k++){int entry=occ[k],face=entry/2;double sign=(entry&1)?-1:1;
  for(int a=0;a<4;a++)sum[a]+=sign*flux[4*face+a];
 }
 for(int a=0;a<4;a++)out[4*i+a]=sum[a];
}
static __global__ void material_direct_energy(const GfsFvPrimitive*c,int n,double*out){int i=blockIdx.x*blockDim.x+threadIdx.x;if(i<n)out[4*i+3]+=c[i].p*out[4*i+2]-c[i].u[0]*out[4*i]-c[i].u[1]*out[4*i+1];}
#endif
static int material_cuda_boundary_rhs2d(const GfsTraceCsr2D*t,const GfsFvPrimitive*c,double gamma,const double*alpha,double*out,GfsMaterialFaceResponse*cache,const double*owner_shift,const double*wall_normal,double*boundary_power){
 if(!t||t->cells<1||t->faces<1||t->entries<1||!c||!out)return 0;
 if(wall_normal&&(!boundary_power||cache||owner_shift))return 0;
 if(wall_normal)for(int k=0;k<t->faces;k++)if(t->face[k].i==t->face[k].j&&!isfinite(wall_normal[k]))return 0;
#if !defined(GFS_MATERIAL_DIRECT_FV) || !defined(GFS_MATERIAL_DIRECT_WALLS)
 if(wall_normal)return 0;
#endif
 if(boundary_power)*boundary_power=0;
 #ifdef GFS_MATERIAL_DIRECT_FV
#ifndef GFS_MATERIAL_DIRECT_WALLS
 for(int j=0;j<t->entries;j++)if(t->entry[j].pad)return 0;
#endif
#else
 if(owner_shift)return 0;
#endif
 if(owner_shift)for(int j=0;j<t->cells;j++)if(!isfinite(owner_shift[j]))return 0;
 if(alpha)for(int i=0;i<3*t->cells;i++)if(!isfinite(alpha[i])||alpha[i]<0||alpha[i]>1)return 0;
 int result=0,*status=(int*)malloc((size_t)t->faces*sizeof(int)),*ds=NULL;GfsFvPrimitive*dc=NULL;GfsTraceFace2D*df=NULL;GfsTraceEntry*de=NULL;double*da=NULL,*dout=NULL,*dshift=NULL;
 GfsMaterialFaceResponse*dr=NULL;
 double*dwall=NULL,*dpower=NULL,*power=NULL;
 double*face_flux=NULL;int *off=NULL,*occ=NULL,*cursor=NULL,*doff=NULL,*docc=NULL;
 if(!status)return 0;
#define MCALL(call) do{cudaError_t err=(call);if(err!=cudaSuccess){fprintf(stderr,"[MATERIAL_CUDA_FAIL] %s\n",cudaGetErrorString(err));goto done;}}while(0)
 MCALL(cudaMalloc((void**)&dc,(size_t)t->cells*sizeof(*dc)));MCALL(cudaMalloc((void**)&df,(size_t)t->faces*sizeof(*df)));MCALL(cudaMalloc((void**)&de,(size_t)t->entries*sizeof(*de)));MCALL(cudaMalloc((void**)&ds,(size_t)t->faces*sizeof(int)));MCALL(cudaMalloc((void**)&dout,(size_t)t->cells*4*sizeof(double)));
 MCALL(cudaMemcpy(dc,c,(size_t)t->cells*sizeof(*dc),cudaMemcpyHostToDevice));MCALL(cudaMemcpy(df,t->face,(size_t)t->faces*sizeof(*df),cudaMemcpyHostToDevice));MCALL(cudaMemcpy(de,t->entry,(size_t)t->entries*sizeof(*de),cudaMemcpyHostToDevice));
 if(alpha){MCALL(cudaMalloc((void**)&da,(size_t)t->cells*3*sizeof(double)));MCALL(cudaMemcpy(da,alpha,(size_t)t->cells*3*sizeof(double),cudaMemcpyHostToDevice));}
 if(owner_shift){MCALL(cudaMalloc((void**)&dshift,(size_t)t->cells*sizeof(double)));MCALL(cudaMemcpy(dshift,owner_shift,(size_t)t->cells*sizeof(double),cudaMemcpyHostToDevice));}
 if(wall_normal){
  power=(double*)malloc((size_t)t->faces*sizeof(double));if(!power)goto done;
  MCALL(cudaMalloc((void**)&dwall,(size_t)t->faces*sizeof(double)));
  MCALL(cudaMalloc((void**)&dpower,(size_t)t->faces*sizeof(double)));
  MCALL(cudaMemcpy(dwall,wall_normal,(size_t)t->faces*sizeof(double),cudaMemcpyHostToDevice));
 }
 MCALL(cudaMemset(dout,0,(size_t)t->cells*4*sizeof(double)));
 if(cache){MCALL(cudaMalloc((void**)&dr,(size_t)t->faces*sizeof(*dr)));}
#ifdef GFS_MATERIAL_DIRECT_FV
 if(t->faces>INT_MAX/2)goto done;
 off=(int*)calloc((size_t)t->cells+1,sizeof(int));occ=(int*)malloc((size_t)t->faces*2*sizeof(int));cursor=(int*)malloc((size_t)t->cells*sizeof(int));
 if(!off||!occ||!cursor)goto done;
 for(int k=0;k<t->faces;k++){int i=t->face[k].i,j=t->face[k].j;if(i<0||i>=t->cells||j<0||j>=t->cells)goto done;off[i+1]++;if(j!=i)off[j+1]++;}
 for(int i=0;i<t->cells;i++){off[i+1]+=off[i];cursor[i]=off[i];}
 for(int k=0;k<t->faces;k++){int i=t->face[k].i,j=t->face[k].j;occ[cursor[i]++]=2*k;if(j!=i)occ[cursor[j]++]=2*k+1;}
 MCALL(cudaMalloc((void**)&face_flux,(size_t)t->faces*4*sizeof(double)));
 MCALL(cudaMalloc((void**)&doff,((size_t)t->cells+1)*sizeof(int)));MCALL(cudaMalloc((void**)&docc,(size_t)off[t->cells]*sizeof(int)));
 MCALL(cudaMemcpy(doff,off,((size_t)t->cells+1)*sizeof(int),cudaMemcpyHostToDevice));MCALL(cudaMemcpy(docc,occ,(size_t)off[t->cells]*sizeof(int),cudaMemcpyHostToDevice));
#endif
 material_faces<<<(t->faces+127)/128,128>>>(dc,df,de,da,gamma,t->faces,dout,ds,dr,dshift,face_flux,dwall,dpower);MCALL(cudaGetLastError());
 MCALL(cudaMemcpy(status,ds,(size_t)t->faces*sizeof(int),cudaMemcpyDeviceToHost));
 for(int i=0;i<t->faces;i++)if(status[i]!=1){
  const GfsTraceFace2D f=t->face[i];GfsFvPrimitive q[2]={};int starts[2]={f.sl,f.sr},lengths[2]={f.nl,f.nr};
  for(int side=0;side<2;side++)for(int j=starts[side];j<starts[side]+lengths[side];j++){
   const GfsTraceEntry e=t->entry[j];int owner=side?f.j:f.i;
   GfsMaterialTraceWeights w=gfs_material_entry_weights(&e,owner,alpha?alpha+3*owner:NULL);
   q[side].rho+=w.rho*c[e.cell].rho;q[side].p+=w.pressure*c[e.cell].p;
   for(int axis=0;axis<2;axis++)q[side].u[axis]+=w.velocity*gfs_trace_velocity_sign(&e,axis)*c[e.cell].u[axis];
  }
#ifdef GFS_MATERIAL_GRADIENT_PROJECT
  for(int side=0;side<2;side++){q[side].rho=f.thermo[2*side];q[side].p=f.thermo[2*side+1];}
#endif
  if(owner_shift){q[0].p-=owner_shift[f.i];q[1].p-=owner_shift[f.j];}
  double normal[3]={f.nx,f.ny,0};GfsMaterialFaceResponse r;
  int host=gfs_material_face_response(2,q[0],q[1],normal,gamma,&r);
  #ifdef GFS_MATERIAL_STAR_ENTROPY
  if(host==1){double theta;if(!gfs_material_star_entropy_limit(2,c[f.i],c[f.j],normal,gamma,&r.pstar,&r.ustar,&theta))host=0;}
  #endif
  fprintf(stderr,"[MATERIAL_FACE_FAILURE] face=%d owners=%d,%d gpu=%d host=%d gamma=%.17g normal=%.17g,%.17g left=%.17g,%.17g,%.17g,%.17g right=%.17g,%.17g,%.17g,%.17g\n",i,f.i,f.j,status[i],host,gamma,f.nx,f.ny,q[0].rho,q[0].p,q[0].u[0],q[0].u[1],q[1].rho,q[1].p,q[1].u[0],q[1].u[1]);
  const char *path=getenv("GFS_MATERIAL_FAILED_FACE_PATH");
  if(path&&*path){FILE *fp=fopen(path,"wx");if(fp){
   fprintf(fp,"side,cell,owner,weight,pad,rho,p,vx,vy,alpha_rho,alpha_p,alpha_v\n");
   for(int side=0;side<2;side++)for(int j=starts[side];j<starts[side]+lengths[side];j++){
    GfsTraceEntry e=t->entry[j];int owner=side?f.j:f.i;
    fprintf(fp,"%d,%d,%d,%.17g,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",side,e.cell,owner,e.weight,e.pad,c[e.cell].rho,c[e.cell].p,c[e.cell].u[0],c[e.cell].u[1],alpha?alpha[3*owner]:1.,alpha?alpha[3*owner+1]:1.,alpha?alpha[3*owner+2]:1.);
   }fclose(fp);
  }}
  result=status[i]==2?2:0;goto done;
 }
#ifdef GFS_MATERIAL_DIRECT_FV
 material_direct_gather<<<(t->cells+127)/128,128>>>(t->cells,doff,docc,face_flux,dout);MCALL(cudaGetLastError());
 material_direct_energy<<<(t->cells+127)/128,128>>>(dc,t->cells,dout);MCALL(cudaGetLastError());
#endif
 MCALL(cudaMemcpy(out,dout,(size_t)t->cells*4*sizeof(double),cudaMemcpyDeviceToHost));
 if(cache){MCALL(cudaMemcpy(cache,dr,(size_t)t->faces*sizeof(*dr),cudaMemcpyDeviceToHost));}
 if(wall_normal){
  MCALL(cudaMemcpy(power,dpower,(size_t)t->faces*sizeof(double),cudaMemcpyDeviceToHost));
  for(int k=0;k<t->faces;k++){if(!isfinite(power[k]))goto done;*boundary_power+=power[k];}
  if(!isfinite(*boundary_power))goto done;
 }
 result=1;
 done:cudaFree(dwall);cudaFree(dpower);free(power);cudaFree(face_flux);cudaFree(doff);cudaFree(docc);free(off);free(occ);free(cursor);cudaFree(dshift);cudaFree(dr);cudaFree(dc);cudaFree(df);cudaFree(de);cudaFree(ds);cudaFree(da);cudaFree(dout);free(status);return result;
#undef MCALL
}

extern "C" int gfs_material_cuda_owner_rhs_cache2d(const GfsTraceCsr2D*t,const GfsFvPrimitive*c,double gamma,const double*alpha,double*out,GfsMaterialFaceResponse*cache,const double*owner_shift){
 return material_cuda_boundary_rhs2d(t,c,gamma,alpha,out,cache,owner_shift,NULL,NULL);
}
extern "C" int gfs_material_cuda_rhs_moving_walls2d(const GfsTraceCsr2D*t,const GfsFvPrimitive*c,double gamma,const double*alpha,const double*wall_normal,double*out,double*boundary_power){
 if(!wall_normal||!boundary_power)return 0;
 return material_cuda_boundary_rhs2d(t,c,gamma,alpha,out,NULL,NULL,wall_normal,boundary_power);
}

extern "C" int gfs_material_cuda_rhs_cache2d(const GfsTraceCsr2D*t,const GfsFvPrimitive*c,double gamma,const double*alpha,double*out,GfsMaterialFaceResponse*cache){return gfs_material_cuda_owner_rhs_cache2d(t,c,gamma,alpha,out,cache,NULL);}

extern "C" int gfs_material_cuda_rhs2d(const GfsTraceCsr2D*t,const GfsFvPrimitive*c,double gamma,const double*alpha,double*out){return gfs_material_cuda_rhs_cache2d(t,c,gamma,alpha,out,NULL);}

/* Independent RHS vectors share one immutable geometry/impedance upload. */
static __global__ void material_actions(const GfsTraceFace2D*faces,const GfsTraceEntry*entry,const GfsMaterialFaceResponse*cache,const double*alpha,int nf,int n,const double*g,double*dg,double*ct,int*valid){
 int k=blockIdx.x*blockDim.x+threadIdx.x,row=blockIdx.y;if(k>=nf)return;
 const GfsTraceFace2D f=faces[k];const GfsMaterialFaceResponse r=cache[k];
 g+=(size_t)row*n;dg+=(size_t)row*n;ct+=(size_t)row*2*n;
 int starts[2]={f.sl,f.sr},counts[2]={f.nl,f.nr};double jump=0;
 for(int s=0;s<2;s++)for(int j=starts[s];j<starts[s]+counts[s];j++){
  const GfsTraceEntry e=entry[j];int owner=s?f.j:f.i;
  GfsMaterialTraceWeights w=gfs_material_entry_weights(&e,owner,alpha?alpha+3*owner:NULL);
  jump+=(s?-1:1)*w.pressure*g[e.cell];
 }
 double z=r.zleft+r.zright,dp=g[f.i]-g[f.j];
 if(!(z>0)||!isfinite(jump)||!isfinite(dp)){atomicExch(valid,0);return;}
 double rate=f.area*jump/z;
 for(int s=0;s<2;s++)for(int j=starts[s];j<starts[s]+counts[s];j++){
  const GfsTraceEntry e=entry[j];int owner=s?f.j:f.i;
  GfsMaterialTraceWeights w=gfs_material_entry_weights(&e,owner,alpha?alpha+3*owner:NULL);
  double force=w.reversible*f.area*(s?r.alpha_right:r.alpha_left)*dp;
  atomicAdd(dg+e.cell,(s?-1:1)*w.pressure*rate);
  atomicAdd(ct+2*e.cell,force*f.nx);
  atomicAdd(ct+2*e.cell+1,force*f.ny*gfs_trace_velocity_sign(&e,1));
 }
}
#ifdef GFS_MATERIAL_CUDA_ACTIONS_GATHER
typedef struct {int face,entry,side;} MaterialOccurrence;
static __global__ void material_action_rates(const GfsTraceFace2D*faces,const GfsTraceEntry*entries,const GfsMaterialFaceResponse*cache,const double*alpha,int nf,int n,const double*g,double*rate,int*valid){
 int k=blockIdx.x*blockDim.x+threadIdx.x,row=blockIdx.y;if(k>=nf)return;
 const GfsTraceFace2D f=faces[k];g+=(size_t)row*n;double jump=0;
 int starts[2]={f.sl,f.sr},counts[2]={f.nl,f.nr};
 for(int side=0;side<2;side++)for(int j=starts[side];j<starts[side]+counts[side];j++){
  const GfsTraceEntry e=entries[j];int owner=side?f.j:f.i;
  GfsMaterialTraceWeights w=gfs_material_entry_weights(&e,owner,alpha?alpha+3*owner:NULL);
  jump+=(side?-1:1)*w.pressure*g[e.cell];
 }
 double z=cache[k].zleft+cache[k].zright;
 if(!(z>0)||!isfinite(jump)||!isfinite(g[f.i]-g[f.j])){atomicExch(valid,0);return;}
 rate[(size_t)row*nf+k]=f.area*jump/z;
}
static __global__ void material_action_gather(const int*off,const MaterialOccurrence*occ,const GfsTraceFace2D*faces,const GfsTraceEntry*entries,const GfsMaterialFaceResponse*cache,const double*alpha,int nf,int n,const double*g,const double*rate,double*dg,double*ct){
 int i=blockIdx.x*blockDim.x+threadIdx.x,row=blockIdx.y;if(i>=n)return;
 g+=(size_t)row*n;rate+=(size_t)row*nf;double d=0,cx=0,cy=0;
 for(int j=off[i];j<off[i+1];j++){
  MaterialOccurrence o=occ[j];GfsTraceFace2D f=faces[o.face];GfsTraceEntry e=entries[o.entry];GfsMaterialFaceResponse r=cache[o.face];
  int owner=o.side?f.j:f.i;
  GfsMaterialTraceWeights w=gfs_material_entry_weights(&e,owner,alpha?alpha+3*owner:NULL);
  double force=w.reversible*f.area*(o.side?r.alpha_right:r.alpha_left)*(g[f.i]-g[f.j]);
  d+=(o.side?-1:1)*w.pressure*rate[o.face];cx+=force*f.nx;cy+=force*f.ny*gfs_trace_velocity_sign(&e,1);
 }
 dg[(size_t)row*n+i]=d;ct[(size_t)row*2*n+2*i]=cx;ct[(size_t)row*2*n+2*i+1]=cy;
}
#endif
extern "C" int gfs_material_cuda_actions2d(const GfsTraceCsr2D*t,const GfsMaterialFaceResponse*cache,const double*alpha,int rows,const double*g,double*dg,double*ct){
 if(!t||t->cells<1||t->faces<1||t->entries<1||!cache||rows<1||rows>65535||!g||!dg||!ct)return 0;
 if(alpha)for(int i=0;i<3*t->cells;i++)if(!isfinite(alpha[i])||alpha[i]<0||alpha[i]>1)return 0;
 int ok=0,valid=1,*dv=NULL;GfsTraceFace2D*df=NULL;GfsTraceEntry*de=NULL;GfsMaterialFaceResponse*dc=NULL;double*da=NULL,*gg=NULL,*dd=NULL,*cc=NULL;
 size_t bytes=(size_t)rows*t->cells*sizeof(double);
#ifdef GFS_MATERIAL_CUDA_ACTIONS_GATHER
 int *off=NULL,*pos=NULL,*d_off=NULL,total=0;MaterialOccurrence *occ=NULL,*d_occ=NULL;double*d_rate=NULL;
#endif
#define ACALL(call) do{cudaError_t err=(call);if(err!=cudaSuccess){fprintf(stderr,"[MATERIAL_ACTION_CUDA_FAIL] %s\n",cudaGetErrorString(err));goto done;}}while(0)
 ACALL(cudaMalloc((void**)&df,(size_t)t->faces*sizeof(*df)));ACALL(cudaMalloc((void**)&de,(size_t)t->entries*sizeof(*de)));ACALL(cudaMalloc((void**)&dc,(size_t)t->faces*sizeof(*dc)));ACALL(cudaMalloc((void**)&dv,sizeof(int)));
 ACALL(cudaMalloc((void**)&gg,bytes));ACALL(cudaMalloc((void**)&dd,bytes));ACALL(cudaMalloc((void**)&cc,2*bytes));
 ACALL(cudaMemcpy(df,t->face,(size_t)t->faces*sizeof(*df),cudaMemcpyHostToDevice));ACALL(cudaMemcpy(de,t->entry,(size_t)t->entries*sizeof(*de),cudaMemcpyHostToDevice));ACALL(cudaMemcpy(dc,cache,(size_t)t->faces*sizeof(*dc),cudaMemcpyHostToDevice));ACALL(cudaMemcpy(gg,g,bytes,cudaMemcpyHostToDevice));ACALL(cudaMemcpy(dv,&valid,sizeof(int),cudaMemcpyHostToDevice));
 if(alpha){ACALL(cudaMalloc((void**)&da,(size_t)t->cells*3*sizeof(double)));ACALL(cudaMemcpy(da,alpha,(size_t)t->cells*3*sizeof(double),cudaMemcpyHostToDevice));}
 ACALL(cudaMemset(dd,0,bytes));ACALL(cudaMemset(cc,0,2*bytes));
#ifdef GFS_MATERIAL_CUDA_ACTIONS_GATHER
 off=(int*)calloc(t->cells+1,sizeof(int));pos=(int*)malloc(t->cells*sizeof(int));if(!off||!pos)goto done;
 for(int k=0;k<t->faces;k++){
  const GfsTraceFace2D*f=t->face+k;int starts[2]={f->sl,f->sr},counts[2]={f->nl,f->nr};
  for(int side=0;side<2;side++)for(int j=starts[side];j<starts[side]+counts[side];j++){
   if(j<0||j>=t->entries||t->entry[j].cell<0||t->entry[j].cell>=t->cells)goto done;
   off[t->entry[j].cell+1]++;
  }
 }
 for(int i=1;i<=t->cells;i++)off[i]+=off[i-1];total=off[t->cells];
 occ=(MaterialOccurrence*)malloc((size_t)total*sizeof(*occ));if(!occ)goto done;memcpy(pos,off,t->cells*sizeof(int));
 for(int k=0;k<t->faces;k++){
  const GfsTraceFace2D*f=t->face+k;int starts[2]={f->sl,f->sr},counts[2]={f->nl,f->nr};
  for(int side=0;side<2;side++)for(int j=starts[side];j<starts[side]+counts[side];j++)occ[pos[t->entry[j].cell]++]={k,j,side};
 }
 ACALL(cudaMalloc((void**)&d_off,(t->cells+1)*sizeof(int)));ACALL(cudaMalloc((void**)&d_occ,(size_t)total*sizeof(*occ)));ACALL(cudaMalloc((void**)&d_rate,(size_t)rows*t->faces*sizeof(double)));
 ACALL(cudaMemcpy(d_off,off,(t->cells+1)*sizeof(int),cudaMemcpyHostToDevice));ACALL(cudaMemcpy(d_occ,occ,(size_t)total*sizeof(*occ),cudaMemcpyHostToDevice));
 material_action_rates<<<dim3((t->faces+127)/128,rows),128>>>(df,de,dc,da,t->faces,t->cells,gg,d_rate,dv);
 material_action_gather<<<dim3((t->cells+127)/128,rows),128>>>(d_off,d_occ,df,de,dc,da,t->faces,t->cells,gg,d_rate,dd,cc);
#else
 material_actions<<<dim3((t->faces+127)/128,rows),128>>>(df,de,dc,da,t->faces,t->cells,gg,dd,cc,dv); 
#endif
 ACALL(cudaGetLastError());ACALL(cudaMemcpy(&valid,dv,sizeof(int),cudaMemcpyDeviceToHost));if(!valid)goto done;
 ACALL(cudaMemcpy(dg,dd,bytes,cudaMemcpyDeviceToHost));ACALL(cudaMemcpy(ct,cc,2*bytes,cudaMemcpyDeviceToHost));ok=1;
 done:
#ifdef GFS_MATERIAL_CUDA_ACTIONS_GATHER
 free(off);free(pos);free(occ);cudaFree(d_off);cudaFree(d_occ);cudaFree(d_rate);
#endif
 cudaFree(df);cudaFree(de);cudaFree(dc);cudaFree(dv);cudaFree(da);cudaFree(gg);cudaFree(dd);cudaFree(cc);return ok;
#undef ACALL
}
