#ifndef GFS_MATERIAL_OPEN_STREAM_H
#define GFS_MATERIAL_OPEN_STREAM_H
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Whole material cells leaving an open domain become exterior samples.
 * Position, weight and primitive state are continuous at transfer. Exterior
 * coordinates subsequently follow the recorded velocity; W and thermodynamic
 * state are prescribed at departure. Thus known Wdot=0 is the derivative of
 * the stored boundary data, not a substitute for copied owned-cell Wdot.
 * This is a ballistic outflow approximation, not an exact transparent BC for
 * arbitrary acoustic waves or backflow. No internal mass exchange or repair.
 * Shared 1/2/3D bookkeeping; dimensional evolution is qualified separately.
 */
typedef struct {
 int64_t id;
 double time,x[3],v[3],mass,U,V,W;
} GfsMaterialOpenSample;
typedef struct {
 int dim;size_t count,capacity;GfsMaterialOpenSample*sample;
} GfsMaterialOpenStream;
static inline int gfs_material_open_sample_valid(int dim,const GfsMaterialOpenSample*a){
 if(!a||dim<1||dim>3||!isfinite(a->time)||!(a->mass>0)||!(a->U>0)||!(a->V>0)||
    !isfinite(a->mass)||!isfinite(a->U)||!isfinite(a->V)||!isfinite(a->W))return 0;
 for(int k=0;k<dim;k++)if(!isfinite(a->x[k])||!isfinite(a->v[k]))return 0;
 return 1;
}
static inline void gfs_material_open_free(GfsMaterialOpenStream*b){
 if(!b)return;free(b->sample);memset(b,0,sizeof(*b));
}
static inline int gfs_material_open_push(GfsMaterialOpenStream*b,const GfsMaterialOpenSample*a){
 if(!b||!gfs_material_open_sample_valid(b->dim,a))return 0;
 for(size_t i=0;i<b->count;i++)if(b->sample[i].id==a->id)return 0;
 if(b->count==b->capacity){
  size_t next=b->capacity?2*b->capacity:32;
  if(next<b->capacity||next>SIZE_MAX/sizeof(*b->sample))return 0;
  void*p=realloc(b->sample,next*sizeof(*b->sample));if(!p)return 0;
  b->sample=p;b->capacity=next;
 }
 b->sample[b->count++]=*a;return 1;
}
static inline int gfs_material_open_position(int dim,const GfsMaterialOpenSample*a,
 double time,const double period[3],double out[3]){
 if(!gfs_material_open_sample_valid(dim,a)||!period||!out||!isfinite(time)||time<a->time)return 0;
 for(int k=0;k<3;k++){
  out[k]=k<dim?a->x[k]+(time-a->time)*a->v[k]:0;
  if(!isfinite(out[k])||!isfinite(period[k])||period[k]<0)return 0;
  if(k<dim&&period[k]>0){out[k]=fmod(out[k],period[k]);if(out[k]<0)out[k]+=period[k];}
 }
 return 1;
}
/* Call only after an accepted step. Sampling during RK stages is immutable. */
static inline int gfs_material_open_prune(GfsMaterialOpenStream*b,double time,
 int axis,double limit,const double period[3]){
 if(!b||axis<0||axis>=b->dim||!period||period[axis]!=0||!isfinite(limit))return 0;
 for(size_t i=0;i<b->count;i++){
  double x[3];if(!gfs_material_open_position(b->dim,b->sample+i,time,period,x)||b->sample[i].v[axis]<=0)return 0;
 }
 size_t n=0;for(size_t i=0;i<b->count;i++){
  double x[3];gfs_material_open_position(b->dim,b->sample+i,time,period,x);
  if(x[axis]<=limit)b->sample[n++]=b->sample[i];
 }
 b->count=n;return 1;
}
/* Explicit records avoid persisting struct padding. The time must match the
 * owned checkpoint on read; absent exterior history is not silently rebuilt.
 * This native-binary format is versioned, not a cross-endian interchange.
 */
static inline int gfs_material_open_write(const char*path,double time,const GfsMaterialOpenStream*b){
 if(!path||!b||b->dim<1||b->dim>3||!isfinite(time))return 0;
 const char magic[8]={'G','F','S','O','P','E','N','1'};uint32_t dim=b->dim;uint64_t n=b->count;
 FILE*f=fopen(path,"wb");if(!f)return 0;int ok=1;
 if(fwrite(magic,1,8,f)!=8||fwrite(&dim,4,1,f)!=1||fwrite(&n,8,1,f)!=1||fwrite(&time,8,1,f)!=1)ok=0;
 for(size_t i=0;i<b->count&&ok;i++){
  const GfsMaterialOpenSample*a=b->sample+i;
  double z[11]={a->time,a->x[0],a->x[1],a->x[2],a->v[0],a->v[1],a->v[2],a->mass,a->U,a->V,a->W};
  if(!gfs_material_open_sample_valid(b->dim,a)||a->time>time||fwrite(&a->id,8,1,f)!=1||fwrite(z,8,11,f)!=11)ok=0;
 }
 if(ferror(f))ok=0;if(fclose(f))ok=0;return ok;
}
static inline int gfs_material_open_read(const char*path,double time,GfsMaterialOpenStream*out){
 if(!path||!out||!isfinite(time))return 0;FILE*f=fopen(path,"rb");if(!f)return 0;
 char magic[8];uint32_t dim=0;uint64_t n=0;double saved=0;GfsMaterialOpenStream b={0};int ok=0;
 if(fread(magic,1,8,f)!=8||memcmp(magic,"GFSOPEN1",8)||fread(&dim,4,1,f)!=1||fread(&n,8,1,f)!=1||fread(&saved,8,1,f)!=1||saved!=time||dim<1||dim>3||n>SIZE_MAX/sizeof(*b.sample))goto done;
 b.dim=dim;
 for(uint64_t i=0;i<n;i++){
  int64_t id;double z[11];if(fread(&id,8,1,f)!=1||fread(z,8,11,f)!=11)goto done;
  GfsMaterialOpenSample a={id,z[0],{z[1],z[2],z[3]},{z[4],z[5],z[6]},z[7],z[8],z[9],z[10]};
  if(a.time>time||!gfs_material_open_push(&b,&a))goto done;
 }
 if(fgetc(f)!=EOF||ferror(f))goto done;
 *out=b;memset(&b,0,sizeof(b));ok=1;
 done:gfs_material_open_free(&b);fclose(f);return ok;
}
#endif
