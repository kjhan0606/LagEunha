#ifndef GFS_MATERIAL_VOLUME_INVERSE_H
#define GFS_MATERIAL_VOLUME_INVERSE_H
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <mkl_pardiso.h>

/* Dimension-independent dV/dW assembly: each face contributes
 * c(dWi-dWowner), c=area/(2*distance). A prescribed exterior uses owner=-1.
 * The geometry caller supplies true face measure/distance and copied owners.
 * No gauge removal is appropriate for this anchored open-domain inverse.
 * This CPU sparse-direct implementation is a qualification path, not a
 * statement about final production cost or multi-rank support.
 */

#include <stdio.h>
#include <time.h>
static double profiler_clock(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
typedef struct {void*pt[64];MKL_INT ip[64];int ready,n,nz;MKL_INT*ia,*ja;double*A;} ProfileCache;
static ProfileCache profile_cache;
static int profile_reuse;
static int profile_clear(void){
 ProfileCache*q=&profile_cache;int ok=1;
 if(q->ready){MKL_INT maxf=1,mnum=1,mt=11,phase=-1,n=q->n,idum=0,nrhs=1,msg=0,err=0;double dummy=0;
  pardiso(q->pt,&maxf,&mnum,&mt,&phase,&n,q->A,q->ia,q->ja,&idum,&nrhs,q->ip,&msg,&dummy,&dummy,&err);ok=err==0;}
 free(q->ia);free(q->ja);free(q->A);memset(q,0,sizeof(*q));return ok;
}

static void gfs_lu_basis_clear(void);
static void profile_release_at_exit(void){if(!profile_clear())fprintf(stderr,"[MATERIAL_VOLUME_FACTOR_RELEASE_FAILED]\n");gfs_lu_basis_clear();}
static int profile_exit_registered;
typedef struct {int col;double value;} GfsVolumeEntry;
typedef struct {double relative_residual;int entries,perturbed_pivots,error;} GfsVolumeLinearStats;
static inline int gfs_volume_entry_order(const void*a,const void*b){
 int x=((const GfsVolumeEntry*)a)->col,y=((const GfsVolumeEntry*)b)->col;
 return (x>y)-(x<y);
}
#include "gfs_volume_lu_precondition.h"
#include "gfs_volume_anchor_graph.h"
static inline int gfs_material_volume_linear_correction(int n,const int*offset,
 const int*owner,const double*coefficient,const double*rhs,double*delta,
 GfsVolumeLinearStats*stats){
 gfs_lu_precondition_enabled=gfs_lu_precondition_requested();
 profile_reuse=getenv("GFS_MATERIAL_VOLUME_SYMBOLIC_REUSE")&&atoi(getenv("GFS_MATERIAL_VOLUME_SYMBOLIC_REUSE"));
 if(profile_reuse&&!profile_exit_registered){if(atexit(profile_release_at_exit))return 0;profile_exit_registered=1;}
 double profile_start=profiler_clock(),assembly_s=0,analysis_s=0,numeric_s=0,solve_s=0,check_s=0,release_s=0,stamp=0;int same=-1,solver_mode=0;
 if(!stats)return 0;memset(stats,0,sizeof(*stats));
 if(n<1||!offset||!owner||!coefficient||!rhs||!delta||offset[0]!=0)return 0;
 int nf=offset[n],largest=0;
 if(nf<0||nf>INT_MAX-n)return 0;
 for(int i=0;i<n;i++){
  if(offset[i+1]<offset[i]||!isfinite(rhs[i]))return 0;
  int count=offset[i+1]-offset[i];if(count>largest)largest=count;
 }
 MKL_INT*ia=malloc(((size_t)n+1)*sizeof(*ia));
 MKL_INT*ja=malloc(((size_t)nf+n)*sizeof(*ja));
 double*A=malloc(((size_t)nf+n)*sizeof(*A));
 GfsVolumeEntry*row=malloc(((size_t)largest+1)*sizeof(*row));
 int ok=0,nz=0;
 MKL_INT maxf=1,mnum=1,mt=11,phase=13,nn=n,idum=0,nrhs=1,msglvl=0,err=0;
 if(!ia||!ja||!A||!row)goto done;
 for(int i=0;i<n;i++){
  ia[i]=nz;double diagonal=0;int count=0;
  for(int k=offset[i];k<offset[i+1];k++){
   int j=owner[k];double c=coefficient[k];
   if(j < -1||j>=n||!(c>=0)||!isfinite(c))goto done;
   /* Cancel self-image dependence symbolically, not as large +/- entries. */
   if(j==i||c==0)continue;
   diagonal+=c;if(j>=0)row[count++]=(GfsVolumeEntry){j,-c};
  }
  if(!(diagonal>0)||!isfinite(diagonal))goto done;
  row[count++]=(GfsVolumeEntry){i,diagonal};
  qsort(row,count,sizeof(*row),gfs_volume_entry_order);
  for(int k=0;k<count;){
   int j=row[k].col;double value=0;
   do{value+=row[k++].value;}while(k<count&&row[k].col==j);
   if(value!=0){ja[nz]=j;A[nz++]=value;}
  }
 }
 ia[n]=nz;stats->entries=nz;

 assembly_s=profiler_clock()-profile_start;
 ProfileCache*q=&profile_cache;
 same=profile_reuse&&q->ready&&q->n==n&&q->nz==nz&&
  !memcmp(q->ia,ia,((size_t)n+1)*sizeof(*ia))&&!memcmp(q->ja,ja,(size_t)nz*sizeof(*ja));
 if(gfs_lu_precondition_enabled){
  int reached=0;stamp=profiler_clock();
  int anchored=gfs_volume_anchor_graph(n,offset,owner,coefficient,&reached);
  check_s+=profiler_clock()-stamp;
  fprintf(stderr,"[MATERIAL_VOLUME_ANCHOR_GRAPH] rows=%d reached=%d certified=%d same_pattern=%d\n",n,reached,anchored,same);
  if(anchored!=1)goto done;
 }
 /* The previous LU remains a preconditioner for M. Its CSR is never
  * overwritten unless a fresh CURRENT J factorization is performed.
  * ID/anchor basis changes have already invalidated M in prepare_basis. */
 if(gfs_lu_precondition_enabled&&profile_reuse&&q->ready&&q->n==n){
  stamp=profiler_clock();int reused=gfs_lu_right_gmres(q,n,offset,owner,coefficient,rhs,delta,stats);
  solve_s+=profiler_clock()-stamp;
  solver_mode=reused==1?1:2;
  if(reused==1){ok=1;goto done;}
  if(reused<0){if(!profile_clear())goto done;same=0;}
 }
 if(!same){
  stamp=profiler_clock();if(!profile_clear())goto done;release_s+=profiler_clock()-stamp;
  q->n=n;q->nz=nz;q->ia=ia;q->ja=ja;q->A=A;ia=ja=NULL;A=NULL;
  pardisoinit(q->pt,&mt,q->ip);q->ip[34]=1;q->ip[7]=2;q->ready=1;
  phase=11;stamp=profiler_clock();pardiso(q->pt,&maxf,&mnum,&mt,&phase,&nn,q->A,q->ia,q->ja,&idum,&nrhs,q->ip,&msglvl,(void*)rhs,delta,&err);analysis_s=profiler_clock()-stamp;
  if(err){stats->error=err;goto done;}
 }else memcpy(q->A,A,(size_t)nz*sizeof(*A));
 phase=22;stamp=profiler_clock();pardiso(q->pt,&maxf,&mnum,&mt,&phase,&nn,q->A,q->ia,q->ja,&idum,&nrhs,q->ip,&msglvl,(void*)rhs,delta,&err);numeric_s=profiler_clock()-stamp;
 if(err){stats->error=err;goto done;}
 phase=33;stamp=profiler_clock();pardiso(q->pt,&maxf,&mnum,&mt,&phase,&nn,q->A,q->ia,q->ja,&idum,&nrhs,q->ip,&msglvl,(void*)rhs,delta,&err);solve_s+=profiler_clock()-stamp;
 stats->error=err;stats->perturbed_pivots=q->ip[13];
 if(err||q->ip[13])goto done;
 stamp=profiler_clock();

 {
  long double scale=0,residual=0;
  for(int i=0;i<n;i++){
   if(!isfinite(delta[i]))goto done;
   long double value=0;
   /* Check the physical difference form, independently of sparse assembly. */
   for(int k=offset[i];k<offset[i+1];k++){
    int j=owner[k];if(j==i)continue;
    value+=(long double)coefficient[k]*((long double)delta[i]-(j>=0?delta[j]:0));
   }
   long double r=fabsl(value-rhs[i]);if(r>residual)residual=r;
   long double b=fabsl(rhs[i]);if(b>scale)scale=b;
  }
  stats->relative_residual=scale>0?(double)(residual/scale):(double)residual;
  ok=isfinite(stats->relative_residual)&&stats->relative_residual<=1e-10;
 }

 check_s+=profiler_clock()-stamp;
done:
 if(!profile_reuse||!ok){stamp=profiler_clock();if(!profile_clear())ok=0;release_s+=profiler_clock()-stamp;}
 fprintf(stderr,"[MATERIAL_VOLUME_FACTOR] reuse=%d same_pattern=%d solver_mode=%d assembly=%.9g analysis=%.9g numeric=%.9g solve=%.9g certificate=%.9g release=%.9g total=%.9g residual=%.17g pivots=%d error=%d\n",profile_reuse,same,solver_mode,assembly_s,analysis_s,numeric_s,solve_s,check_s,release_s,profiler_clock()-profile_start,stats->relative_residual,stats->perturbed_pivots,stats->error);
 free(row);free(A);free(ja);free(ia);return ok;
}
#endif

