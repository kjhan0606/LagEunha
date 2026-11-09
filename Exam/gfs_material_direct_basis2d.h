#ifndef GFS_MATERIAL_DIRECT_BASIS2D_H
#define GFS_MATERIAL_DIRECT_BASIS2D_H
#include <mkl_pardiso.h>
#include <mkl_service.h>
#include <time.h>
#if defined(GFS_MATERIAL_DIRECT_THREADED) && defined(_OPENMP)
#include <omp.h>
#endif
static inline double gfs_direct_clock2d(void){struct timespec t;clock_gettime(CLOCK_MONOTONIC,&t);return t.tv_sec+1e-9*t.tv_nsec;}
/* Optional native CPU backend. One immutable geometry, one batch, no retained
 * factor across stages. No regularization; certify on original directed L. */
typedef struct {int i,j;double a;} GfsDirectEntry;
static int gfs_direct_entry_cmp(const void*a,const void*b){const GfsDirectEntry*x=a,*y=b;return x->i!=y->i?(x->i>y->i?1:-1):(x->j>y->j)-(x->j<y->j);}
typedef struct {void*pt[64];MKL_INT ip[64],n,nz;int started,ready;double*A;MKL_INT*ia,*ja;} GfsMaterialDirectFactor2D;
static inline void gfs_material_direct_factor_clear2d(GfsMaterialDirectFactor2D*f){
 if(!f)return;
 if(f->started){MKL_INT maxf=1,mnum=1,mt=2,phase=-1,idum=0,nrhs=1,msglvl=0,err=0;double dummy=0;
  pardiso(f->pt,&maxf,&mnum,&mt,&phase,&f->n,f->A,f->ia,f->ja,&idum,&nrhs,f->ip,&msglvl,&dummy,&dummy,&err);}
 free(f->A);free(f->ia);free(f->ja);memset(f,0,sizeof(*f));
}
static inline int gfs_material_direct_basis_cached2d(int n,int nf,const GfsTraceDirected2D*e,const unsigned char*wall,int nr,const GfsMaterialActiveFace2D*active,double*out,GfsMaterialDirectFactor2D*retained){
 if(n<2||nf<1||nr<1||!e||!active||!out)return 0;
 int profile=getenv("GFS_MATERIAL_ACTIVE_DETAIL")!=NULL;double ticks[4]={0},last=profile?gfs_direct_clock2d():0;int phase_index=0;
 int ok=0,m=n-1,nz=0;GfsMaterialDirectFactor2D local={0},*factor=retained?retained:&local;MKL_INT mt=2,err=0,maxf=1,mnum=1,phase=12,nn=m,nrhs=nr,msglvl=0,idum=0;
 /* Bounded one-matrix cache, opt-in. Numerical factors require exact CSR/value
  * equality; symbolic analysis alone may be reused for identical CSR structure.
  * Changed values are always numerically re-factored, and every solution is
  * certified on the current directed operator. Called serially outside OpenMP. */
 static GfsMaterialDirectFactor2D across_rhs={0};
 const char *cross_env=getenv("GFS_MATERIAL_FACTOR_CROSS_RHS");
 int cross_rhs=cross_env&&!strcmp(cross_env,"1");
 if(cross_rhs)factor=&across_rhs;
 double*coefficient=calloc((size_t)nf,sizeof(double)),*diag=calloc(n,sizeof(double)),*A=malloc((size_t)(nf+n)*sizeof(double)),*B=calloc((size_t)nr*m,sizeof(double)),*X=calloc((size_t)nr*m,sizeof(double)),*res=calloc(n,sizeof(double));
 GfsDirectEntry*en=malloc((size_t)(nf+n)*sizeof(*en));MKL_INT*ia=calloc(n,sizeof(*ia)),*ja=malloc((size_t)(nf+n)*sizeof(*ja));
 int solve_threads=1;
#if defined(GFS_MATERIAL_DIRECT_THREADED) && defined(_OPENMP)
 solve_threads=omp_get_max_threads();if(solve_threads<1)solve_threads=1;
#endif
 int old_threads=mkl_set_num_threads_local(solve_threads);
 if(!coefficient||!diag||!A||!B||!X||!res||!en||!ia||!ja)goto done;
 for(int k=0;k<nf;k++){
  int i=e[k].i,j=e[k].j;if(i<0||i>=n||j<0||j>=n)goto done;
  double d=hypot(e[k].dx,e[k].dy);if(!(d>0)||!isfinite(d))goto done;
  if(wall&&wall[k]){if(wall[k]!=1||!gfs_material_stationary_wall2d(e+k))goto done;continue;}
  double c=e[k].area/(2*d);if(!(c>0)||!isfinite(c)||i==j)goto done;coefficient[k]=c;diag[i]+=c;
  if(i>0&&j>i)en[nz++]=(GfsDirectEntry){i-1,j-1,-c};
 }
 for(int i=1;i<n;i++){if(!(diag[i]>0))goto done;en[nz++]=(GfsDirectEntry){i-1,i-1,diag[i]};}
 qsort(en,nz,sizeof(*en),gfs_direct_entry_cmp);int unique=0;
 for(int k=0;k<nz;k++){if(unique&&en[k].i==en[unique-1].i&&en[k].j==en[unique-1].j)en[unique-1].a+=en[k].a;else en[unique++]=en[k];}
 for(int k=0;k<unique;k++){ia[en[k].i+1]++;ja[k]=en[k].j;A[k]=en[k].a;}for(int i=1;i<=m;i++)ia[i]+=ia[i-1];
 for(int k=0;k<nr;k++){int i=active[k].i,j=active[k].j;double s=active[k].sign;if(i<0||i>=n||j<0||j>=n||i==j||fabs(s)!=1)goto done;if(i)B[k*m+i-1]=s;if(j)B[k*m+j-1]=-s;}
 int same=factor->ready&&factor->n==m&&factor->nz==unique&&
  memcmp(factor->ia,ia,(size_t)(m+1)*sizeof(*ia))==0&&memcmp(factor->ja,ja,(size_t)unique*sizeof(*ja))==0&&memcmp(factor->A,A,(size_t)unique*sizeof(*A))==0;
 if(profile){double now=gfs_direct_clock2d();ticks[0]=now-last;last=now;}phase_index=1;
 int same_pattern=factor->ready&&factor->n==m&&factor->nz==unique&&
  memcmp(factor->ia,ia,(size_t)(m+1)*sizeof(*ia))==0&&memcmp(factor->ja,ja,(size_t)unique*sizeof(*ja))==0;
 if(!same){
  if(same_pattern){
   /* Reuse symbolic analysis only. Re-factor every changed value on the
    * current matrix; retain the original residual/pivot certification. */
   memcpy(factor->A,A,(size_t)unique*sizeof(*A));phase=22;
  }else{
   gfs_material_direct_factor_clear2d(factor);factor->n=m;factor->nz=unique;factor->A=A;factor->ia=ia;factor->ja=ja;A=NULL;ia=ja=NULL;
   pardisoinit(factor->pt,&mt,factor->ip);factor->ip[34]=1;factor->ip[7]=2;factor->started=1;
  }
  factor->ready=0;
  pardiso(factor->pt,&maxf,&mnum,&mt,&phase,&nn,factor->A,factor->ia,factor->ja,&idum,&nrhs,factor->ip,&msglvl,B,X,&err);if(err||factor->ip[13])goto done;
  factor->ready=1;
 }
 if(profile)fprintf(stderr,"[MATERIAL_SYMBOLIC_REUSE] reused=%d numeric_refactor=%d\n",same_pattern,!same);
 if(getenv("GFS_MATERIAL_ACTIVE_DETAIL"))fprintf(stderr,"[MATERIAL_FACTOR_REUSE] reused=%d rows=%d rhs=%d\n",same,m,nr);
 if(profile){double now=gfs_direct_clock2d();ticks[1]=now-last;last=now;}phase_index=2;
 phase=33;pardiso(factor->pt,&maxf,&mnum,&mt,&phase,&nn,factor->A,factor->ia,factor->ja,&idum,&nrhs,factor->ip,&msglvl,B,X,&err);if(err)goto done;
 if(profile){double now=gfs_direct_clock2d();ticks[2]=now-last;last=now;}phase_index=3;
 for(int k=0;k<nr;k++){
  double*q=out+(size_t)k*n,mean=0;q[0]=0;for(int i=1;i<n;i++){q[i]=X[k*m+i-1];if(!isfinite(q[i]))goto done;mean+=q[i]/n;}for(int i=0;i<n;i++)q[i]-=mean;
  memset(res,0,(size_t)n*sizeof(double));res[active[k].i]=-active[k].sign;res[active[k].j]=active[k].sign;
  for(int z=0;z<nf;z++)if(!wall||!wall[z])res[e[z].i]+=coefficient[z]*(q[e[z].i]-q[e[z].j]);
  double r2=0;for(int i=0;i<n;i++)r2+=res[i]*res[i];if(!isfinite(r2)||sqrt(r2/2)>1e-10)goto done;
 }
 ok=1;
 done:
 if(profile){ticks[phase_index]+=gfs_direct_clock2d()-last;fprintf(stderr,"[MATERIAL_FACTOR_DETAIL] rhs=%d ok=%d assembly_ms=%.9g factor_ms=%.9g triangular_ms=%.9g certificate_ms=%.9g\n",nr,ok,1000*ticks[0],1000*ticks[1],1000*ticks[2],1000*ticks[3]);}
 if(!ok||(!retained&&!cross_rhs))gfs_material_direct_factor_clear2d(factor);
 mkl_set_num_threads_local(old_threads);free(coefficient);free(diag);free(A);free(B);free(X);free(res);free(en);free(ia);free(ja);return ok;
}
static inline int gfs_material_direct_basis2d(int n,int nf,const GfsTraceDirected2D*e,const unsigned char*wall,int nr,const GfsMaterialActiveFace2D*active,double*out){return gfs_material_direct_basis_cached2d(n,nf,e,wall,nr,active,out,NULL);}
#endif
