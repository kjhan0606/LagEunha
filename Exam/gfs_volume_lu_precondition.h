#ifndef GFS_LU_PRECONDITIONER_H
#define GFS_LU_PRECONDITIONER_H
/* Offline qualification. Old LU is only a right preconditioner. The
 * current raw-face J always supplies matvecs and the acceptance residual. */
static int gfs_lu_precondition_enabled;
static double gfs_lu_norm(int n,const double*x){
 long double v=0;for(int i=0;i<n;i++)v+=(long double)x[i]*x[i];return(double)sqrtl(v);
}
static double gfs_lu_dot(int n,const double*x,const double*y){
 long double v=0;for(int i=0;i<n;i++)v+=(long double)x[i]*y[i];return(double)v;
}
static int gfs_lu_current_apply(int n,const int*offset,const int*owner,const double*c,
 const double*x,double*y){
 for(int i=0;i<n;i++){
  long double v=0;
  for(int k=offset[i];k<offset[i+1];k++){
   int j=owner[k];if(j==i)continue;
   v+=(long double)c[k]*((long double)x[i]-(j>=0?x[j]:0));
  }
  y[i]=(double)v;if(!isfinite(y[i]))return 0;
 }
 return 1;
}
static double gfs_lu_true_residual(int n,const int*offset,const int*owner,const double*c,
 const double*b,const double*x){
 long double scale=0,residual=0;
 for(int i=0;i<n;i++){
  if(!isfinite(x[i]))return HUGE_VAL;
  long double v=0;
  for(int k=offset[i];k<offset[i+1];k++){
   int j=owner[k];if(j==i)continue;
   v+=(long double)c[k]*((long double)x[i]-(j>=0?x[j]:0));
  }
  long double r=fabsl(v-b[i]);if(r>residual)residual=r;
  if(fabsl(b[i])>scale)scale=fabsl(b[i]);
 }
 return scale>0?(double)(residual/scale):(double)residual;
}
static int gfs_lu_old_apply(ProfileCache*q,const double*b,double*x){
 MKL_INT maxf=1,mnum=1,mt=11,phase=33,n=q->n,idum=0,nrhs=1,msg=0,err=0;
 /* q->A, ia, ja, permutations and factors all remain those of OLD M.
  * Internal iterative refinement therefore also acts on M, not current J. */
 pardiso(q->pt,&maxf,&mnum,&mt,&phase,&n,q->A,q->ia,q->ja,&idum,&nrhs,q->ip,&msg,
         (void*)b,x,&err);
 return !err&&!q->ip[13];
}
static int gfs_lu_right_gmres(ProfileCache*q,int n,const int*offset,const int*owner,
 const double*c,const double*b,double*x,GfsVolumeLinearStats*stats){
 const int limit=n<8?n:8; /* Computational work cap; never an acceptance tolerance. */
 double *v=calloc((size_t)(limit+1)*n,sizeof(double));
 double *z=calloc((size_t)limit*n,sizeof(double)),*w=malloc((size_t)n*sizeof(double));
 double h[72]={0},cs[8]={0},sn[8]={0},g[9]={0},y[8]={0};
 int answer=0,iterations=0;double last=HUGE_VAL;
 if(!v||!z||!w)goto done;
 double beta=gfs_lu_norm(n,b);if(!isfinite(beta))goto done;
 if(beta==0){memset(x,0,(size_t)n*sizeof(double));last=gfs_lu_true_residual(n,offset,owner,c,b,x);answer=last<=1e-10;goto done;}
 g[0]=beta;for(int i=0;i<n;i++)v[i]=b[i]/beta;
 for(int j=0;j<limit;j++){
  double*zj=z+(size_t)j*n;
  if(!gfs_lu_old_apply(q,v+(size_t)j*n,zj)){answer=-1;goto done;}
  if(!gfs_lu_current_apply(n,offset,owner,c,zj,w))goto done;
  /* Reorthogonalize; do not use a preconditioned residual as the test. */
  for(int pass=0;pass<2;pass++)for(int k=0;k<=j;k++){
   double*vk=v+(size_t)k*n,dot=gfs_lu_dot(n,vk,w);if(!isfinite(dot))goto done;
   h[k*8+j]+=dot;for(int i=0;i<n;i++)w[i]-=dot*vk[i];
  }
  double next=gfs_lu_norm(n,w);if(!isfinite(next))goto done;h[(j+1)*8+j]=next;
  if(next>0)for(int i=0;i<n;i++)v[(size_t)(j+1)*n+i]=w[i]/next;
  for(int k=0;k<j;k++){
   double a=h[k*8+j],d=h[(k+1)*8+j];h[k*8+j]=cs[k]*a+sn[k]*d;h[(k+1)*8+j]=-sn[k]*a+cs[k]*d;
  }
  double a=h[j*8+j],d=h[(j+1)*8+j],norm=hypot(a,d);if(!(norm>0)||!isfinite(norm))goto done;
  cs[j]=a/norm;sn[j]=d/norm;h[j*8+j]=norm;h[(j+1)*8+j]=0;
  g[j+1]=-sn[j]*g[j];g[j]=cs[j]*g[j];
  for(int k=j;k>=0;k--){
   long double value=g[k];for(int l=k+1;l<=j;l++)value-=(long double)h[k*8+l]*y[l];
   y[k]=(double)(value/h[k*8+k]);if(!isfinite(y[k]))goto done;
  }
  memset(x,0,(size_t)n*sizeof(double));
  for(int k=0;k<=j;k++)for(int i=0;i<n;i++)x[i]+=y[k]*z[(size_t)k*n+i];
  iterations=j+1;last=gfs_lu_true_residual(n,offset,owner,c,b,x);
  if(isfinite(last)&&last<=1e-10){answer=1;break;}
  if(next==0)break;
 }
done:
 if(answer==1){stats->relative_residual=last;stats->error=0;stats->perturbed_pivots=0;}
 fprintf(stderr,"[VOLUME_LU_PRECONDITION] iterations=%d converged=%d current_raw_residual=%.17g\n",iterations,answer,last);
 free(w);free(z);free(v);return answer;
}
/* Row/anchor identity belongs to the caller's physical basis. The numerical
 * cache has its own lifetime; factorization may clear it after preparation. */
typedef struct {
 int n,valid; long long *id; unsigned char *anchor;
} GfsVolumeLuBasis;
static GfsVolumeLuBasis gfs_volume_lu_basis;
static void gfs_lu_basis_clear(void) {
 free(gfs_volume_lu_basis.id);free(gfs_volume_lu_basis.anchor);
 memset(&gfs_volume_lu_basis,0,sizeof(gfs_volume_lu_basis));
}
static int gfs_lu_precondition_requested(void) {
 const char *v=getenv("GFS_MATERIAL_VOLUME_LU_PRECONDITION");
 return v && atoi(v)!=0;
}
static int gfs_lu_prepare_basis(int n,const long long *id,const int *off,
 const int *owner,const double *coefficient) {
 if(!gfs_lu_precondition_requested())return 1;
 if(n<1||!id||!off||!owner||!coefficient||off[0]!=0||off[n]<0)return 0;
 for(int i=0;i<n;i++)if(off[i+1]<off[i])return 0;
 unsigned char *anchor=calloc((size_t)n,1);if(!anchor)return 0;
 for(int i=0;i<n;i++)for(int k=off[i];k<off[i+1];k++) {
  if(owner[k]<-1||owner[k]>=n||!(coefficient[k]>=0)||!isfinite(coefficient[k])) {
   free(anchor);return 0;
  }
  if(owner[k]==-1&&coefficient[k]>0)anchor[i]=1;
 }
 GfsVolumeLuBasis *b=&gfs_volume_lu_basis;
 int same=b->valid&&b->n==n&&
  !memcmp(b->id,id,(size_t)n*sizeof(*id))&&!memcmp(b->anchor,anchor,(size_t)n);
 if(same){free(anchor);return 1;}
 long long *copy=malloc((size_t)n*sizeof(*copy));
 if(!copy){free(anchor);return 0;}
 memcpy(copy,id,(size_t)n*sizeof(*copy));
 int had_factors=profile_cache.ready;
 if(!profile_clear()){free(copy);free(anchor);return 0;}
 gfs_lu_basis_clear();b->n=n;b->id=copy;b->anchor=anchor;b->valid=1;
 fprintf(stderr,"[MATERIAL_VOLUME_LU_BASIS] rows=%d invalidated_factors=%d\n",n,had_factors);
 return 1;
}

#endif
