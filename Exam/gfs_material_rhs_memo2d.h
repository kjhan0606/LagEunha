#ifndef GFS_MATERIAL_RHS_MEMO2D_H
#define GFS_MATERIAL_RHS_MEMO2D_H
/* Exact-byte pure RHS memo. Padding can cause misses, never false hits.
   Caller must skip this when response exports or other side effects are needed.
   Single-thread adapter only. No geometry, constraint, or weight-solve reuse. */
static unsigned char *gfs_rhs_memo_key;
static double *gfs_rhs_memo_value;
static size_t gfs_rhs_memo_size;
static int gfs_rhs_memo_valid;
static int gfs_rhs_memo_prepare(const GfsTraceCsr2D*t,const GfsFvPrimitive*s,
 double gamma,const double*a,double*out,unsigned char**key,size_t*size){
 int counts[4]={t->cells,t->faces,t->entries,a!=NULL};
 size_t f=(size_t)t->faces*sizeof(*t->face),e=(size_t)t->entries*sizeof(*t->entry),
 c=(size_t)t->cells*sizeof(*s),al=a?(size_t)t->cells*3*sizeof(double):0;
 *size=sizeof(counts)+sizeof(gamma)+f+e+c+al;*key=malloc(*size);if(!*key)return 0;
 unsigned char*p=*key;
 #define GFS_MEMO_COPY(x,n) do{memcpy(p,x,n);p+=n;}while(0)
 GFS_MEMO_COPY(counts,sizeof(counts));GFS_MEMO_COPY(&gamma,sizeof(gamma));
 GFS_MEMO_COPY(t->face,f);GFS_MEMO_COPY(t->entry,e);GFS_MEMO_COPY(s,c);if(al)GFS_MEMO_COPY(a,al);
 #undef GFS_MEMO_COPY
 if(gfs_rhs_memo_valid&&*size==gfs_rhs_memo_size&&!memcmp(*key,gfs_rhs_memo_key,*size)){
 memcpy(out,gfs_rhs_memo_value,(size_t)t->cells*4*sizeof(double));return 1;}
 return 0;
}
static void gfs_rhs_memo_store(unsigned char*key,size_t size,int n,const double*out,int ok){
 if(!ok||!key){gfs_rhs_memo_valid=0;free(key);return;}
 double*v=malloc((size_t)n*4*sizeof(double));if(!v){free(key);gfs_rhs_memo_valid=0;return;}
 memcpy(v,out,(size_t)n*4*sizeof(double));free(gfs_rhs_memo_key);free(gfs_rhs_memo_value);
 gfs_rhs_memo_key=key;gfs_rhs_memo_value=v;gfs_rhs_memo_size=size;gfs_rhs_memo_valid=1;
}
#endif
