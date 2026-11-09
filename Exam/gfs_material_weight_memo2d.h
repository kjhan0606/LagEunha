#ifndef GFS_MATERIAL_WEIGHT_MEMO2D_H
#define GFS_MATERIAL_WEIGHT_MEMO2D_H
/* Exact pure-solve reuse, serialized adapter only. Full input bytes, no
   tolerance relaxation or stage-history assumption. */
static int gfs_material_weight_memo2d(int n,int nf,const GfsTraceDirected2D*e,
 const GfsFvPrimitive*c,const double*rhs,double*q,double*v,int*it,double*rel,const unsigned char*w){
 static unsigned char*saved;static double*values;static size_t saved_size;static int saved_it;static double saved_rel;
 int dims[3]={n,nf,w!=NULL};size_t ne=(size_t)nf*sizeof(*e),nc=(size_t)n*sizeof(*c),nr=(size_t)n*4*sizeof(double),nw=w?(size_t)nf:0;
 size_t size=sizeof(dims)+ne+nc+nr+nw;unsigned char*key=malloc(size),*p=key;
 if(!key)return gfs_material_weight_rate_boundaries2d(n,nf,e,c,rhs,q,v,it,rel,w);
 #define WM_COPY(x,z) do{memcpy(p,x,z);p+=z;}while(0)
 WM_COPY(dims,sizeof(dims));WM_COPY(e,ne);WM_COPY(c,nc);WM_COPY(rhs,nr);if(nw)WM_COPY(w,nw);
 #undef WM_COPY
 if(saved&&size==saved_size&&!memcmp(saved,key,size)){memcpy(q,values,(size_t)n*sizeof(double));memcpy(v,values+n,(size_t)n*sizeof(double));*it=saved_it;*rel=saved_rel;free(key);fprintf(stderr,"[WEIGHT_MEMO] hit=1\n");return 1;}
 int ok=gfs_material_weight_rate_boundaries2d(n,nf,e,c,rhs,q,v,it,rel,w);
 double*copy=ok?malloc((size_t)2*n*sizeof(double)):NULL;
 free(saved);free(values);saved=NULL;values=NULL;
 if(copy){memcpy(copy,q,(size_t)n*sizeof(double));memcpy(copy+n,v,(size_t)n*sizeof(double));saved=key;saved_size=size;values=copy;saved_it=*it;saved_rel=*rel;}else free(key);
 fprintf(stderr,"[WEIGHT_MEMO] hit=0\n");return ok;
}
#endif
