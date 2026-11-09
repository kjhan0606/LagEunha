/* MPI ownership adapter for the existing globally coupled material solve.
 * Persistent particles are decomposed by the user's k-section DD. A
 * replicated operator view supplies the global constraint graph; expensive
 * GPU cell tessellation is partitioned by those owners and assembled without
 * floating-point reduction. Rank zero solves the global constraints and
 * broadcasts the resulting rates; memory scaling remains O(N) per rank.
 * No migration occurs inside a rejected SSP trial. Accepted particles are
 * returned to their previous owners and pmigrate applies the DD pivots once.
 */
#include <limits.h>
#include "gfs_material_mpi.h"
typedef struct {unsigned long long id;int rank;} GfsMaterialMpiOwner;
static MPI_Comm material_mpi_comm=MPI_COMM_NULL;
static int material_mpi_id=-1,material_mpi_count;
static GfsMaterialMpiOwner *material_mpi_owners;
static int material_mpi_owner_order(const void*a,const void*b){
 unsigned long long x=((const GfsMaterialMpiOwner*)a)->id,y=((const GfsMaterialMpiOwner*)b)->id;return(x>y)-(x<y);
}
static int material_mpi_particle_order(const void*a,const void*b){
 unsigned long long x=PINDX((const treevorostressrk4particletype*)a),y=PINDX((const treevorostressrk4particletype*)b);return(x>y)-(x<y);
}
int gfs_material_mpi_rank(void){return material_mpi_id;}
int gfs_material_mpi_rows(int n,const long long*ids,unsigned char*owned){
 if(material_mpi_id<0)return 0;
 for(int i=0;i<n;i++){
  GfsMaterialMpiOwner key={(unsigned long long)ids[i],0};
  GfsMaterialMpiOwner *o=bsearch(&key,material_mpi_owners,material_mpi_count,sizeof(key),material_mpi_owner_order);
  if(!o)MPI_Abort(material_mpi_comm,99);
  owned[i]=o->rank==material_mpi_id;
 }
 return 1;
}
void gfs_material_mpi_bcast(void*data,size_t bytes){
 if(material_mpi_id<0)return;
 for(size_t at=0;at<bytes;){int count=bytes-at>INT_MAX?INT_MAX:(int)(bytes-at);
  MPI_Bcast((unsigned char*)data+at,count,MPI_BYTE,0,material_mpi_comm);at+=count;
 }
}
void gfs_material_mpi_merge(void*rows,size_t stride,int n,const unsigned char*owned){
 if(material_mpi_id<0)return;
 unsigned char *p=rows;size_t bytes=(size_t)n*stride;
 for(int i=0;i<n;i++)if(!owned[i])memset(p+(size_t)i*stride,0,stride);
 /* Exactly one nonzero byte contributor per row, including signed zero.
  * Chunk the byte count rather than narrowing a size_t to an MPI int. */
 for(size_t at=0;at<bytes;){int count=bytes-at>INT_MAX?INT_MAX:(int)(bytes-at);
  MPI_Allreduce(MPI_IN_PLACE,p+at,count,MPI_UNSIGNED_CHAR,MPI_BOR,material_mpi_comm);at+=count;
 }
}
void gfs_material_mpi_counters(unsigned long long*v,int n){
 if(material_mpi_id>=0)MPI_Allreduce(MPI_IN_PLACE,v,n,MPI_UNSIGNED_LONG_LONG,MPI_SUM,material_mpi_comm);
}
void gfs_material_mpi_maxima(double*v,int n){
 if(material_mpi_id>=0)MPI_Allreduce(MPI_IN_PLACE,v,n,MPI_DOUBLE,MPI_MAX,material_mpi_comm);
}
static void material_mpi_self_dd(DoDeInfo*d,const DoDeInfo*old,SimParameters*s){
 *d=*old;d->com=MPI_COMM_SELF;d->myid=0;d->nid=1;d->npivot=0;d->pivot=NULL;
 d->nsubgroup=d->subgroupsize=1;d->subgroupid=0;
 d->lgroup.r.rmin[0]=SIMBOX(s).x.min;d->lgroup.r.rmax[0]=SIMBOX(s).x.max;
 d->lgroup.r.rmin[1]=SIMBOX(s).y.min;d->lgroup.r.rmax[1]=SIMBOX(s).y.max;
}
typedef struct {
 SimParameters *original,view;
 DoDeInfo base_dd,tree_dd;
} GfsMaterialMpiView;
void *gfs_material_mpi_view_begin(void*opaque,void**ticket){
 SimParameters*s=opaque;*ticket=NULL;
 if(NID(s)==1)return s;
 if(material_mpi_id>=0||!GAS_GPU_ENABLED(s)||TVORORK4_DDINFO(s)[0].n_size!=sizeof(treevorostressrk4particletype))MPI_Abort(MPI_COMM(s),99);
 /* Existing RK4 example: mpirks builds pivots for this particle payload;
  * pmigrate below uses those very pivots. Do not replace this by stripes. */
 static SimParameters *decomposed=NULL;
 if(decomposed!=s){
  ptrdiff_t minimum=VORO_NP(s);MPI_Allreduce(MPI_IN_PLACE,&minimum,1,MPI_PTRDIFF_T,MPI_MIN,MPI_COMM(s));
  /* The legacy kth-pivot estimator requires a sample on each rank. Empty
   * initial domains retain the valid initial k-section pivots; no particle
   * is fabricated or transferred by a different decomposition scheme. */
  if(minimum>0)mpirks((void**)&VORORK4_TBP(s),&VORO_NP(s),sizeof(treevorostressrk4particletype),
    &TVORORK4_DDFUNC(s),TVORORK4_DDINFO(s),MPI_COMM(s),&GRIDINFO(s),oldway,OldCom);
  decomposed=s;
 }
 int ranks=NID(s),rank=MYID(s);ptrdiff_t local=VORO_NP(s);
 if(local<0||local>INT_MAX/(int)sizeof(treevorostressrk4particletype))MPI_Abort(MPI_COMM(s),99);
 int *counts=malloc((size_t)ranks*sizeof(int)),*offset=malloc((size_t)ranks*sizeof(int));
 if(!counts||!offset)MPI_Abort(MPI_COMM(s),99);
 int bytes=(int)(local*sizeof(treevorostressrk4particletype));MPI_Allgather(&bytes,1,MPI_INT,counts,1,MPI_INT,MPI_COMM(s));
 size_t totalbytes=0;for(int r=0;r<ranks;r++){offset[r]=(int)totalbytes;totalbytes+=counts[r];if(totalbytes>INT_MAX)MPI_Abort(MPI_COMM(s),99);}
 int n=totalbytes/sizeof(treevorostressrk4particletype);if(n<1)MPI_Abort(MPI_COMM(s),99);
 treevorostressrk4particletype *all=my_malloc(totalbytes);
 MPI_Allgatherv(VORORK4_TBP(s),bytes,MPI_BYTE,all,counts,offset,MPI_BYTE,MPI_COMM(s));
 material_mpi_owners=malloc((size_t)n*sizeof(*material_mpi_owners));if(!material_mpi_owners)MPI_Abort(MPI_COMM(s),99);
 for(int r=0;r<ranks;r++)for(size_t k=offset[r]/sizeof(*all);k<(offset[r]+counts[r])/sizeof(*all);k++)material_mpi_owners[k]=(GfsMaterialMpiOwner){PINDX(all+k),r};
 qsort(material_mpi_owners,n,sizeof(*material_mpi_owners),material_mpi_owner_order);
 qsort(all,n,sizeof(*all),material_mpi_particle_order);
 for(int i=1;i<n;i++)if(PINDX(all+i)==PINDX(all+i-1))MPI_Abort(MPI_COMM(s),99);
 GfsMaterialMpiView *v=calloc(1,sizeof(*v));if(!v)MPI_Abort(MPI_COMM(s),99);
 v->original=s;v->view=*s;SimParameters*view=&v->view;
 material_mpi_self_dd(&v->base_dd,VORORK4_DDINFO(s),s);material_mpi_self_dd(&v->tree_dd,TVORORK4_DDINFO(s),s);
 MYID(view)=0;NID(view)=1;MPI_COMM(view)=MPI_COMM_SELF;NDDINFO(view)=1;
 VORORK4_DDINFO(view)=&v->base_dd;TVORORK4_DDINFO(view)=&v->tree_dd;
 VORO_NP(view)=VORO_TNP(view)=NPSUM(view)=n;VORORK4_TBP(view)=(void*)all;
 VORORK4_TBPP(view)=NULL;VORO_NPAD(view)=0;VORO_BASICCELL(view)=NULL;
 Xmin_HydroExam(view)=SIMBOX(s).x.min;Xmax_HydroExam(view)=SIMBOX(s).x.max;
 Ymin_HydroExam(view)=SIMBOX(s).y.min;Ymax_HydroExam(view)=SIMBOX(s).y.max;
 material_mpi_comm=MPI_COMM(s);material_mpi_id=rank;material_mpi_count=n;
 fprintf(stderr,"[MATERIAL_KSECTION] rank=%d ranks=%d owned=%td total=%d global_graph_replica=1 constraint_root=0\n",rank,ranks,local,n);
 free(counts);free(offset);*ticket=v;return view;
}
void gfs_material_mpi_view_end(void*ticket,int allow_boundary_events){
 if(!ticket)return;
 GfsMaterialMpiView*v=ticket;SimParameters*s=v->original,*view=&v->view;
 int n=VORO_NP(view),rank=MYID(s);
 if(n<1||(!allow_boundary_events&&n!=material_mpi_count))MPI_Abort(MPI_COMM(s),99);
 treevorostressrk4particletype*all=(void*)VORORK4_TBP(view);
 qsort(all,n,sizeof(*all),material_mpi_particle_order);
 /* Different GPU/solver results must never become silently different copies
  * of one accepted global state. Check physical values, not padding/pointers. */
 uint64_t hash=UINT64_C(14695981039346656037);
 for(int i=0;i<n;i++){
  if(i&&PINDX(all+i)==PINDX(all+i-1))MPI_Abort(MPI_COMM(s),99);
  uint64_t id=PINDX(all+i);
  const unsigned char*ib=(const unsigned char*)&id;
  for(size_t k=0;k<sizeof(id);k++){hash^=ib[k];hash*=UINT64_C(1099511628211);}
  double q[8]={all[i].x,all[i].y,all[i].vx,all[i].vy,all[i].mass,
   all[i].weight_thermal.U,all[i].weight_thermal.V,all[i].weight.W};
  const unsigned char*b=(const unsigned char*)q;
  for(size_t k=0;k<sizeof(q);k++){hash^=b[k];hash*=UINT64_C(1099511628211);}
 }
 uint64_t lo=hash,hi=hash;MPI_Allreduce(MPI_IN_PLACE,&lo,1,MPI_UINT64_T,MPI_MIN,MPI_COMM(s));MPI_Allreduce(MPI_IN_PLACE,&hi,1,MPI_UINT64_T,MPI_MAX,MPI_COMM(s));
 if(lo!=hi){fprintf(stderr,"[MATERIAL_MPI] replicated constraint state mismatch\n");MPI_Abort(MPI_COMM(s),99);}
 ptrdiff_t selected=0;
 for(int i=0;i<n;i++){
  GfsMaterialMpiOwner key={PINDX(all+i),0};GfsMaterialMpiOwner*o=bsearch(&key,material_mpi_owners,material_mpi_count,sizeof(key),material_mpi_owner_order);
  if(!o&&!allow_boundary_events)MPI_Abort(MPI_COMM(s),99);
  if((o?o->rank:0)==rank)selected++;
 }
 treevorostressrk4particletype *destination=my_malloc((size_t)(selected?selected:1)*sizeof(*destination));ptrdiff_t k=0;
 for(int i=0;i<n;i++){
  GfsMaterialMpiOwner key={PINDX(all+i),0};GfsMaterialMpiOwner*o=bsearch(&key,material_mpi_owners,material_mpi_count,sizeof(key),material_mpi_owner_order);
  if((o?o->rank:0)==rank)destination[k++]=all[i];
 }
 my_free(VORORK4_TBP(s));VORORK4_TBP(s)=(void*)destination;VORO_NP(s)=selected;
 GAS_dtold(s)=GAS_dtold(view);BASICCELL_CELLWIDTH(s)=BASICCELL_CELLWIDTH(view);
 if(VORORK4_TBPP(view))my_free(VORORK4_TBPP(view));my_free(all);
 gpu_get_context()->tess_faces_on_device=0;
 material_mpi_id=-1;material_mpi_comm=MPI_COMM_NULL;material_mpi_count=0;
 free(material_mpi_owners);material_mpi_owners=NULL;free(v);
 migrateTreeVorork4Particles(s);
}
static double gfs_material_ssprk2_mpi(SimParameters*s,
 void(*padding)(SimParameters*,postype),Voro2D_point*(*neighbors)(SimParameters*,int,int,int*),
 treevorork4particletype*(*cells)(SimParameters*,int,int,int*),
 void(*linked)(SimParameters*,postype,postype,postype,postype,postype,void(*)(SimParameters*,postype)),
 void(*post)(SimParameters*)){
 void*ticket=NULL;SimParameters*view=gfs_material_mpi_view_begin(s,&ticket);
 double h=gfs_material_ssprk2_native(view,padding,neighbors,cells,linked,post);
 if(ticket){double bounds[2]={h,-h};MPI_Allreduce(MPI_IN_PLACE,bounds,2,MPI_DOUBLE,MPI_MIN,MPI_COMM(s));if(bounds[0]!=-bounds[1])MPI_Abort(MPI_COMM(s),99);}
 gfs_material_mpi_view_end(ticket,0);return h;
}
