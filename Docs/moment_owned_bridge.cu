/* Native per-owner MPI stage; canonical flux is computed exactly once upstream.
 * Reuses stage_cell/cell and the validated factor ABI without ghost updates.
 * Build separately from the legacy fast-math translation unit.
 */
#include "moment_backend_bridge.cu"
#include "moment_owned_backend.h"

extern "C" int lagmoment_select_visible_device(int ordinal,int *count,char *name,unsigned char *uuid){
 if(!count||!name||!uuid)return 1;
 int n=0;if(cudaGetDeviceCount(&n)!=cudaSuccess||ordinal<0||ordinal>=n)return 2;
 if(cudaSetDevice(ordinal)!=cudaSuccess)return 3;
 cudaDeviceProp prop;if(cudaGetDeviceProperties(&prop,ordinal)!=cudaSuccess)return 4;
 *count=n;std::memcpy(name,prop.name,256);std::memcpy(uuid,prop.uuid.bytes,16);return 0;
}

/* Preserve the canonical addition order while transferring only ledger terms.
 * The old path copied/scanned every face node on the GPU in one thread even
 * though only physical boundary nodes contribute to the first four entries. */
__global__ void bridge_compact_ordered_ledger_kernel(int wall_nodes,int frozen_cells,int active_cells,
    const double *wall,const double *frozen,const double *active,double *result){
 if(threadIdx.x==0&&blockIdx.x==0){
  double value[12]={0.};
  for(int k=0;k<wall_nodes;k++){
   const double *f=wall+5*k;
   value[0]-=f[0];value[1]-=f[1];value[2]-=f[2];value[3]-=f[3]-f[4];
  }
  for(int i=0;i<frozen_cells;i++){
   const double *r=frozen+5*i;
   value[4]-=r[0];value[5]-=r[1];value[6]-=r[2];value[7]-=r[3]-r[4];
  }
  for(int i=0;i<active_cells;i++){
   const double *g=active+5*i;
   value[8]+=g[0];value[9]+=g[1];value[10]+=g[2];value[11]+=g[3]-g[4];
  }
  for(int k=0;k<12;k++)result[k]=value[k];
 }
}

extern "C" int lagmoment_canonical_nodes(int backend,int threads,int nn,
    const double *a,double *result){
 if((backend!=0&&backend!=1)||threads<1||nn<1||!a||!result)return 1;
 try{
  std::vector<double> nf(19*size_t(nn));
  if(backend==0){
   omp_set_num_threads(threads);
   #pragma omp parallel for
   for(int k=0;k<nn;k++)node(a+33*k,nf.data()+19*k);
  }else{
   Device<double> da,dnf;
   if(!da.set(33*size_t(nn),a)||!dnf.set(nf.size()))return 2;
   nodes_kernel<<<(nn+127)/128,128>>>(da.p,dnf.p,nn);
   if(cudaGetLastError()!=cudaSuccess||cudaDeviceSynchronize()!=cudaSuccess)return 3;
   if(!dnf.get(nf.data()))return 4;
  }
  for(int k=0;k<nn;k++)if(nf[19*k+18]!=0.)return 100+int(nf[19*k+18]);
  for(double v:nf)if(!isfinite(v))return 5;
  std::copy(nf.begin(),nf.end(),result);return 0;
 }catch(...){return 6;}
}

extern "C" int lagmoment_ordered_ledgers(int backend,int threads,int n,int nn,
 const double *nf,const int64_t *owners,const int64_t *fr,const double *stage,
 const double *hydro,double *result){
 if((backend!=0&&backend!=1)||threads<1||n<1||nn<1||!nf||!owners||!fr||!stage||!hydro||!result)return 1;
 double profile_start=omp_get_wtime();
 try{
  std::vector<double> out(12,0.),wall_terms,frozen_terms,active_terms;
  if(backend==1){wall_terms.reserve(5*size_t(nn)/32);frozen_terms.reserve(5*size_t(n));active_terms.reserve(5*size_t(n));}
  for(int k=0;k<nn;k++){
   const int64_t left=owners[2*k],right=owners[2*k+1];const double *f=nf+19*k;
   if(left<0||left>=n||right< -1||right>=n||left==right||f[18]!=0.)return 1;
   for(int q=0;q<19;q++)if(!isfinite(f[q]))return 5;
   if(backend==1&&right<0)for(int q:{7,8,17,10,11})wall_terms.push_back(f[q]);
  }
  for(int i=0;i<n;i++){
   const double *s=stage+SWIDTH*i,*h=hydro+WIDTH*i;
   if((fr[i]!=0&&fr[i]!=1)||s[93]!=0.||h[105]!=0.)return 1;
   for(int q=0;q<SWIDTH;q++)if(!isfinite(s[q]))return 5;
   for(int q=0;q<WIDTH;q++)if(!isfinite(h[q]))return 5;
   if(backend==1){
    if(fr[i])for(int q:{7,8,17,10,11})frozen_terms.push_back(h[q]);
    else {const double *g=s+30;for(int q:{7,8,17,10,11})active_terms.push_back(g[q]);}
   }
  }
  if(backend==0){
   ledger(n,nn,nf,owners,fr,hydro,out.data());gravity_ledger(n,fr,stage,out.data()+8);
  }else{
   Device<double> dw,df,da,dr;
   const int nw=int(wall_terms.size()/5),nfr=int(frozen_terms.size()/5),na=int(active_terms.size()/5);
   if(!dw.set(wall_terms.size(),wall_terms.data())||!df.set(frozen_terms.size(),frozen_terms.data())
    ||!da.set(active_terms.size(),active_terms.data())||!dr.set(12))return 2;
   bridge_compact_ordered_ledger_kernel<<<1,1>>>(nw,nfr,na,dw.p,df.p,da.p,dr.p);
   if(cudaGetLastError()!=cudaSuccess||cudaDeviceSynchronize()!=cudaSuccess)return 3;
   if(!dr.get(out.data()))return 4;
  }
  if(std::getenv("LAGFORCE_PROFILE_LEDGER")){
   std::fprintf(stderr,"[LAGMOMENT_ORDERED_LEDGER] backend=%s nodes=%d wall_nodes=%zu cells=%d frozen=%zu active=%zu h2d_bytes=%zu seconds=%.9g\n",
    backend==1?"cuda":"cpu",nn,wall_terms.size()/5,n,frozen_terms.size()/5,active_terms.size()/5,
    (wall_terms.size()+frozen_terms.size()+active_terms.size())*sizeof(double),omp_get_wtime()-profile_start);
   std::fflush(stderr);
  }
  for(double v:out)if(!isfinite(v))return 5;
  std::copy(out.begin(),out.end(),result);return 0;
 }catch(...){return 6;}
}

extern "C" int lagmoment_owned_stage(int backend,int threads,uint64_t epoch,uint64_t expected_epoch,
 int nr,int n,int nn,int64_t ni,int64_t ng,const int64_t *owned,const int64_t *node_ids,
 const double *a,const double *nf,const double *raw,const int64_t *owners,const int64_t *off,
 const int64_t *inc,const int64_t *fr,const double *extra,const double *L,
 const double *gn,const int64_t *goff,double *stage_out,double *q_out,double *hydro_out){
 if((backend!=0&&backend!=1)||threads<1||epoch!=expected_epoch||nr<1||n<1||n>nr||nn<1||ni<1||ng<1
  ||!owned||!node_ids||!a||!nf||!raw||!owners||!off||!inc||!fr||!extra||!L||!gn||!goff
  ||!stage_out||!q_out||!hydro_out)return 1;
 if(off[0]!=0||off[n]!=ni||goff[0]!=0||goff[n]!=ng)return 1;
 try{
  std::vector<int> local(nr,-1);
  for(int i=0;i<n;i++){
   if(owned[i]<0||owned[i]>=nr||(i&&owned[i]<=owned[i-1]))return 1;
   if(off[i]>off[i+1]||goff[i]>=goff[i+1]||(fr[i]!=0&&fr[i]!=1))return 1;
   local[owned[i]]=i;
  }
  int64_t expected=0;
  for(int k=0;k<nn;k++){
   if(node_ids[k]<0||(k&&node_ids[k]<=node_ids[k-1]))return 1;
   if(owners[2*k]<0||owners[2*k]>=nr||owners[2*k+1]<-1||owners[2*k+1]>=nr
    ||owners[2*k]==owners[2*k+1])return 1;
   expected+=local[owners[2*k]]>=0;
   expected+=owners[2*k+1]>=0&&local[owners[2*k+1]]>=0;
   if(nf[19*k+18]!=0.)return 100+int(nf[19*k+18]);
  }
  if(expected!=ni)return 1;
  for(int i=0;i<n;i++){
   int64_t previous=-1;
   for(int64_t t=off[i];t<off[i+1];t++){
    int64_t code=inc[t];if(code==0||code>nn||code<-int64_t(nn))return 1;
    int64_t k=code>0?code-1:-code-1;
    if(node_ids[k]<=previous||owners[2*k+(code>0?1:0)]!=owned[i])return 1;
    previous=node_ids[k];
   }
  }
  for(int64_t k=0;k<19*int64_t(nn);k++)if(!isfinite(nf[k]))return 5;
  std::vector<double> out(SWIDTH*size_t(n)),q(16*size_t(ni)),hydro(WIDTH*size_t(n));
  if(backend==0){
   omp_set_num_threads(threads);
   #pragma omp parallel for
   for(int i=0;i<n;i++)stage_cell(i,a,nf,raw,off,inc,fr,extra,L,gn,goff,out.data()+SWIDTH*i,q.data(),hydro.data()+WIDTH*i);
  }else{
   Device<double> da,dnf,dr,de,dl,dgn,dout,dq,dh;
   Device<int64_t> doff,dinc,dfr,dgoff;
   if(!da.set(33*size_t(nn),a)||!dnf.set(19*size_t(nn),nf)||!dr.set(19*size_t(n),raw)
    ||!de.set(20*size_t(nn),extra)||!dl.set(4*size_t(n),L)||!dgn.set(7*size_t(ng),gn)
    ||!dout.set(out.size())||!dq.set(q.size())||!dh.set(hydro.size())
    ||!doff.set(n+1,off)||!dinc.set(ni,inc)||!dfr.set(n,fr)||!dgoff.set(n+1,goff))return 2;
   bridge_stage_kernel<<<(n+63)/64,64>>>(n,da.p,dnf.p,dr.p,doff.p,dinc.p,dfr.p,de.p,dl.p,dgn.p,dgoff.p,dout.p,dq.p,dh.p);
   if(cudaGetLastError()!=cudaSuccess||cudaDeviceSynchronize()!=cudaSuccess)return 3;
   if(!dout.get(out.data())||!dq.get(q.data())||!dh.get(hydro.data()))return 4;
  }
  for(int i=0;i<n;i++)if(out[SWIDTH*i+93]!=0.)return 100+int(out[SWIDTH*i+93]);
  for(double v:out)if(!isfinite(v))return 5;
  for(double v:q)if(!isfinite(v))return 5;
  for(double v:hydro)if(!isfinite(v))return 5;
  std::copy(out.begin(),out.end(),stage_out);std::copy(q.begin(),q.end(),q_out);
  std::copy(hydro.begin(),hydro.end(),hydro_out);return 0;
 }catch(...){return 6;}
}
