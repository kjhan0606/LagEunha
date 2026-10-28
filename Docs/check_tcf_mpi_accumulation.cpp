#include <mpi.h>
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <vector>
#include "tcf_full_face_probe.h"
template<class T> std::vector<T> read(const char*path){
 FILE*f=fopen(path,"rb"); if(!f) MPI_Abort(MPI_COMM_WORLD,2);
 fseek(f,0,SEEK_END);long bytes=ftell(f);rewind(f);
 if(bytes%sizeof(T))MPI_Abort(MPI_COMM_WORLD,3);
 std::vector<T> a(bytes/sizeof(T));
 if(fread(a.data(),sizeof(T),a.size(),f)!=a.size())MPI_Abort(MPI_COMM_WORLD,4);
 fclose(f);return a;
}
int main(int argc,char**argv){
 int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
 int rank,size;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&size);
 if(provided<MPI_THREAD_FUNNELED)MPI_Abort(MPI_COMM_WORLD,5);
 auto input=read<double>("inputs.bin");auto gpu=read<double>("gpu.bin");
 auto owners=read<int64_t>("owners.bin");size_t n=owners.size();
 if(input.size()!=13*n||gpu.size()!=6*n)MPI_Abort(MPI_COMM_WORLD,6);
 size_t cells=*std::max_element(owners.begin(),owners.end())+1,L=cells*6;
 std::vector<double> reference(L,0),gpu_sum(L,0);
 for(size_t k=0;k<n;k++){double o[6];tcf_full_face_probe(&input[13*k],o);
  for(int c=0;c<6;c++){reference[6*owners[k]+c]+=o[c];gpu_sum[6*owners[k]+c]+=gpu[6*k+c];}}
 int requested=argc>1?atoi(argv[1]):1;omp_set_dynamic(0);omp_set_num_threads(requested);
 bool pass=true;
 for(int mode=0;mode<2;mode++){
  std::vector<double> threads(requested*L,0),local(L,0),total(L,0);
  std::vector<long long> counts(requested,0);int team=0;
#pragma omp parallel
  {
   int tid=omp_get_thread_num();double*sum=&threads[tid*L];
#pragma omp single
   team=omp_get_num_threads();
#pragma omp for schedule(static)
   for(size_t k=0;k<n;k++){
    if((mode? k%size:owners[k]%size)!=(size_t)rank)continue;
    double o[6];tcf_full_face_probe(&input[13*k],o);
    for(int c=0;c<6;c++)sum[6*owners[k]+c]+=o[c];counts[tid]++;
   }
  }
  long long local_count=0,global_count=0;
  for(int t=0;t<team;t++){local_count+=counts[t];for(size_t k=0;k<L;k++)local[k]+=threads[t*L+k];}
  MPI_Allreduce(local.data(),total.data(),L,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(&local_count,&global_count,1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
  double err=0,gerr=0,scale=0;int finite=1,minteam,maxteam;
  for(size_t k=0;k<L;k++){finite &= std::isfinite(total[k]);err=std::max(err,fabs(total[k]-reference[k]));
   gerr=std::max(gerr,fabs(total[k]-gpu_sum[k]));scale=std::max(scale,fabs(reference[k]));}
  MPI_Allreduce(&team,&minteam,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);
  MPI_Allreduce(&team,&maxteam,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
  pass &= finite && global_count==(long long)n && minteam==requested && maxteam==requested && err<1e-12*scale && gerr<1e-12*scale;
  if(rank==0)printf("{\"ranks\":%d,\"teams_min_max\":[%d,%d],\"partition\":\"%s\",\"faces\":%lld,\"finite\":%d,\"max_abs_cpu\":%.17g,\"max_abs_gpu_sum\":%.17g,\"reference_max_abs\":%.17g}\n",size,minteam,maxteam,mode?"face_index":"owner_id",global_count,finite,err,gerr,scale);
 }
 MPI_Finalize();return pass?0:1;
}
