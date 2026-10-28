#include "full_moment_provider_mpi.h"
#include <string.h>
static int collective(MPI_Comm c,int status){int out;MPI_Allreduce(&status,&out,1,MPI_INT,MPI_MAX,c);return out;}
int full_moment_collective_release(MPI_Comm c,FullMomentProvider *provider,FullMomentProviderRequest *request,const char *geometry){
 int rank;MPI_Comm_rank(c,&rank);FullMomentProviderResponse response;
 request->command=FMP_RELEASE;request->outer=1;
 int status=full_moment_provider_call(provider,request,&response);
 if(collective(c,status))return 1;
 /* No service can still be hashing/reading the old shared mesh now. */
 if(!rank){
  if(!geometry||strlen(geometry)>=sizeof request->geometry_path)status=1;
  else{strcpy(request->geometry_path,geometry);request->command=FMP_DISPOSE_GEOMETRY;status=full_moment_provider_call(provider,request,&response);}
 }
 return collective(c,status);
}
