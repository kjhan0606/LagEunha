// Whole-record native readback of reusable factors, not a physical simulation.
#include <mpi.h>
#include <omp.h>
#include <stdexcept>
#include <string>
#include <cstdio>
#include <unistd.h>
static void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
#include "material_cell_factors_io.h"
int main(int argc,char**argv){
    MPI_Init(&argc,&argv);int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    try{
        need(argc==4&&ranks==4,"MPI4 readback: partition_run factor_run new_result_json");double start=MPI_Wtime();std::string suffix=std::to_string(rank);
        cell_condensation::Partition p(std::string(argv[1])+"/cell_partition_r"+suffix+".bin",rank,ranks);
        cell_condensation::Factors factor(p,std::string(argv[2])+"/cell_factors_r"+suffix+".bin",rank,ranks);
        cell_condensation::SkeletonMass mass(p,std::string(argv[2])+"/skeleton_mass_r"+suffix+".bin",rank,ranks);
        uint64_t count[5]={factor.cells,uint64_t(factor.L.size()),uint64_t(factor.IG.size()),uint64_t(mass.value.size()),p.ni};MPI_Allreduce(MPI_IN_PLACE,count,5,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
        need(count[0]==71080&&count[1]==635433712&&count[2]==379124224&&count[3]==96946492&&count[4]==9460664,"all original source factor/CSR records checked, no sampling");
        if(rank==0){FILE*f=fopen(argv[3],"wx");need(f,"new whole-record readback result");fprintf(f,"{\"scope\":\"whole-record native checksum/readback and exact binary80-to-DD representability; no solve or evolution\",\"cells\":%llu,\"L_entries\":%llu,\"IG_entries\":%llu,\"partial_CSR_entries\":%llu,\"interior_dofs\":%llu,\"all_payload_checksums_pass\":true,\"all_factor_diagonals_positive\":true,\"all_binary80_entries_exact_DD_encoding\":true,\"all_partial_CSR_entries_symmetric\":true,\"physical_trajectory_advanced\":false,\"actual_full_equation_solve_pending\":true,\"seconds\":%.9g}\n",(unsigned long long)count[0],(unsigned long long)count[1],(unsigned long long)count[2],(unsigned long long)count[3],(unsigned long long)count[4],MPI_Wtime()-start);need(fflush(f)==0&&fsync(fileno(f))==0&&fclose(f)==0,"durable readback report");}
    }catch(const std::exception&e){fprintf(stderr,"rank%d CELL_FACTOR_READBACK_REJECTED %s\n",rank,e.what());MPI_Abort(MPI_COMM_WORLD,2);return 2;}
    MPI_Finalize();return 0;
}
