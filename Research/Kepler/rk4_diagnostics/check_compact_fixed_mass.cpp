// Synthetic fixed-M reader and finite-mass law arithmetic regression.
// This does not read accepted simulation inputs or run a physical step.
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include <unistd.h>
#include "material_compact_checkpoint.h"
using R=long double;
static void need(bool b,const char*s){compact_p3::require(b,s);}
#include "material_full_owned.h"
#include "material_distributed_mass.h"
#include "material_compact_mass.h"
namespace st=shared_trace;
#include "material_compact_source.h"
#include "material_compact_source_geometry.h"
#include "material_compact_source_context.h"
#include "material_compact_source_target.h"
#include "material_compact_target_registry.h"
#include "material_compact_mass_law.h"

static void write_payload(const std::string&path,const std::array<R,3>&mass){
    FILE*f=fopen(path.c_str(),"wb");need(f,"new isolated synthetic fixed-M fixture");
    uint64_t h[4]={0x4d464d4153533031ULL,1,uint64_t(mass.size())*ST_MAX*ST_MAX,sizeof(R)};
    need(fwrite(h,sizeof(uint64_t),4,f)==4,"synthetic fixed-M header");
    std::array<R,size_t(ST_MAX)*ST_MAX>block{};
    for(R total:mass){block.fill(0);for(int a=0;a<10;a++)block[a*ST_MAX+a]=total/10;
        need(fwrite(block.data(),sizeof(R),block.size(),f)==block.size(),"synthetic fixed-M triangle block");}
    need(fclose(f)==0,"close synthetic fixed-M fixture");
}
static bool rejects(const std::string&path,const std::vector<uint64_t>&triangle_cell,size_t cells){
    try{(void)compact_mass_law::read_fixed_mass(path,triangle_cell,cells);return false;}
    catch(const std::exception&){return true;}
}
int main(int argc,char**argv){
    MPI_Init(&argc,&argv);int rank=0,ranks=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    try{
        need(ranks==1,"one-rank finite fixed-M payload fixture");
        std::string path="/tmp/compact_fixed_mass_reader_"+std::to_string(getpid())+".bin";
        std::array<R,3>triangle_mass={2,3,4};write_payload(path,triangle_mass);
        auto payload=compact_mass_law::read_fixed_mass(path,{0,0,1},2);
        need(payload.triangles==3&&payload.cell_mass[0]==5&&payload.cell_mass[1]==4&&payload.total_mass==9,"first-ten P3 fixed-M contraction summed to original source cells");
        auto pair=compact_mass_law::lossless_pair(R(1)/R(7));
        double pair_wire[2]={pair.hi,pair.lo};
        need(compact_mass_law::from_pair(pair_wire)==R(1)/R(7),"binary80 mass survives fieldwise two-double transport");
        auto law_entry=compact_mass_law::make_entry(7,0,5,2,4.8L);
        need(law_entry.mu==2.5L&&fabsl(law_entry.point_to_fixed_relative-0.04L)<1e-18L,"mu=mref/Vsource with coefficient one, while point/M discrepancy remains explicit");

        std::string bad_count=path+".bad_count";write_payload(bad_count,triangle_mass);
        FILE*f=fopen(bad_count.c_str(),"r+b");need(f,"open isolated malformed fixture");uint64_t wrong=uint64_t(2)*ST_MAX*ST_MAX;
        need(fseeko(f,2*sizeof(uint64_t),SEEK_SET)==0&&fwrite(&wrong,sizeof(wrong),1,f)==1&&fclose(f)==0,"write malformed fixed-M count");
        bool rejected_count=rejects(bad_count,{0,0,1},2);
        std::string truncated=path+".truncated";write_payload(truncated,triangle_mass);
        need(truncate(truncated.c_str(),off_t(4*sizeof(uint64_t)+2*ST_MAX*ST_MAX*sizeof(R)))==0,"truncate isolated fixed-M fixture");
        bool rejected_truncated=rejects(truncated,{0,0,1},2);
        unlink(path.c_str());unlink(bad_count.c_str());unlink(truncated.c_str());
        need(rejected_count&&rejected_truncated,"wrong block count and incomplete payload rejected");
        if(rank==0)std::cout<<std::setprecision(18)
            <<"{\"scope\":\"synthetic original fixed-M binary80 reader; no accepted source or physical step\","
            <<"\"triangles\":3,\"cells\":2,\"cell_mref\":["<<(double)payload.cell_mass[0]<<","<<(double)payload.cell_mass[1]<<"],"
            <<"\"total_mref\":"<<(double)payload.total_mass<<",\"mu_example\":"<<(double)law_entry.mu<<",\"point_mass_discrepancy_example\":"<<(double)law_entry.point_to_fixed_relative<<",\"lossless_pair_pass\":true,"
            <<"\"wrong_count_rejected\":"<<(rejected_count?"true":"false")<<",\"truncated_payload_rejected\":"<<(rejected_truncated?"true":"false")<<","
            <<"\"actual_fixed_mass_read\":false,\"weight_law_solved\":false,\"new_physical_step\":false}\n";
        MPI_Finalize();return 0;
    }catch(const std::exception&e){std::cerr<<"FIXED_MASS_FIXTURE_REJECTED "<<e.what()<<"\n";MPI_Abort(MPI_COMM_WORLD,2);return 2;}
}
