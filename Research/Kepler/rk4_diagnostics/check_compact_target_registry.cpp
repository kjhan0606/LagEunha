// Finite MPI serialization/connectivity checks. No accepted state, target
// power cells, geometry certification, transfer or physical step is executed.
#include <mpi.h>
#include <iostream>
#include <iomanip>
#include "material_compact_checkpoint.h"
using R=long double;
static void need(bool b,const char*s){compact_p3::require(b,s);}
#include "material_full_owned.h"
#include "material_compact_mass.h"
namespace st=shared_trace;
#include "material_compact_source.h"
#include "material_compact_source_geometry.h"
#include "material_compact_source_context.h"
#include "material_compact_source_target.h"
#include "material_compact_target_registry.h"
struct Record {uint64_t cell;struct {GfsMaterialTriangle geometry;}triangle;};
struct Packet {double parameters[4]={0,0,0,.047};std::vector<Record>record;};
int main(int argc,char**argv){MPI_Init(&argc,&argv);try{{
    int rank=0,ranks=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    need(ranks==1||ranks==2||ranks==4,"finite MPI1/2/4 fixture");
    Packet packet;double outer[4][2]={{4,4},{5,4},{5,5},{4,5}};
    for(int j=0;j<4;j++){
        Record r={};r.cell=10*(j+1);for(int k=0;k<2;k++){r.triangle.geometry.x[0][k]=outer[j][k];r.triangle.geometry.x[1][k]=outer[(j+1)%4][k];r.triangle.geometry.x[2][k]=4.5;}packet.record.push_back(r);
    }
    auto f=st::topology(packet,ranks);
    for(size_t j=0;j<f.polynomial.size();j++){
        double x=f.polynomial[j][0],y=f.polynomial[j][1];
        f.coefficient[j].x[0]=1.02*x+.03*y+.2;f.coefficient[j].x[1]=-.01*x+.98*y-.1;
        f.coefficient[j].v[0]=-.2*y+.01;f.coefficient[j].v[1]=.2*x-.03;
    }
    for(size_t j=f.polynomial.size();j<f.coefficient.size();j++)for(int k=0;k<2;k++){
        f.coefficient[j].x[k]=ldexp(double((j+1)*(k+1)),-32);
        f.coefficient[j].v[k]=-ldexp(double((j+2)*(k+1)),-28);
    }
    compact_source::TargetSource local;local.source_steps=5;local.source_epoch=13;local.source_time=.1;local.global_sites=4;local.global_exterior=4;
    for(int j=0;j<4;j++)local.source_shard_hash[j]=111+j;
    for(int j=0;j<4;j++){
        int owner=j%std::max(1,ranks-1);if(owner!=rank)continue; // last MPI rank intentionally empty
        compact_source::Site s;s.cell=10*(j+1);s.source_triangle=50+j;s.owner=owner;s.chart=j%3;s.source_mass=R(j+1)+ldexpl(1.L,-62);
        for(int k=0;k<2;k++){
            s.reference[k]=(outer[j][k]+outer[(j+1)%4][k]+4.5)/3;s.position[k]=s.reference[k]+.1;s.velocity[k]=.125*(j+k);s.position_error[k]=ldexp(1.,-45);s.velocity_error[k]=ldexp(1.,-47);
            if(j%2==0){s.canonical_bound=true;s.position_low[k]=ldexp(double(j+1),-56);s.velocity_low[k]=-ldexp(double(j+1),-57);s.canonical_position_error[k]=s.canonical_velocity_error[k]=ldexp(1.,-105);}
        }local.site.push_back(s);
        const auto&t=f.triangle[j];compact_source::Exterior e;e.cell=s.cell;e.source_triangle=s.source_triangle;e.owner=owner;e.edge=900+j;e.local_edge=0;e.orientation=t.endpoint[0][0]==0?1:-1;
        const auto&q=t.edge[0];e.vertex[0]=100+q.poly[0];e.vertex[1]=100+q.poly[3];
        e.frozen=compact_geometry::freeze_trace(q,f.coefficient,e.edge,local.source_shard_hash[rank],local.source_epoch);
        for(int a=0;a<2;a++)for(int k=0;k<2;k++)e.reference[a][k]=t.vertex[t.endpoint[0][a]][k];
        for(int k=0;k<2;k++){e.lower[k]=3;e.upper[k]=7;}local.exterior.push_back(e);
    }
    auto registry=compact_target::gather(local,MPI_COMM_WORLD);
    need(registry.site.size()==4&&registry.exterior.size()==4&&registry.boundary_cycle.size()==4,"complete small geometry registry");
    need(registry.source_steps==5&&registry.source_epoch==13&&registry.source_time==.1&&registry.source_shard_hash==local.source_shard_hash,"exact registry source provenance");
    for(int j=0;j<4;j++)need(registry.site[j].source_mass==R(j+1)+ldexpl(1.L,-62),"binary80 mass low bits retained");
    for(int j=0;j<4;j++){
        const auto&s=registry.site[j];need(s.canonical_bound==(j%2==0),"explicit canonical precision provenance retained");
        if(s.canonical_bound)for(int k=0;k<2;k++)need(s.position_low[k]==ldexp(double(j+1),-56)&&s.velocity_low[k]==-ldexp(double(j+1),-57)&&s.canonical_position_error[k]==ldexp(1.,-105)&&s.canonical_velocity_error[k]==ldexp(1.,-105),"canonical high/low/tail wire exact");
    }
    uint64_t minimum=registry.wire_hash,maximum=minimum;MPI_Allreduce(MPI_IN_PLACE,&minimum,1,MPI_UINT64_T,MPI_MIN,MPI_COMM_WORLD);MPI_Allreduce(MPI_IN_PLACE,&maximum,1,MPI_UINT64_T,MPI_MAX,MPI_COMM_WORLD);need(minimum==maximum,"same global fieldwise wire digest");
    R curve_difference=0;
    for(const auto&e:registry.exterior){
        size_t j=e.edge-900;auto original=compact_geometry::freeze_trace(f.triangle[j].edge[0],f.coefficient,e.edge,local.source_shard_hash[e.owner],13);
        auto before=compact_target::pack_edges({e},e.owner,local.source_shard_hash[e.owner],13);
        auto round=compact_target::unpack_edge(before.tag.data(),before.value.data());
        auto after=compact_target::pack_edges({round},e.owner,local.source_shard_hash[e.owner],13);
        need(before.tag==after.tag&&before.value.size()==after.value.size(),"full fieldwise curve tag round trip");
        for(size_t k=0;k<before.value.size();k++)need(compact_target::bits(before.value[k])==compact_target::bits(after.value[k]),"every frozen curve double round trips bitwise");
        for(double u:{0.,.125,.5,.875,1.}){
            STJet a[2],b[2];need(compact_geometry::trace(original,STJet(u,1,0),a)&&compact_geometry::trace(e.frozen,STJet(u,1,0),b),"same frozen canonical curve function after MPI");
            for(int k=0;k<2;k++)curve_difference=std::max({curve_difference,fabsl(R(a[k].v)-b[k].v),fabsl(R(a[k].x)-b[k].x)});
        }
    }
    need(curve_difference==0,"frozen exterior values/derivatives unchanged");
    unsigned rejected=0;
    for(int test=0;test<8;test++){
        auto bad=local;
        if(rank==0){
            if(test==0)bad.source_epoch++;
            if(test==1)bad.site[0].owner=1;
            if(test==2)bad.site[0].cell=40;
            if(test==3)bad.exterior.pop_back();
            if(test==4)bad.exterior[0].frozen.coefficient[0].x[0]=nextafter(bad.exterior[0].frozen.coefficient[0].x[0],INFINITY);
            if(test==5)bad.exterior[0].orientation*=-1;
            if(test==6)bad.site[0].position_error[0]=-1;
            if(test==7)bad.exterior[0].frozen.function.poly[1]=0;
        }
        int caught=0;try{auto unused=compact_target::gather(bad,MPI_COMM_WORLD);(void)unused;}catch(const std::exception&){caught=1;}
        int all=0;MPI_Allreduce(&caught,&all,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);need(all==1,"malformed registry rejected collectively, including empty rank");rejected++;
    }
    // A pair of closed triangles has correct vertex degrees but is NOT the
    // original square's one boundary component. Exercise the cycle check.
    auto disconnected=registry;disconnected.exterior.clear();disconnected.boundary_cycle.clear();
    for(int component=0;component<2;component++)for(int j=0;j<3;j++){
        auto e=registry.exterior[0];e.vertex[0]=1000+3*component+j;e.vertex[1]=1000+3*component+(j+1)%3;e.orientation=1;
        for(int a=0;a<2;a++)for(int k=0;k<2;k++){
            double x=4+component+double((j+a)%3)/8;e.reference[a][k]=x;e.frozen.coefficient[a?3:0].x[k]=x;e.frozen.coefficient[a?3:0].v[k]=0;
        }disconnected.exterior.push_back(e);
    }
    int caught=0;try{compact_target::verify_boundary(disconnected);}catch(const std::exception&){caught=1;}need(caught==1,"disconnected extra closed boundary cycle rejected");rejected++;
    if(rank==0)std::cout<<std::setprecision(17)<<"{\"status\":\"SMALL_SOURCE_GEOMETRY_REGISTRY_MPI_PASS_SYNTHETIC_ONLY\",\"ranks\":"<<ranks<<",\"sites\":4,\"curves\":4,\"boundary_components\":1,\"empty_rank_exercised\":"<<(ranks>1?"true":"false")<<",\"binary80_mass_exact\":true,\"curve_round_trip_difference\":"<<double(curve_difference)<<",\"negative_cases_rejected\":"<<rejected<<",\"actual_source_imported\":false,\"actual_partition_proven\":false,\"new_physical_step\":false}\n";
}MPI_Finalize();return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);return 1;}}
