// Finite synthetic bounded-plane geometry checks; not curved Ω* or a remap.
#include <mpi.h>
#include <omp.h>
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
#include "material_compact_power_cell.h"
namespace cp=compact_power;
int main(int argc,char**argv){MPI_Init(&argc,&argv);try{{
    int ranks=0;MPI_Comm_size(MPI_COMM_WORLD,&ranks);need(ranks==1,"one rank finite convex power geometry fixture");
    compact_target::Registry registry;registry.exterior.resize(1);
    for(int k=0;k<2;k++){registry.exterior[0].lower[k]=-1.25;registry.exterior[0].upper[k]=1.25;}
    for(int j=0;j<25;j++){
        compact_source::Site s;s.cell=10+j;s.position[0]=-.75+.375*(j%5)+double((j*37)%29-14)/8192;s.position[1]=-.75+.375*(j/5)+double((j*17)%31-15)/8192;s.canonical_bound=true;
        for(int k=0;k<2;k++){s.position_low[k]=ldexp(double(j%3-1),-60);s.position_error[k]=1e-12;s.canonical_position_error[k]=1e-31;}registry.site.push_back(s);
    }
    cp::CellBudget budget;budget.maximum_nodes=4096;budget.maximum_planes=512;budget.maximum_corners=128;
    size_t planes=0,nodes=0,excluded=0,shared_vertices=0;R maximum_coverage_width=0;int used_threads=0;double construction_seconds=0;
    for(int pass=0;pass<2;pass++){
        std::vector<R>weight(25);if(pass)for(size_t j=0;j<weight.size();j++)weight[j]=R(int((j*7)%13)-6)/4096;
        cp::Spatial tree(registry,weight);cp::S area;
        auto external_weights=weight;cp::Spatial sealed(registry,external_weights);external_weights[0]+=1;
        need(sealed.weight[0]==weight[0]&&sealed.weight_hash==tree.weight_hash,"each geometry trial owns an immutable weight copy");
        struct Seen {std::array<cp::S,2>x;size_t count=0;};std::map<cp::CornerKey,Seen>corner;
        std::vector<cp::BoundedCell>cells(registry.site.size());bool valid=true;double start=MPI_Wtime();
        #pragma omp parallel reduction(&&:valid) reduction(max:used_threads)
        {
            used_threads=std::max(used_threads,omp_get_num_threads());
            #pragma omp for schedule(static)
            for(size_t i=0;i<registry.site.size();i++)try{cells[i]=cp::bounded_cell(tree,i,budget);}catch(const std::exception&){valid=false;}
        }
        construction_seconds+=MPI_Wtime()-start;need(valid,"all parallel bounded power cells constructed");
        for(size_t i=0;i<registry.site.size();i++){
            const auto&cell=cells[i];need(!cell.corner.empty(),"all small synthetic power cells positive");planes+=cell.tested_sites.size();nodes+=cell.nodes;excluded+=cell.excluded.size();
            for(const auto&c:cell.corner){
                // Independent exhaustive half-space verification includes
                // every site omitted by the spatial node exclusion proof.
                for(size_t j=0;j<registry.site.size();j++)if(j!=i){
                    cp::Constraint constraint;constraint.plane=cp::plane(registry.site[i],weight[i],registry.site[j],weight[j]);constraint.key={1,constraint.plane.cell[0],constraint.plane.cell[1]};constraint.side=cell.cell==constraint.plane.cell[0]?1:-1;
                    need(cp::side(constraint,c)<=0,"exhaustive all-site power inequality, including every pruned node");
                }
                auto at=corner.emplace(c.key,Seen{c.x,0});auto&prior=at.first->second;
                for(int k=0;k<2;k++)need(prior.x[k].n==c.x[k].n&&prior.x[k].e==c.x[k].e,"shared corner uses bit-identical canonical expression across cells");prior.count++;
            }
            const auto&origin=cell.corner[0].x;
            for(size_t j=1;j+1<cell.corner.size();j++){
                const auto&a=cell.corner[j].x;const auto&b=cell.corner[j+1].x;
                area=area+((a[0]-origin[0])*(b[1]-origin[1])-(a[1]-origin[1])*(b[0]-origin[0]))/cp::S(2);
            }
        }
        for(const auto&entry:corner){const auto&key=entry.first;size_t expected=key.kind==0?1:(key.kind==1?2:3);need(entry.second.count==expected,"bounding corners/edge events/interior triples have complete cell incidence");shared_vertices+=expected>1;}
        need(cp::sf::lower(area)<=cp::Q(6.25L)&&cp::sf::upper(area)>=cp::Q(6.25L),"sum of all bounded power areas encloses the full source-enclosing rectangle");
        maximum_coverage_width=std::max(maximum_coverage_width,R(cp::sf::upper(area)-cp::sf::lower(area)));
        need(cp::sf::upper(area)-cp::sf::lower(area)<cp::Q(1e-24L),"bounded power coverage arithmetic width");
    }
    need(excluded>0&&planes<2*25*24,"spatial exclusion exercised, fewer planes than exhaustive construction");
    unsigned rejected=0;
    auto bad=registry;bad.site.resize(4);double square[4][2]={{-.5,-.5},{.5,-.5},{.5,.5},{-.5,.5}};
    for(int j=0;j<4;j++){bad.site[j].canonical_bound=false;for(int k=0;k<2;k++){bad.site[j].position[k]=square[j][k];bad.site[j].position_error[k]=bad.site[j].position_low[k]=bad.site[j].canonical_position_error[k]=0;}}
    std::vector<R>weight(4);cp::Spatial degenerate(bad,weight);bool caught=false;
    try{auto unused=cp::bounded_cell(degenerate,0,budget);(void)unused;}catch(const std::exception&){caught=true;}need(caught,"unresolved four-site concurrence not broken by tolerance or jitter");rejected++;
    weight[1]=100;cp::Spatial empty(bad,weight);auto missing=cp::bounded_cell(empty,0,budget);need(missing.corner.empty(),"a weight-swallowed cell is reported empty, never filled/repaired");rejected++;
    std::vector<R>regular_weights(registry.site.size());cp::Spatial regular(registry,regular_weights);auto short_budget=budget;short_budget.maximum_planes=1;caught=false;
    try{auto unused=cp::bounded_cell(regular,12,short_budget);(void)unused;}catch(const std::exception&){caught=true;}need(caught,"finite per-cell work budget enforced");rejected++;
    std::cout<<std::setprecision(17)<<"{\"status\":\"BOUNDED_PLANE_POWER_CELLS_PASS_SYNTHETIC_ONLY\",\"passes\":2,\"cells_per_pass\":25,\"openmp_threads\":"<<used_threads<<",\"construction_seconds\":"<<construction_seconds<<",\"tested_planes\":"<<planes<<",\"exhaustive_plane_count\":1200,\"visited_nodes\":"<<nodes<<",\"excluded_nodes\":"<<excluded<<",\"shared_canonical_vertices\":"<<shared_vertices<<",\"maximum_area_enclosure_width\":"<<double(maximum_coverage_width)<<",\"negative_cases\":"<<rejected<<",\"actual_curved_domain_clipped\":false,\"actual_source_imported\":false,\"actual_partition_proven\":false,\"new_physical_step\":false}\n";
}MPI_Finalize();return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);return 1;}}
