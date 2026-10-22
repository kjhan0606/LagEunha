// Canonical finite geometry fixture only: no actual source guard/partition.
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
#include "material_compact_power_curve.h"
namespace cp=compact_power;
template<class F> static compact_source::Exterior exterior(F x){
    compact_source::Exterior e;e.edge=77;e.frozen.edge=77;e.frozen.source_hash=11;e.frozen.source_epoch=5;e.frozen.function.rank=0;e.frozen.function.eps=.047;
    for(int j=0;j<4;j++){
        e.frozen.function.poly[j]=j;e.frozen.coefficient[j].x[0]=x(j);
        e.frozen.coefficient[j].x[1]=4+double(j*j)/16;
    }
    e.reference[0][0]=3;e.reference[0][1]=4;e.reference[1][0]=6;e.reference[1][1]=4;
    return e;
}
static cp::Q polynomial(const compact_source::Exterior&e,cp::Q t){
    // Independent Q128 direct Lagrange value, not production decomposition.
    cp::Q sum=0;for(int j=0;j<4;j++){
        cp::Q N=1;for(int k=0;k<4;k++)if(j!=k)N*=((3*t)-k)/cp::Q(j-k);
        sum+=cp::Q(e.frozen.coefficient[j].x[0])*N;
    }return sum;
}
int main(int argc,char**argv){MPI_Init(&argc,&argv);try{{
    int ranks=0;MPI_Comm_size(MPI_COMM_WORLD,&ranks);need(ranks==1,"one rank finite canonical geometry fixture");
    compact_source::Site a,b;a.cell=10;b.cell=20;a.position[0]=4;b.position[0]=5;a.position[1]=b.position[1]=4;
    cp::RootBudget budget;budget.maximum_boxes=2048;budget.maximum_depth=32;budget.maximum_newton=32;budget.parameter_width=ldexpq(1,-44);
    unsigned cases=0;size_t boxes=0;R maximum_width=0;
    auto check=[&](const cp::Roots&r,size_t n){
        need(r.complete&&r.unresolved==cp::Unresolved::none&&r.event.size()==n,"complete canonical isolated event count");cases++;boxes+=r.boxes;
        for(size_t j=0;j<r.event.size();j++){const auto&e=r.event[j];need(e.cell[0]==10&&e.cell[1]==20&&e.edge==77&&e.ordinal==j,"canonical plane/curve/event identity");need(e.upper-e.lower<=budget.parameter_width,"explicit caller geometry width");maximum_width=std::max(maximum_width,R(e.upper-e.lower));}
    };
    auto line=exterior([](int j){return double(3+j);});auto p=cp::plane(a,0,b,0);auto roots=cp::intersections(p,cp::curve(line),budget);check(roots,1);
    need(roots.event[0].lower<=cp::Q(.5L)&&roots.event[0].upper>=cp::Q(.5L)&&roots.event[0].derivative_sign==1,"bracketed exact affine crossing retained");
    auto reversed=cp::intersections(cp::plane(b,0,a,0),cp::curve(line),budget);check(reversed,1);
    need(roots.event[0].lower==reversed.event[0].lower&&roots.event[0].upper==reversed.event[0].upper,"pair reversal preserves canonical event bits");
    auto weighted=cp::intersections(cp::plane(a,.125L,b,.25L),cp::curve(line),budget);check(weighted,1);
    auto shifted=cp::intersections(cp::plane(a,.625L,b,.75L),cp::curve(line),budget);check(shifted,1);
    need(weighted.event[0].lower==shifted.event[0].lower&&weighted.event[0].upper==shifted.event[0].upper,"common exactly representable weight shift leaves power plane unchanged");
    cp::Q weighted_root=(cp::Q(4.4375L)-3)/3;need(weighted.event[0].lower<=weighted_root&&weighted.event[0].upper>=weighted_root,"weight sign agrees with |x-q|²-W inequality");
    auto outside=exterior([](int j){return double(6+j);});check(cp::intersections(p,cp::curve(outside),budget),0);
    auto triple=exterior([](int j){R t=R(j)/3;return double(4.5L+(t-.21L)*(t-.47L)*(t-.81L));});
    auto three=cp::intersections(p,cp::curve(triple),budget);check(three,3);
    cp::Q expected[3]={cp::Q(.21L),cp::Q(.47L),cp::Q(.81L)};
    for(int j=0;j<3;j++){
        cp::Q t=expected[j];for(int it=0;it<12;it++){
            cp::Q h=ldexpq(1,-40),derivative=(polynomial(triple,t+h)-polynomial(triple,t-h))/(2*h);
            t-=(polynomial(triple,t)-cp::Q(4.5L))/derivative;
        }
        need(t>=three.event[j].lower&&t<=three.event[j].upper,"independent binary128 cubic root lies inside canonical event enclosure");
    }
    auto unresolved_budget=budget;unresolved_budget.maximum_boxes=64;unresolved_budget.maximum_depth=8;
    auto coincident=exterior([](int){return 4.5;});auto unresolved=cp::intersections(p,cp::curve(coincident),unresolved_budget);
    need(!unresolved.complete&&unresolved.unresolved!=cp::Unresolved::none,"coincident curve is not silently flattened/tied or declared empty");cases++;
    auto endpoint=exterior([](int j){return 4.5+double(j);});unresolved=cp::intersections(p,cp::curve(endpoint),budget);
    need(!unresolved.complete&&unresolved.unresolved==cp::Unresolved::endpoint_or_coincidence,"endpoint event requires explicit algebraic adapter, never snapped");cases++;
    auto uncertain_a=a;uncertain_a.position_error[0]=.001;
    unresolved=cp::intersections(cp::plane(uncertain_a,0,b,0),cp::curve(line),budget);
    need(!unresolved.complete&&unresolved.unresolved==cp::Unresolved::root_width,"source q uncertainty cannot be discarded to certify narrow root");cases++;
    uncertain_a.canonical_bound=true;uncertain_a.canonical_position_error[0]=1e-30;
    auto refined=cp::intersections(cp::plane(uncertain_a,0,b,0),cp::curve(line),budget);check(refined,1);
    need(refined.event[0].lower<=cp::Q(.5L)&&refined.event[0].upper>=cp::Q(.5L),"explicit independently refined SAME q enclosure resolves cache precision floor");
    std::cout<<std::setprecision(17)<<"{\"status\":\"CANONICAL_POWER_CURVE_EVENTS_PASS_SYNTHETIC_ONLY\",\"cases\":"<<cases<<",\"certified_boxes\":"<<boxes<<",\"largest_root_parameter_width\":"<<double(maximum_width)<<",\"triple_crossing_roots\":3,\"weight_sign_and_gauge_pass\":true,\"pair_reversal_pass\":true,\"coincident_endpoint_uncertain_cases_unresolved\":true,\"actual_source_imported\":false,\"actual_partition_proven\":false,\"new_physical_step\":false}\n";
}MPI_Finalize();return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);return 1;}}
