// Synthetic transport/classification fixture only. No accepted source is
// imported and the guarded full target_source_inputs entrypoint is not run.
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
#include "material_compact_source_fields.h"
#include "material_compact_power_curve.h"
namespace sf=shared_function;
struct Record {uint64_t cell;struct {GfsMaterialTriangle geometry;}triangle;};
struct Packet {double parameters[4]={0,0,0,.047};std::vector<Record>record;};
int main(int argc,char**argv){MPI_Init(&argc,&argv);try{{
    int ranks=0;MPI_Comm_size(MPI_COMM_WORLD,&ranks);need(ranks==1,"one-rank synthetic source fixture");
    Packet packet;double vertices[2][3][2]={{{-1,0},{1,0},{0,1}},{{-1,0},{0,-1},{1,0}}};
    // Keep these deliberately coarse synthetic triangles within the original
    // radial-Taylor basis domain |L+Q|<.75. A coarse triangle crossing r=0
    // is not a valid source element of that existing native implementation.
    for(auto&triangle:vertices)for(auto&point:triangle)for(double&x:point)x+=4;
    for(int s=0;s<2;s++){Record r={};r.cell=7;for(int j=0;j<3;j++)for(int k=0;k<2;k++)r.triangle.geometry.x[j][k]=vertices[s][j][k];packet.record.push_back(r);}
    compact_source::Model m;m.owned.f=st::topology(packet,1);m.owned.source_triangle={42,12};m.owned.source_cell={7,7};
    auto&f=m.owned.f;std::vector<uint64_t>ids(f.coefficient.size());std::vector<int>owner(ids.size());for(size_t j=0;j<ids.size();j++)ids[j]=j;
    m.layout=std::make_unique<material_distributed::Ownership>(ids,owner,MPI_COMM_WORLD);
    // The evolved field is affine translation/shear and linear velocity.
    // Original analytic residual functions have ZERO state coefficients.
    // In particular this is not original Kepler orbit re-initialization.
    for(size_t j=0;j<f.polynomial.size();j++){
        double x=f.polynomial[j][0],y=f.polynomial[j][1];f.coefficient[j].x[0]=1.2*x+.1*y+.3;f.coefficient[j].x[1]=-.2*x+1.1*y-.4;
        f.coefficient[j].v[0]=-.7*y+.2;f.coefficient[j].v[1]=.7*x-.1;
    }
    m.generating=f.coefficient;
    for(size_t j=0;j<f.polynomial.size();j++){
        m.generating[j].x[0]=f.polynomial[j][0];m.generating[j].x[1]=f.polynomial[j][1];
        m.generating[j].v[0]=m.generating[j].v[1]=0;
    }
    double points[6][2]={{-1,0},{0,0},{0,.25},{0,-.25},{-.5,.5},{0,ldexp(1.,-40)}};
    for(auto&point:points)for(double&x:point)x+=4;
    int chart[6]={2,1,0,0,1,0};uint64_t source[6]={12,12,42,12,42,42};R maximum=0,maximum_enclosure=0;
    for(int j=0;j<6;j++){
        auto site=compact_source::transported_site(m,{0,1},points[j],3);
        need(site.canonical_bound,"actual represented source supplies canonical high/low/tail geometry");
        for(int k=0;k<2;k++)need(site.canonical_position_error[k]<1e-25&&site.canonical_velocity_error[k]<1e-25,"source evaluation keeps sub-double canonical precision for later geometry");
        need(site.cell==7&&site.chart==chart[j]&&site.source_triangle==source[j]&&site.source_mass==3,"exact source chart and deterministic minimum-ID tie");
        double x=points[j][0],y=points[j][1],position[2]={1.2*x+.1*y+.3,-.2*x+1.1*y-.4},velocity[2]={-.7*y+.2,.7*x-.1};
        for(int k=0;k<2;k++){
            maximum=std::max({maximum,fabsl(R(site.position[k])-position[k]),fabsl(R(site.velocity[k])-velocity[k])});
            maximum_enclosure=std::max({maximum_enclosure,R(site.position_error[k]),R(site.velocity_error[k])});
        }
    }
    auto current=sf::compile(f,0,f.coefficient),generating=sf::compile(f,0,m.generating);
    compact_source_fields::Law law;law.gamma=5./3;law.P0=.01;law.eps=.047;
    law.current_delta=.1;law.generating_delta=.05;
    sf::S X(sf::Q(4)),Y(sf::Q(4.25));
    auto field=compact_source_fields::evaluate(current,generating,X,Y,law);
    GfsEnrichConfig cfg={};cfg.field.ic.kind=2;cfg.field.ic.cx=cfg.field.ic.cy=12;
    cfg.field.ic.rin=2;cfg.field.ic.rout=8;cfg.field.ic.eps=law.eps;
    cfg.field.ic.gamma=law.gamma;cfg.field.ic.pressure=law.P0;
    double native[4];need(gfs_enrich_reference(&cfg,16,16.25,native,0,0),"native accepted radial density law");
    double j=1.34+law.current_delta,j0=1+law.generating_delta;
    double rho=native[0]/j,heat=law.P0/(law.gamma-1)*pow(j0/j,law.gamma);
    need(field.strictly_inside_source&&!field.native_field_error_certified&&!field.quadrature_error_certified&&
        fabsl(R(field.current_j.n)-j)<1e-12L&&fabsl(R(field.generating_j.n)-j0)<1e-12L&&
        fabsl(R(field.rho_physical.n)-rho)<1e-12L&&fabsl(R(field.heat_physical.n)-heat)<1e-12L,
        "selected generating-map coefficients and corrected physical rho/heat match accepted native formulas");
    Packet tiny;Record tr={};tr.cell=8;
    const double tiny_vertex[3][2]={{4,4},{4.125,4},{4,4.125}};
    for(int a=0;a<3;a++)for(int k=0;k<2;k++)tr.triangle.geometry.x[a][k]=tiny_vertex[a][k];
    tiny.record.push_back(tr);auto tiny_f=st::topology(tiny,1);
    for(size_t a=0;a<tiny_f.polynomial.size();a++){
        double x=tiny_f.polynomial[a][0],y=tiny_f.polynomial[a][1];
        tiny_f.coefficient[a].x[0]=1.2*x+.1*y+.3;
        tiny_f.coefficient[a].x[1]=-.2*x+1.1*y-.4;
        tiny_f.coefficient[a].v[0]=-.7*y+.2;
        tiny_f.coefficient[a].v[1]=.7*x-.1;
    }
    auto tiny_generating=tiny_f.coefficient;
    for(size_t a=0;a<tiny_f.polynomial.size();a++){
        tiny_generating[a].x[0]=tiny_f.polynomial[a][0];
        tiny_generating[a].x[1]=tiny_f.polynomial[a][1];
        tiny_generating[a].v[0]=tiny_generating[a].v[1]=0;
    }
    auto tiny_current_map=sf::compile(tiny_f,0,tiny_f.coefficient);
    auto tiny_initial_map=sf::compile(tiny_f,0,tiny_generating);
    auto tiny_bound=sf::tight_real_bounds(tiny_current_map);
    need(tiny_bound.L<1,"small synthetic affine source has a certified whole-triangle contraction bound");
    sf::S tiny_X(sf::Q(4.03125)),tiny_Y(sf::Q(4.03125));
    auto tiny_reference_field=compact_source_fields::evaluate(tiny_current_map,tiny_initial_map,tiny_X,tiny_Y,law);
    auto tiny_physical_field=compact_source_fields::evaluate_at_physical(tiny_current_map,tiny_initial_map,
        tiny_bound,tiny_reference_field.physical_x[0],tiny_reference_field.physical_x[1],law);
    auto inside=compact_source_fields::enclose_physical_box(tiny_current_map,tiny_bound,
        tiny_reference_field.physical_x[0],tiny_reference_field.physical_x[1]);
    sf::S wide_x(tiny_reference_field.physical_x[0].n,sf::Q(.5L));
    auto unresolved=compact_source_fields::enclose_physical_box(tiny_current_map,tiny_bound,
        wide_x,tiny_reference_field.physical_x[1]);
    need(inside.entire_physical_box_inside&&inside.inverse_radius>=0&&
        !unresolved.entire_physical_box_inside,
        "strict contained source box passes while a source-edge-crossing box remains unresolved");
    need(fabsl(R(tiny_physical_field.rho_physical.n)-R(tiny_reference_field.rho_physical.n))<1e-12L&&
        fabsl(R(tiny_physical_field.heat_physical.n)-R(tiny_reference_field.heat_physical.n))<1e-12L,
        "known source triangle physical inverse evaluates the same current physical field");
    const auto&first_edge=f.triangle[0].edge[0];
    int other_edge=-1;
    for(int e=0;e<3;e++)if(f.triangle[1].edge[e].poly[1]==first_edge.poly[1])other_edge=e;
    need(other_edge>=0,"two synthetic source triangles share the original horizontal edge");
    auto edge0=compact_power::source_edge_curve(f,0,0,41,11,3);
    auto edge1=compact_power::source_edge_curve(f,1,other_edge,41,11,3);
    for(double u:{.2,.8}){
        double native_position[2],native_velocity[2];
        need(st_trace(first_edge,f.coefficient.data(),u,native_position,native_velocity),
            "native current shared source edge trace");
        auto p0=compact_power::position(edge0,sf::exact(u));
        auto p1=compact_power::position(edge1,sf::exact(u));
        for(int k=0;k<2;k++)need(fabsl(R(p0[k].v.n)-native_position[k])<1e-12L&&
            fabsl(R(p1[k].v.n)-R(p0[k].v.n))<1e-24L,
            "both incident source triangles provide the same canonical current edge curve");
    }
    compact_power::Plane crossing={};crossing.cell[0]=10;crossing.cell[1]=20;
    crossing.origin[0]=sf::exact(5.55);crossing.normal[0]=sf::S(1);
    compact_power::RootBudget root_budget;root_budget.maximum_boxes=4096;
    root_budget.maximum_depth=48;root_budget.maximum_newton=48;
    root_budget.parameter_width=sf::Q(1e-12);
    auto root0=compact_power::intersections(crossing,edge0,root_budget);
    auto root1=compact_power::intersections(crossing,edge1,root_budget);
    need(root0.complete&&root1.complete&&root0.event.size()==1&&root1.event.size()==1&&
        root0.event[0].edge==41&&root1.event[0].edge==41&&
        root0.event[0].lower<sf::Q(.5208333333333333L)&&
        root0.event[0].upper>sf::Q(.5208333333333333L)&&
        root0.event[0].lower==root1.event[0].lower&&root0.event[0].upper==root1.event[0].upper,
        "one shared bounded power-plane root on the same current source edge");
    need(maximum<1e-12L&&maximum_enclosure<1e-10L,"actual represented source transport and canonical enclosure");
    std::cout<<std::setprecision(17)<<"{\"status\":\"SOURCE_VERTEX_EDGE_INTERIOR_TRANSPORT_PASS_SYNTHETIC_ONLY\",\"queries\":6,\"minimum_source_triangle_tie\":true,\"evolved_affine_field_used\":true,\"maximum_field_difference\":"<<double(maximum)<<",\"maximum_native_canonical_enclosure_error\":"<<double(maximum_enclosure)<<",\"actual_source_imported\":false,\"actual_target_source_inputs_executed\":false,\"new_physical_step\":false}\n";
}MPI_Finalize();return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);return 1;}}
