// Synthetic whole-boundary graph tests, not an actual-source transaction.
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
#include "material_compact_domain.h"
#include "material_compact_chart_jacobian.h"
#include "material_compact_chart_area.h"
#include "material_compact_chart_degree.h"
namespace cp=compact_power;namespace cd=compact_domain;
using Q=cp::Q;
static compact_target::Registry curved_square(){
    compact_target::Registry r;double site[5][2]={{-.55,-.43},{.47,-.51},{.55,.61},{-.51,.43},{.03,.02}};
    for(int j=0;j<5;j++){compact_source::Site s;s.cell=10+j;s.canonical_bound=true;for(int k=0;k<2;k++){s.reference[k]=s.position[k]=site[j][k];s.position_low[k]=ldexp(double(j%3-1),-60);s.position_error[k]=1e-12;s.canonical_position_error[k]=1e-31;}r.site.push_back(s);}
    double vertex[4][2]={{-1,-1},{1,-1},{1,1},{-1,1}};
    for(int j=0;j<4;j++){
        compact_source::Exterior e;e.edge=50+j;e.vertex[0]=100+j;e.vertex[1]=100+(j+1)%4;e.cell=r.site[j].cell;e.owner=0;e.orientation=j==3?-1:1;e.frozen.edge=e.edge;e.frozen.source_hash=11;e.frozen.source_epoch=5;e.frozen.function.rank=0;e.frozen.function.eps=.047;
        if(e.orientation==-1)std::swap(e.vertex[0],e.vertex[1]);
        for(int a=0;a<4;a++){
            e.frozen.function.poly[a]=a;R t=R(a)/3;if(e.orientation==-1)t=1-t;
            R x=(1-t)*vertex[j][0]+t*vertex[(j+1)%4][0],y=(1-t)*vertex[j][1]+t*vertex[(j+1)%4][1],bend=.125L*t*(1-t);
            if(j==0)y-=bend;if(j==1)x+=bend;if(j==2)y+=bend;if(j==3)x-=bend;
            e.frozen.coefficient[a].x[0]=double(x);e.frozen.coefficient[a].x[1]=double(y);
        }
        for(int k=0;k<2;k++){
            e.reference[0][k]=vertex[e.orientation==1?j:(j+1)%4][k];e.reference[1][k]=vertex[e.orientation==1?(j+1)%4:j][k];e.frozen.function.A[k]=e.reference[0][k];e.frozen.function.D[k]=e.reference[1][k]-e.reference[0][k];
            e.lower[k]=std::min(vertex[j][k],vertex[(j+1)%4][k]);e.upper[k]=std::max(vertex[j][k],vertex[(j+1)%4][k]);
        }
        if(j==0)e.lower[1]-=.05;if(j==1)e.upper[0]+=.05;if(j==2)e.upper[1]+=.05;if(j==3)e.lower[0]-=.05;r.exterior.push_back(e);
    }return r;
}
static compact_target::Registry u_domain(){
    compact_target::Registry r;for(int j=0;j<2;j++){compact_source::Site s;s.cell=10+j;s.position[0]=j?.625:0;s.position[1]=j?.5:-.5;r.site.push_back(s);}
    double v[8][2]={{-1,-1},{1,-1},{1,1},{.4,1},{.4,-.3},{-.4,-.3},{-.4,1},{-1,1}};
    for(int j=0;j<8;j++){
        compact_source::Exterior e;e.edge=50+j;e.vertex[0]=100+j;e.vertex[1]=100+(j+1)%8;e.cell=10;e.orientation=1;e.frozen.edge=e.edge;e.frozen.function.eps=.047;
        for(int a=0;a<4;a++){e.frozen.function.poly[a]=a;for(int k=0;k<2;k++)e.frozen.coefficient[a].x[k]=double(((R)(3-a)*v[j][k]+R(a)*v[(j+1)%8][k])/3);}
        for(int k=0;k<2;k++){e.reference[0][k]=v[j][k];e.reference[1][k]=v[(j+1)%8][k];e.lower[k]=std::min(v[j][k],v[(j+1)%8][k])-.001;e.upper[k]=std::max(v[j][k],v[(j+1)%8][k])+.001;}r.exterior.push_back(e);
    }return r;
}
int main(int argc,char**argv){MPI_Init(&argc,&argv);try{{
    int ranks=0;MPI_Comm_size(MPI_COMM_WORLD,&ranks);need(ranks==1,"one-rank finite whole-boundary graph fixture");
    cp::CellBudget cell_budget;cell_budget.maximum_nodes=4096;cell_budget.maximum_planes=512;cell_budget.maximum_corners=128;
    cd::Budget budget;budget.maximum_root_calls=4096;budget.maximum_winding_calls=256;budget.maximum_arcs=128;budget.root.maximum_boxes=2048;budget.root.maximum_depth=32;budget.root.maximum_newton=32;budget.root.parameter_width=ldexpq(1,-40);
    auto registry=curved_square();compact_target::verify_boundary(registry);std::vector<R>weight(5);cp::Spatial tree(registry,weight);cd::Domain domain(registry);
    struct Incidence {size_t count=0;int orientation=0;};std::map<std::pair<cd::Key,cd::Key>,Incidence>interior;std::vector<cd::Graph>graphs;
    size_t curves=0,arcs=0,root_calls=0,root_boxes=0,winding=0;
    size_t chart_patches=0,chart_curved_patches=0,chart_midpoints=0;std::vector<compact_chart::Candidate>chart_candidate;
    compact_chart::Budget chart_budget;chart_budget.maximum_vertices=128;chart_budget.maximum_ears=128;chart_budget.maximum_patches=512;
    auto check_chart=[&](const cd::Domain&d,const cd::Graph&g){
        auto candidate=compact_chart::connect(d,g,chart_budget);chart_patches+=candidate.patch.size();chart_midpoints+=candidate.curved_midpoints;
        struct Edge {size_t count=0,curved=0,boundary=0;int orientation=0;};std::map<std::pair<compact_chart::Key,compact_chart::Key>,Edge>edges;cp::S area;
        for(const auto&p:candidate.patch){
            int curved=0;area=area+compact_chart::orientation(p.vertex[0],p.vertex[1],p.vertex[2]);
            for(int k=0;k<3;k++){
                curved+=p.edge[k].curved;auto a=p.vertex[k].key,b=p.vertex[(k+1)%3].key;int direction=1;if(b<a){std::swap(a,b);direction=-1;}auto&e=edges[{a,b}];e.count++;e.orientation+=direction;e.curved+=p.edge[k].curved;e.boundary+=p.edge[k].boundary;
                if(p.edge[k].curved)need(p.edge[k].endpoint[0]==p.vertex[k].key&&p.edge[k].endpoint[1]==p.vertex[(k+1)%3].key,"same implicit endpoint function on curved patch edge");
            }
            need(curved<=1,"at most one same-source curved edge per collar fan patch");chart_curved_patches+=curved;
        }
        for(const auto&e:edges){need((e.second.count==1&&e.second.orientation!=0&&e.second.boundary==1)||(e.second.count==2&&e.second.orientation==0&&e.second.curved==0&&e.second.boundary==0),"physical straight/curved boundary tags retained; internal fan edges occur in opposite pairs");}
        need(cp::sf::lower(area)>0,"positive total candidate reference area, not a curved physical-area proof");
        return candidate;
    };
    for(size_t i=0;i<5;i++){
        auto cell=cp::bounded_cell(tree,i,cell_budget);auto graph=cd::intersect(domain,tree,i,cell,budget);need(graph.component.size()==1&&!graph.generator_anchor_used,"one clipped synthetic cell, no unbound actual-source shortcut");
        for(const auto&a:graph.arc){
            if(a.curved)curves++;
            else{
                need(a.line.kind==1,"no artificial enclosing-box edge survives as physical boundary");cd::Key from=graph.vertex[a.from].key,to=graph.vertex[a.to].key;int orientation=1;if(to<from){std::swap(from,to);orientation=-1;}auto&entry=interior[{from,to}];entry.count++;entry.orientation+=orientation;
            }
        }
        arcs+=graph.arc.size();root_calls+=graph.root_calls;root_boxes+=graph.root_boxes;winding+=graph.winding_calls;graphs.push_back(std::move(graph));
        chart_candidate.push_back(check_chart(domain,graphs.back()));
    }
    for(const auto&entry:interior)need(entry.second.count==2&&entry.second.orientation==0,"each physical interior straight arc appears exactly twice with opposite orientation");
    // Join exterior pieces by exact endpoint/root identities, not by a gap
    // tolerance on rounded parameters. Keep interval widths and both traces.
    for(size_t j=0;j<registry.exterior.size();j++){
        struct Piece {Q lo,hi;cd::Key from,to;Q alo,ahi,blo,bhi;};std::vector<Piece>pieces;
        for(const auto&graph:graphs)for(const auto&a:graph.arc)if(a.curved&&a.source_curve==j){
            auto from=graph.vertex[a.from].key,to=graph.vertex[a.to].key;Q alo=a.parameter_lower[0],ahi=a.parameter_upper[0],blo=a.parameter_lower[1],bhi=a.parameter_upper[1];
            if(registry.exterior[j].orientation==-1){std::swap(from,to);std::swap(alo,blo);std::swap(ahi,bhi);}pieces.push_back({alo,bhi,from,to,alo,ahi,blo,bhi});
        }
        need(!pieces.empty(),"every original curved exterior retained");std::sort(pieces.begin(),pieces.end(),[](const Piece&a,const Piece&b){return a.lo<b.lo;});
        need(pieces.front().from==cd::endpoint_key(registry.exterior[j],0)&&pieces.back().to==cd::endpoint_key(registry.exterior[j],1)&&pieces.front().alo==0&&pieces.back().bhi==1,"whole canonical curve covered from original endpoint to endpoint");
        for(size_t k=1;k<pieces.size();k++)need(pieces[k-1].to==pieces[k].from&&pieces[k-1].blo==pieces[k].alo&&pieces[k-1].bhi==pieces[k].ahi,"same implicit root joins neighboring curved pieces without discarded strip or duplicate interval");
    }
    auto chart_certificate=compact_chart::certify_global_connectivity(registry,chart_candidate,0,512);
    need(chart_certificate.global_components==1&&chart_certificate.exact_C0_connectivity&&chart_certificate.oriented_source_boundary_complete&&
        chart_certificate.manifold_vertex_stars&&chart_certificate.disk_topology&&chart_certificate.euler_characteristic==1&&
        !chart_certificate.curved_J_positive&&!chart_certificate.degree_one_coverage&&chart_certificate.geometry_modes==2,
        "global chart edge pairing, connected disk topology and complete oriented source boundary; Jacobian/degree remain unclaimed");
    need(chart_certificate.chart_hash==compact_chart::chart_collection_hash(chart_candidate),"connectivity certificate seals the complete canonical chart candidate set");
    auto altered_chart=chart_candidate;altered_chart[0].patch[0].vertex[0].x[0].n+=ldexpq(1,-80);bool stale_chart_rejected=false;
    try{compact_chart::require_same_chart_complex(registry,altered_chart,chart_certificate);}catch(const std::exception&){stale_chart_rejected=true;}
    need(stale_chart_rejected,"Jacobian stage rejects a chart set changed after global connectivity certification");
    compact_chart::JacobianBudget jacobian_budget;jacobian_budget.maximum_boxes_per_patch=4096;jacobian_budget.maximum_depth=32;
    auto chart_jacobian=compact_chart::certify_all_chart_jacobians(registry,chart_candidate,jacobian_budget);
    need(chart_jacobian.all_reference_J_positive&&chart_jacobian.patch_lower.size()==chart_certificate.patch.size()&&chart_jacobian.curved_patches==chart_curved_patches,
        "continuous interval Jacobian lower bound on every synthetic whole-reference chart; no sampled-J acceptance");
    compact_chart::AreaBudget area_budget;area_budget.maximum_interior_boxes_per_patch=65536;
    area_budget.maximum_boundary_boxes_per_patch=65536;area_budget.maximum_depth=32;
    area_budget.maximum_absolute_area_error_per_patch=2e-2Q;
    compact_chart::AreaBuilder area_builder(registry,chart_candidate,chart_certificate,chart_jacobian,area_budget);
    compact_chart::S interior_area,boundary_area,gauss8_area,gauss12_area;
    size_t area_nodes=0,area_boxes=0,boundary_boxes=0,gauss8_nodes=0,gauss12_nodes=0;
    for(size_t p=0;p<chart_certificate.patch.size();p++){
        auto rule=area_builder.patch_rule(p);
        need(rule.positive_nodes&&rule.independent_area_agreement&&!rule.continuous_physical_coverage&&
            rule.chart_hash==chart_certificate.chart_hash&&rule.source_epoch==registry.source_epoch,
            "positive target-area nodes and independent boundary integral bind to the same synthetic chart");
        interior_area=interior_area+rule.interior_area;boundary_area=boundary_area+rule.boundary_area;
        area_nodes+=rule.node.size();area_boxes+=rule.interior_boxes;boundary_boxes+=rule.boundary_boxes;
        auto g8=area_builder.positive_gauss_rule(p,8),g12=area_builder.positive_gauss_rule(p,12);
        need(g8.positive_reference_and_physical_weights&&g12.positive_reference_and_physical_weights&&
            !g8.field_integration_error_certified&&!g12.field_integration_error_certified&&
            g8.node.size()==64&&g12.node.size()==144&&g8.chart_hash==chart_certificate.chart_hash&&
            g12.chart_hash==chart_certificate.chart_hash,
            "positive bounded-order Gauss reference/physical rule follows the same sealed synthetic curved chart");
        gauss8_area=gauss8_area+g8.physical_area;gauss12_area=gauss12_area+g12.physical_area;
        gauss8_nodes+=g8.node.size();gauss12_nodes+=g12.node.size();
    }
    Q analytic_curved_square_area=Q(25)/Q(6);
    need(compact_chart::sf::lower(interior_area)<=analytic_curved_square_area&&
        compact_chart::sf::upper(interior_area)>=analytic_curved_square_area&&
        compact_chart::sf::lower(boundary_area)<=analytic_curved_square_area&&
        compact_chart::sf::upper(boundary_area)>=analytic_curved_square_area,
        "independent positive interior and Green-boundary area enclosures contain the analytic curved-square area");
    need(fabsq(gauss8_area.n-analytic_curved_square_area)<1e-10Q&&
        fabsq(gauss12_area.n-analytic_curved_square_area)<1e-10Q&&
        fabsq(gauss8_area.n-gauss12_area.n)<1e-11Q,
        "positive 8x8 and 12x12 curved-chart Gauss rules recover the analytic synthetic area and agree");
    compact_source::GeometryBound unproven_source_geometry;bool unbound_degree_rejected=false;
    try{compact_chart::certify_degree_one(registry,chart_candidate,chart_certificate,chart_jacobian,unproven_source_geometry);}catch(const std::exception&){unbound_degree_rejected=true;}
    need(unbound_degree_rejected,"degree-one coverage rejects charts without the accepted source global-injectivity certificate");
    auto cap_registry=curved_square();cap_registry.site.resize(2);cap_registry.site[0].position[0]=.75;cap_registry.site[1].position[0]=1;for(auto&s:cap_registry.site)s.position[1]=0;
    std::vector<R>cap_weight={.0703125L,0};cp::Spatial cap_tree(cap_registry,cap_weight);cd::Domain cap_domain(cap_registry);auto cap_cell=cp::bounded_cell(cap_tree,1,cell_budget);auto cap=cd::intersect(cap_domain,cap_tree,1,cap_cell,budget);
    need(cap.component.size()==1&&cap.arc.size()==2&&cap.vertex.size()==2,"valid two-arc curved cap is retained, not rejected for lacking polygon star shape");
    auto cap_candidate=check_chart(cap_domain,cap);need(cap_candidate.curved_midpoints==1&&cap_candidate.patch.size()==3,"one cap midpoint produces a positive three-patch reference fan");
    auto nonconvex=u_domain();std::vector<R>uw(2);cp::Spatial utree(nonconvex,uw);cd::Domain udomain(nonconvex);auto ucell=cp::bounded_cell(utree,1,cell_budget);auto disjoint=cd::intersect(udomain,utree,1,ucell,budget);
    need(disjoint.component.size()==2,"one convex power cell intersect nonconvex domain yields two retained components");
    bool caught=false;auto broken=graphs[0];broken.arc.pop_back();try{cd::close(broken);}catch(const std::exception&){caught=true;}need(caught,"open physical boundary graph rejected");
    caught=false;try{compact_source::GeometryBound g;domain.bind_actual_source(g);}catch(const std::exception&){caught=true;}need(caught,"unproven source cannot enable transported-generator anchor shortcut");
    caught=false;auto stale=graphs[0];stale.source_epoch++;try{compact_chart::connect(domain,stale,chart_budget);}catch(const std::exception&){caught=true;}need(caught,"stale source graph rejected before chart connectivity");
    caught=false;auto stale_cell=cp::bounded_cell(tree,0,cell_budget);stale_cell.weight_hash++;try{cd::intersect(domain,tree,0,stale_cell,budget);}catch(const std::exception&){caught=true;}need(caught,"stale weight-trial cell rejected before curved clipping");
    std::vector<compact_chart::Node>bowtie(4);double xy[4][2]={{0,0},{1,1},{0,1},{1,0}};for(int j=0;j<4;j++){bowtie[j].key={uint64_t(j)};for(int k=0;k<2;k++)bowtie[j].x[k]=cp::S(cp::Q(xy[j][k]));}
    caught=false;try{compact_chart::polygon_guard(bowtie);}catch(const std::exception&){caught=true;}need(caught,"self-crossing chord polygon rejected");
    std::cout<<std::setprecision(17)<<"{\"status\":\"POSITIVE_CURVED_GAUSS_AREA_PASS_SYNTHETIC_ONLY\",\"gauss8_nodes\":"<<gauss8_nodes
        <<",\"gauss12_nodes\":"<<gauss12_nodes<<",\"gauss8_area_center\":"<<(double)gauss8_area.n
        <<",\"gauss12_area_center\":"<<(double)gauss12_area.n<<",\"analytic_curved_square_area\":"<<(double)analytic_curved_square_area
        <<",\"field_integration_error_certified\":false,\"actual_source_imported\":false,\"new_physical_step\":false}\n";
    std::cout<<std::setprecision(17)<<"{\"status\":\"CURVED_CHART_AREA_AND_JACOBIAN_PASS_SYNTHETIC_ONLY\",\"partition_cells\":5,\"retained_arcs\":"<<arcs<<",\"curved_arcs\":"<<curves<<",\"shared_straight_arcs\":"<<interior.size()<<",\"root_calls\":"<<root_calls<<",\"root_boxes\":"<<root_boxes<<",\"winding_calls\":"<<winding<<",\"chart_patches_including_cap\":"<<chart_patches<<",\"curved_chart_patches\":"<<chart_curved_patches<<",\"chart_midpoints\":"<<chart_midpoints<<",\"global_chart_patches\":"<<chart_certificate.patch.size()<<",\"global_internal_edges\":"<<chart_certificate.internal_edges<<",\"global_shared_power_edges\":"<<chart_certificate.straight_shared_edges<<",\"global_components\":"<<chart_certificate.global_components<<",\"global_unique_nodes\":"<<chart_certificate.unique_nodes<<",\"global_unique_edges\":"<<chart_certificate.unique_edges<<",\"global_boundary_edges\":"<<chart_certificate.boundary_edges<<",\"global_vertex_stars\":"<<chart_certificate.vertex_stars<<",\"global_euler_characteristic\":"<<chart_certificate.euler_characteristic<<",\"global_disk_topology_pass\":true,\"global_geometry_modes\":"<<chart_certificate.geometry_modes<<",\"global_exact_C0_connectivity_pass\":true,\"chart_candidate_hash_replay_rejected\":true,\"whole_oriented_source_boundary_pass\":true,\"whole_reference_J_positive\":true,\"minimum_reference_J_lower_bound\":"<<(double)chart_jacobian.minimum_lower<<",\"Jacobian_interval_boxes\":"<<chart_jacobian.boxes<<",\"Jacobian_interval_subdivisions\":"<<chart_jacobian.subdivisions<<",\"positive_target_area_nodes\":"<<area_nodes<<",\"target_area_interval_boxes\":"<<area_boxes<<",\"target_boundary_interval_boxes\":"<<boundary_boxes<<",\"target_area_lower\":"<<(double)compact_chart::sf::lower(interior_area)<<",\"target_area_upper\":"<<(double)compact_chart::sf::upper(interior_area)<<",\"independent_boundary_area_lower\":"<<(double)compact_chart::sf::lower(boundary_area)<<",\"independent_boundary_area_upper\":"<<(double)compact_chart::sf::upper(boundary_area)<<",\"analytic_curved_square_area\":"<<(double)analytic_curved_square_area<<",\"unbound_source_geometry_degree_gate_rejected\":true,\"degree_one_physical_coverage\":false,\"whole_exterior_implicit_root_joins_pass\":true,\"two_arc_cap_pass\":true,\"disconnected_components_retained\":2,\"open_graph_unbound_anchor_stale_graph_and_weights_rejected\":true,\"actual_source_imported\":false,\"actual_partition_proven\":false,\"actual_projection_executed\":false,\"new_physical_step\":false}\n";
}MPI_Finalize();return 0;}catch(const std::exception&e){std::cerr<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);return 1;}}
