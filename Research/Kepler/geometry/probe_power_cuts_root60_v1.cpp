// Actual step5 target chart diagnostic. A deterministic, very small weight
// perturbation tests a nondegenerate trial geometry. This is not a solution
// of the screened fixed-mass weight law.
#include "../material_compact_rezone_parallel_v2.cpp"
#include <atomic>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <quadmath.h>

static std::string exact_quad_text(__float128 value){
    char text[128];
    compact_p3::require(quadmath_snprintf(text,sizeof(text),"%.36Qg",value)>0,
        "serialize exact binary128 event endpoint in roundtrip decimal");
    return text;
}
static std::vector<std::string> cut_csv_fields(const std::string&line){
    std::vector<std::string>out;
    std::stringstream stream(line);std::string field;
    while(std::getline(stream,field,','))out.push_back(field);
    return out;
}
static __float128 cut_parse_quad(const std::string&word){
    char*end=nullptr;__float128 value=strtoflt128(word.c_str(),&end);
    compact_p3::require(end==word.c_str()+word.size()&&finiteq(value)&&
        exact_quad_text(value)==word,
        "lossless finite binary128 canonical cut parameter");
    return value;
}

static R diagnostic_weight(uint64_t id){
    // SplitMix64 supplies a stable nonlocal tie break for the nearly
    // cocircular transported lattice. The entire trial range is 1e-18;
    // these values are sealed exactly as binary80 law inputs downstream.
    uint64_t z=id+0x9e3779b97f4a7c15ULL;
    z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;
    z=(z^(z>>27))*0x94d049bb133111ebULL;
    z^=z>>31;
    return (R(z)/R(UINT64_MAX)-R(.5L))*R(1e-18L);
}

static compact_chart::S diagnostic_green_gauss(
        const compact_chart::Patch&patch,const compact_power::Curve*curve,
        const std::vector<double>&nodes,const std::vector<double>&weights){
    namespace cc=compact_chart;
    const auto&A=patch.vertex[0].x,&B=patch.vertex[1].x,&C=patch.vertex[2].x;
    cc::S area=(cc::area_cross(A,B)+cc::area_cross(B,C)+cc::area_cross(C,A))/cc::S(2);
    if(curve)for(size_t k=0;k<nodes.size();k++)
        area=area+cc::sf::exact(weights[k])*
            cc::boundary_integrand(patch,*curve,cc::Q(nodes[k]),cc::Q(nodes[k]));
    return area;
}

static compact_source::GeometryBound load_geometry(
        const compact_source::Model&source,const std::string&directory){
    int rank=source.layout->rank;
    std::string path=directory+"/geometry_r"+std::to_string(rank)+".bin";
    FILE*file=fopen(path.c_str(),"rb");
    compact_p3::require(file,"previous actual step5 source geometry witness readable");
    auto read=[&](void*data,size_t size){
        compact_p3::require(fread(data,1,size,file)==size,
            "complete previous actual source geometry witness");
    };
    uint64_t header[9];read(header,sizeof(header));
    const auto&o=source.owned;
    compact_p3::require(header[0]==0x4c46414354534731ULL&&header[1]==1&&
        header[2]==uint64_t(rank)&&header[3]==4&&
        header[4]==o.f.triangle.size()&&
        header[5]==sizeof(shared_function::RealBounds)&&
        header[6]==source.accepted.coefficient_hash&&
        header[7]==source.accepted.physical.epoch&&
        header[8]==source.accepted.physical.steps,
        "same current coefficient/ABI/owner/epoch for retained geometry proof");
    compact_source::GeometryBound g;
    g.source_coefficient_hash=header[6];g.source_epoch=header[7];g.source_steps=header[8];
    read(&g.source_time,sizeof(g.source_time));
    read(&g.maximum_Dchi_minus_I,sizeof(g.maximum_Dchi_minus_I));
    read(&g.minimum_canonical_J,sizeof(g.minimum_canonical_J));
    compact_p3::require(g.source_time==source.accepted.physical.time&&
        g.maximum_Dchi_minus_I>=0&&g.maximum_Dchi_minus_I<1&&
        g.minimum_canonical_J>0,
        "same accepted time and strict previous global source geometry bounds");
    g.triangle.resize(header[4]);g.complete.resize(header[4]);
    read(g.triangle.data(),g.triangle.size()*sizeof(shared_function::RealBounds));
    read(g.complete.data(),g.complete.size());
    int trailing=fgetc(file),closed=fclose(file);
    compact_p3::require(trailing==EOF&&closed==0,
        "no trailing or truncated saved source geometry witness");
    R local_max=0,local_min=std::numeric_limits<R>::infinity();
    for(size_t s=0;s<g.triangle.size();s++){
        const auto&b=g.triangle[s];
        compact_p3::require(g.complete[s]==1&&finiteq(b.L)&&b.L>=0&&b.L<1,
            "every retained source triangle has a complete contraction bound");
        for(auto displacement:b.displacement)
            compact_p3::require(finiteq(displacement)&&displacement>=0,
                "finite nonnegative retained source displacement");
        using S=shared_function::S;
        auto determinant=(S(1)-S(b.FminusI[0]))*(S(1)-S(b.FminusI[3]))-
            S(b.FminusI[1])*S(b.FminusI[2]);
        auto lower=shared_function::lower(determinant);
        compact_p3::require(finiteq(lower)&&lower>0,
            "strict positive retained source triangle Jacobian lower bound");
        local_max=std::max(local_max,nextafterl(R(b.L),INFINITY));
        local_min=std::min(local_min,nextafterl(R(lower),-INFINITY));
    }
    MPI_Allreduce(MPI_IN_PLACE,&local_max,1,MPI_LONG_DOUBLE,MPI_MAX,source.layout->comm);
    MPI_Allreduce(MPI_IN_PLACE,&local_min,1,MPI_LONG_DOUBLE,MPI_MIN,source.layout->comm);
    compact_p3::require(local_max==g.maximum_Dchi_minus_I&&
        local_min==g.minimum_canonical_J,
        "recomputed global extrema match retained same-state geometry proof");
    g.external_edges=compact_source::check_shared_edges(source);
    g.exact_shared_trace=true;g.canonical_global_injective=true;
    return g;
}

int main(int argc,char**argv){
    int provided=0;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
    int rank=0,ranks=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    try{
        compact_p3::require((argc==8||argc==9||argc==10||argc==11||argc==13||argc==14)&&ranks==4&&provided>=MPI_THREAD_FUNNELED,
            "actual weighted trial target diagnostic needs source, geometry and optional site/weight input");
        bool external_weights=argc>=9&&std::filesystem::is_regular_file(argv[8]);
        bool certified_full=argc==11&&std::string(argv[10])=="all_certified";
        bool whole_patch_test=argc==11&&std::string(argv[10])=="whole_patch";
        bool boundary_index_test=argc==11&&std::string(argv[10])=="boundary_index";
        bool boundary_edge_test=argc==11&&std::string(argv[10])=="boundary_edges";
        bool boundary_scale_test=argc==11&&std::string(argv[10])=="boundary_edges_1024";
        bool boundary_failed_test=argc==11&&std::string(argv[10])=="boundary_failed_six";
        bool boundary_parallel_prefix=argc==11&&std::string(argv[10])=="boundary_parallel_prefix";
        bool boundary_parallel_range=argc==13&&std::string(argv[10])=="boundary_parallel_range";
        bool boundary_cut_range=argc==14&&std::string(argv[10])=="boundary_cut_range";
        bool boundary_witness_test=argc==11&&std::string(argv[10])=="boundary_witness";
        size_t single_site=argc==9&&!external_weights?std::stoull(argv[8]):
            argc==11&&!certified_full&&!whole_patch_test&&
                !boundary_index_test&&!boundary_edge_test&&
                !boundary_scale_test&&!boundary_failed_test&&
                !boundary_parallel_prefix&&!boundary_witness_test?
                std::stoull(argv[10]):SIZE_MAX;
        compact_p3::require(argc<10||external_weights,
            "damped weight direction needs a complete external weight ledger");
        R damping=argc>=10?std::stold(argv[9]):R(1);
        compact_p3::require(std::isfinite(damping)&&damping>0&&damping<=1,
            "finite positive Newton damping no greater than one");
        compact_source::Inputs input;
        input.initial=argv[1];input.immutable_manifest=argv[2];
        input.positive_directory=argv[3];input.whole_audit=argv[4];
        input.reference_receipt=argv[5];input.accepted_prefix=argv[6];
        input.accepted_steps=5;input.accepted_time=0.00011263377801805954;
        const double begun=MPI_Wtime();
        {
            auto source=compact_source::import_accepted_state(input,MPI_COMM_WORLD);
            auto geometry=load_geometry(source,argv[7]);
            if(rank==0)std::cout<<"ACTUAL_STEP5_CACHED_GEOMETRY_PASS seconds="
                <<MPI_Wtime()-begun<<std::endl;
            auto context=compact_source::derive_accepted_context(source,input,geometry);
            auto target=compact_source::extract_target_inputs(source,input,geometry,context);
            auto registry=compact_target::gather(target,MPI_COMM_WORLD);
            auto law=compact_mass_law::bind_actual(input,source,geometry,context,
                registry,MPI_COMM_WORLD);
            if(rank==0)std::cout<<"ACTUAL_STEP5_PRETARGET_PASS sites="
                <<registry.site.size()<<" seconds="<<MPI_Wtime()-begun<<std::endl;
            if(rank==0||whole_patch_test||boundary_parallel_range||boundary_cut_range){
                namespace cp=compact_power;namespace cd=compact_domain;
                namespace cc=compact_chart;
                std::vector<R>weights(registry.site.size());
                if(external_weights){
                    std::ifstream file(argv[8]);compact_p3::require(bool(file),
                        "read complete external Newton direction ledger");
                    std::string line;compact_p3::require(bool(std::getline(file,line)),
                        "external Newton direction header");
                    for(size_t i=0;i<weights.size();i++){
                        compact_p3::require(bool(std::getline(file,line)),
                            "one external Newton direction per actual cell");
                        std::stringstream stream(line);std::vector<std::string>field;
                        std::string value;
                        while(std::getline(stream,value,','))field.push_back(value);
                        compact_p3::require(field.size()==6&&
                            std::stoull(field[0])==registry.site[i].cell,
                            "sorted complete original cell weight direction");
                        R initial=std::stold(field[1]),delta=std::stold(field[2]);
                        compact_p3::require(initial==diagnostic_weight(registry.site[i].cell)&&
                            std::isfinite(delta),
                            "external Newton direction starts at same sealed trial weights");
                        weights[i]=initial+damping*delta;
                        compact_p3::require(std::isfinite(weights[i]),
                            "finite damped power weight for every actual cell");
                    }
                    compact_p3::require(!std::getline(file,line),
                        "no extra external weight entries");
                    std::cout<<std::setprecision(21)
                        <<"ACTUAL_STEP5_EXTERNAL_NEWTON_WEIGHTS damping="<<damping
                        <<" fixed_mass_law_solution=0"<<std::endl;
                }else{
                    for(size_t i=0;i<weights.size();i++)
                        weights[i]=diagnostic_weight(registry.site[i].cell);
                    std::cout<<"ACTUAL_STEP5_DIAGNOSTIC_WEIGHT_RANGE [-5e-19,5e-19]"
                        <<" fixed_mass_law_solution=0"<<std::endl;
                }
                cp::Spatial sites(registry,weights);
                cd::Domain domain(registry);domain.bind_actual_source(geometry);
                cp::CellBudget cell_budget;
                cell_budget.maximum_nodes=65536;
                cell_budget.maximum_planes=single_site==SIZE_MAX?8192:71080;
                cell_budget.maximum_corners=512;
                cd::Budget domain_budget;
                domain_budget.maximum_root_calls=8192;
                domain_budget.maximum_winding_calls=4096;
                domain_budget.maximum_arcs=4096;
                domain_budget.root.maximum_boxes=8192;
                domain_budget.root.maximum_depth=48;
                domain_budget.root.maximum_newton=48;
                domain_budget.root.parameter_width=ldexpq(1,-40);
                cc::Budget chart_budget;
                chart_budget.maximum_vertices=4096;
                chart_budget.maximum_ears=4096;
                chart_budget.maximum_patches=16384;
                std::vector<cc::Candidate>candidate(registry.site.size());
                std::atomic<bool>failed(false);
                std::atomic<size_t>progress(0);
                size_t completed=0,patches=0,root_calls=0,first=SIZE_MAX;
                std::string reason;
                #pragma omp parallel for schedule(dynamic,1) reduction(+:completed,patches,root_calls)
                for(size_t i=0;i<registry.site.size();i++){
                    if(single_site!=SIZE_MAX&&i!=single_site)continue;
                    if(failed.load(std::memory_order_relaxed))continue;
                    try{
                        auto cell=cp::bounded_cell(sites,i,cell_budget);
                        if(i==single_site)std::cout
                            <<"ACTUAL_WEIGHTED_SELECTED_POWER index="<<i
                            <<" cell="<<cell.cell
                            <<" nodes="<<cell.nodes
                            <<" planes="<<cell.tested_sites.size()
                            <<" excluded="<<cell.excluded.size()
                            <<" corners="<<cell.corner.size()
                            <<" seconds="<<MPI_Wtime()-begun<<std::endl;
                        if(i==single_site){
                            const auto&q=registry.site[i];
                            std::cout<<std::setprecision(21)
                                <<"ACTUAL_WEIGHTED_SELECTED_SITE q0="<<q.position[0]
                                <<" q1="<<q.position[1]<<std::endl;
                            for(size_t k=0;k<cell.corner.size();k++){
                                const auto&v=cell.corner[k];
                                std::cout<<"ACTUAL_WEIGHTED_SELECTED_CORNER k="<<k
                                    <<" x0="<<R(v.x[0].n)
                                    <<" e0="<<R(v.x[0].e)
                                    <<" x1="<<R(v.x[1].n)
                                    <<" e1="<<R(v.x[1].e)<<std::endl;
                            }
                        }
                        auto graph=cd::intersect_actual_domain(domain,sites,i,cell,
                            domain_budget);
                        if(i==0){
                            std::cout<<std::setprecision(21)
                                <<"ACTUAL_WEIGHTED_CELL0_INPUT corners="<<cell.corner.size()
                                <<" vertices="<<graph.vertex.size()
                                <<" arcs="<<graph.arc.size()
                                <<" components="<<graph.component.size()
                                <<" candidates="<<graph.curve_candidates
                                <<" roots="<<graph.root_calls<<std::endl;
                            for(size_t j=0;j<graph.vertex.size();j++){
                                const auto&v=graph.vertex[j];
                                std::cout<<"ACTUAL_WEIGHTED_CELL0_VERTEX j="<<j
                                    <<" key_kind="<<v.key.kind;
                                for(auto id:v.key.id)std::cout<<" key="<<id;
                                for(int d=0;d<2;d++)std::cout<<" x"<<d
                                    <<"="<<R(v.x[d].n)<<" radius"<<d
                                    <<"="<<R(v.x[d].e);
                                std::cout<<std::endl;
                            }
                            for(size_t j=0;j<graph.arc.size();j++){
                                const auto&a=graph.arc[j];
                                std::cout<<"ACTUAL_WEIGHTED_CELL0_ARC j="<<j
                                    <<" from="<<a.from<<" to="<<a.to
                                    <<" curved="<<int(a.curved)
                                    <<" line_kind="<<a.line.kind
                                    <<" line_a="<<a.line.a<<" line_b="<<a.line.b
                                    <<" source_curve="<<a.source_curve<<std::endl;
                            }
                        }
                        auto chart=cc::connect_actual_domain(domain,sites,graph,chart_budget);
                        compact_p3::require(!chart.patch.empty(),
                            "actual source generator has a nonempty target chart");
                        root_calls+=graph.root_calls;patches+=chart.patch.size();
                        candidate[i]=std::move(chart);completed++;
                        size_t done=progress.fetch_add(1)+1;
                        if(done%5000==0)std::cout<<"ACTUAL_STEP5_WEIGHTED_PROGRESS completed="
                            <<done<<" seconds="<<MPI_Wtime()-begun<<std::endl;
                    }catch(const std::exception&e){
                        failed.store(true,std::memory_order_relaxed);
                        #pragma omp critical
                        {if(i<first){first=i;reason=e.what();}}
                    }
                }
                std::cout<<"ACTUAL_STEP5_WEIGHTED_CHART_BUILD completed="<<completed
                    <<" total="<<registry.site.size()<<" patches="<<patches
                    <<" root_calls="<<root_calls<<" first_failure_index="<<first
                    <<" first_failure_cell="
                    <<(first<registry.site.size()?registry.site[first].cell:UINT64_MAX)
                    <<" reason="<<reason<<" seconds="<<MPI_Wtime()-begun
                    <<std::endl;
                compact_p3::require(!failed.load()&&
                    completed==(single_site==SIZE_MAX?registry.site.size():size_t(1)),
                    "all actual weighted trial target chart candidates constructed");
                if(single_site!=SIZE_MAX){
                    if(single_site==0&&candidate[0].patch.size()>4){
                        const auto&p=candidate[0].patch[4];
                        const auto&edge=p.edge[1];
                        auto curve=cp::curve(registry.exterior.at(edge.source));
                        std::cout<<"ACTUAL_SELECTED_PATCH4_GEOMETRY rank="
                            <<curve.scalar.native.rank<<" source="<<edge.source
                            <<" edge="<<curve.edge<<" u0="<<R(edge.parameter[0].n)
                            <<" u1="<<R(edge.parameter[1].n)<<std::endl;
                        for(int vi=0;vi<3;vi++)std::cout<<std::setprecision(21)
                            <<"ACTUAL_SELECTED_PATCH4_VERTEX i="<<vi
                            <<" x="<<R(p.vertex[vi].x[0].n)
                            <<" y="<<R(p.vertex[vi].x[1].n)
                            <<" ex="<<R(p.vertex[vi].x[0].e)
                            <<" ey="<<R(p.vertex[vi].x[1].e)<<std::endl;
                        for(int si=0;si<=4;si++)for(int ti=0;ti<=8;ti++){
                            cc::Q s=cc::Q(si)/4,t=cc::Q(ti)/8;
                            auto J=cc::chart_J_interval(p,curve,edge.source,s,s,t,t);
                            std::cout<<std::setprecision(15)
                                <<"ACTUAL_SELECTED_PATCH4_J s="<<R(s)
                                <<" t="<<R(t)
                                <<" lower="<<R(cc::sf::lower(J))
                                <<" upper="<<R(cc::sf::upper(J))<<std::endl;
                        }
                    }
                    cc::JacobianBudget selected_budget;
                    selected_budget.maximum_boxes_per_patch=4096;
                    selected_budget.maximum_depth=32;
                    auto selected_jacobian=cc::certify_chart_jacobians(
                        registry,candidate[single_site],selected_budget);
                    std::cout<<"ACTUAL_WEIGHTED_SINGLE_SITE_PASS index="<<single_site
                        <<" minimum_J="<<R(selected_jacobian.minimum_lower)
                        <<" seconds="<<MPI_Wtime()-begun<<std::endl;
                    cc::AreaBudget selected_area_budget;
                    selected_area_budget.maximum_boundary_boxes_per_patch=65536;
                    selected_area_budget.maximum_absolute_area_error_per_patch=1e-11Q;
                    cc::S selected_area;
                    size_t area_boxes=0,area_curved=0;
                    for(size_t p=0;p<candidate[single_site].patch.size();p++){
                        const auto&patch=candidate[single_site].patch[p];
                        cp::Curve descriptor;const cp::Curve*curve=nullptr;
                        if(patch.edge[1].curved){
                            descriptor=cp::curve(registry.exterior.at(patch.edge[1].source));
                            curve=&descriptor;area_curved++;
                        }
                        size_t boxes=0,subdivisions=0;
                        try{
                            auto patch_area=cc::boundary_area(patch,curve,selected_area_budget,
                                boxes,subdivisions);
                            selected_area=selected_area+patch_area;area_boxes+=boxes;
                            if(curve)std::cout<<std::setprecision(21)
                                <<"ACTUAL_SELECTED_CERTIFIED_AREA_PATCH index="<<p
                                <<" boxes="<<boxes
                                <<" lower="<<R(cc::sf::lower(patch_area))
                                <<" upper="<<R(cc::sf::upper(patch_area))<<std::endl;
                        }catch(const std::exception&e){
                            std::cout<<"ACTUAL_SELECTED_CERTIFIED_AREA_FAIL index="<<p
                                <<" boxes="<<boxes<<" reason="<<e.what()<<std::endl;
                            throw;
                        }
                    }
                    std::cout<<std::setprecision(21)
                        <<"ACTUAL_SELECTED_CERTIFIED_AREA_CELL cell="
                        <<candidate[single_site].cell<<" curved="<<area_curved
                        <<" boxes="<<area_boxes
                        <<" lower="<<R(cc::sf::lower(selected_area))
                        <<" upper="<<R(cc::sf::upper(selected_area))
                        <<" source="<<law.entries()[single_site].Vsource
                        <<" seconds="<<MPI_Wtime()-begun<<std::endl;
                }else{
                auto connectivity=cc::certify_global_connectivity(registry,candidate,
                    law.law_hash(),2000000);
                std::cout<<"ACTUAL_STEP5_WEIGHTED_TOPOLOGY_PASS patches="
                    <<connectivity.patch.size()<<" seconds="<<MPI_Wtime()-begun
                    <<std::endl;
                cc::JacobianBudget jacobian_budget;
                jacobian_budget.maximum_boxes_per_patch=4096;
                jacobian_budget.maximum_depth=32;
                auto jacobian=cc::certify_all_chart_jacobians(registry,candidate,
                    jacobian_budget);
                auto degree=cc::certify_degree_one(registry,candidate,
                    connectivity,jacobian,geometry);
                std::cout<<"ACTUAL_STEP5_WEIGHTED_DEGREE_PASS covered="
                    <<int(degree.degree_one_physical_coverage)
                    <<" minimum_J="<<R(jacobian.minimum_lower)
                    <<" seconds="<<MPI_Wtime()-begun<<std::endl;
                if(boundary_index_test||boundary_edge_test||
                    boundary_scale_test||boundary_failed_test||
                    boundary_parallel_prefix||boundary_parallel_range||
                    boundary_cut_range||boundary_witness_test){
                    compact_overlap_boundary::Budget boundary_budget;
                    boundary_budget.maximum_segments=1000000;
                    boundary_budget.maximum_nodes=500000;
                    boundary_budget.maximum_candidates_per_edge=8192;
                    boundary_budget.maximum_root_calls_per_edge=8192;
                    boundary_budget.maximum_events_per_edge=8192;
                    boundary_budget.leaf_size=8;
                    boundary_budget.maximum_depth=64;
                    compact_overlap_boundary::Index boundary(sites,candidate,
                        connectivity,boundary_budget);
                    std::cout<<"ACTUAL_STEP5_TARGET_BOUNDARY_INDEX_PASS segments="
                        <<boundary.segments()<<" physical="<<boundary.physical_segments()
                        <<" internal="<<boundary.internal_segments()
                        <<" nodes="<<boundary.nodes()
                        <<" trace_ties="<<boundary.trace_ties()
                        <<" mass_law_solved=0 seconds="<<MPI_Wtime()-begun
                        <<std::endl;
                    if(boundary_cut_range){
                        namespace ce=compact_overlap_events;
                        size_t first_edge=std::stoull(argv[11]);
                        size_t end_edge=std::stoull(argv[12]);
                        cc::require(first_edge<end_edge&&
                            end_edge<=source.owned.global_edges&&
                            end_edge-first_edge<=50000,
                            "finite original edge range for sealed canonical cut ledger");
                        std::filesystem::path chunk=argv[13];
                        std::filesystem::path canonical=chunk.parent_path()/
                            ("power_edge_reconcile_"+std::to_string(first_edge)+
                             "_"+std::to_string(end_edge)+"_root60_v1");
                        std::ifstream manifest(chunk/
                            ("source_edge_manifest_rank"+std::to_string(rank)+".csv"));
                        std::string header,line;
                        cc::require(bool(manifest)&&bool(std::getline(manifest,header))&&
                            bool(std::getline(manifest,line)),
                            "read sealed source-edge range manifest for native cuts");
                        auto m=cut_csv_fields(line);
                        cc::require(m.size()==18&&
                            std::stoull(m[0])==size_t(rank)&&
                            std::stoull(m[1])==first_edge&&
                            std::stoull(m[2])==end_edge&&
                            std::stoull(m[3])==source.owned.global_edges&&
                            std::stoull(m[4])==registry.wire_hash&&
                            std::stoull(m[5])==registry.source_shard_hash[rank]&&
                            std::stoull(m[6])==registry.source_epoch&&
                            std::stoull(m[9])==law.law_hash()&&
                            std::stoull(m[10])==boundary.chart_hash()&&
                            std::stoull(m[11])==boundary.weight_hash(),
                            "same source, target, mass law, weight and original edge range for canonical cuts");
                        std::ifstream csv(canonical/
                            ("canonical_owner"+std::to_string(rank)+".csv"));
                        cc::require(bool(csv)&&bool(std::getline(csv,header))&&
                            header.find("event_lower,event_upper")!=std::string::npos,
                            "read owner canonical event CSV for native cut ordering");
                        std::vector<ce::Canonical>events;
                        while(std::getline(csv,line)){
                            auto f=cut_csv_fields(line);
                            cc::require(f.size()==20&&events.size()<200000,
                                "complete finite canonical source-edge cut input row");
                            ce::Canonical c;
                            c.event.edge=std::stoull(f[0]);
                            c.event.cell[0]=std::stoull(f[1]);
                            c.event.cell[1]=std::stoull(f[2]);
                            c.event.ordinal=std::stoull(f[3]);
                            c.event.lower=cut_parse_quad(f[4]);
                            c.event.upper=cut_parse_quad(f[5]);
                            c.event.derivative_sign=std::stoi(f[6]);
                            c.target_segment=std::stoull(f[7]);
                            c.target_endpoint=std::stoi(f[8]);
                            c.source_vertex=std::stoull(f[9]);
                            c.edge_owner=std::stoi(f[10]);
                            c.incidences=std::stoi(f[11]);
                            c.source_registry_hash=registry.wire_hash;
                            c.source_epoch=registry.source_epoch;
                            c.target_chart_hash=boundary.chart_hash();
                            c.weight_hash=boundary.weight_hash();
                            cc::require(c.event.edge>=first_edge&&
                                c.event.edge<end_edge&&c.edge_owner==rank&&
                                (c.incidences==1||c.incidences==2),
                                "canonical owner event in sealed original edge range");
                            for(int k=0;k<c.incidences;k++){
                                int sr=std::stoi(f[12+k]);
                                cc::require(sr>=0&&sr<ranks,
                                    "canonical incidence from accepted source shard");
                                c.source_rank[k]=sr;
                                c.source_hash[k]=registry.source_shard_hash[sr];
                                c.incident_triangle[k]=std::stoull(f[14+k]);
                                c.local_edge[k]=std::stoi(f[16+k]);
                                c.orientation[k]=std::stoi(f[18+k]);
                            }
                            events.push_back(c);
                        }
                        cc::require(csv.eof(),
                            "complete canonical source-edge cut input");
                        compact_overlap_cuts::Budget cut_budget;
                        cut_budget.maximum_events=200000;
                        cut_budget.maximum_cuts=200000;
                        cut_budget.maximum_segments_per_cut=1024;
                        auto ordered=compact_overlap_cuts::order(boundary,events,cut_budget);
                        cc::require(ordered.complete&&
                            ordered.cut.size()+ordered.grouped_endpoint_events==events.size(),
                            "all canonical events represented once in ordered source-edge cuts");
                        std::cout<<"NATIVE_SOURCE_EDGE_CUT_ORDER_PASS rank="<<rank
                            <<" edge_range=["<<first_edge<<","<<end_edge<<")"
                            <<" events="<<events.size()
                            <<" cuts="<<ordered.cut.size()
                            <<" grouped_target_node_events="<<ordered.grouped_endpoint_events
                            <<" mass_law_solved=0 seconds="<<MPI_Wtime()-begun
                            <<std::endl;
                    }
                    if(boundary_parallel_prefix||boundary_parallel_range){
                        size_t first_edge=boundary_parallel_range?std::stoull(argv[11]):0;
                        size_t end_edge=boundary_parallel_range?std::stoull(argv[12]):10000;
                        cc::require(first_edge<end_edge&&
                            end_edge<=source.owned.global_edges&&
                            end_edge-first_edge<=50000,
                            "finite explicit checkpointable source-edge range");
                        cp::RootBudget root_budget;
                        root_budget.maximum_boxes=8192;
                        root_budget.maximum_depth=48;
                        root_budget.maximum_newton=48;
                        root_budget.parameter_width=ldexpq(1,-40);
                        compact_overlap_boundary::ProductionBudget production_budget;
                        production_budget.maximum_global_edges=4000000;
                        production_budget.maximum_local_incidences=4000000;
                        production_budget.maximum_local_unique_edges=end_edge-first_edge;
                        production_budget.maximum_total_candidate_segments=1000000;
                        production_budget.maximum_total_root_calls=1000000;
                        production_budget.maximum_total_root_boxes=50000000;
                        production_budget.maximum_local_wires=500000;
                        auto batch=cp::produce_actual_source_edge_boundary_wires_range(
                            source,geometry,registry,boundary,degree,root_budget,
                            production_budget,first_edge,end_edge);
                        std::cout<<"ACTUAL_STEP5_PARALLEL_SOURCE_EDGE_PREFIX"
                            <<" rank="<<rank
                            <<" edge_range=["<<first_edge<<","<<end_edge<<")"
                            <<" unique="<<batch.unique_edges
                            <<" candidates="<<batch.candidate_segments
                            <<" root_calls="<<batch.root_calls
                            <<" root_boxes="<<batch.root_boxes
                            <<" wires="<<batch.wire.size()
                            <<" complete_global="<<int(batch.complete)
                            <<" seconds="<<MPI_Wtime()-begun<<std::endl;
                        if(boundary_parallel_range){
                            auto path=std::filesystem::path(argv[0]).parent_path()/
                                ("source_edge_wires_rank"+std::to_string(rank)+".csv");
                            std::ofstream csv(path);
                            cc::require(bool(csv),
                                "create sealed per-rank source-edge wire range ledger");
                            csv<<"source_registry_hash,source_hash,source_epoch,target_chart_hash,weight_hash,edge_range_first,edge_range_end,event_cell_low,event_cell_high,event_edge,event_ordinal,event_lower,event_upper,event_derivative_sign,target_segment,target_endpoint,source_vertex,source_rank,incident_triangle,local_edge,orientation,edge_owner\n";
                            for(const auto&w:batch.wire)
                                csv<<w.source_registry_hash<<","<<w.source_hash<<","
                                    <<w.source_epoch<<","<<w.target_chart_hash<<","
                                    <<w.weight_hash<<","<<first_edge<<","<<end_edge<<","
                                    <<w.event.cell[0]<<","<<w.event.cell[1]<<","
                                    <<w.event.edge<<","<<w.event.ordinal<<","
                                    <<exact_quad_text(w.event.lower)<<","
                                    <<exact_quad_text(w.event.upper)<<","
                                    <<w.event.derivative_sign<<","
                                    <<w.target_segment<<","<<w.target_endpoint<<","
                                    <<w.source_vertex<<","<<w.source_rank<<","
                                    <<w.incident_triangle<<","<<w.local_edge<<","
                                    <<w.orientation<<","<<w.edge_owner<<"\n";
                            csv.close();cc::require(bool(csv),
                                "complete per-rank source-edge wire range ledger flush");
                            auto manifest_path=std::filesystem::path(argv[0]).parent_path()/
                                ("source_edge_manifest_rank"+std::to_string(rank)+".csv");
                            std::ofstream manifest(manifest_path);
                            cc::require(bool(manifest),
                                "create per-rank source-edge range manifest");
                            manifest<<"rank,edge_range_first,edge_range_end,global_edges,source_registry_hash,source_shard_hash,source_epoch,source_steps,source_time_bits,mass_law_hash,target_chart_hash,weight_hash,unique_edges,candidate_segments,root_calls,root_boxes,wire_rows,complete_global\n";
                            manifest<<rank<<","<<first_edge<<","<<end_edge<<","
                                <<source.owned.global_edges<<","<<registry.wire_hash<<","
                                <<registry.source_shard_hash[rank]<<","<<registry.source_epoch<<","
                                <<registry.source_steps<<","<<compact_target::bits(registry.source_time)<<","
                                <<law.law_hash()<<","<<boundary.chart_hash()<<","
                                <<boundary.weight_hash()<<","<<batch.unique_edges<<","
                                <<batch.candidate_segments<<","<<batch.root_calls<<","
                                <<batch.root_boxes<<","<<batch.wire.size()<<","
                                <<int(batch.complete)<<"\n";
                            manifest.close();cc::require(bool(manifest),
                                "complete per-rank source-edge range manifest flush");
                        }
                    }
                    if(boundary_witness_test){
                        struct Witness {size_t local;int edge;size_t segment;cc::Q lo,hi;};
                        const Witness witness[]={
                            {144562,2,163399,.204956Q,.205078Q},
                            {162632,0,183991,.955078Q,.956055Q},
                            {176185,0,199718,.775208Q,.775269Q},
                            {185220,2,209435,.984619Q,.984741Q},
                            {225878,2,255466,.849442Q,.849457Q},
                            {266536,2,596797,.903320Q,.903809Q}};
                        auto path=std::filesystem::path(argv[0]).parent_path()/
                            "source_edge_witness.csv";
                        std::ofstream csv(path);cc::require(bool(csv),
                            "create centered source-edge failure witness ledger");
                        csv<<std::setprecision(21)
                            <<"edge,segment,physical,u_lower,u_upper,fmid_lower,fmid_upper,dmid_lower,dmid_upper,fbox_lower,fbox_upper,dbox_lower,dbox_upper,second_lower,second_upper\n";
                        for(const auto&w:witness){
                            auto curve=cp::canonical_actual_source_edge(
                                source,geometry,w.local,w.edge);
                            const auto&segment=boundary.segment(w.segment);
                            cc::Q mid=(w.lo+w.hi)/2;
                            auto fm=cp::value(segment.plane,
                                cp::position_centered(curve,cc::S(mid)));
                            auto fb=cp::value(segment.plane,
                                cp::position_centered(curve,cc::sf::range(w.lo,w.hi)));
                            auto jets=cp::position_second(curve,cc::sf::range(w.lo,w.hi));
                            cc::S second;
                            for(int d=0;d<2;d++)
                                second=second+segment.plane.normal[d]*jets[d].dd;
                            csv<<curve.edge<<","<<w.segment<<","
                                <<int(segment.physical_boundary)<<","
                                <<R(w.lo)<<","<<R(w.hi)<<","
                                <<R(cc::sf::lower(fm.v))<<","<<R(cc::sf::upper(fm.v))<<","
                                <<R(cc::sf::lower(fm.x))<<","<<R(cc::sf::upper(fm.x))<<","
                                <<R(cc::sf::lower(fb.v))<<","<<R(cc::sf::upper(fb.v))<<","
                                <<R(cc::sf::lower(fb.x))<<","<<R(cc::sf::upper(fb.x))<<","
                                <<R(cc::sf::lower(second))<<","<<R(cc::sf::upper(second))
                                <<"\n";
                        }
                        csv.close();cc::require(bool(csv),
                            "complete centered source-edge failure witness ledger flush");
                        std::cout<<"ACTUAL_STEP5_SOURCE_EDGE_FAILURE_WITNESSES count="
                            <<sizeof(witness)/sizeof(*witness)
                            <<" seconds="<<MPI_Wtime()-begun
                            <<std::endl;
                    }
                    if(boundary_edge_test||boundary_scale_test||boundary_failed_test){
                        cp::RootBudget root_budget;
                        root_budget.maximum_boxes=8192;
                        root_budget.maximum_depth=48;
                        root_budget.maximum_newton=48;
                        root_budget.parameter_width=ldexpq(1,-40);
                        auto path=std::filesystem::path(argv[0]).parent_path()/
                            "source_edge_boundary_sample.csv";
                        std::ofstream csv(path);cc::require(bool(csv),
                            "create actual source-edge boundary diagnostic ledger");
                        csv<<"edge,local_triangle,local_edge,candidate_segments,root_calls,root_boxes,crossings,status\n";
                        const auto&o=source.owned;
                        size_t passed=0,unresolved=0;
                        std::set<uint64_t>sampled;
                        auto inspect=[&](size_t s,int k){
                            uint64_t id=o.triangle_edge[3*s+k];
                            if(!sampled.insert(id).second)return;
                            size_t candidates=0,calls=0,boxes=0,crossings=0;
                            std::string status="PASS";
                            try{
                                auto curve=cp::canonical_actual_source_edge(
                                    source,geometry,s,k);
                                auto roots=boundary.intersections(curve,root_budget);
                                cc::require(roots.complete,
                                    "complete selected actual source-edge intersection");
                                candidates=roots.candidate_segments;
                                calls=roots.root_calls;boxes=roots.root_boxes;
                                crossings=roots.crossing.size();passed++;
                            }catch(const std::exception&e){
                                status=e.what();unresolved++;
                            }
                            csv<<id<<","<<s<<","<<k<<","<<candidates<<","
                                <<calls<<","<<boxes<<","<<crossings<<","
                                <<std::quoted(status)<<"\n";
                        };
                        if(boundary_failed_test){
                            const std::pair<size_t,int>six[]={
                                {33106,2},{59814,2},{74281,0},
                                {169150,2},{181113,0},{183617,0}};
                            for(const auto&item:six)inspect(item.first,item.second);
                        }else{
                            for(size_t s=0;s<8;s++)for(int k=0;k<3;k++)inspect(s,k);
                            const size_t spread=boundary_scale_test?1024:64;
                            for(size_t j=0;j<spread;j++){
                                size_t s=j*(o.f.triangle.size()-1)/(spread-1);
                                inspect(s,int(j%3));
                            }
                        }
                        csv.close();cc::require(bool(csv),
                            "complete selected actual source-edge boundary ledger flush");
                        std::cout<<"ACTUAL_STEP5_SOURCE_EDGE_BOUNDARY_SAMPLE passed="
                            <<passed<<" unresolved="<<unresolved
                            <<" sampled="<<sampled.size()
                            <<" seconds="<<MPI_Wtime()-begun<<std::endl;
                    }
                }else if(whole_patch_test){
                    compact_source_spatial::Budget source_budget;
                    source_budget.maximum_entries=1200000;
                    source_budget.maximum_nodes=500000;
                    source_budget.maximum_results_per_query=8192;
                    source_budget.leaf_size=8;
                    source_budget.maximum_depth=64;
                    source_budget.ranks=4;
                    compact_source_spatial::WholePatchBudget containment_budget;
                    containment_budget.maximum_patches=500000;
                    containment_budget.maximum_total_candidates_per_rank=100000000;
                    auto containment=compact_source_spatial::certify_whole_patch_source_containment(
                        source,geometry,registry,source_budget,candidate,connectivity,
                        degree,containment_budget);
                    if(rank==0){
                        auto path=std::filesystem::path(argv[0]).parent_path()/
                            "whole_patch_containment.csv";
                        std::ofstream csv(path);cc::require(bool(csv),
                            "create actual whole-patch source containment ledger");
                        csv<<"flat_patch,target_cell,source_triangle,source_owner,certified\n";
                        for(size_t p=0;p<containment.patch.size();p++){
                            const auto&entry=containment.patch[p];
                            csv<<p<<","<<connectivity.patch[p].cell<<","
                                <<entry.triangle<<","<<entry.owner<<","
                                <<int(entry.certified)<<"\n";
                        }
                        csv.close();cc::require(bool(csv),
                            "complete actual whole-patch containment ledger flush");
                        std::cout<<"ACTUAL_STEP5_WHOLE_PATCH_CONTAINMENT_RESULT certified="
                            <<containment.certified<<" unresolved="<<containment.unresolved
                            <<" patches="<<containment.patch.size()
                            <<" candidates_rank0="<<containment.candidates
                            <<" mass_law_solved=0 seconds="<<MPI_Wtime()-begun
                            <<std::endl;
                    }
                }else if(certified_full){
                    struct CellAreaResult {cc::S area;size_t boxes=0,curved=0;R relative_bound=0;std::string reason;};
                    std::vector<CellAreaResult> result(candidate.size());
                    std::atomic<size_t> certified_progress(0);
                    #pragma omp parallel for schedule(dynamic,16)
                    for(size_t i=0;i<candidate.size();i++){
                        try{
                            cc::AreaBudget budget;
                            budget.maximum_boundary_boxes_per_patch=65536;
                            budget.maximum_absolute_area_error_per_patch=1e-11Q;
                            for(const auto&patch:candidate[i].patch){
                                cp::Curve descriptor;const cp::Curve*curve=nullptr;
                                if(patch.edge[1].curved){
                                    descriptor=cp::curve(registry.exterior.at(patch.edge[1].source));
                                    curve=&descriptor;result[i].curved++;
                                }
                                size_t boxes=0,subdivisions=0;
                                result[i].area=result[i].area+cc::boundary_area(
                                    patch,curve,budget,boxes,subdivisions);
                                result[i].boxes+=boxes;
                            }
                            R source_area=law.entries()[i].Vsource;
                            R lower=R(cc::sf::lower(result[i].area));
                            R upper=R(cc::sf::upper(result[i].area));
                            cc::require(lower>0&&lower<=upper&&source_area>0,
                                "positive certified target and source cell area");
                            result[i].relative_bound=std::max(
                                fabsl(lower-source_area),fabsl(upper-source_area))/source_area;
                            if(result[i].relative_bound>1e-8L)
                                result[i].reason="certified area residual exceeds 1e-8";
                        }catch(const std::exception&e){result[i].reason=e.what();}
                        size_t done=certified_progress.fetch_add(1)+1;
                        if(done%10000==0)std::cout
                            <<"ACTUAL_STEP5_CERTIFIED_AREA_PROGRESS completed="<<done
                            <<" seconds="<<MPI_Wtime()-begun<<std::endl;
                    }
                    auto path=std::filesystem::path(argv[0]).parent_path()/
                        "certified_area_residual.csv";
                    std::ofstream csv(path);cc::require(bool(csv),
                        "create complete certified area residual ledger");
                    csv<<std::setprecision(21)
                        <<"cell,source_area,target_lower,target_upper,relative_bound,curved,boxes,status\n";
                    size_t certified=0,unresolved=0,total_boxes=0;
                    R maximum_relative_bound=0;size_t worst=SIZE_MAX;
                    for(size_t i=0;i<result.size();i++){
                        const auto&r=result[i];total_boxes+=r.boxes;
                        R source_area=law.entries()[i].Vsource;
                        if(r.reason.empty())certified++;else unresolved++;
                        if(r.relative_bound>maximum_relative_bound){
                            maximum_relative_bound=r.relative_bound;worst=i;
                        }
                        csv<<registry.site[i].cell<<","<<source_area<<","
                            <<R(cc::sf::lower(r.area))<<","<<R(cc::sf::upper(r.area))
                            <<","<<r.relative_bound<<","<<r.curved<<","<<r.boxes
                            <<","<<(r.reason.empty()?"PASS":r.reason)<<"\n";
                    }
                    csv.close();cc::require(bool(csv),"complete certified area ledger flush");
                    std::cout<<std::setprecision(21)
                        <<"ACTUAL_STEP5_CERTIFIED_AREA_RESULT certified="<<certified
                        <<" unresolved="<<unresolved<<" total="<<result.size()
                        <<" boxes="<<total_boxes<<" maximum_relative_bound="
                        <<maximum_relative_bound<<" worst_index="<<worst
                        <<" weight_law_solved=0 seconds="<<MPI_Wtime()-begun
                        <<std::endl;
                }else{
                // Paired-order Green quadrature gives only a diagnostic
                // geometric area residual; its integration error is NOT certified.
                std::vector<double>node8,weight8,node16,weight16;
                shared_trace::gauss(8,node8,weight8);
                shared_trace::gauss(16,node16,weight16);
                std::vector<cc::S> area8(candidate.size()),area16(candidate.size());
                std::vector<size_t> curved_per_cell(candidate.size());
                std::atomic<bool> area_failed(false);
                std::atomic<size_t> area_progress(0);
                size_t first_area_failure=SIZE_MAX;std::string area_reason;
                #pragma omp parallel for schedule(dynamic,64)
                for(size_t i=0;i<candidate.size();i++){
                    if(area_failed.load(std::memory_order_relaxed))continue;
                    try{
                        cc::S cell_area,cell_area8;size_t curved=0;
                        for(const auto&patch:candidate[i].patch){
                            const cp::Curve*curve=nullptr;cp::Curve descriptor;
                            if(patch.edge[1].curved){
                                descriptor=cp::curve(registry.exterior.at(patch.edge[1].source));
                                curve=&descriptor;curved++;
                            }
                            cell_area8=cell_area8+diagnostic_green_gauss(
                                patch,curve,node8,weight8);
                            cell_area=cell_area+diagnostic_green_gauss(
                                patch,curve,node16,weight16);
                        }
                        cc::require(cc::sf::lower(cell_area)>0,
                            "positive complete target-cell Green area diagnostic");
                        area8[i]=cell_area8;area16[i]=cell_area;
                        curved_per_cell[i]=curved;
                        size_t done=area_progress.fetch_add(1)+1;
                        if(done%10000==0)std::cout
                            <<"ACTUAL_STEP5_WEIGHTED_AREA_PROGRESS completed="<<done
                            <<" seconds="<<MPI_Wtime()-begun<<std::endl;
                    }catch(const std::exception&e){
                        area_failed.store(true,std::memory_order_relaxed);
                        #pragma omp critical
                        {if(i<first_area_failure){first_area_failure=i;area_reason=e.what();}}
                    }
                }
                if(area_failed.load())throw std::runtime_error(
                    "Green area index="+std::to_string(first_area_failure)+
                    " cell="+std::to_string(registry.site[first_area_failure].cell)+
                    " "+area_reason);
                size_t curved_area_patches=0;
                R maximum_relative_residual=0,maximum_order_disagreement=0;
                size_t worst_cell=SIZE_MAX;R worst_signed_residual=0;
                cc::S total_target_area;
                auto csv_path=std::filesystem::path(argv[0]).parent_path()/
                    "gauss_area_residual.csv";
                std::ofstream csv(csv_path);
                cc::require(bool(csv),"create complete per-cell diagnostic area ledger");
                csv<<std::setprecision(21)
                    <<"cell,source_area,target_gauss8,target_gauss16,relative_residual,curved_patches\n";
                for(size_t i=0;i<candidate.size();i++){
                    curved_area_patches+=curved_per_cell[i];
                    total_target_area=total_target_area+area16[i];
                    R source_area=law.entries()[i].Vsource;
                    R lo=R(cc::sf::lower(area16[i])),hi=R(cc::sf::upper(area16[i]));
                    maximum_order_disagreement=std::max(maximum_order_disagreement,
                        fabsl(R(area16[i].n-area8[i].n))/source_area);
                    R signed_residual=((lo+hi)/2-source_area)/source_area;
                    csv<<registry.site[i].cell<<","<<source_area<<","
                        <<R(area8[i].n)<<","<<R(area16[i].n)<<","
                        <<signed_residual<<","<<curved_per_cell[i]<<"\n";
                    if(fabsl(signed_residual)>maximum_relative_residual){
                        maximum_relative_residual=fabsl(signed_residual);
                        worst_signed_residual=signed_residual;worst_cell=i;
                    }
                }
                csv.close();cc::require(bool(csv),"complete per-cell diagnostic area ledger flush");
                std::cout<<std::setprecision(21)
                    <<"ACTUAL_STEP5_WEIGHTED_GAUSS_AREA_DIAGNOSTIC cells="<<candidate.size()
                    <<" curved_patches="<<curved_area_patches
                    <<" target_gauss16_sum=["<<R(cc::sf::lower(total_target_area))
                    <<","<<R(cc::sf::upper(total_target_area))<<"]"
                    <<" worst_cell="<<worst_cell
                    <<" worst_signed_relative_residual="<<worst_signed_residual
                    <<" maximum_relative_gauss8_16_difference="<<maximum_order_disagreement
                    <<" integration_error_certified=0"
                    <<" weight_law_solved=0"
                    <<" seconds="<<MPI_Wtime()-begun<<std::endl;
                // Linearized power-area geometry for a subsequent mass-law
                // solve. Every internal edge is counted from both cells and
                // must have agreeing interval lengths before serialization.
                struct Interface {cc::S side[2];size_t pieces[2]={};};
                std::map<std::pair<uint64_t,uint64_t>,Interface>interface;
                for(size_t i=0;i<candidate.size();i++)for(const auto&patch:candidate[i].patch)
                    for(int k=0;k<3;k++){
                        const auto&e=patch.edge[k];
                        if(!e.boundary||e.curved||e.line.kind!=1)continue;
                        cc::require(candidate[i].cell==e.line.a||candidate[i].cell==e.line.b,
                            "internal power boundary belongs to one canonical side");
                        int side=candidate[i].cell==e.line.a?0:1;
                        const auto&A=patch.vertex[k].x,&B=patch.vertex[(k+1)%3].x;
                        cc::S dx=B[0]-A[0],dy=B[1]-A[1];
                        cc::S length=cc::sf::mf::root(dx*dx+dy*dy);
                        cc::require(cc::sf::lower(length)>0,
                            "strictly positive actual internal power-edge length");
                        auto&entry=interface[{e.line.a,e.line.b}];
                        entry.side[side]=entry.side[side]+length;entry.pieces[side]++;
                    }
                auto adjacency_path=std::filesystem::path(argv[0]).parent_path()/
                    "power_area_jacobian_edges.csv";
                std::ofstream adjacency(adjacency_path);
                cc::require(bool(adjacency),"create complete power-area Jacobian edge ledger");
                adjacency<<std::setprecision(21)
                    <<"cell_low,cell_high,boundary_length,site_distance,conductance,pieces_low,pieces_high\n";
                R max_interface_relative_disagreement=0;
                for(const auto&pair:interface){
                    uint64_t a=pair.first.first,b=pair.first.second;
                    const auto&entry=pair.second;
                    cc::require(entry.pieces[0]>0&&entry.pieces[1]>0,
                        "both actual cells contribute each internal power edge");
                    cc::Q lo=std::max(cc::sf::lower(entry.side[0]),
                        cc::sf::lower(entry.side[1]));
                    cc::Q hi=std::min(cc::sf::upper(entry.side[0]),
                        cc::sf::upper(entry.side[1]));
                    cc::require(lo>0&&lo<=hi,
                        "two actual cells enclose the same internal power-edge length");
                    R length=(R(lo)+R(hi))/2;
                    max_interface_relative_disagreement=std::max(
                        max_interface_relative_disagreement,
                        fabsl(R(entry.side[0].n-entry.side[1].n))/length);
                    auto ia=std::lower_bound(registry.site.begin(),registry.site.end(),a,
                        [](const compact_source::Site&s,uint64_t id){return s.cell<id;});
                    auto ib=std::lower_bound(registry.site.begin(),registry.site.end(),b,
                        [](const compact_source::Site&s,uint64_t id){return s.cell<id;});
                    cc::require(ia!=registry.site.end()&&ib!=registry.site.end()&&
                        ia->cell==a&&ib->cell==b,"power-edge generator identities exist");
                    cc::S separation;
                    for(int d=0;d<2;d++){
                        cc::S q=compact_source::canonical_position(*ib,d)-
                            compact_source::canonical_position(*ia,d);
                        separation=separation+q*q;
                    }
                    cc::S distance=cc::sf::mf::root(separation);
                    cc::require(cc::sf::lower(distance)>0,
                        "positive actual power-neighbor site separation");
                    R d=(R(cc::sf::lower(distance))+R(cc::sf::upper(distance)))/2;
                    adjacency<<a<<","<<b<<","<<length<<","<<d<<","
                        <<length/(2*d)<<","<<entry.pieces[0]<<","
                        <<entry.pieces[1]<<"\n";
                }
                adjacency.close();cc::require(bool(adjacency),
                    "complete actual power-area Jacobian edge ledger flush");
                std::cout<<"ACTUAL_STEP5_WEIGHTED_LINEARIZED_POWER_AREA_EDGES"
                    <<" pairs="<<interface.size()
                    <<" max_relative_interface_length_disagreement="
                    <<max_interface_relative_disagreement
                    <<" integration_error_certified=0"
                    <<" seconds="<<MPI_Wtime()-begun<<std::endl;
                }
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
        MPI_Finalize();return 0;
    }catch(const std::exception&e){
        std::cerr<<"rank"<<rank<<" ACTUAL_STEP5_WEIGHTED_FAIL "<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD,2);return 2;
    }
}
