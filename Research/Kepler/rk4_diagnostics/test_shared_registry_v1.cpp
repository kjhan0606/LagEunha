#define GFS_SHARED_FUNCTION_LIBRARY
#include "test_shared_function_v1.cpp"
#include "material_shared_retained.h"
#include <tuple>

namespace target=shared_target;
namespace {
// These are registry slots, NOT newly evaluated overlap integrals.
struct Slot {
    uint64_t id,source,cell,target;
    uint32_t target_index,owner;
    int32_t coarse,resolution; // 0 excluded,1 pending,2 retained empty,3 whole-group member,4 clipped
    int32_t component,reserved;
    uint64_t points,panels;
};
static_assert(sizeof(Slot)==72,"pointer-free registry slot layout");
struct EdgeEntry {const STEdge* function=nullptr;std::vector<std::array<int,3>> side;};
struct Known {int classification=-1,group=0;uint64_t points=0,panels=0;};
struct RootWork {int edge;uint64_t first,second;bool retained=false;};
static FILE* output(const std::string&dest,const char*name,const char*header){
    FILE*f=fopen((dest+"/"+name).c_str(),"w");need(f,"new registry output");fprintf(f,"%s\n",header);return f;
}
static int classify(const S*image,const std::vector<target::Plane>&planes,int&outside,Q&margin){
    bool whole=true;outside=-1;margin=-__builtin_huge_valq();
    for(size_t j=0;j<planes.size();j++){
        auto value=target::value(planes[j],image[0],image[1]);
        if(sf::lower(value)>0){outside=j;margin=sf::lower(value);return 0;}
        whole=whole&&sf::upper(value)<0;margin=std::max(margin,sf::upper(value));
    }
    return whole?1:2;
}
}

int main(int argc,char**argv){
    int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);int code=2;
    try{
        need(argc==11,"packet saved sites weights originalcall bounds published whole4 roots NEWoutput");
        auto begin=std::chrono::steady_clock::now();auto p=curved::read(argv[1]);
        const auto f=load_saved(std::string(argv[2])+"/shared_state.bin",p);auto coefficients=f.coefficient;
        auto site=sites(argv[3],argv[4]);std::sort(site.begin(),site.end(),[](const Site&a,const Site&b){return a.id<b.id;});
        const size_t ns=f.triangle.size(),nt=site.size();need(ns==384&&nt==24&&f.cell.size()==24,"original sparse 384-source patch and 24 target sites");
        auto retained=target::retained_sources(argv[5],argv[6],f);std::string published=argv[7],whole=argv[8],roots=argv[9],dest=argv[10];
        std::ifstream terminal(published+"/native_exit_code.txt");int exit=-1;need(bool(terminal>>exit)&&exit==0,"retained target publication terminal success");
        FILE*packet=fopen((published+"/caller_target_packet.bin").c_str(),"rb");need(packet,"retained caller packet");uint64_t header[2];need(fread(header,sizeof(header),1,packet)==1&&header[0]==1&&header[1]==68526,"actual epoch1 target68526 publication");fclose(packet);
        std::vector<Known>known(ns);std::ifstream completion(published+"/all384_source_completion.csv");need(bool(completion),"published all384 completion");std::string line;std::getline(completion,line);std::set<int>seen;
        while(std::getline(completion,line)){
            auto c=target::columns(line);int s=std::stoi(c[0]);need(s>=0&&size_t(s)<ns&&seen.insert(s).second&&c[3]=="1"&&c[4]=="1"&&c[5]=="1"&&c[6]=="1","retained all384 coverage and integration complete");
            need(std::stoull(c[1])==f.cell[f.triangle[s].cell],"retained original physical source");
            known[s].classification=std::stoi(c[2]);known[s].group=std::stoi(c[7]);known[s].panels=std::stoull(c[8]);known[s].points=std::stoull(c[10]);
        }
        need(seen.size()==ns,"all known target source slots");
        std::ifstream wholes(whole+"/whole_inside_source_integrals.csv");need(bool(wholes),"retained group members");std::getline(wholes,line);std::set<int>group;
        while(std::getline(wholes,line)){auto c=target::columns(line);int s=std::stoi(c[0]);need(known[s].classification==1&&known[s].group==1&&c[1]=="1"&&group.insert(s).second,"unique retained group member");known[s].panels=std::stoull(c[2]);known[s].points=std::stoull(c[4]);}
        need(group==std::set<int>({240,244,248,252}),"four-source group is retained ONCE, not four copied integrals");

        std::vector<std::vector<target::Plane>>planes(nt);
        std::vector<S>image(2*ns);std::vector<Slot>slot(ns*nt);
        for(size_t s=0;s<ns;s++)for(int k=0;k<2;k++){
            Q lo=Q(f.triangle[s].vertex[0][k]),hi=lo;for(int a=1;a<3;a++){lo=std::min(lo,Q(f.triangle[s].vertex[a][k]));hi=std::max(hi,Q(f.triangle[s].vertex[a][k]));}
            image[2*s+k]=sf::range(lo,hi);image[2*s+k].e=sf::add(image[2*s+k].e,retained[s].real.displacement[k]);
        }
        FILE*targets=rank==0?output(dest,"target_registry.csv","target_index,target,x,y,weight,bounded,strict_site_interior,partition_tie_rule,physical_domain,normalizer_state,all12_tau_available"):nullptr;
        FILE*pairplanes=rank==0?output(dest,"power_plane_registry.csv","first_target,second_target,nx,ny,b,nx_radius,ny_radius,b_radius,first_owns_equality"):nullptr;
        for(size_t t=0;t<nt;t++)planes[t]=target::planes(site,site[t].id);
        for(size_t t=0;t<nt;t++){
            bool inside=true;
            for(const auto&plane:planes[t])inside=inside&&target::sign(target::exact_value(plane,site[t].p[0].n,site[t].p[1].n))<0;
            need(inside,"current W0 generators have strict nonempty target cells");
            if(rank==0){fprintf(targets,"%zu,%llu,%s,%s,%s,%d,1,smallest_target_id_owns_exact_equality,source_patch_intersection_only_no_fill,%s,%d\n",t,(unsigned long long)site[t].id,qs(site[t].p[0].n).c_str(),qs(site[t].p[1].n).c_str(),qs(site[t].w.n).c_str(),target::bounded_target_cell(planes[t]),site[t].id==68526?"retained_original_physical_tau":"pending_positive_physical_lower_bounds",site[t].id==68526);
                for(size_t j=t+1;j<nt;j++){
                    auto it=std::find_if(planes[t].begin(),planes[t].end(),[&](const target::Plane&p){return p.competitor==site[j].id;});need(it!=planes[t].end(),"canonical pair plane");
                    auto other=std::find_if(planes[j].begin(),planes[j].end(),[&](const target::Plane&p){return p.competitor==site[t].id;});
                    need(other!=planes[j].end()&&target::sign(it->eb+other->eb)==0,"opposite exact pair offsets");
                    for(int k=0;k<2;k++)need(target::sign(it->en[k]+other->en[k])==0,"opposite exact pair normals");
                    fprintf(pairplanes,"%llu,%llu,%s,%s,%s,%s,%s,%s,1\n",(unsigned long long)site[t].id,(unsigned long long)site[j].id,qs(it->n[0].n).c_str(),qs(it->n[1].n).c_str(),qs(it->b.n).c_str(),qs(it->n[0].e).c_str(),qs(it->n[1].e).c_str(),qs(it->b.e).c_str());
                }
            }
        }
        if(rank==0){fclose(targets);fclose(pairplanes);}
        #pragma omp parallel for schedule(static)
        for(size_t id=0;id<slot.size();id++)if(id%ranks==size_t(rank)){
            size_t s=id/nt,t=id%nt;int outside;Q margin;
            auto&r=slot[id];r={id,s,f.cell[f.triangle[s].cell],site[t].id,uint32_t(t),uint32_t(id%ranks),0,0,p.record[s].component,0,0,0};
            r.coarse=classify(image.data()+2*s,planes[t],outside,margin);r.resolution=r.coarse?1:0;
            if(site[t].id==68526){const auto&k=known[s];r.resolution=k.classification==0?0:k.classification==1?3:k.points?4:2;r.points=k.points;r.panels=k.panels;}
        }
        std::vector<Slot>all(rank==0?slot.size()*ranks:0);MPI_Gather(slot.data(),slot.size()*sizeof(Slot),MPI_BYTE,all.data(),slot.size()*sizeof(Slot),MPI_BYTE,0,MPI_COMM_WORLD);
        if(rank==0){
            for(size_t id=0;id<slot.size();id++)slot[id]=all[(id%ranks)*slot.size()+id];
            FILE*images=output(dest,"source_image_registry.csv","source,cell,x,x_radius,y,y_radius,L,displacement_x,displacement_y");
            for(size_t s=0;s<ns;s++)fprintf(images,"%zu,%llu,%s,%s,%s,%s,%s,%s,%s\n",s,(unsigned long long)f.cell[f.triangle[s].cell],qs(image[2*s].n).c_str(),qs(image[2*s].e).c_str(),qs(image[2*s+1].n).c_str(),qs(image[2*s+1].e).c_str(),qs(retained[s].real.L).c_str(),qs(retained[s].real.displacement[0]).c_str(),qs(retained[s].real.displacement[1]).c_str());fclose(images);
            FILE*pairs=output(dest,"source_target_registry.csv","slot_id,source,cell,target,target_index,owner_rank,coarse_classification,resolution,component,retained_points,retained_panels,actual_overlap_integral_available,value_binding");
            FILE*outgoing=output(dest,"source_outgoing_incidence.csv","slot_id,source,target,coefficient,slot_byte_offset,numerical_value_available");
            FILE*incoming=output(dest,"target_incoming_incidence.csv","slot_id,source,target,coefficient,slot_byte_offset,numerical_value_available");
            uint64_t counts[5]={},active=0,coarse_counts[3]={},known_points=0,pending=0;
            for(size_t s=0;s<ns;s++){
                int wholes=0,candidates=0;for(size_t t=0;t<nt;t++){wholes+=slot[s*nt+t].coarse==1;candidates+=slot[s*nt+t].coarse!=0;}need(wholes<=1&&candidates>0,"power partition coarse source ownership not contradictory");
            }
            for(const auto&r:slot){
                counts[r.resolution]++;coarse_counts[r.coarse]++;bool numerical=r.resolution==4;const char*binding=r.resolution==0?"certified_zero":r.resolution==2?"retained_certified_empty":r.resolution==3?"retained_four_source_group_0_not_individual_value":r.resolution==4?"retained_per_source_clipped_all12":"pending_actual_atlas_and_quadrature";
                fprintf(pairs,"%llu,%llu,%llu,%llu,%u,%u,%d,%d,%d,%llu,%llu,%d,%s\n",(unsigned long long)r.id,(unsigned long long)r.source,(unsigned long long)r.cell,(unsigned long long)r.target,r.target_index,r.owner,r.coarse,r.resolution,r.component,(unsigned long long)r.points,(unsigned long long)r.panels,numerical,binding);
                if(r.resolution==1||r.resolution==3||r.resolution==4){active++;fprintf(outgoing,"%llu,%llu,%llu,-1,%llu,%d\n",(unsigned long long)r.id,(unsigned long long)r.source,(unsigned long long)r.target,(unsigned long long)(64+r.id*sizeof(Slot)),numerical);fprintf(incoming,"%llu,%llu,%llu,1,%llu,%d\n",(unsigned long long)r.id,(unsigned long long)r.source,(unsigned long long)r.target,(unsigned long long)(64+r.id*sizeof(Slot)),numerical);}
                if(r.resolution==1)pending++;if(r.target==68526)known_points+=r.points;
            }
            fclose(pairs);fclose(outgoing);fclose(incoming);need(known_points==268032&&counts[2]==10&&counts[3]==4&&counts[4]==26,"known target retained without numerical replay");
            FILE*binary=fopen((dest+"/pair_registry.bin").c_str(),"wb");uint64_t h[8]={0x4753465245473031ULL,1,1,sizeof(Slot),slot.size(),nt,ns,uint64_t(ranks)};need(binary&&fwrite(h,sizeof(h),1,binary)==1&&fwrite(slot.data(),sizeof(Slot),slot.size(),binary)==slot.size(),"pointer-free source-target registry");fclose(binary);
            FILE*groups=output(dest,"retained_reduction_groups.csv","group_id,target,members,positive_points,arithmetic_value_instances,individual_source_closure_proved,policy");fprintf(groups,"0,68526,240;244;248;252,45056,1,0,reuse_group_once_individual_whole_source_numerical_closure_pending\n");fclose(groups);

            std::map<int,EdgeEntry>edges;
            for(size_t s=0;s<ns;s++)for(int e=0;e<3;e++){
                const auto&t=f.triangle[s];const auto&function=t.edge[e];auto&edge=edges[function.poly[1]];
                if(edge.function)need(!memcmp(edge.function,&function,sizeof(STEdge)),"identical canonical shared edge function");else edge.function=&function;
                int a=t.endpoint[e][0],b=t.endpoint[e][1];need((a+1)%3==b||(b+1)%3==a,"oriented triangle edge");edge.side.push_back({int(s),e,(a+1)%3==b?1:-1});
            }
            std::set<std::tuple<int,uint64_t,uint64_t>>old_roots;
            std::ifstream input_roots(roots+"/canonical_event_root_packets.csv");need(bool(input_roots),"known complete canonical event records");std::getline(input_roots,line);
            while(std::getline(input_roots,line)){auto c=target::columns(line);need(c[2]=="1","retained canonical root complete");uint64_t other=std::stoull(c[1]);old_roots.emplace(std::stoi(c[0]),std::min(uint64_t(68526),other),std::max(uint64_t(68526),other));}
            std::map<std::tuple<int,uint64_t,uint64_t>,RootWork>jobs;
            FILE*edgecsv=output(dest,"canonical_edge_registry.csv","edge,first_vertex,last_vertex,first_source,first_local,first_orientation,second_source,second_local,second_orientation,physical_external_boundary,interval_lo,interval_hi,fragment_state,policy");
            FILE*external=output(dest,"external_edge_target_candidates.csv","edge,source,target,slot_id,fragment_state");
            uint64_t external_edges=0,external_candidates=0,paired=0;
            for(const auto&entry:edges){
                const auto&e=entry.second;need(e.side.size()==1||e.side.size()==2,"one boundary or two interior edge owners");bool boundary=e.side.size()==1;
                if(!boundary)need(e.side[0][2]+e.side[1][2]==0,"paired source edges have opposite orientations");
                external_edges+=boundary;paired+=!boundary;
                fprintf(edgecsv,"%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,0,1,%s,original_patch_boundary_no_fill_no_periodic_invention\n",entry.first,e.function->poly[0],e.function->poly[3],e.side[0][0],e.side[0][1],e.side[0][2],boundary?-1:e.side[1][0],boundary?-1:e.side[1][1],boundary?0:e.side[1][2],boundary,boundary?"unsubdivided_candidate_only":"paired_same_function");
                for(size_t t=0;t<nt;t++){
                    bool candidate=false;for(const auto&side:e.side){const auto&r=slot[side[0]*nt+t];candidate=candidate||(r.resolution==1||r.resolution==3||r.resolution==4);}
                    if(boundary&&candidate){external_candidates++;fprintf(external,"%d,%d,%llu,%llu,pending_actual_edge_partition\n",entry.first,e.side[0][0],(unsigned long long)site[t].id,(unsigned long long)(e.side[0][0]*nt+t));}
                    if(!candidate)continue;
                    for(const auto&plane:planes[t]){
                        bool possible=false;for(const auto&side:e.side){auto value=target::value(plane,image[2*side[0]],image[2*side[0]+1]);possible=possible||!(sf::lower(value)>0||sf::upper(value)<0);}
                        if(!possible)continue;
                        uint64_t a=std::min(site[t].id,plane.competitor),b=std::max(site[t].id,plane.competitor);auto key=std::make_tuple(entry.first,a,b);
                        jobs.emplace(key,RootWork{entry.first,a,b,old_roots.count(key)>0});
                    }
                }
            }
            fclose(edgecsv);fclose(external);need(edges.size()==616&&2*paired+external_edges==3*ns,"full original edge incidence");
            FILE*work=output(dest,"canonical_root_work_registry.csv","work_id,edge,first_target,second_target,retained_complete_roots,new_root_evaluation_cap,reason");uint64_t reused=0,fresh=0;
            for(const auto&entry:jobs){const auto&r=entry.second;fprintf(work,"%llu,%d,%llu,%llu,%d,%d,%s\n",(unsigned long long)(reused+fresh),r.edge,(unsigned long long)r.first,(unsigned long long)r.second,r.retained,r.retained?0:256,r.retained?"same_unordered_exact_plane_and_canonical_edge":"pending_fixed256_root_budget_no_evaluation_yet");reused+=r.retained;fresh+=!r.retained;}fclose(work);
            FILE*cost=output(dest,"cost_preflight.csv","item,count,interpretation");
            fprintf(cost,"source_target_registry_slots,%zu,exact\nactive_nonzero_or_pending_slots,%llu,unique_incoming_outgoing_identity\npending_actual_atlas_slots,%llu,not_proved_nonempty\nnew_deduplicated_root_jobs,%llu,no_roots_evaluated\nnew_root_window_evaluation_cap,%llu,fixed256_per_job_only_not_wall_time_bound\nretained_root_jobs,%llu,no_numeric_replay\nretained_target_positive_queries,%llu,no_field_replay\nregistry_binary_bytes,%zu,includes64byte_header\nall24_final_caller_packet_bytes,9600,24times400_only_not_quadrature_storage\nnew_panel_count,-1,UNKNOWN_until_actual_atlas\nnew_positive_query_count,-1,UNKNOWN_until_actual_atlas\nnew_positive_query_bytes,-1,UNKNOWN_until_actual_atlas\nfull_quadrature_launch_allowed,0,cost_and_actual_coverage_not_yet_closed\n",slot.size(),(unsigned long long)active,(unsigned long long)pending,(unsigned long long)fresh,(unsigned long long)(fresh*256),(unsigned long long)reused,(unsigned long long)known_points,64+sizeof(Slot)*slot.size());fclose(cost);
            checkpoint(dest+"/caller_unchanged_state.bin",f,1);
            printf("{\"scope\":\"all24-target registry/incidence/cost PREFLIGHT only on original sparse shared patch; no new roots or quadrature or evolution\",\"MPI_ranks\":%d,\"OMP_threads\":%d,\"sources\":384,\"targets\":24,\"pair_slots\":%zu,\"coarse_excluded\":%llu,\"coarse_whole\":%llu,\"coarse_candidates\":%llu,\"pending_actual_atlas_slots\":%llu,\"unique_active_incidence_slots\":%llu,\"canonical_edges\":%zu,\"paired_edges\":%llu,\"external_edges\":%llu,\"external_edge_candidate_incidences\":%llu,\"new_deduplicated_root_jobs\":%llu,\"retained_root_jobs\":%llu,\"target68526_reused_positive_queries\":%llu,\"retained_whole_group_value_count\":1,\"registry_incidence_cost_preflight_complete\":true,\"actual_all_target_overlap_geometry_complete\":false,\"source_partition_boundary_cancellation_complete\":false,\"global_source_image_disjointness_complete\":false,\"all_source_numerical_conservation_complete\":false,\"external_actual_fragments_complete\":false,\"other23_physical_normalizers_complete\":false,\"full_quadrature_launch_allowed\":false,\"new_field_queries\":0,\"accepted_RK4_evolution\":false,\"long_Kepler_validated\":false,\"seconds\":%.17g}\n",ranks,omp_get_max_threads(),slot.size(),(unsigned long long)coarse_counts[0],(unsigned long long)coarse_counts[1],(unsigned long long)coarse_counts[2],(unsigned long long)pending,(unsigned long long)active,edges.size(),(unsigned long long)paired,(unsigned long long)external_edges,(unsigned long long)external_candidates,(unsigned long long)fresh,(unsigned long long)reused,(unsigned long long)known_points,std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count());
            code=0;
        }
        need(!memcmp(coefficients.data(),f.coefficient.data(),coefficients.size()*sizeof(STCoefficient)),"registry never changes native coefficients");MPI_Bcast(&code,1,MPI_INT,0,MPI_COMM_WORLD);
    }catch(const std::exception&e){fprintf(stderr,"rank%d FAIL %s\n",rank,e.what());MPI_Abort(MPI_COMM_WORLD,2);}MPI_Finalize();return code;
}
