#define GFS_SHARED_FUNCTION_LIBRARY
#include "test_shared_function_v1.cpp"
#include "material_shared_retained.h"
#include <tuple>
namespace target=shared_target;namespace sd=shared_duffy;
namespace {
using Key=std::tuple<int,uint64_t,uint64_t>;
struct Slot {int coarse=0,resolution=0;};
struct Query {S x,y;target::Membership membership;};
struct Result {int complete=0,empty=0,reused=0;uint64_t windows=0;target::Domain domain;std::vector<Query>query;std::vector<std::string>event;std::string reason;};
static std::string quote(const std::string&x){std::string out="\"";for(char c:x){if(c=='\"')out+='\"';out+=c;}return out+'\"';}
static std::string events(const target::Node&node,const std::vector<target::Plane>&planes,uint64_t own){std::set<uint64_t>id={own};for(int j:node.active){need(j>=0&&size_t(j)<planes.size(),"active plane identity");id.insert(planes[j].competitor);}std::string s="p";for(auto i:id)s+=';'+std::to_string(i);return s;}
static std::vector<std::string> node_events(const target::Domain&d,const std::map<int,target::Edge>&edges,const std::vector<target::Plane>&planes,uint64_t own){
    std::vector<std::string>out(d.node.size());std::vector<bool>used(d.node.size(),false);for(const auto&b:d.boundary){used[b.first]=true;used[b.last]=true;}
    for(size_t i=0;i<d.node.size();i++){
        const auto&n=d.node[i];if(n.vertex>=0){out[i]="v;"+std::to_string(n.vertex);continue;}if(!used[i]){out[i]="unused;"+std::to_string(i);continue;}
        if(n.active.size()>=2){out[i]=events(n,planes,own);continue;}
        need(n.active.size()==1,"used nonvertex event has explicit power plane");int found=0;
        for(const auto&b:d.boundary)if(b.edge>=0&&(b.first==int(i)||b.last==int(i))){
            const S&t=b.first==int(i)?b.lo:b.hi;const auto&roots=edges.at(b.edge).events[n.active[0]].root;
            for(size_t k=0;k<roots.size();k++)if(roots[k].t.n==t.n&&roots[k].t.e==t.e){uint64_t a=std::min(own,planes[n.active[0]].competitor),z=std::max(own,planes[n.active[0]].competitor);std::string id="e;"+std::to_string(b.edge)+';'+std::to_string(a)+';'+std::to_string(z)+';'+std::to_string(k);need(out[i].empty()||out[i]==id,"single shared root event identity");out[i]=id;found++;}
        }
        need(found>0,"used edge crossing binds the actual canonical root index");
    }return out;
}
}
int main(int argc,char**argv){
    int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);int code=2;
    try{
        need(argc==12,"packet saved sites weights originalcall bounds registry roots knowncoverage NEWoutput originalTargetPublication");auto begin=std::chrono::steady_clock::now();auto p=curved::read(argv[1]);const auto f=load_saved(std::string(argv[2])+"/shared_state.bin",p);auto before=f.coefficient;
        auto site=sites(argv[3],argv[4]);std::sort(site.begin(),site.end(),[](const Site&a,const Site&b){return a.id<b.id;});auto source=target::retained_sources(argv[5],argv[6],f);target::coverage_overlay(argv[9],source);auto known=target::retained_domains(argv[9],source);
        std::string registry=argv[7],rootdir=argv[8],dest=argv[10];std::string line;std::vector<Slot>slot;std::ifstream pairs(registry+"/source_target_registry.csv");need(bool(pairs),"audited pair registry");std::getline(pairs,line);
        std::ifstream publication(std::string(argv[11])+"/native_exit_code.txt");int published=-1;need(bool(publication>>published)&&published==0,"original target publication is retained");
        while(std::getline(pairs,line)){auto c=target::columns(line);size_t id=std::stoull(c[0]);need(id==slot.size()&&std::stoull(c[3])==site[id%24].id,"original registry slot identity");slot.push_back({std::stoi(c[6]),std::stoi(c[7])});}need(slot.size()==9216,"all source-target slots");
        std::vector<S>image(768);std::ifstream images(registry+"/source_image_registry.csv");need(bool(images),"retained enclosure registry");std::getline(images,line);int count=0;
        while(std::getline(images,line)){auto c=target::columns(line);int s=std::stoi(c[0]);need(s==count++,"ordered image boxes");image[2*s]=S(target::number(c[2]),target::number(c[3]));image[2*s+1]=S(target::number(c[4]),target::number(c[5]));}need(count==384,"all image bounds");
        std::ifstream rootexit(rootdir+"/native_exit_code.txt");int exit=-1;need(bool(rootexit>>exit)&&exit==0,"all planned roots complete");std::map<Key,target::Roots>roots;
        for(int owner=0;owner<4;owner++){
            std::ifstream packets(rootdir+"/canonical_event_root_packets_rank"+std::to_string(owner)+".csv");need(bool(packets),"all canonical unordered pair root packets");std::getline(packets,line);
            while(std::getline(packets,line)){auto c=target::columns(line);need(c[5]=="1","actual roots complete");Key key(std::stoi(c[1]),std::stoull(c[2]),std::stoull(c[3]));auto&r=roots[key];r.complete=true;if(c[6]!="-1"){need(std::stoull(c[6])==r.root.size(),"ordered canonical event roots");r.root.push_back({S(target::number(c[8]),target::number(c[9])),std::stoi(c[10])});}snprintf(r.reason,sizeof(r.reason),"same_canonical_edge_unordered_exact_plane_roots");}
        }
        need(roots.size()==2452,"all registered roots recovered without replay");
        std::vector<sd::Compiled>maps(384);int bad=0;
        #pragma omp parallel for schedule(static) reduction(+:bad)
        for(int s=0;s<384;s++)try{maps[s]=sd::compile(f,s);}catch(const std::exception&){bad++;}
        need(!bad,"same canonical maps compile");auto topology=target::topology(maps);
        std::vector<std::map<int,target::Edge>>edge(24);std::vector<std::vector<target::Plane>>planes(24);
        for(int t=0;t<24;t++)if(t%ranks==rank){
            planes[t]=target::planes(site,site[t].id);edge[t]=topology;
            for(auto&entry:edge[t]){
                auto&e=entry.second;e.events.resize(23);for(int j=0;j<23;j++){
                    auto&event=e.events[j];uint64_t a=std::min(site[t].id,planes[t][j].competitor),b=std::max(site[t].id,planes[t][j].competitor);auto it=roots.find({e.key,a,b});
                    if(it!=roots.end())event=it->second;
                    else{auto value=target::value(planes[t][j],image[2*e.source],image[2*e.source+1]);bool separated=sf::lower(value)>0||sf::upper(value)<0;event.complete=separated;snprintf(event.reason,sizeof(event.reason),"%s",separated?"retained_whole_image_strict_plane_sign_excludes_all_edge_roots":"unplanned_unproved_edge_plane_not_needed_or_registry_defect");}
                }
            }
        }
        std::vector<Result>result(9216);
        #pragma omp parallel for schedule(dynamic,1)
        for(size_t id=0;id<slot.size();id++)if(int(id%24)%ranks==rank){
            int s=id/24,t=id%24;auto&r=result[id];const auto&row=slot[id];
            try{
                if(row.resolution==0){r.complete=1;r.empty=1;r.reused=site[t].id==68526;r.reason="certified_whole_image_exclusion";continue;}
                if(row.resolution==2){r.complete=1;r.empty=1;r.reused=1;r.reason="retained_original_actual_empty_domain";continue;}
                if(row.resolution==3||row.coarse==1){r.complete=1;r.reused=site[t].id==68526;r.reason="strict_whole_source_target_inclusion_implicit_original_triangle_boundary";continue;}
                if(row.resolution==4){r.domain=known[s];r.reused=1;r.reason="retained_original_target_actual_domain";}
                else{
                    auto contains=[&](S x,S y){auto membership=target::winding(maps[s],source[s].real,edge[t],x,y);r.windows+=membership.windows;r.query.push_back({x,y,membership});return target::inside_or_outside(membership);};
                    r.domain=target::build(maps[s],source[s].real,edge[t],planes[t],contains,&site[t],false);r.reason=r.domain.empty?"actual_empty_arrangement_winding":"actual_oriented_domain_arrangement_winding";
                }
                r.event=node_events(r.domain,edge[t],planes[t],site[t].id);
                for(const auto&b:r.domain.boundary)if(b.edge<0){int common=0;for(int j:r.domain.node[b.first].active)common+=std::find(r.domain.node[b.last].active.begin(),r.domain.node[b.last].active.end(),j)!=r.domain.node[b.last].active.end();need(common==1,"one exact plane for every internal target interface");}
                r.empty=r.domain.empty;r.complete=1;
            }catch(const std::exception&e){r.reason=e.what();}
        }
        FILE*ledger=fopen((dest+"/actual_domain_registry_rank"+std::to_string(rank)+".csv").c_str(),"w");FILE*nodes=fopen((dest+"/domain_nodes_rank"+std::to_string(rank)+".csv").c_str(),"w");FILE*boundaries=fopen((dest+"/domain_boundaries_rank"+std::to_string(rank)+".csv").c_str(),"w");FILE*queries=fopen((dest+"/actual_membership_queries_rank"+std::to_string(rank)+".csv").c_str(),"w");need(ledger&&nodes&&boundaries&&queries,"actual domain packets and failed arguments");
        fprintf(ledger,"slot_id,source,component,target,coarse_classification,retained_resolution,geometry_complete,empty,reused,loops,boundaries,winding_queries,winding_windows,reason\n");
        fprintf(nodes,"slot_id,node,global_vertex,x,x_radius,y,y_radius,event_id,active_competitors\n");
        fprintf(boundaries,"slot_id,source,target,boundary,first_node,last_node,first_event,last_event,canonical_edge,first_plane_target,second_plane_target,parameter_lo,lo_radius,parameter_hi,hi_radius,external_support_boundary,implicit_whole_source_boundary\n");
        fprintf(queries,"slot_id,query,x,x_radius,y,y_radius,location,ray,failed_edge,root_windows,reason\n");
        uint64_t totals[7]={};
        for(size_t id=0;id<slot.size();id++)if(int(id%24)%ranks==rank){
            int s=id/24,t=id%24;auto&r=result[id];totals[0]++;totals[1]+=!r.complete;totals[2]+=r.complete&&r.empty;totals[3]+=r.complete&&!r.empty&&(slot[id].coarse==1||slot[id].resolution==3);totals[4]+=r.complete&&!r.empty&&!(slot[id].coarse==1||slot[id].resolution==3);totals[5]+=r.query.size();totals[6]+=r.windows;
            fprintf(ledger,"%zu,%d,%d,%llu,%d,%d,%d,%d,%d,%zu,%zu,%zu,%llu,%s\n",id,s,p.record[s].component,(unsigned long long)site[t].id,slot[id].coarse,slot[id].resolution,r.complete,r.empty,r.reused,r.domain.loop.size(),r.domain.boundary.size(),r.query.size(),(unsigned long long)r.windows,quote(r.reason).c_str());
            for(size_t j=0;j<r.query.size();j++){const auto&q=r.query[j];fprintf(queries,"%zu,%zu,%s,%s,%s,%s,%d,%d,%d,%llu,%s\n",id,j,qs(q.x.n).c_str(),qs(q.x.e).c_str(),qs(q.y.n).c_str(),qs(q.y.e).c_str(),int(q.membership.location),q.membership.ray,q.membership.failed_edge,(unsigned long long)q.membership.windows,quote(q.membership.reason).c_str());}
            if(!r.complete||r.empty)continue;
            if(slot[id].coarse==1||slot[id].resolution==3){const auto&triangle=f.triangle[s];for(int e=0;e<3;e++){int a=e,b=(e+1)%3,orientation=triangle.endpoint[e][0]==a?1:-1;int key=triangle.edge[e].poly[1];std::string first="v;"+std::to_string(triangle.id[a]),last="v;"+std::to_string(triangle.id[b]);fprintf(boundaries,"%zu,%d,%llu,%d,%d,%d,%s,%s,%d,0,0,%d,0,%d,0,%d,1\n",id,s,(unsigned long long)site[t].id,e,a,b,first.c_str(),last.c_str(),key,orientation>0?0:1,orientation>0?1:0,topology.at(key).sides.size()==1);}continue;}
            for(size_t j=0;j<r.domain.node.size();j++){const auto&n=r.domain.node[j];std::string active;for(int a:n.active)active+=std::to_string(planes[t][a].competitor)+';';fprintf(nodes,"%zu,%zu,%d,%s,%s,%s,%s,%s,%s\n",id,j,n.vertex,qs(n.p[0].n).c_str(),qs(n.p[0].e).c_str(),qs(n.p[1].n).c_str(),qs(n.p[1].e).c_str(),r.event[j].c_str(),active.c_str());}
            for(size_t j=0;j<r.domain.boundary.size();j++){
                const auto&b=r.domain.boundary[j];uint64_t a=0,z=0;if(b.edge<0){std::vector<int>common;for(int k:r.domain.node[b.first].active)if(std::find(r.domain.node[b.last].active.begin(),r.domain.node[b.last].active.end(),k)!=r.domain.node[b.last].active.end())common.push_back(k);need(common.size()==1,"one exact plane for every internal target interface");a=std::min(site[t].id,planes[t][common[0]].competitor);z=std::max(site[t].id,planes[t][common[0]].competitor);}
                fprintf(boundaries,"%zu,%d,%llu,%zu,%d,%d,%s,%s,%d,%llu,%llu,%s,%s,%s,%s,%d,0\n",id,s,(unsigned long long)site[t].id,j,b.first,b.last,r.event[b.first].c_str(),r.event[b.last].c_str(),b.edge,(unsigned long long)a,(unsigned long long)z,qs(b.lo.n).c_str(),qs(b.lo.e).c_str(),qs(b.hi.n).c_str(),qs(b.hi.e).c_str(),b.edge>=0&&topology.at(b.edge).sides.size()==1);
            }
        }
        fclose(ledger);fclose(nodes);fclose(boundaries);fclose(queries);MPI_Allreduce(MPI_IN_PLACE,totals,7,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);need(totals[0]==9216,"all original recipients recorded");
        need(!memcmp(before.data(),f.coefficient.data(),before.size()*sizeof(STCoefficient)),"actual domains preserve state");
        if(rank==0){checkpoint(dest+"/caller_unchanged_state.bin",f,1);printf("{\"scope\":\"actual all24x384 domain arrangements only, canonical root reuse, immutable original support; no quadrature/evolution\",\"MPI_ranks\":%d,\"OMP_threads\":%d,\"pair_slots\":%llu,\"incomplete_domains\":%llu,\"complete_empty_or_excluded\":%llu,\"complete_whole_source\":%llu,\"complete_clipped_domains\":%llu,\"actual_membership_queries\":%llu,\"actual_membership_root_windows\":%llu,\"all_domain_arrangements_complete\":%s,\"source_partition_boundary_cancellation_proved\":false,\"global_source_image_disjointness_proved\":false,\"full_quadrature_allowed\":false,\"new_physical_quadrature_queries\":0,\"accepted_RK4_evolution\":false,\"long_Kepler_validated\":false,\"seconds\":%.17g}\n",ranks,omp_get_max_threads(),(unsigned long long)totals[0],(unsigned long long)totals[1],(unsigned long long)totals[2],(unsigned long long)totals[3],(unsigned long long)totals[4],(unsigned long long)totals[5],(unsigned long long)totals[6],totals[1]?"false":"true",std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count());}code=totals[1]?2:0;
    }catch(const std::exception&e){fprintf(stderr,"rank%d FAIL %s\n",rank,e.what());MPI_Abort(MPI_COMM_WORLD,2);}MPI_Finalize();return code;
}
