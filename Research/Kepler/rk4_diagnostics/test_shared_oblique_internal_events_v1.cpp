#define GFS_PENDING29_PROJECTION_LIBRARY
#include "test_shared_pending29_projection_v1.cpp"
namespace {
struct InternalEvent29 {uint64_t slot,source,target,chart,face,side,index,count,complete,windows;S parameter,s,t,reference[2],projection[2];char reason[192];};
static sd::Jet side_function(const target::Face&face,const Frame29&f,int side,S u){
    S s,t,ds,dt;
    if(side==0){s=S(f.cut)+(S(1)-S(f.cut))*u;t=S(0);ds=S(1)-S(f.cut);}
    else if(side==1){s=S(1);t=u;dt=S(1);}
    else{need(side==2,"one original internal radial edge");s=S(1)-(S(1)-S(f.cut))*u;t=S(1);ds=-(S(1)-S(f.cut));}
    auto a=target::face(face,f.original,s,t);sd::Jet out(a.v,a.g[0]*ds+a.g[1]*dt);out.H[0][0]=a.H[0][0]*ds*ds+(a.H[0][1]+a.H[1][0])*ds*dt+a.H[1][1]*dt*dt;return out;
}
}
int main(int argc,char**argv){int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);try{
    need(argc==7&&ranks==1,"packet saved sites weights complete_referencepacket NEWdest");auto begin=std::chrono::steady_clock::now();auto p=curved::read(argv[1]);const auto f=load_saved(std::string(argv[2])+"/shared_state.bin",p);auto before=f.coefficient;auto site=sites(argv[3],argv[4]);std::sort(site.begin(),site.end(),[](const Site&a,const Site&b){return a.id<b.id;});std::string previous=argv[5],dest=argv[6];auto frames=read29<Frame29>(previous+"/oblique_coordinate_frames_rank0.bin");need(frames.size()==14,"only14 failed original oblique charts");std::map<int,sd::Compiled>maps;for(const auto&z:frames)if(!maps.count(z.source))maps.emplace(z.source,sd::compile(f,z.source));
    std::vector<std::vector<InternalEvent29>>records(frames.size());std::vector<std::string>reasons(frames.size());std::vector<uint64_t>completed(frames.size()),jobs(frames.size()),windows(frames.size());
    #pragma omp parallel for schedule(static)
    for(size_t j=0;j<frames.size();j++){const auto&z=frames[j];try{need(z.passed&&z.cut>0,"retained outer chart proof");const auto&m=maps.at(z.source);auto ps=target::planes(site,z.target);std::vector<target::Face>faces;for(const auto&plane:ps)faces.push_back(target::pullback(m,plane));
        // Standard chart part0 has canonical source-edge side0; part1 side2.
        // Neither is replayed. Side3 is strictly excluded by the saved proof.
        std::vector<int>sides=z.original.part==0?std::vector<int>{1,2}:std::vector<int>{0,1};
        for(int side:sides)for(size_t face=0;face<ps.size();face++)if(z.active_mask&(uint64_t(1)<<face)){
            auto roots=target::scalar_roots([&](S u){return side_function(faces[face],z,side,u);});jobs[j]++;windows[j]+=roots.windows;
            for(size_t k=0;k<std::max(size_t(1),roots.root.size());k++){InternalEvent29 e={};e.slot=z.slot;e.source=z.source;e.target=z.target;e.chart=z.chart;e.face=face;e.side=side;e.index=roots.root.empty()?UINT64_MAX:k;e.count=roots.root.size();e.complete=roots.complete;e.windows=roots.windows;snprintf(e.reason,sizeof(e.reason),"%s",roots.reason);
                if(!roots.root.empty()){e.parameter=roots.root[k];if(side==0){e.s=S(z.cut)+(S(1)-S(z.cut))*e.parameter;e.t=S(0);}else if(side==1){e.s=S(1);e.t=e.parameter;}else{e.s=S(1)-(S(1)-S(z.cut))*e.parameter;e.t=S(1);}auto r=target::frame(m,z.original,e.t);for(int q=0;q<2;q++)e.reference[q]=r.A[q]+e.s*r.D[q];project29(z,e.reference,e.projection);}records[j].push_back(e);
            }
            need(roots.complete,roots.reason);
        }
        completed[j]=1;reasons[j]="complete_only_new_internal_chart_boundary_root_sets";
    }catch(const std::exception&e){reasons[j]=e.what();}}
    std::ofstream summary(dest+"/internal_chart_boundary_status.csv"),events(dest+"/internal_chart_boundary_events.csv");summary<<"slot_id,source,target,chart,complete,new_internal_root_jobs,root_windows,reason\n";events<<"slot_id,source,target,chart,face,side,index,count,complete,windows,parameter,parameter_radius,s,s_radius,t,t_radius,u,u_radius,v,v_radius,reason\n";std::vector<InternalEvent29>all;uint64_t done=0,totaljobs=0,totalwindows=0,rootcount=0;
    for(size_t j=0;j<frames.size();j++){const auto&z=frames[j];done+=completed[j];totaljobs+=jobs[j];totalwindows+=windows[j];summary<<z.slot<<','<<z.source<<','<<z.target<<','<<z.chart<<','<<completed[j]<<','<<jobs[j]<<','<<windows[j]<<','<<reasons[j]<<'\n';append(all,records[j]);for(const auto&e:records[j]){rootcount+=e.index!=UINT64_MAX;events<<e.slot<<','<<e.source<<','<<e.target<<','<<e.chart<<','<<e.face<<','<<e.side<<','<<e.index<<','<<e.count<<','<<e.complete<<','<<e.windows;for(S a:{e.parameter,e.s,e.t,e.projection[0],e.projection[1]})events<<','<<qs(a.n)<<','<<qs(a.e);events<<','<<e.reason<<'\n';}}
    summary.close();events.close();packet(dest+"/internal_chart_boundary_events_rank0.bin",all,0);packet(dest+"/unchanged_oblique_frames_rank0.bin",frames,0);need(!memcmp(before.data(),f.coefficient.data(),before.size()*sizeof(STCoefficient)),"same original state and physical heat");printf("{\"scope\":\"ONLY14 failed oblique charts, two genuinely new internal boundary edges each; canonical edge roots and completed fields/panels never replayed\",\"completed_internal_boundary_root_charts\":%llu,\"pending_internal_boundary_root_charts\":%llu,\"new_internal_boundary_root_jobs\":%llu,\"new_internal_root_work_windows\":%llu,\"stored_rootset_rows\":%zu,\"stored_internal_roots\":%llu,\"OMP_threads\":%d,\"new_inverse_solves\":0,\"new_canonical_root_jobs\":0,\"normal_field_or_panel_replays\":0,\"new_physical_quadrature_queries\":0,\"remaining_incomplete_atlas_slots\":29,\"actual_event_slab_atlas_still_pending\":true,\"full_quadrature_allowed\":false,\"state_unchanged\":true,\"seconds\":%.17g}\n",(unsigned long long)done,(unsigned long long)(14-done),(unsigned long long)totaljobs,(unsigned long long)totalwindows,all.size(),(unsigned long long)rootcount,omp_get_max_threads(),std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count());
}catch(const std::exception&e){fprintf(stderr,"oblique internal boundary events FAIL %s\n",e.what());MPI_Abort(MPI_COMM_WORLD,2);}MPI_Finalize();return 0;}
