#define GFS_BASELINE_FIELDS_LIBRARY
#include "test_shared_baseline_fields_v1.cpp"
#include "material_shared_reference80_stream.h"
#include "material_shared_runtime_reference.h"
namespace rr=runtime_reference;
namespace {
struct RegularRealInput {uint64_t group,slot,source,target,chart,dependent,node,passed,lower_iterations,upper_iterations;S free;sd::AffineAngular angular_lo,angular_hi;S radial_lo,radial_hi,derivative[2],argument;Q argument_limit;int32_t phase,edge;char reason[192];};
struct ObliqueRealInput {uint64_t group,slot,source,target,chart,first_panel,node,passed;S u;op::Branches branches;S argument;Q argument_limit;int32_t phase,edge;char reason[192];};
static_assert(sizeof(RegularRealInput)==1456&&sizeof(ObliqueRealInput)==800,"saved real cache ABI");
struct InputWitness {uint64_t slot,source,target,kind,panel,iu,iv,old_bad_edges,new_bad_edges;S X[2],lambda[3],old_lambda[3],new_lambda[3],old_sum[3],new_sum[3];};
static_assert(sizeof(InputWitness)==624,"native reference input witness ABI");
struct Item {uint64_t slot,source,target,kind,panel,group;TypedPanel regular;Panel oblique;};
struct Task {Item item;int node,iv;const RegularRealInput*regular=nullptr;const ObliqueRealInput*oblique=nullptr;};
}
int main(int argc,char**argv){int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);int ranks;MPI_Comm_size(MPI_COMM_WORLD,&ranks);try{
 need(argc==7&&ranks==1,"packet saved mixed_inputs regular_real_cache oblique_real_cache NEWdest");auto start=std::chrono::steady_clock::now();auto p=curved::read(argv[1]);const auto f=load_saved(std::string(argv[2])+"/shared_state.bin",p);const auto before=f.coefficient;std::string input=argv[3],dest=argv[6],line;
 auto regular=read29<RegularRealInput>(argv[4]);auto oblique=read29<ObliqueRealInput>(argv[5]);need(regular.size()==104336&&oblique.size()==4384,"all completed actual real node inputs");
 std::map<uint64_t,std::string>paths;std::ifstream registry(input+"/../shared_atlas_adoption_v7/typed_packet_registry.csv");need(bool(registry),"same v7 typed packet registry");std::getline(registry,line);while(std::getline(registry,line)){auto c=target::columns(line);paths.emplace(std::stoull(c[0]),c[1]);}
 std::map<std::string,std::vector<TypedPanel>>rplans;std::map<std::string,std::vector<Panel>>oplans;std::vector<Item>whole;std::map<uint64_t,std::vector<Item>>grouped;std::ifstream jobs(input+"/panel_jobs.csv");need(bool(jobs),"fixed remaining panels");std::getline(jobs,line);
 while(std::getline(jobs,line)){auto c=target::columns(line);Item x={};x.slot=std::stoull(c[0]);x.source=std::stoull(c[1]);x.target=std::stoull(c[2]);x.panel=std::stoull(c[4]);auto path=paths.at(std::stoull(c[5]));uint64_t index=std::stoull(c[6]);int64_t group=std::stoll(c[7]);x.group=group;
  if(c[3]=="regular"){if(!rplans.count(path))rplans.emplace(path,baseline_read(path));x.regular=rplans.at(path).at(index);x.kind=x.regular.geometry.whole?0:1;if(x.kind==0)whole.push_back(x);else if(group>=0)grouped[group].push_back(x);}
  else{need(c[3]=="oblique"&&group>=0,"same oblique runtime job");if(!oplans.count(path))oplans.emplace(path,read29<Panel>(path));x.oblique=oplans.at(path).at(index);x.kind=2;grouped[group].push_back(x);}
 }
 std::map<uint64_t,Frame29>frame;std::map<uint64_t,op::Linear>linear;std::map<uint64_t,sd::Compiled>maps;std::map<std::string,std::vector<Frame29>>frames;
 std::ifstream groups(input+"/real_graph_groups.csv");need(bool(groups),"same group frame references");std::getline(groups,line);while(std::getline(groups,line)){auto c=target::columns(line);if(c[4]!="oblique_same_frame_real_branches")continue;auto path=c[11];if(!frames.count(path))frames.emplace(path,read29<Frame29>(path));uint64_t g=std::stoull(c[0]);frame.emplace(g,frames.at(path).at(std::stoull(c[12])));const auto&fr=frame.at(g);if(!maps.count(fr.source))maps.emplace(fr.source,sd::compile(f,fr.source));linear.emplace(g,op::linear(maps.at(fr.source),fr));}
 std::vector<Task>tasks;for(const auto&x:whole)for(int iu:{0,15})for(int iv:{0,15})tasks.push_back({x,iu,iv,nullptr,nullptr});
 for(const auto&r:regular){need(r.passed,"only completed real branch nodes");const auto&v=grouped.at(r.group);for(size_t which:{size_t(0),v.size()-1}){for(int iv:{0,15})tasks.push_back({v[which],int(r.node),iv,&r,nullptr});if(v.size()==1)break;}}
 for(const auto&o:oblique){need(o.passed,"only completed oblique nodes");const auto&v=grouped.at(o.group);for(size_t which:{size_t(0),v.size()-1}){for(int iv:{0,15})tasks.push_back({v[which],int(o.node),iv,nullptr,&o});if(v.size()==1)break;}}
 need(tasks.size()<700000,"fixed native input-only witness budget");auto gauss=target::gauss16();std::vector<std::vector<InputWitness>>witness(omp_get_max_threads());uint64_t old_bad=0,new_bad=0,failures=0;
 #pragma omp parallel for schedule(dynamic,256) reduction(+:old_bad,new_bad,failures)
 for(size_t i=0;i<tasks.size();i++){const auto&task=tasks[i];const auto&x=task.item;try{const auto&t=f.triangle[x.source];rr::Point q;
  if(x.kind==0)q=rr::whole(t,x.regular,gauss.x[task.node],gauss.x[task.iv]);
  else if(x.kind==1){const auto&panel=x.regular.geometry;S eta=panel.eta_lo+(panel.eta_hi-panel.eta_lo)*gauss.x[task.iv];q=rr::regular(t,x.regular,task.regular,task.regular->free,eta);}
  else{const auto&b=x.oblique;S eta=b.eta0+(b.eta1-b.eta0)*gauss.x[task.iv];q=rr::oblique(t,frame.at(x.group),linear.at(x.group),b,task.oblique->u,eta,task.oblique->branches);}
  using T=long double;namespace r80=shared_reference80_ball;T X=T(q.X[0].n),Y=T(q.X[1].n),dx=X-T(t.vertex[0][0]),dy=Y-T(t.vertex[0][1]),ux=T(t.vertex[1][0])-T(t.vertex[0][0]),uy=T(t.vertex[1][1])-T(t.vertex[0][1]),vx=T(t.vertex[2][0])-T(t.vertex[0][0]),vy=T(t.vertex[2][1])-T(t.vertex[0][1]),det=ux*vy-uy*vx;T old[3];old[1]=(dx*vy-dy*vx)/det;old[2]=(ux*dy-uy*dx)/det;old[0]=1-old[1]-old[2];mf::B<T>o[3],n[3];InputWitness w={};w.slot=x.slot;w.source=x.source;w.target=x.target;w.kind=x.kind;w.panel=x.panel;w.iu=task.node;w.iv=task.iv;
  for(int k=0;k<2;k++)w.X[k]=q.X[k];for(int j=0;j<3;j++){w.lambda[j]=q.lambda[j];o[j]=r80::cache(old[j],q.lambda[j]);n[j]=r80::cache(T(q.lambda[j].n),q.lambda[j]);w.old_lambda[j]=S(Q(o[j].n),Q(o[j].e));w.new_lambda[j]=S(Q(n[j].n),Q(n[j].e));}
  for(int e=0;e<3;e++){int a=t.endpoint[e][0],b=t.endpoint[e][1];auto osum=o[a]+o[b],nsum=n[a]+n[b];w.old_sum[e]=S(Q(osum.n),Q(osum.e));w.new_sum[e]=S(Q(nsum.n),Q(nsum.e));if(!(osum.n>osum.e))w.old_bad_edges|=uint64_t(1)<<e;if(!(nsum.n>nsum.e))w.new_bad_edges|=uint64_t(1)<<e;}
  old_bad+=w.old_bad_edges!=0;new_bad+=w.new_bad_edges!=0;if(w.old_bad_edges||w.new_bad_edges)witness[omp_get_thread_num()].push_back(w);
 }catch(const std::exception&e){failures++;
  #pragma omp critical(gfs_native_input_failure)
  {fprintf(stderr,"reference_input_failure slot=%llu panel=%llu kind=%llu reason=%s\n",(unsigned long long)x.slot,(unsigned long long)x.panel,(unsigned long long)x.kind,e.what());}
 }}
 std::vector<InputWitness>all;for(const auto&v:witness)append(all,v);std::sort(all.begin(),all.end(),[](const InputWitness&a,const InputWitness&b){return std::tie(a.slot,a.kind,a.panel,a.iu,a.iv)<std::tie(b.slot,b.kind,b.panel,b.iu,b.iv);});packet(dest+"/native_barycentric_input_witnesses_rank0.bin",all,0);need(!memcmp(before.data(),f.coefficient.data(),before.size()*sizeof(STCoefficient)),"reference input audit preserves state and heat");
 printf("{\"scope\":\"fixed actual Gauss input-only endpoints of remaining whole panels and first/last eta panels of every prepared regular/oblique free group; no full field or physical query; not a full transfer or a global field-error proof\",\"checked_actual_reference_inputs\":%zu,\"legacy_Cartesian_inverse_inputs_with_nonseparated_edge_sum\":%llu,\"direct_reference_lambda_inputs_with_nonseparated_edge_sum\":%llu,\"geometry_input_failures\":%llu,\"saved_input_witnesses\":%zu,\"new_full_field_queries\":0,\"new_physical_quadrature_queries\":0,\"accepted_RK4_evolution\":false,\"long_Kepler_validated\":false,\"state_unchanged\":true,\"OMP_threads\":%d,\"seconds\":%.17g}\n",tasks.size(),(unsigned long long)old_bad,(unsigned long long)new_bad,(unsigned long long)failures,all.size(),omp_get_max_threads(),std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());
 }catch(const std::exception&e){fprintf(stderr,"native reference input preflight FAIL %s\n",e.what());MPI_Abort(MPI_COMM_WORLD,2);}MPI_Finalize();return 0;}
