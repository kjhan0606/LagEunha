#define GFS_PENDING_CLEARANCE_LIBRARY
#include "test_shared_pending_clearance_v1.cpp"
namespace {
struct SelectedDerivative {uint64_t slot,source,target,chart,face,dependent;sd::Chart reference;S gradient[2],direction[2],derivative;};
static_assert(sizeof(SelectedDerivative)==544,"selected derivative proof ABI");
static std::vector<TypedPanel> baseline_read(const std::string&path){
 FILE*fp=fopen(path.c_str(),"rb");need(fp,"baseline saved panel input");uint64_t h[8];
 need(fread(h,sizeof(h),1,fp)==1&&(h[0]==0x4753415749544e31ULL||h[0]==0x475341544c415331ULL)&&h[1]==1&&h[2]==1&&h[3]==704&&h[6]==24&&h[7]==384,"baseline frozen typed header");
 std::vector<TypedPanel>v(h[4]);need(v.empty()||fread(v.data(),sizeof(TypedPanel),v.size(),fp)==v.size(),"baseline exact panel records");need(fgetc(fp)==EOF,"baseline packet exact end");fclose(fp);return v;
}
static void baseline_charts(const sd::Compiled&m,const std::vector<target::Face>&faces,Clearance&w,std::vector<SelectedDerivative>&proof){
 auto shape=target::shape_bounds(m);std::map<std::pair<int,int>,bool>seen;
 for(const auto&b:w.geometry.panel){const auto&p=b.geometry;int ci=p.strip.chart;need(ci>=0&&ci<6,"saved baseline reference chart index");auto&c=w.geometry.chart[ci];
  if(!w.geometry.seen[ci]){c.reference=p.reference;c.complete=true;c.coordinate.derivative.resize(faces.size());c.coordinate.dependent=p.whole?-1:p.strip.dependent;w.geometry.seen[ci]=true;}
  else need(!memcmp(&c.reference,&p.reference,sizeof(sd::Chart)),"same saved reference for every baseline chart panel");
  if(p.whole)continue;
  if(c.coordinate.dependent<0)c.coordinate.dependent=p.strip.dependent;
  need(c.coordinate.dependent==p.strip.dependent&&b.dependent==p.strip.dependent,"same saved dependent coordinate");
  auto frame=target::frame(m,p.reference,sf::range(0,1));
  for(int face:{p.strip.lower,p.strip.upper})if(face>=0&&!seen[{ci,face}]){
   need(size_t(face)<faces.size(),"selected saved boundary face in range");seen[{ci,face}]=true;
   SelectedDerivative z={};z.slot=b.slot;z.source=b.source;z.target=b.target;z.chart=ci;z.face=face;z.dependent=p.strip.dependent;z.reference=p.reference;
   auto cone=target::gradient_cone(faces[face],shape);
   for(int k=0;k<2;k++){z.gradient[k]=cone.gradient[k];z.direction[k]=p.strip.dependent?frame.E[k]:frame.D[k];z.derivative=z.derivative+z.gradient[k]*z.direction[k];}
   need(sf::lower(z.derivative)>0||sf::upper(z.derivative)<0,"selected saved graph derivative remains separated");
   c.coordinate.derivative[face]=z.derivative;c.coordinate.active.push_back(face);proof.push_back(z);
  }
 }
}
}
#ifndef GFS_BASELINE_FIELDS_LIBRARY
int main(int argc,char**argv){int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);try{
 need(argc==7&&ranks==1,"packet saved sites weights exact_slot_jobs NEWdest");auto start=std::chrono::steady_clock::now();auto p=curved::read(argv[1]);const auto f=load_saved(std::string(argv[2])+"/shared_state.bin",p);auto before=f.coefficient;auto site=sites(argv[3],argv[4]);std::sort(site.begin(),site.end(),[](const Site&a,const Site&b){return a.id<b.id;});std::string dest=argv[6],line;
 std::map<std::string,std::vector<TypedPanel>>cache;std::map<uint64_t,Clearance>selected;std::map<uint64_t,std::string>kind;std::ifstream input(argv[5]);need(bool(input),"fixed missing-field input table");std::getline(input,line);uint64_t total=0,historical=0;
 while(std::getline(input,line)){if(!line.empty()&&line.back()=='\r')line.pop_back();auto c=target::columns(line);need(c.size()==7,"exact seven input columns");uint64_t slot=std::stoull(c[0]),source=std::stoull(c[1]),target=std::stoull(c[2]),first=std::stoull(c[4]),count=std::stoull(c[5]);need(!selected.count(slot),"unique missing slot");auto&w=selected[slot].geometry;w.slot=slot;w.source=source;w.target=target;kind[slot]=c[6];need(c[6]=="new_complex_evidence_and_full_field"||c[6]=="recover_historical_pass_data","exact work classification");if(!cache.count(c[3]))cache.emplace(c[3],baseline_read(c[3]));const auto&v=cache.at(c[3]);need(count>0&&first+count<=v.size(),"fixed input range");
  for(uint64_t j=0;j<count;j++){auto b=v[first+j];need(b.slot==slot&&b.source==source&&b.target==target&&b.index==j,"same original plan identity/index");w.panel.push_back(b);}total+=count;if(c[6]=="recover_historical_pass_data")historical+=2*count;
 }
 need(total==34955&&historical==7952&&selected.size()==500,"exact missing34955panels/69910axes, none of saved110332replayed");
 std::map<int,sd::Compiled>maps;std::vector<Clearance*>work;for(auto&e:selected){work.push_back(&e.second);if(!maps.count(e.second.geometry.source))maps.emplace(e.second.geometry.source,sd::compile(f,e.second.geometry.source));}std::vector<std::vector<SelectedDerivative>>derivatives(work.size());uint64_t finished=0;
 #pragma omp parallel for schedule(dynamic,1)
 for(size_t j=0;j<work.size();j++){auto&w=*work[j];try{std::vector<target::Face>faces;const auto&m=maps.at(w.geometry.source);for(const auto&plane:target::planes(site,w.geometry.target))faces.push_back(target::pullback(m,plane));baseline_charts(m,faces,w,derivatives[j]);check(p,m,faces,w);}catch(const std::exception&e){w.reason=e.what();}
  #pragma omp critical(gfs_baseline_progress)
  {finished++;if(!w.complete||finished%25==0||finished==work.size()){fprintf(stderr,"finished=%llu/500 slot=%llu fields=%zu/%zu complete=%d reason=%s\n",(unsigned long long)finished,(unsigned long long)w.geometry.slot,w.fields.size(),2*w.geometry.panel.size(),w.complete,w.reason.c_str());fflush(stderr);}}
 }
 std::vector<RegularAxis>axes;std::vector<FullField>fields;std::vector<pr::Boundary>boundaries;std::vector<SelectedDerivative>proof;std::vector<TypedPanel>unchanged;uint64_t complete=0;
 std::ofstream csv(dest+"/regular_geometry_density_clearance.csv");csv<<"slot_id,source,target,regular_panels,checked_axes,passed_axes,clearance_complete,cached_boundary_records,work_kind,reason\n";
 for(size_t j=0;j<work.size();j++){const auto&w=*work[j];append(axes,w.axes);append(fields,w.fields);append(boundaries,w.boundaries);append(proof,derivatives[j]);append(unchanged,w.geometry.panel);complete+=w.complete;csv<<w.geometry.slot<<','<<w.geometry.source<<','<<w.geometry.target<<','<<w.geometry.panel.size()<<','<<w.axes.size()<<','<<w.fields.size()<<','<<w.complete<<','<<w.boundaries.size()<<','<<kind.at(w.geometry.slot)<<','<<w.reason<<'\n';}csv.close();
 packet(dest+"/regular_axis_proofs_rank0.bin",axes,0);packet(dest+"/passed_regular_full_fields_rank0.bin",fields,0);packet(dest+"/prepared_regular_boundaries_rank0.bin",boundaries,0);packet(dest+"/selected_derivative_proofs_rank0.bin",proof,0);packet(dest+"/unchanged_baseline_panels_rank0.bin",unchanged,0);
 need(!memcmp(before.data(),f.coefficient.data(),before.size()*sizeof(STCoefficient)),"baseline recovery preserves all physical coefficients/heat");
 printf("{\"scope\":\"only34955baseline panels missing full fields; historical source133 clearance retained, same frozen CONTROL state and exact plans; save all passing fields\",\"regular_panels\":34955,\"slot_jobs\":500,\"completed_regular_slot_clearances\":%llu,\"new_complex_axis_evaluations\":%zu,\"new_complex_axis_passes\":%zu,\"saved_full_fields\":%zu,\"cached_boundary_records\":%zu,\"selected_derivative_records\":%zu,\"historical_pass_axes_scheduled_for_data_recovery\":7952,\"previous_successful_field_replays\":0,\"new_arrangements\":0,\"new_event_root_jobs\":0,\"new_inverse_solves\":0,\"new_physical_quadrature_queries\":0,\"all12_analytic_remainders_complete\":false,\"full_quadrature_allowed\":false,\"state_unchanged\":true,\"OMP_threads\":%d,\"seconds\":%.17g}\n",(unsigned long long)complete,axes.size(),fields.size(),fields.size(),boundaries.size(),proof.size(),omp_get_max_threads(),std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());
 }catch(const std::exception&e){fprintf(stderr,"baseline full-field recovery FAIL %s\n",e.what());MPI_Abort(MPI_COMM_WORLD,2);}MPI_Finalize();return 0;}
#endif
