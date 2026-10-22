#define GFS_PENDING_CLEARANCE_LIBRARY
#include "test_shared_pending_clearance_v1.cpp"
#include "material_shared_logistic_density.h"
namespace {
struct DensityWitness {uint64_t slot,source,target,panel,axis,passed,old_failed;S radius;target::DensityDisk density;char reason[192];};
static_assert(sizeof(DensityWitness)==672,"saved density witness ABI");
}
int main(int argc,char**argv){int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);int ranks;MPI_Comm_size(MPI_COMM_WORLD,&ranks);try{
 need(argc==4&&ranks==1,"original packet saved_baseline_axis_packet NEWdest");auto start=std::chrono::steady_clock::now();const auto p=curved::read(argv[1]);auto axes=read29<RegularAxis>(argv[2]);std::string dest=argv[3];std::vector<DensityWitness>witness;
 std::ofstream params(dest+"/original_density_parameters.csv");params<<"parameter,value\n";for(int i=0;i<10;i++)params<<i<<','<<qs(Q(p.parameters[i]))<<'\n';params.close();
 uint64_t passed=0;for(const auto&a:axes)if(!a.passed){DensityWitness z={};z.slot=a.slot;z.source=a.source;z.target=a.target;z.panel=a.panel;z.axis=a.axis;z.radius=a.reference_radius;z.old_failed=1;
  try{need(std::string(a.reason)=="original logistic density poles not excluded on complex domain","only saved logistic representation failures");bool reproduced=false;try{(void)target::complex_density_bound(p,z.radius);}catch(const std::exception&e){reproduced=std::string(e.what())==a.reason;}need(reproduced,"same saved radius reproduces old density-bound failure");z.density=target::stable_reference_density(p,z.radius);z.passed=1;passed++;snprintf(z.reason,sizeof(z.reason),"same_density_disk_two_sided_exponential_denominators_exclude_zero");}
  catch(const std::exception&e){snprintf(z.reason,sizeof(z.reason),"%s",e.what());}witness.push_back(z);
 }
 need(witness.size()==64,"exact64 stored first failure witnesses, no tail evaluation");packet(dest+"/same_radius_density_witnesses_rank0.bin",witness,0);
 std::ofstream csv(dest+"/density_witnesses.csv");csv<<"slot_id,source,target,panel,axis,passed,reason\n";for(const auto&z:witness)csv<<z.slot<<','<<z.source<<','<<z.target<<','<<z.panel<<','<<z.axis<<','<<z.passed<<','<<z.reason<<'\n';csv.close();
 printf("{\"scope\":\"same64 saved failed radius disks only; equivalent two-sided logistic identities and unchanged density profile; no field or suffix replay\",\"saved_old_failures\":64,\"same_disk_density_passes\":%llu,\"new_complex_field_queries\":0,\"new_physical_quadrature_queries\":0,\"unattempted_suffix_axes_not_certified\":11264,\"all12_complete\":false,\"full_quadrature_allowed\":false,\"seconds\":%.17g}\n",(unsigned long long)passed,std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());
 }catch(const std::exception&e){fprintf(stderr,"density witness FAIL %s\n",e.what());MPI_Abort(MPI_COMM_WORLD,2);}MPI_Finalize();return 0;}
