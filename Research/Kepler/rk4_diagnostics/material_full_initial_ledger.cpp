// Consume the completed full initial operator. No point query, factorization,
// new mass equation or trajectory is performed here.
#include <mpi.h>
#include <omp.h>
#include <algorithm>
#include <array>
#include <vector>
#include <map>
#include <string>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <stdexcept>
static void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
using R=long double;
#include "material_full_owned.h"
#include "material_full_initial_batch.h"
#include "material_distributed_mass.h"
namespace shared_trace_rhs {using Vec=std::vector<R>;}
namespace st=shared_trace;
#include "material_geometric_pressure_geometry.h"
namespace md=material_distributed;namespace gp=geometric_pressure;
static R val(MFDoubleDouble a){return R(a.hi)+a.lo;}
template<size_t N>static void header(FILE*f,const std::array<uint64_t,N>&expected){
    need(f,"saved initial ledger input");std::array<uint64_t,N>h;material_full::read(f,h.data(),N);need(h==expected,"saved initial ledger input identity/ABI");
}
static void complete(FILE*f){need(fgetc(f)==EOF&&fclose(f)==0,"complete saved ledger input");}
static int run(const std::string&saved,const std::string&curv,MPI_Comm comm){
    int rank,ranks;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&ranks);need(ranks==4,"saved native four-rank ownership");double begin=MPI_Wtime();std::string suffix=std::to_string(rank);
    auto m=material_full::load_basis(saved+"/full_conditioned_basis_r"+suffix+".bin");size_t nd=m.global.size(),nt=m.f.triangle.size(),nc=m.f.cell.size();md::Ownership layout(m.global,m.owner,comm);need(layout.global_owned==m.global_dofs,"full unique owner count");
    // Unlike a floating-point equality check, this includes signed-zero bits.
    for(int d=0;d<4;d++){md::Vec v(nd);for(size_t i=0;i<nd;i++)v[i]=d<2?m.f.coefficient[i].x[d]:m.f.coefficient[i].v[d-2];layout.pull(v);for(size_t i=0;i<nd;i++){double got=double(v[i]),want=d<2?m.f.coefficient[i].x[d]:m.f.coefficient[i].v[d-2];need(memcmp(&got,&want,sizeof(double))==0,"literal owner/ghost coefficient bits including signed zero");}}
    std::array<md::Vec,2>a={md::Vec(nd),md::Vec(nd)};FILE*f=fopen((saved+"/full_initial_acceleration_r"+suffix+".bin").c_str(),"rb");header(f,std::array<uint64_t,4>{0x4d46414343454c31ULL,1,nd,sizeof(R)});for(auto&v:a){material_full::read(f,v.data(),nd);auto before=v;layout.pull(v);need(v==before,"saved acceleration owner/ghost values");}complete(f);
    f=fopen((saved+"/full_fixed_mass_r"+suffix+".bin").c_str(),"rb");header(f,std::array<uint64_t,4>{0x4d464d4153533031ULL,1,nt*ST_MAX*ST_MAX,sizeof(R)});md::Vec mass(nt*ST_MAX*ST_MAX);material_full::read(f,mass.data(),mass.size());complete(f);
    std::vector<std::array<R,2*ST_MAX>>force(nt);std::vector<std::array<R,MF_INITIAL_METRICS>>initial(nc);std::vector<std::array<R,MF_CURV_METRICS>>point(nc);
    std::vector<MFInitialBlock>batch(MF_BATCH_TRIANGLES);f=fopen((saved+"/full_raw_initial_blocks_r"+suffix+".bin").c_str(),"rb");header(f,std::array<uint64_t,6>{0x4d46494e49545231ULL,1,nt,sizeof(MFInitialBlock),uint64_t(rank),m.global_dofs});
    for(size_t at=0;at<nt;at+=batch.size()){size_t n=std::min(batch.size(),nt-at);material_full::read(f,batch.data(),n);
        #pragma omp parallel for schedule(static)
        for(size_t j=0;j<n;j++){size_t s=at+j;const auto&t=m.f.triangle[s];R T[ST_MAX][ST_MAX]={};for(int k=0;k<t.n;k++)T[k][k]=1;for(int k=0;k<2;k++){int b=t.bubble_offset+k;for(int z=0;z<t.n;z++)T[b][z]=0;for(int z=0;z<2;z++){T[b][t.bubble_offset+z]=t.bubble_encode[k][z];T[b][9]-=R(t.bubble_encode[k][z])*t.bubble_shift[z];}}for(int d=0;d<2;d++)for(int k=0;k<t.n;k++)for(int z=0;z<t.n;z++)force[s][d*ST_MAX+k]+=T[k][z]*val(batch[j].force[d][z]);}
        for(size_t j=0;j<n;j++){need(batch[j].bad_point<0,"retained initial batch valid");for(int k=0;k<MF_INITIAL_METRICS;k++)initial[m.f.triangle[at+j].cell][k]+=val(batch[j].metric[k]);}
    }complete(f);
    std::vector<MFCurvatureBlock>cb(MF_BATCH_TRIANGLES);f=fopen((curv+"/full_curvature_blocks_r"+suffix+".bin").c_str(),"rb");header(f,std::array<uint64_t,5>{0x4d46435552563031ULL,1,nt,sizeof(MFCurvatureBlock),uint64_t(rank)});
    for(size_t at=0;at<nt;at+=cb.size()){size_t n=std::min(cb.size(),nt-at);material_full::read(f,cb.data(),n);for(size_t j=0;j<n;j++){need(cb[j].bad_point<0,"retained initial curvature valid");for(int k=0;k<MF_CURV_METRICS;k++)point[m.f.triangle[at+j].cell][k]+=val(cb[j].metric[k]);}}complete(f);
    f=fopen((saved+"/full_generating_thermal_r"+suffix+".bin").c_str(),"rb");header(f,std::array<uint64_t,6>{0x4d46544845524d31ULL,1,nc,sizeof(R),uint64_t(rank),1});double law[4];material_full::read(f,law,4);need(law[0]>1&&law[1]>0&&law[2]==.047&&law[3]==0,"explicit fresh constant-pressure generating-context import");std::vector<std::array<R,2>>delta(nc);material_full::read(f,delta.data(),nc);complete(f);R gamma=law[0],P=law[1];
    auto geometry=gp::geometry(m.f);std::vector<std::array<R,2>>ext(nc);
    for(const auto&e:geometry.edge){uint64_t id=(m.global[e.reference.poly[1]]-m.global_vertices)/2;bool external=m.external_edge.at(id);if(external)need(e.side.size()==1,"true global external side");for(const auto&side:e.side){const auto&t=m.f.triangle[side.triangle];for(int j=0;j<e.n;j++){int local=std::find(t.id,t.id+t.n,e.id[j])-t.id;need(local<t.n,"saved ledger canonical element identity");for(int d=0;d<2;d++){R fp=R(side.sign)*P*e.H[j][d];force[side.triangle][d*ST_MAX+local]+=fp;if(external){R traction=-fp;force[side.triangle][d*ST_MAX+local]+=traction;const auto&c=m.f.coefficient[e.id[j]];ext[side.cell][0]+=R(c.v[d])*traction;ext[side.cell][1]+=(d==1?R(c.x[0]):-R(c.x[1]))*traction;}}}}}
    // Each element owns its residual; scatter below is serial and deterministic.
    // The cell ledger stores work, torque, CM force, Iraw-dot, actual M work/Ldot.
    std::vector<std::array<R,8>>element(nt);R load=0,bubble=0;
    #pragma omp parallel for schedule(static) reduction(max:load,bubble)
    for(size_t s=0;s<nt;s++){const auto&t=m.f.triangle[s];const R*M=mass.data()+s*ST_MAX*ST_MAX;auto&v=element[s];for(int j=0;j<t.n;j++){const auto&c=m.f.coefficient[t.id[j]];for(int d=0;d<2;d++){R ma=0,mv=0;for(int k=0;k<t.n;k++){ma+=M[j*ST_MAX+k]*a[d][t.id[k]];mv+=M[j*ST_MAX+k]*m.f.coefficient[t.id[k]].v[d];}R r=ma-force[s][d*ST_MAX+j];load=std::max(load,fabsl(force[s][d*ST_MAX+j]));force[s][d*ST_MAX+j]=r;if(j>=t.bubble_offset)bubble=std::max(bubble,fabsl(r));v[0]+=R(c.v[d])*r;v[1]+=(d==1?R(c.x[0]):-R(c.x[1]))*r;if(m.f.dof[t.id[j]].kind<3)v[2+d]+=ma;v[4]+=2*R(c.x[d])*mv;v[5]+=R(c.v[d])*ma;v[6]+=(d==1?R(c.x[0]):-R(c.x[1]))*ma;}}}
    std::vector<std::array<R,8>>cell(nc);std::array<md::Vec,2>reaction={md::Vec(nd),md::Vec(nd)};
    for(size_t s=0;s<nt;s++){const auto&t=m.f.triangle[s];for(int k=0;k<8;k++)cell[t.cell][k]+=element[s][k];for(int d=0;d<2;d++)for(int j=0;j<t.n;j++)reaction[d][t.id[j]]+=force[s][d*ST_MAX+j];}
    R shared=0;for(auto&v:reaction){layout.sum(v);for(auto i:layout.owned)shared=std::max(shared,fabsl(v[i]));}R scales[3]={load,bubble,shared};MPI_Allreduce(MPI_IN_PLACE,scales,3,MPI_LONG_DOUBLE,MPI_MAX,comm);need(scales[0]>0,"nonzero initial load scale");
    FILE*csv=fopen(("full_initial_local_ledger_r"+suffix+".csv").c_str(),"wx");need(csv,"new full initial cell/rotation ledger");fprintf(csv,"cell,U,K,reaction_work,reaction_torque,external_work,external_torque,energy_error_over_U,angular_error_over_2K,resolved_spin,inertia,omega,Erot,spin_dot,Erot_dot,remaining_Kdot,mass_Kdot_minus_point,mass_Ldot_minus_point\n");
    R high[5]={0,0,0,0,0},total[5]={};
    for(size_t c=0;c<nc;c++){const auto&i=initial[c];const auto&p=point[c];const auto&r=cell[c];R U=P*geometry.volume[c]/(gamma-1),K=i[5],udot=-P*geometry.vdot[c],e=p[5]+udot-p[7]-ext[c][0]-r[0],l=p[6]-p[8]-ext[c][1]-r[1];need(U>0&&K>0&&i[0]>0,"positive source energies/mass");R rx=i[10]/i[0],ry=i[11]/i[0],vx=i[12]/i[0],vy=i[13]/i[0],spin=i[6]-(rx*i[13]-ry*i[12]),I=i[14]-i[0]*(rx*rx+ry*ry),Idot=r[4]-2*(rx*i[12]+ry*i[13]),Sdot=p[6]-(rx*r[3]-ry*r[2]);need(I>0&&std::isfinite(I),"positive resolved cell internal inertia");R omega=spin/I,Er=spin*spin/(2*I),Erdot=omega*Sdot-.5L*omega*omega*Idot,remaining=p[5]-vx*r[2]-vy*r[3]-Erdot;
        high[0]=std::max(high[0],fabsl(e)/U);high[1]=std::max(high[1],fabsl(l)/(2*K));high[2]=std::max(high[2],fabsl((geometry.volume[c]-i[2])/i[1]-delta[c][0]));high[3]=std::max(high[3],fabsl(r[5]-p[5])/U);high[4]=std::max(high[4],fabsl(r[6]-p[6])/(2*K));total[0]+=e;total[1]+=l;total[2]+=U;total[3]+=K;total[4]+=Er;
        fprintf(csv,"%llu,%.21Lg,%.21Lg,%.21Lg,%.21Lg,%.21Lg,%.21Lg,%.21Lg,%.21Lg,%.21Lg,%.21Lg,%.21Lg,%.21Lg,%.21Lg,%.21Lg,%.21Lg,%.21Lg,%.21Lg\n",(unsigned long long)m.f.cell[c],U,K,r[0],r[1],ext[c][0],ext[c][1],fabsl(e)/U,fabsl(l)/(2*K),spin,I,omega,Er,Sdot,Erdot,remaining,r[5]-p[5],r[6]-p[6]);
    }material_full::durable(csv);MPI_Allreduce(MPI_IN_PLACE,high,5,MPI_LONG_DOUBLE,MPI_MAX,comm);MPI_Allreduce(MPI_IN_PLACE,total,5,MPI_LONG_DOUBLE,MPI_SUM,comm);bool pass=high[0]<1e-8L&&high[1]<1e-11L&&scales[1]/scales[0]<1e-11L&&scales[2]/scales[0]<1e-11L&&high[2]==0;
    if(rank==0){FILE*out=fopen("full_initial_local_ledger_result.json","wx");need(out,"new initial ledger report");fprintf(out,"{\"scope\":\"all full initial cells; saved fixed M/element gravity/actual acceleration/point-work, canonical pressure and true external traction; no queries or mass solve replay\",\"global_scalar_dofs\":%llu,\"all_owner_ghost_coefficient_bits_match_including_signed_zero\":true,\"global_ID_order_checked_without_unsigned_difference\":true,\"every_cell_energy_rate_error_over_U\":%.21Lg,\"every_cell_angular_rate_error_over_2K\":%.21Lg,\"zero_trace_bubble_reaction_relative\":%.21Lg,\"global_shared_reaction_residual_relative\":%.21Lg,\"generating_delta0_restore_difference\":%.21Lg,\"every_cell_mass_vs_point_Kdot_error_over_U\":%.21Lg,\"every_cell_mass_vs_point_Ldot_error_over_2K\":%.21Lg,\"physical_U\":%.21Lg,\"physical_K\":%.21Lg,\"resolved_rotational_kinetic_energy\":%.21Lg,\"rotation_is_diagnostic_decomposition_not_energy_reservoir\":true,\"imported_thermal_law\":\"FRESH_CONSTANT_PRESSURE_GENERATING_CONTEXT\",\"seconds\":%.9g,\"cell_local_hybrid_rate_audit_complete\":%s,\"accepted_steps\":0,\"long_Kepler_validated\":false}\n",(unsigned long long)m.global_dofs,high[0],high[1],scales[1]/scales[0],scales[2]/scales[0],high[2],high[3],high[4],total[2],total[3],total[4],MPI_Wtime()-begin,pass?"true":"false");material_full::durable(out);}return pass?0:2;
}
int main(int argc,char**argv){int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);int rc=2;try{need(argc==3&&provided>=MPI_THREAD_FUNNELED,"initial ledger saved-state and curvature directories");rc=run(argv[1],argv[2],MPI_COMM_WORLD);}catch(const std::exception&e){int rank;MPI_Comm_rank(MPI_COMM_WORLD,&rank);fprintf(stderr,"rank%d FULL_INITIAL_LOCAL_LEDGER_FATAL %s\n",rank,e.what());MPI_Abort(MPI_COMM_WORLD,2);}MPI_Finalize();return rc;}
