// Assemble T^T M_old T from retained native binary80 blocks. No quadrature,
// heat repair, physical force evaluation, or time evolution is performed.
#include <mpi.h>
#include <omp.h>
#include <algorithm>
#include <array>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>
#include <unistd.h>
static void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
using R=long double;using U=uint64_t;
#include "material_full_owned.h"
#include "material_distributed_mass.h"
#include "material_compact_injection_io.h"
namespace md=material_distributed;
struct Writer {
    FILE*f;U checksum=14695981039346656037ULL;
    explicit Writer(const std::string&path):f(fopen(path.c_str(),"wbx")){need(f,"new restricted mass output");}
    template<class T>void put(const T*p,size_t n){need(fwrite(p,sizeof(T),n,f)==n,"restricted mass write");checksum=compact_io::hash(p,n*sizeof(T),checksum);}
    void finish(){need(fwrite(&checksum,sizeof(checksum),1,f)==1,"restricted mass checksum");material_full::durable(f);f=nullptr;}
    ~Writer(){if(f)fclose(f);}
};
struct Group {size_t begin=0,count=0;std::vector<uint32_t>column;};
struct Record {
    std::array<R,ST_MAX*ST_MAX>matrix{};
    std::array<R,2*ST_MAX>force{};
    // old K, new K, old L, new L, old work, new work, old torque, new torque,
    // old-measure squared velocity and position differences.
    std::array<R,10>metric{};
};
static Group group_for(const material_full::Owned&m,const compact_io::Injection&t,size_t begin,size_t count){
    Group out;out.begin=begin;out.count=count;
    for(size_t s=begin;s<begin+count;s++)for(int a=0;a<m.f.triangle[s].n;a++){
        size_t i=m.f.triangle[s].id[a];for(U k=t.row[i];k<t.row[i+1];k++)if(t.weight[k]!=0)out.column.push_back(t.column[k]);
    }
    std::sort(out.column.begin(),out.column.end());out.column.erase(std::unique(out.column.begin(),out.column.end()),out.column.end());return out;
}
static Record restrict_group(const material_full::Owned&m,const compact_io::Injection&t,const std::vector<R>&mass,
        const std::array<std::vector<R>,2>&acceleration,const Group&g){
    Record out;int np=g.column.size(),n=np+2;need(np>=10&&np<=13&&n<=ST_MAX,"bounded actual compact element support");
    for(size_t s=g.begin;s<g.begin+g.count;s++){
        const auto&triangle=m.f.triangle[s];const R*M=mass.data()+s*ST_MAX*ST_MAX;
        R S[ST_MAX][ST_MAX]={},MS[ST_MAX][ST_MAX]={},load[2][ST_MAX]={};
        R velocity_error[2][ST_MAX]={},position_error[2][ST_MAX]={};
        for(int a=0;a<triangle.n;a++){
            size_t i=triangle.id[a];
            for(U k=t.row[i];k<t.row[i+1];k++){
                auto at=std::lower_bound(g.column.begin(),g.column.end(),t.column[k]);need(at!=g.column.end()&&*at==t.column[k],"actual group injection support");S[a][at-g.column.begin()]+=t.weight[k];
            }
            for(int d=0;d<2;d++)S[a][np+d]=t.residual[d][i];
            for(int d=0;d<2;d++){
                R x=0,v=0;
                for(int k=0;k<n;k++){size_t c=k<np?g.column[k]:t.border[k-np];x+=S[a][k]*t.initial[c].x[d];v+=S[a][k]*t.initial[c].v[d];}
                velocity_error[d][a]=v-m.f.coefficient[i].v[d];position_error[d][a]=x-m.f.coefficient[i].x[d];
            }
        }
        for(int a=0;a<triangle.n;a++)for(int b=0;b<triangle.n;b++){
            R z=M[a*ST_MAX+b];
            for(int k=0;k<n;k++)if(S[b][k]!=0)MS[a][k]+=z*S[b][k];
            size_t i=triangle.id[a],j=triangle.id[b];
            for(int d=0;d<2;d++){
                load[d][a]+=z*acceleration[d][j];
                out.metric[0]+=.5L*R(m.f.coefficient[i].v[d])*z*m.f.coefficient[j].v[d];
                out.metric[2]+=(d==1?R(m.f.coefficient[i].x[0]):-R(m.f.coefficient[i].x[1]))*z*m.f.coefficient[j].v[d];
                out.metric[8]+=velocity_error[d][a]*z*velocity_error[d][b];
                out.metric[9]+=position_error[d][a]*z*position_error[d][b];
            }
        }
        for(int a=0;a<triangle.n;a++){
            size_t i=triangle.id[a];
            for(int d=0;d<2;d++){
                out.metric[4]+=R(m.f.coefficient[i].v[d])*load[d][a];
                out.metric[6]+=(d==1?R(m.f.coefficient[i].x[0]):-R(m.f.coefficient[i].x[1]))*load[d][a];
            }
            for(int k=0;k<n;k++){
                for(int j=k;j<n;j++)out.matrix[k*ST_MAX+j]+=S[a][k]*MS[a][j];
                for(int d=0;d<2;d++)out.force[d*ST_MAX+k]+=S[a][k]*load[d][a];
            }
        }
    }
    for(int k=0;k<n;k++)for(int j=0;j<k;j++)out.matrix[k*ST_MAX+j]=out.matrix[j*ST_MAX+k];
    for(int k=0;k<n;k++){
        need(out.matrix[k*ST_MAX+k]>0&&std::isfinite(out.matrix[k*ST_MAX+k]),"positive restricted element diagonal (not an SPD proof)");
        size_t i=k<np?g.column[k]:t.border[k-np];
        for(int d=0;d<2;d++){
            out.metric[5]+=R(t.initial[i].v[d])*out.force[d*ST_MAX+k];
            out.metric[7]+=(d==1?R(t.initial[i].x[0]):-R(t.initial[i].x[1]))*out.force[d*ST_MAX+k];
        }
        for(int j=0;j<n;j++){
            size_t c=j<np?g.column[j]:t.border[j-np];R z=out.matrix[k*ST_MAX+j];
            for(int d=0;d<2;d++){
                out.metric[1]+=.5L*R(t.initial[i].v[d])*z*t.initial[c].v[d];
                out.metric[3]+=(d==1?R(t.initial[i].x[0]):-R(t.initial[i].x[1]))*z*t.initial[c].v[d];
            }
        }
    }
    return out;
}
int main(int argc,char**argv){
    MPI_Init(&argc,&argv);int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    try{
        need(argc==4&&ranks==4,"usage MPI4: assembler saved_initial_dir injection_run new_output_dir");double begin=MPI_Wtime();std::string suffix=std::to_string(rank);
        compact_io::Injection t(std::string(argv[2])+"/compact_injection_r"+suffix+".bin",rank,ranks);
        auto m=material_full::load_basis(std::string(argv[1])+"/full_conditioned_basis_r"+suffix+".bin");
        need(t.fine_global==m.global&&t.old_total==m.global_dofs,"same retained conditioned field");
        md::Ownership old_layout(m.global,m.owner,MPI_COMM_WORLD),layout(t.global,t.owner,MPI_COMM_WORLD);
        need(layout.global_owned==1581183,"same approved compact ownership");
        FILE*f=fopen((std::string(argv[1])+"/full_fixed_mass_r"+suffix+".bin").c_str(),"rb");need(f,"authoritative old mass");U h[4];material_full::read(f,h,4);
        need(h[0]==0x4d464d4153533031ULL&&h[1]==1&&h[2]==m.f.triangle.size()*ST_MAX*ST_MAX&&h[3]==sizeof(R),"original binary80 mass ABI");
        std::vector<R>mass(h[2]);material_full::read(f,mass.data(),mass.size());need(fgetc(f)==EOF&&fclose(f)==0,"complete original mass");
        f=fopen((std::string(argv[1])+"/full_initial_acceleration_r"+suffix+".bin").c_str(),"rb");need(f,"retained original solved acceleration");material_full::read(f,h,4);
        need(h[0]==0x4d46414343454c31ULL&&h[1]==1&&h[2]==m.global.size()&&h[3]==sizeof(R),"retained acceleration ABI");
        std::array<std::vector<R>,2>acceleration;for(auto&a:acceleration){a.resize(m.global.size());material_full::read(f,a.data(),a.size());auto before=a;old_layout.pull(a);need(a==before,"original acceleration owner/ghost bitmatch");}
        need(fgetc(f)==EOF&&fclose(f)==0,"complete retained acceleration");
        need(m.f.triangle.size()%4==0,"owned whole parent groups");std::vector<Group>groups;U coarse=0,fine=0,maxpoly=0;
        for(size_t s=0;s<m.f.triangle.size();s+=4){
            need(m.source_triangle[s]%4==0,"parent child start");for(int k=1;k<4;k++)need(m.source_triangle[s+k]==m.source_triangle[s]+U(k),"parent siblings share owner");
            Group parent=group_for(m,t,s,4);
            if(parent.column.size()==10){coarse++;groups.push_back(std::move(parent));}
            else {need(parent.column.size()>=25&&parent.column.size()<=28,"fixed protected-core polynomial footprint");for(int k=0;k<4;k++){Group child=group_for(m,t,s+k,1);maxpoly=std::max<U>(maxpoly,child.column.size());fine++;groups.push_back(std::move(child));}}
        }
        Writer out(std::string(argv[3])+"/compact_mass_elements_r"+suffix+".bin");
        U header[8]={0x4d434d4153533031ULL,1,U(rank),U(ranks),U(groups.size()),t.total,sizeof(R),t.old_total};out.put(header,8);
        f=fopen((std::string(argv[3])+"/progress_r"+suffix+".tsv").c_str(),"wx");need(f,"new native assembly progress");fprintf(f,"groups_done\tgroups_total\tseconds\n");
        std::array<R,10>metric{};R border[4]={};std::array<md::Vec,2>load={md::Vec(t.global.size()),md::Vec(t.global.size())};
        const size_t batch=1024;std::vector<Record>records(batch);
        for(size_t at=0;at<groups.size();at+=batch){
            size_t count=std::min(batch,groups.size()-at);int failures=0;
            #pragma omp parallel for schedule(static) reduction(+:failures)
            for(size_t j=0;j<count;j++)try{records[j]=restrict_group(m,t,mass,acceleration,groups[at+j]);}catch(...){failures++;}
            need(!failures,"restricted mass element construction; no clipping or factor repair");
            for(size_t j=0;j<count;j++){
                const auto&g=groups[at+j];const auto&r=records[j];U np=g.column.size(),n=np+2;
                U description[3]={m.source_triangle[g.begin],g.count,np};out.put(description,3);out.put(g.column.data(),g.column.size());
                for(U k=0;k<n;k++)out.put(r.matrix.data()+k*ST_MAX,n);
                for(int d=0;d<2;d++)out.put(r.force.data()+d*ST_MAX,n);
                for(int k=0;k<10;k++)metric[k]+=r.metric[k];
                for(int a=0;a<2;a++)for(int b=0;b<2;b++)border[2*a+b]+=r.matrix[(np+a)*ST_MAX+np+b];
                for(U k=0;k<n;k++){size_t i=k<np?g.column[k]:t.border[k-np];for(int d=0;d<2;d++)load[d][i]+=r.force[d*ST_MAX+k];}
            }
            fprintf(f,"%zu\t%zu\t%.9g\n",at+count,groups.size(),MPI_Wtime()-begin);fflush(f);
        }
        out.finish();material_full::durable(f);for(auto&v:load)layout.sum(v);
        Writer rhs(std::string(argv[3])+"/compact_retained_rhs_r"+suffix+".bin");U rh[6]={0x4d43524853303031ULL,1,U(rank),U(ranks),U(t.global.size()),sizeof(R)};rhs.put(rh,6);for(auto&v:load)rhs.put(v.data(),v.size());rhs.finish();
        MPI_Allreduce(MPI_IN_PLACE,metric.data(),10,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);MPI_Allreduce(MPI_IN_PLACE,border,4,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
        U counts[4]={U(groups.size()),coarse,fine,maxpoly};MPI_Allreduce(MPI_IN_PLACE,counts,3,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);MPI_Allreduce(MPI_IN_PLACE,&counts[3],1,MPI_UINT64_T,MPI_MAX,MPI_COMM_WORLD);
        need(counts[0]==351112&&counts[1]==262440&&counts[2]==88672&&counts[3]<=13,"actual fixed candidate element count");
        const R U0=.000863999588012891662344L;R kerror=(metric[1]-metric[0])/U0,lerror=(metric[3]-metric[2])/(2*metric[0]);
        R work=(metric[5]-metric[4])/U0,torque=(metric[7]-metric[6])/(2*metric[0]);
        need(fabsl(kerror)<1e-8L&&fabsl(lerror)<1e-11L&&fabsl(work)<1e-8L&&fabsl(torque)<1e-11L,"original cold K/L/work/torque assembly comparison");
        need(border[0]>0&&border[3]>0&&border[0]*border[3]-border[1]*border[2]>0,"positive raw global residual Gram (not full Schur proof)");
        if(rank==0){
            FILE*result=fopen((std::string(argv[3])+"/mass_restriction_result.json").c_str(),"wx");need(result,"new mass congruence report");
            fprintf(result,"{\"scope\":\"actual saved binary80 M congruence, grouped compact element blocks and weak restriction of M_old times retained acceleration; no new quadrature, solve or evolution\",\"compact_scalar_dofs\":%llu,\"compact_elements\":%llu,\"coarse_parent_elements\":%llu,\"protected_fine_elements\":%llu,\"maximum_polynomial_element_support\":%llu,\"K_old\":%.21Lg,\"K_new\":%.21Lg,\"K_difference_over_U\":%.21Lg,\"L_difference_over_2K\":%.21Lg,\"work_difference_over_U\":%.21Lg,\"torque_difference_over_2K\":%.21Lg,\"old_measure_velocity_difference_squared\":%.21Lg,\"old_measure_position_difference_squared\":%.21Lg,\"raw_residual_Gram\":[%.21Lg,%.21Lg,%.21Lg,%.21Lg],\"suggested_border_norms\":[%.21Lg,%.21Lg],\"force_reference_inherits_original_acceleration_solve_residual\":true,\"full_Schur_positivity_and_true_residual_solve_pending\":true,\"physical_pressure_curvature_and_heat_gates_pending\":true,\"assembly_cold_comparison_pass\":true,\"seconds\":%.9g,\"new_accepted_step\":false,\"long_Kepler_validated\":false}\n",(unsigned long long)t.total,(unsigned long long)counts[0],(unsigned long long)counts[1],(unsigned long long)counts[2],(unsigned long long)counts[3],metric[0],metric[1],kerror,lerror,work,torque,metric[8],metric[9],border[0],border[1],border[2],border[3],sqrtl(border[0]),sqrtl(border[3]),MPI_Wtime()-begin);material_full::durable(result);
        }
    }catch(const std::exception&e){fprintf(stderr,"rank%d COMPACT_MASS_RESTRICTION_REJECTED %s\n",rank,e.what());MPI_Abort(MPI_COMM_WORLD,2);return 2;}
    MPI_Finalize();return 0;
}
