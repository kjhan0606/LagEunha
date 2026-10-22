// Exact algebra on retained original binary80 M. No new quadrature or evolution.
#include <mpi.h>
#include <omp.h>
#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>
#include <unistd.h>
using U=uint64_t;using R=long double;
static void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
#include "material_full_owned.h"
#include "material_distributed_mass.h"
#include "material_cell_partition_io.h"
namespace cc=cell_condensation;namespace md=material_distributed;
struct Writer {
    FILE*f;U hash=14695981039346656037ULL;
    explicit Writer(const std::string&p):f(fopen(p.c_str(),"wbx")){need(f,"new original cell factors");}
    template<class T>void put(const T*p,size_t n){need(fwrite(p,sizeof(T),n,f)==n,"cell factor write");auto*b=reinterpret_cast<const unsigned char*>(p);for(size_t j=0;j<n*sizeof(T);j++){hash^=b[j];hash*=1099511628211ULL;}}
    void finish(){need(fwrite(&hash,sizeof(U),1,f)==1,"cell factor checksum");material_full::durable(f);f=nullptr;}
    ~Writer(){if(f)fclose(f);}
};
static size_t packed(size_t i,size_t j){return i*(i+1)/2+j;}
struct Result {std::vector<R>L,IG,S;R residual=0,min_relative_pivot=1;};
static Result factor_cell(const material_full::Owned&m,const cc::Partition&p,const cc::Cell&c,const std::vector<R>&mass){
    size_t ni=c.I.size(),ng=c.G.size(),n=ni+ng;std::vector<R>A(n*n);Result out;out.L.resize(ni*(ni+1)/2);out.IG.resize(ni*ng);out.S.resize(ng*ng);
    for(U s=c.begin;s<c.begin+c.count;s++){
        const auto&t=m.f.triangle[s];need(m.source_cell[s]==c.source,"same original physical cell during M assembly");int local[ST_MAX];
        for(int a=0;a<t.n;a++){
            uint32_t id=t.id[a];auto at=std::lower_bound(c.I.begin(),c.I.end(),id);
            if(at!=c.I.end()&&*at==id)local[a]=at-c.I.begin();
            else{auto old=std::lower_bound(p.old_G.begin(),p.old_G.end(),id);need(old!=p.old_G.end()&&*old==id,"trace coefficient in original skeleton");uint32_t g=old-p.old_G.begin();auto pos=std::lower_bound(c.G.begin(),c.G.end(),g);need(pos!=c.G.end()&&*pos==g,"source trace coefficient support");local[a]=ni+(pos-c.G.begin());}
        }
        const R*M=mass.data()+s*ST_MAX*ST_MAX;
        for(int a=0;a<t.n;a++)for(int b=0;b<t.n;b++){R z=M[a*ST_MAX+b];need(std::isfinite(z)&&z==M[b*ST_MAX+a],"finite symmetric authoritative original block");A[size_t(local[a])*n+local[b]]+=z;}
    }
    R scale=0;for(size_t i=0;i<ni;i++)scale=std::max(scale,A[i*n+i]);need(scale>0,"positive source interior scale");
    for(size_t i=0;i<ni;i++)for(size_t j=0;j<=i;j++){
        R value=A[i*n+j];need(value==A[j*n+i],"assembled exact interior symmetry");for(size_t k=0;k<j;k++)value-=out.L[packed(i,k)]*out.L[packed(j,k)];
        if(i==j){if(!(value>0&&std::isfinite(value))){fprintf(stderr,"CELL_II_NONPOSITIVE source=%llu pivot=%zu value=%.21Lg scale=%.21Lg\n",(unsigned long long)c.source,i,value,scale);need(false,"source interior Cholesky; no repair");}out.min_relative_pivot=std::min(out.min_relative_pivot,value/scale);out.L[packed(i,j)]=sqrtl(value);}
        else out.L[packed(i,j)]=value/out.L[packed(j,j)];
    }
    R error=0,norm=0;
    for(size_t i=0;i<ni;i++)for(size_t j=0;j<=i;j++){
        R product=0;for(size_t k=0;k<=j;k++)product+=out.L[packed(i,k)]*out.L[packed(j,k)];R z=A[i*n+j],d=product-z,weight=i==j?1:2;error+=weight*d*d;norm+=weight*z*z;
    }
    out.residual=sqrtl(error/norm);need(std::isfinite(out.residual)&&out.residual<1e-16L,"all original cell II factor residuals");
    // Y=L^-1 M_IG. Store M_IG, not an explicit M_II inverse.
    std::vector<R>Y(ni*ng);
    for(size_t i=0;i<ni;i++)for(size_t g=0;g<ng;g++){
        R z=A[i*n+ni+g];need(z==A[(ni+g)*n+i],"same original IG/GI block");out.IG[i*ng+g]=z;
        for(size_t k=0;k<i;k++)z-=out.L[packed(i,k)]*Y[k*ng+g];Y[i*ng+g]=z/out.L[packed(i,i)];
    }
    for(size_t a=0;a<ng;a++)for(size_t b=0;b<=a;b++){
        R z=A[(ni+a)*n+ni+b];need(z==A[(ni+b)*n+ni+a],"same original GG symmetry");for(size_t i=0;i<ni;i++)z-=Y[i*ng+a]*Y[i*ng+b];need(std::isfinite(z),"finite exact cell Schur complement");out.S[a*ng+b]=out.S[b*ng+a]=z;
    }
    std::vector<R>Ls(ng*(ng+1)/2);
    for(size_t i=0;i<ng;i++)for(size_t j=0;j<=i;j++){
        R z=out.S[i*ng+j];for(size_t k=0;k<j;k++)z-=Ls[packed(i,k)]*Ls[packed(j,k)];
        if(i==j){if(!(z>0&&std::isfinite(z))){fprintf(stderr,"CELL_SCHUR_NONPOSITIVE source=%llu pivot=%zu value=%.21Lg\n",(unsigned long long)c.source,i,z);need(false,"original cell Schur Cholesky; no repair");}Ls[packed(i,j)]=sqrtl(z);}else Ls[packed(i,j)]=z/Ls[packed(j,j)];
    }
    return out;
}
int main(int argc,char**argv){
    MPI_Init(&argc,&argv);int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    try{
        need(argc==5&&ranks==4,"MPI4 usage: factor original_initial structure_run verified_capacity_json new_output");
        double begin=MPI_Wtime();std::ifstream input(argv[3]);need(bool(input),"measured H200 hardware capacity evidence for CPU-only setup");std::string capacity((std::istreambuf_iterator<char>(input)),{});
        bool allocated=capacity.find("\"planned_allocation_pass\":true")!=std::string::npos;
        bool hardware_bound=capacity.find("\"measured_hardware_capacity_bound_pass\":true")!=std::string::npos;
        need((allocated||hardware_bound)&&capacity.find("\"twenty_percent_headroom\":true")!=std::string::npos&&capacity.find("\"planned_all_ranks_bytes\":106707056260")!=std::string::npos,"measured fixed layout hardware capacity with twenty percent headroom; fresh allocation mandatory before GPU execution");
        std::string suffix=std::to_string(rank),dest=argv[4];cc::Partition p(std::string(argv[2])+"/cell_partition_r"+suffix+".bin",rank,ranks);
        auto m=material_full::load_basis(std::string(argv[1])+"/full_conditioned_basis_r"+suffix+".bin");need(m.global.size()==p.nf&&m.f.triangle.size()==p.nt,"same unchanged original basis");
        for(size_t i=0;i<p.ng;i++)need(m.global[p.old_G[i]]==p.global_G[i]&&m.owner[p.old_G[i]]==p.owner_G[i],"original skeleton provenance and ownership");
        md::Ownership layout(p.global_G,p.owner_G,MPI_COMM_WORLD);need(layout.global_owned==p.total,"unchanged unique skeleton owners");
        FILE*f=fopen((std::string(argv[1])+"/full_fixed_mass_r"+suffix+".bin").c_str(),"rb");need(f,"original retained binary80 mass");U h[4];material_full::read(f,h,4);need(h[0]==0x4d464d4153533031ULL&&h[1]==1&&h[2]==p.nt*ST_MAX*ST_MAX&&h[3]==sizeof(R),"authoritative original mass ABI");
        std::vector<R>mass(h[2]);material_full::read(f,mass.data(),mass.size());need(fgetc(f)==EOF&&fclose(f)==0,"complete authoritative original mass");
        Writer factors(dest+"/cell_factors_r"+suffix+".bin");U header[9]={0x4d43464143543031ULL,1,U(rank),U(ranks),U(p.cell.size()),p.nf,p.ng,p.total,sizeof(R)};factors.put(header,9);
        std::vector<R>CSR(p.column.size());R worst=0,minpivot=1;U L_count=0,IG_count=0;
        f=fopen((dest+"/progress_r"+suffix+".tsv").c_str(),"wx");need(f,"new bounded factor progress");fprintf(f,"cells_done\tcells_total\tseconds\tmax_factor_relative_residual\n");
        const size_t batch=128;std::vector<Result>result(batch);
        for(size_t at=0;at<p.cell.size();at+=batch){
            size_t count=std::min(batch,p.cell.size()-at);int failed=0;
            #pragma omp parallel for schedule(dynamic,1) reduction(+:failed)
            for(size_t j=0;j<count;j++)try{result[j]=factor_cell(m,p,p.cell[at+j],mass);}catch(const std::exception&e){fprintf(stderr,"rank%d source=%llu FACTOR_FAIL %s\n",rank,(unsigned long long)p.cell[at+j].source,e.what());failed++;}
            need(!failed,"original cell factors positive and accurate; no repairs");
            for(size_t j=0;j<count;j++){
                const auto&c=p.cell[at+j];const auto&r=result[j];U desc[5]={c.source,c.begin,c.count,U(c.I.size()),U(c.G.size())};factors.put(desc,5);factors.put(r.L.data(),r.L.size());factors.put(r.IG.data(),r.IG.size());
                L_count+=r.L.size();IG_count+=r.IG.size();worst=std::max(worst,r.residual);minpivot=std::min(minpivot,r.min_relative_pivot);
                for(size_t a=0;a<c.G.size();a++){uint32_t i=c.G[a];for(size_t b=0;b<c.G.size();b++){uint32_t k=c.G[b];auto it=std::lower_bound(p.column.begin()+p.row[i],p.column.begin()+p.row[i+1],k);need(it!=p.column.begin()+p.row[i+1]&&*it==k,"preallocated exact Schur support");CSR[it-p.column.begin()]+=r.S[a*c.G.size()+b];}}
            }
            fprintf(f,"%zu\t%zu\t%.9g\t%.21Lg\n",at+count,p.cell.size(),MPI_Wtime()-begin,worst);fflush(f);need(MPI_Wtime()-begin<8*3600,"fixed eight-hour cell factor setup budget");
        }
        factors.finish();material_full::durable(f);
        md::Vec diagonal(p.ng);
        for(size_t i=0;i<p.ng;i++){
            for(U j=p.row[i];j<p.row[i+1];j++){
                uint32_t k=p.column[j];need(std::isfinite(CSR[j]),"finite assembled original Schur entry");auto other=std::lower_bound(p.column.begin()+p.row[k],p.column.begin()+p.row[k+1],uint32_t(i));need(other!=p.column.begin()+p.row[k+1]&&*other==i&&CSR[other-p.column.begin()]==CSR[j],"exact assembled partial Schur symmetry");if(k==i)diagonal[i]=CSR[j];
            }
        }
        layout.sum(diagonal);for(R z:diagonal)need(z>0&&std::isfinite(z),"positive assembled original skeleton Jacobi");
        Writer matrix(dest+"/skeleton_mass_r"+suffix+".bin");U mh[8]={0x4d43534353523031ULL,1,U(rank),U(ranks),p.ng,U(p.column.size()),p.total,sizeof(R)};matrix.put(mh,8);matrix.put(p.row.data(),p.row.size());matrix.put(p.column.data(),p.column.size());matrix.put(CSR.data(),CSR.size());matrix.put(diagonal.data(),diagonal.size());matrix.finish();
        MPI_Allreduce(MPI_IN_PLACE,&worst,1,MPI_LONG_DOUBLE,MPI_MAX,MPI_COMM_WORLD);MPI_Allreduce(MPI_IN_PLACE,&minpivot,1,MPI_LONG_DOUBLE,MPI_MIN,MPI_COMM_WORLD);U counts[2]={L_count,IG_count};MPI_Allreduce(MPI_IN_PLACE,counts,2,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
        if(rank==0){need(counts[0]==635433712&&counts[1]==379124224,"same exact original source factors and IG counts");f=fopen((dest+"/factor_result.json").c_str(),"wx");need(f,"new original exact condensation report");
            fprintf(f,"{\"scope\":\"all original binary80 source II Cholesky and IG blocks, positive source Schur complements and assembled partial skeleton CSR; no new quadrature or physical evolution\",\"original_scalar_dofs\":%llu,\"skeleton_scalar_dofs\":%llu,\"cells\":71080,\"II_packed_entries\":%llu,\"IG_entries\":%llu,\"max_cell_II_relative_factor_residual\":%.21Lg,\"min_cell_II_pivot_over_cell_diagonal_scale\":%.21Lg,\"all_cell_II_and_Schur_positive\":true,\"assembled_symmetric\":true,\"fresh_GPU_allocation_pending\":%s,\"full_original_equation_residual_pending\":true,\"DD_GPU_recovery_pending\":true,\"cache_implementation_pending\":true,\"physical_trajectory_advanced\":false,\"long_Kepler_validated\":false,\"seconds\":%.9g}\n",(unsigned long long)p.old_total,(unsigned long long)p.total,(unsigned long long)counts[0],(unsigned long long)counts[1],worst,minpivot,allocated?"false":"true",MPI_Wtime()-begin);material_full::durable(f);}
    }catch(const std::exception&e){fprintf(stderr,"rank%d ORIGINAL_CELL_CONDENSATION_REJECTED %s\n",rank,e.what());MPI_Abort(MPI_COMM_WORLD,2);return 2;}
    MPI_Finalize();return 0;
}
