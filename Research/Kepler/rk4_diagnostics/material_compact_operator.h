#ifndef MATERIAL_COMPACT_OPERATOR_H
#define MATERIAL_COMPACT_OPERATOR_H
// The finite-functional RHS component of the ONE compact transaction. It has
// no commit/restart/trajectory API and does not certify source/target overlap.
#include "material_compact_mass.h"
#include "material_compact_volume.h"
#include "material_compact_thermo_cuda.h"
#include <cstring>
namespace compact_operator {
using R=long double;namespace md=material_distributed;
struct Seed {
    compact_thermo::Config law;
    std::vector<compact_geometry::Patch>patch;
    std::vector<compact_p3::ThermalDensity>field;
    std::vector<STCoefficient>generating;
    std::vector<double>generating_delta;
    std::vector<uint8_t>constant;
    std::vector<compact_volume::CanonicalEdge>edge;
    compact_volume::Rule polynomial_edge_rule,curved_edge_rule;
    std::vector<uint64_t>offset;
    std::vector<compact_thermo::QuadraturePoint>quadrature;
    // All arrays are immutable after accepted projection and belong to one
    // representation epoch. Quadrature/rank/domain/error proofs are external
    // mandatory transaction gates; a Seed object is not an acceptance receipt.
};
template<class T>struct Device {
    T*ptr=nullptr;size_t n=0;
    Device()=default;Device(const Device&)=delete;Device&operator=(const Device&)=delete;
    ~Device(){cudaFree(ptr);}
    void allocate(size_t count){md::require(!ptr&&count>0,"single nonempty compact device allocation");n=count;md::require(cudaMalloc(reinterpret_cast<void**>(&ptr),n*sizeof(T))==cudaSuccess,"compact device allocation");}
    void upload(const T*source){md::require(ptr&&cudaMemcpy(ptr,source,n*sizeof(T),cudaMemcpyHostToDevice)==cudaSuccess,"compact fixed/stage device upload");}
};
struct Stage {
    std::vector<std::array<R,MF_EVOL_METRICS>>cell;
    std::array<R,MF_EVOL_METRICS>global={};
    std::vector<MFEvolvedCell>context;
    std::array<md::Vec,2>force;
    compact_volume::Geometry geometry;
    R force_work=0,force_torque=0,external_work=0,external_torque=0;
    R energy_identity=0,torque_identity=0,energy_scale=0,torque_scale=0;
    R minimum_J=0,minimum_jhat=0,maximum_face_derivative_error=0;
    uint64_t points_per_pass=0;
    double wall_seconds=0;
};
inline void check_ghosts(const md::Ownership&layout,const std::vector<STCoefficient>&c){
    compact_mass::collective_require(c.size()==layout.global.size(),layout.comm,"same current compact state/layout extent");
    bool valid=true;for(const auto&q:c)for(int k=0;k<2;k++)valid&=std::isfinite(q.x[k])&&std::isfinite(q.v[k]);
    compact_mass::collective_require(valid,layout.comm,"finite compact coefficients before MPI exchange");md::Vec values(c.size());bool same=true;
    for(int d=0;d<4;d++){
        for(size_t j=0;j<c.size();j++)values[j]=d<2?c[j].x[d]:c[j].v[d-2];layout.pull(values);
        for(size_t j=0;j<c.size();j++){double actual=double(values[j]),stored=d<2?c[j].x[d]:c[j].v[d-2];same&=memcmp(&actual,&stored,sizeof(double))==0;}
    }
    compact_mass::collective_require(same,layout.comm,"current compact owned/ghost bits agree; no hidden state replacement");
}
class Operator {
    const Seed&seed;const md::Ownership&layout;
    Device<compact_geometry::Patch>patch;
    Device<STCoefficient>current,generating;
    Device<compact_p3::ThermalDensity>field;
    Device<MFEvolvedCell>context;
    Device<uint64_t>offset;
    Device<compact_thermo::QuadraturePoint>quadrature;
    Device<MFEvolvedBlock>output;
    std::vector<MFEvolvedBlock>batch;
public:
    Operator(const Seed&s,const md::Ownership&l):seed(s),layout(l){
        bool valid=!s.patch.empty()&&s.patch.size()==s.field.size()&&s.generating.size()==l.global.size()&&s.constant.size()==l.global.size();
        valid&=!s.generating_delta.empty()&&s.offset.size()==s.patch.size()+1&&!s.quadrature.empty();
        valid&=s.law.gamma>1&&s.law.background_pressure>0&&s.law.eps>0&&std::isfinite(s.law.gamma)&&std::isfinite(s.law.background_pressure)&&std::isfinite(s.law.eps);
        if(valid){
            valid&=s.offset.front()==0&&s.offset.back()==s.quadrature.size();
            for(double x:s.generating_delta)valid&=std::isfinite(x);
            for(auto x:s.constant)valid&=x==0||x==1;
            for(size_t k=0;k<s.patch.size();k++){
                const auto&p=s.patch[k];valid&=p.polynomial.cell<s.generating_delta.size()&&p.curved<=1&&s.offset[k]<s.offset[k+1]&&s.offset[k+1]-s.offset[k]<=INT_MAX;
                valid&=compact_p3::nonnegative_controls(s.field[k].mass)&&compact_p3::nonnegative_controls(s.field[k].internal_energy);
                for(auto id:p.polynomial.id)valid&=id<l.global.size()&&(id<s.constant.size()&&s.constant[id]==1);
                if(p.curved)for(auto id:p.modes.id)valid&=id<l.global.size()&&(id<s.constant.size()&&s.constant[id]==0);
            }
            for(const auto&q:s.quadrature){double total=0;valid&=q.weight>0&&std::isfinite(q.weight);for(double x:q.l){valid&=x>=0&&std::isfinite(x);total+=x;}valid&=fabs(total-1)<=32*DBL_EPSILON;}
        }
        compact_mass::collective_require(valid,l.comm,"immutable positive compact fields and explicit finite functional");
        double low[3]={s.law.gamma,s.law.background_pressure,s.law.eps},high[3]={low[0],low[1],low[2]};
        MPI_Allreduce(MPI_IN_PLACE,low,3,MPI_DOUBLE,MPI_MIN,l.comm);MPI_Allreduce(MPI_IN_PLACE,high,3,MPI_DOUBLE,MPI_MAX,l.comm);
        compact_mass::collective_require(low[0]==high[0]&&low[1]==high[1]&&low[2]==high[2],l.comm,"same actual EOS/gravity law on all compact ranks");
        check_ghosts(l,s.generating);valid=true;
        try{
            patch.allocate(s.patch.size());patch.upload(s.patch.data());current.allocate(s.generating.size());
            generating.allocate(s.generating.size());generating.upload(s.generating.data());
            field.allocate(s.field.size());field.upload(s.field.data());context.allocate(s.generating_delta.size());
            offset.allocate(s.offset.size());offset.upload(s.offset.data());quadrature.allocate(s.quadrature.size());quadrature.upload(s.quadrature.data());
            output.allocate(std::min(size_t(MF_BATCH_TRIANGLES),s.patch.size()));batch.resize(output.n);
        }catch(const std::exception&){valid=false;}
        compact_mass::collective_require(valid,l.comm,"persistent compact CUDA operator allocation/upload");
    }
    Operator(const Operator&)=delete;Operator&operator=(const Operator&)=delete;
    Stage evaluate(const std::vector<STCoefficient>&c){
        double begin=MPI_Wtime();check_ghosts(layout,c);Stage result;bool valid=true;
        try{result.geometry=compact_volume::assemble(seed.edge,c,seed.generating_delta.size(),seed.polynomial_edge_rule,seed.curved_edge_rule,seed.constant);}catch(const std::exception&){valid=false;}
        compact_mass::collective_require(valid,layout.comm,"same-stage actual curved cell boundary functional");
        result.context.resize(seed.generating_delta.size());for(size_t k=0;k<result.context.size();k++)result.context[k].generating_delta=seed.generating_delta[k];
        result.cell.resize(result.context.size());for(auto&f:result.force)f.resize(c.size());
        valid=cudaMemcpy(current.ptr,c.data(),c.size()*sizeof(STCoefficient),cudaMemcpyHostToDevice)==cudaSuccess;
        compact_mass::collective_require(valid,layout.comm,"current compact RK state device upload");
        R minJ=std::numeric_limits<R>::infinity(),minjhat=minJ;
        for(int pass=0;pass<3;pass++){
            valid=cudaMemcpy(context.ptr,result.context.data(),context.n*sizeof(MFEvolvedCell),cudaMemcpyHostToDevice)==cudaSuccess;
            compact_mass::collective_require(valid,layout.comm,"same-stage compact cell context upload");
            for(auto&cell:result.cell)cell={};uint64_t points=0;
            for(size_t at=0;at<seed.patch.size();at+=batch.size()){
                size_t count=std::min(batch.size(),seed.patch.size()-at);
                int code=compact_material_batch_device(patch.ptr,patch.n,current.ptr,generating.ptr,current.n,field.ptr,context.ptr,context.n,offset.ptr,quadrature.ptr,quadrature.n,at,count,pass,seed.law,output.ptr,nullptr);
                if(code||cudaMemcpy(batch.data(),output.ptr,count*sizeof(MFEvolvedBlock),cudaMemcpyDeviceToHost)!=cudaSuccess){valid=false;break;}
                for(size_t j=0;j<count;j++){
                    const auto&p=seed.patch[at+j];const auto&b=batch[j];valid&=b.bad_point<0&&b.minJ>0;points+=b.points;minJ=std::min(minJ,R(b.minJ));
                    if(pass)minjhat=std::min(minjhat,R(b.min_j));
                    auto&cell=result.cell[p.polynomial.cell];
                    for(int k=0;k<MF_EVOL_METRICS;k++){R x=R(b.metric[k].hi)+b.metric[k].lo;valid&=std::isfinite(x);cell[k]+=x;}
                    if(pass==2){
                        int n=p.curved?12:10;for(int a=0;a<n;a++){
                            uint64_t id=a<10?p.polynomial.id[a]:p.modes.id[a-10];
                            for(int d=0;d<2;d++){R f=R(b.force[d][a].hi)+b.force[d][a].lo;valid&=std::isfinite(f);result.force[d][id]+=f;}
                        }
                    }
                }
            }
            valid&=points==seed.quadrature.size();
            compact_mass::collective_require(valid,layout.comm,"complete finite compact CUDA pass, positive geometry and thermal state");
            result.points_per_pass=points;
            for(size_t k=0;k<result.cell.size();k++){
                const auto&m=result.cell[k];auto&C=result.context[k];valid&=m[0]>0&&m[3]>0;
                if(pass==0&&m[0]>0){C.delta=double((result.geometry.cell[k].volume-m[1])/m[0]);C.delta_dot=double((result.geometry.cell[k].volume_dot-m[2])/m[0]);}
                if(pass==1&&m[0]>0)C.mean_deltaP=double(m[4]/m[0]);
                valid&=std::isfinite(C.delta)&&std::isfinite(C.delta_dot)&&std::isfinite(C.mean_deltaP);
            }
            compact_mass::collective_require(valid,layout.comm,"finite current compact geometric-pressure context; no heat fit");
        }
        std::vector<R>pressure(result.context.size());for(size_t k=0;k<pressure.size();k++)pressure[k]=R(seed.law.background_pressure)+result.context[k].mean_deltaP;
        std::vector<compact_volume::V2>face_force;valid=true;
        try{face_force=compact_volume::mean_pressure_force(seed.edge,result.geometry,pressure,c.size(),seed.law.background_pressure);}catch(const std::exception&){valid=false;}
        compact_mass::collective_require(valid,layout.comm,"same canonical mean-pressure and physical-exterior force");
        for(size_t j=0;j<c.size();j++)for(int d=0;d<2;d++)result.force[d][j]+=face_force[j][d];
        for(auto&f:result.force)layout.sum(f);
        for(const auto&cell:result.cell)for(int k=0;k<MF_EVOL_METRICS;k++)result.global[k]+=cell[k];
        MPI_Allreduce(MPI_IN_PLACE,result.global.data(),MF_EVOL_METRICS,MPI_LONG_DOUBLE,MPI_SUM,layout.comm);
        for(size_t j:layout.owned){
            R work=R(c[j].v[0])*result.force[0][j]+R(c[j].v[1])*result.force[1][j];
            R torque=R(c[j].x[0])*result.force[1][j]-R(c[j].x[1])*result.force[0][j];
            result.force_work+=work;result.force_torque+=torque;result.energy_scale+=fabsl(work);result.torque_scale+=fabsl(torque);
        }
        for(size_t j=0;j<seed.edge.size();j++)if(seed.edge[j].physical_exterior){
            const auto&e=seed.edge[j];R p=-R(e.side[0].sign)*seed.law.background_pressure;
            for(int a=0;a<e.basis.n;a++){
                auto id=e.basis.id[a];auto f=compact_volume::times(p,result.geometry.edge[j].H[a]);
                result.external_work+=R(c[id].v[0])*f[0]+R(c[id].v[1])*f[1];result.external_torque+=R(c[id].x[0])*f[1]-R(c[id].x[1])*f[0];
            }
        }
        R totals[6]={result.force_work,result.force_torque,result.external_work,result.external_torque,result.energy_scale,result.torque_scale};
        MPI_Allreduce(MPI_IN_PLACE,totals,6,MPI_LONG_DOUBLE,MPI_SUM,layout.comm);
        result.force_work=totals[0];result.force_torque=totals[1];result.external_work=totals[2];result.external_torque=totals[3];
        result.energy_scale=totals[4]+fabsl(result.global[6])+fabsl(result.global[11])+fabsl(result.external_work);
        result.torque_scale=totals[5]+fabsl(result.global[12])+fabsl(result.external_torque);
        result.energy_identity=result.force_work+result.global[6]-result.global[11]-result.external_work;
        result.torque_identity=result.force_torque-result.global[12]-result.external_torque;
        MPI_Allreduce(&minJ,&result.minimum_J,1,MPI_LONG_DOUBLE,MPI_MIN,layout.comm);MPI_Allreduce(&minjhat,&result.minimum_jhat,1,MPI_LONG_DOUBLE,MPI_MIN,layout.comm);
        MPI_Allreduce(&result.geometry.derivative_error,&result.maximum_face_derivative_error,1,MPI_LONG_DOUBLE,MPI_MAX,layout.comm);
        uint64_t local_points=result.points_per_pass;MPI_Allreduce(&local_points,&result.points_per_pass,1,MPI_UINT64_T,MPI_SUM,layout.comm);
        double elapsed=MPI_Wtime()-begin;MPI_Allreduce(&elapsed,&result.wall_seconds,1,MPI_DOUBLE,MPI_MAX,layout.comm);
        return result;
    }
};
} // namespace compact_operator
#endif
