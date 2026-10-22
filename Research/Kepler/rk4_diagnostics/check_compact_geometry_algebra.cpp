// Finite algebra/API fixture, NOT an accepted-source guard or physical run.
// No checkpoint, original trajectory, force oracle, or remap is read here.
#include "material_compact_geometry_cuda.h"
#include <iostream>
#include <iomanip>
#include <vector>
#include <stdexcept>
#include <algorithm>
namespace cg=compact_geometry;
void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
void cuda_need(cudaError_t e){need(e==cudaSuccess,cudaGetErrorString(e));}
template<class T>struct Device {
    T*ptr=nullptr;size_t n;
    explicit Device(size_t count):n(count){cuda_need(cudaMalloc(reinterpret_cast<void**>(&ptr),n*sizeof(T)));}
    ~Device(){cudaFree(ptr);}
    void put(const T*x){cuda_need(cudaMemcpy(ptr,x,n*sizeof(T),cudaMemcpyHostToDevice));}
    std::vector<T>get(){std::vector<T>v(n);cuda_need(cudaMemcpy(v.data(),ptr,n*sizeof(T),cudaMemcpyDeviceToHost));return v;}
};
double difference(double a,double b){return fabs(a-b)/std::max({1.,fabs(a),fabs(b)});}
int main(){try{
    cg::Patch p;double vertices[3][2]={{0,0},{1,0},{0,1}};
    for(int j=0;j<3;j++)for(int k=0;k<2;k++)p.polynomial.vertex[j][k]=vertices[j][k];
    compact_p3::geometry(p.polynomial);for(int a=0;a<10;a++)p.polynomial.id[a]=a;
    p.curved=1;p.modes.component=0;p.modes.id[0]=10;p.modes.id[1]=11;p.modes.scale[0]=.02;p.modes.scale[1]=.03;
    STEdge e={};e.rank=0;e.A[0]=1;e.D[0]=-1;e.D[1]=1;e.eps=.047;
    std::vector<STCoefficient>old(4);
    for(int a=0;a<4;a++){
        e.poly[a]=a;double t=a/3.;old[a].x[0]=1-t+.02*t*(1-t);old[a].x[1]=t+.03*t*(1-t)*(t-.5);
    }
    p.trace=cg::freeze_trace(e,old,0,1,0);cg::interpolate_chart(p);cg::check_modes(p,12);
    std::vector<STCoefficient>c(12);
    for(int a=0;a<10;a++){auto x=cg::initial_polynomial_position(p,a);for(int k=0;k<2;k++)c[a].x[k]=x[k];}
    c[10]=cg::initial_mode(p.modes,0);c[11]=cg::initial_mode(p.modes,1);
    double omega=.7;
    for(int a=0;a<12;a++){
        c[a].v[0]=-omega*c[a].x[1]+(a<10?.13:0);
        c[a].v[1]=omega*c[a].x[0]+(a<10?-.07:0);
    }
    compact_p3::ThermalDensity f;for(int a=0;a<10;a++){f.mass[a]=1.25;f.internal_energy[a]=1e-6;}
    auto q0=compact_p3::mass_rule(p.polynomial);
    std::vector<compact_p3::QuadraturePoint>quadrature(q0.begin(),q0.end());
    auto mass=cg::mass(p,f,quadrature);long double M[12][12]={};
    for(int a=0;a<12;a++)for(int b=0;b<12;b++)M[a][b]=mass[a*12+b];
    // Global rank is not inferred from diagonal entries. This synthetic
    // polynomial fixture has two independent quartic completion functions.
    long double chol[12][12]={},minpivot=1e300L;
    for(int a=0;a<12;a++)for(int b=0;b<=a;b++){
        long double v=M[a][b];for(int k=0;k<b;k++)v-=chol[a][k]*chol[b][k];
        if(a==b){need(v>0&&std::isfinite(v),"synthetic completed mass positive");minpivot=std::min(minpivot,v);chol[a][b]=sqrtl(v);}
        else chol[a][b]=v/chol[b][b];
    }
    long double rigid_error=0,point_K=0,point_L=0,point_mass=0;
    for(const auto&q:quadrature){
        cg::Basis b;GfsMaterialPoint12 x;need(cg::basis(p,q.l,&b)&&cg::point(b,c.data(),&x),"completed rigid map");
        double expected[4]={0,-omega,omega,0};for(int k=0;k<4;k++)rigid_error=std::max(rigid_error,fabsl(x.grad_v[k]-expected[k]));
        point_mass+=q.weight*1.25L;point_K+=.5L*q.weight*1.25L*(x.velocity[0]*x.velocity[0]+x.velocity[1]*x.velocity[1]);
        point_L+=q.weight*1.25L*(x.position[0]*x.velocity[1]-x.position[1]*x.velocity[0]);
    }
    long double metric_K=0,metric_L=0,metric_mass=0;
    for(int a=0;a<12;a++)for(int b=0;b<12;b++){
        metric_K+=.5L*M[a][b]*(static_cast<long double>(c[a].v[0])*c[b].v[0]+static_cast<long double>(c[a].v[1])*c[b].v[1]);
        metric_L+=M[a][b]*(static_cast<long double>(c[a].x[0])*c[b].v[1]-static_cast<long double>(c[a].x[1])*c[b].v[0]);
        if(a<10&&b<10)metric_mass+=M[a][b];
    }
    need(rigid_error<2e-13L&&fabsl(metric_K-point_K)<2e-15L&&fabsl(metric_L-point_L)<2e-15L&&fabsl(metric_mass-point_mass)<2e-15L,"mass/rigid-rotation represented-field closure");
    // Independent dynamic geometry velocities, including off-axis motion.
    // Freezing the two added mode velocities makes the work check fail.
    c[10].v[0]=.043;c[10].v[1]=-.027;c[11].v[0]=-.031;c[11].v[1]=.019;
    std::vector<cg::Sample>samples;for(const auto&q:quadrature)samples.push_back({0,{q.l[0],q.l[1],q.l[2]}});
    samples.push_back({0,{1,0,0}});samples.push_back({0,{0,.63,.37}});samples.push_back({0,{0,1,0}});samples.push_back({0,{0,0,1}});
    size_t n=samples.size();std::vector<cg::StressInput>si(n);std::vector<cg::BoundaryInput>bi(n);
    std::vector<cg::SampleValue>cpu(n);std::vector<cg::StressValue>sv(n);std::vector<cg::BoundaryValue>bv(n);
    double stress_energy=0,stress_torque=0,boundary_work=0,boundary_torque=0,geometry_error=0;
    for(size_t i=0;i<n;i++){
        cg::Basis b;need(cg::basis(p,samples[i].l,&b)&&cg::point(b,c.data(),&cpu[i].map),"completed arbitrary velocity map");
        double B[10];compact_p3::bernstein(samples[i].l,B);cpu[i].rho_reference=compact_p3::field(f.mass,B);cpu[i].energy_reference=compact_p3::field(f.internal_energy,B);cpu[i].valid=1;
        STJet E[2];need(cg::extension(p,samples[i].l,E),"synthetic original curve chart");
        for(int k=0;k<2;k++){
            double ref=0;for(int a=0;a<3;a++)ref+=samples[i].l[a]*p.polynomial.vertex[a][k];
            geometry_error=std::max(geometry_error,fabs(cpu[i].map.position[k]-ref-E[k].v));
        }
        si[i].point=samples[i];si[i].weight=.03;si[i].rho_reference=1.25;
        si[i].sigma[0]=-.4;si[i].sigma[1]=.08;si[i].sigma[2]=-.7;si[i].gravity[0]=-.21;si[i].gravity[1]=.17;
        sv[i]=cg::stress(b,c.data(),cpu[i].map,si[i].rho_reference,si[i].weight,si[i].sigma,si[i].gravity[0],si[i].gravity[1]);need(sv[i].valid,"CPU stress valid");
        stress_energy=std::max(stress_energy,fabs(sv[i].energy_residual)/std::max(1.,sv[i].energy_scale));
        stress_torque=std::max(stress_torque,fabs(sv[i].torque_residual)/std::max(1.,sv[i].torque_scale));
        // Boundary API algebra is checked at each sample; only the explicitly
        // exterior samples below are used for the direct trace/work check.
        bi[i].point=samples[i];bi[i].normal[0]=bi[i].normal[1]=sqrt(.5);bi[i].weight=.03;bi[i].pressure=1e-6;
        bv[i]=cg::exterior(b,c.data(),cpu[i].map,bi[i].normal[0],bi[i].normal[1],bi[i].weight,bi[i].pressure);need(bv[i].valid,"CPU exterior valid");
        boundary_work=std::max(boundary_work,fabs(bv[i].work_residual));boundary_torque=std::max(boundary_torque,fabs(bv[i].torque_residual));
        if(samples[i].l[0]==0){
            STJet gamma[2];need(cg::trace(p.trace,STJet(samples[i].l[2],1,0),gamma),"same exterior trace derivative");
            double ds=sqrt(2.),actual_work=-bi[i].weight/ds*bi[i].pressure*(cpu[i].map.velocity[0]*gamma[1].x-cpu[i].map.velocity[1]*gamma[0].x);
            need(fabs(actual_work-bv[i].work)<1e-20,"actual curved trace normal/work");
        }
    }
    need(stress_energy<2e-13&&stress_torque<2e-13&&boundary_work<1e-20&&boundary_torque<1e-20&&geometry_error<2e-14,"same-basis stress/boundary/geometry identities");
    cuda_need(cudaSetDevice(0));Device<cg::Patch>dp(1);dp.put(&p);Device<STCoefficient>dc(c.size());dc.put(c.data());Device<compact_p3::ThermalDensity>df(1);df.put(&f);
    Device<cg::Sample>dq(n);dq.put(samples.data());Device<cg::SampleValue>dv(n);Device<cg::StressInput>ds(n);ds.put(si.data());Device<cg::StressValue>dsv(n);
    Device<cg::BoundaryInput>db(n);db.put(bi.data());Device<cg::BoundaryValue>dbv(n);
    cuda_need(cudaError_t(compact_geometry_sample_device(dp.ptr,1,dc.ptr,c.size(),df.ptr,dq.ptr,n,dv.ptr,nullptr)));
    cuda_need(cudaError_t(compact_geometry_stress_device(dp.ptr,1,dc.ptr,c.size(),ds.ptr,n,dsv.ptr,nullptr)));
    cuda_need(cudaError_t(compact_geometry_exterior_device(dp.ptr,1,dc.ptr,c.size(),db.ptr,n,dbv.ptr,nullptr)));cuda_need(cudaDeviceSynchronize());
    auto gpu=dv.get();auto gs=dsv.get();auto gb=dbv.get();double gpu_error=0;
    for(size_t i=0;i<n;i++){
        need(gpu[i].valid&&gs[i].valid&&gb[i].valid,"GPU completed samples/work valid");
        const auto&a=cpu[i].map;const auto&b=gpu[i].map;
        for(int k=0;k<2;k++){gpu_error=std::max(gpu_error,difference(a.position[k],b.position[k]));gpu_error=std::max(gpu_error,difference(a.velocity[k],b.velocity[k]));}
        for(int k=0;k<4;k++){gpu_error=std::max(gpu_error,difference(a.F[k],b.F[k]));gpu_error=std::max(gpu_error,difference(a.grad_v[k],b.grad_v[k]));}
        gpu_error=std::max(gpu_error,difference(a.J,b.J));gpu_error=std::max(gpu_error,difference(cpu[i].rho_reference,gpu[i].rho_reference));gpu_error=std::max(gpu_error,difference(cpu[i].energy_reference,gpu[i].energy_reference));
        for(int k=0;k<2;k++)for(int a0=0;a0<12;a0++){gpu_error=std::max(gpu_error,difference(sv[i].force[k][a0],gs[i].force[k][a0]));gpu_error=std::max(gpu_error,difference(bv[i].force[k][a0],gb[i].force[k][a0]));}
        gpu_error=std::max({gpu_error,difference(sv[i].Udot,gs[i].Udot),difference(sv[i].gravity_power,gs[i].gravity_power),difference(bv[i].work,gb[i].work),difference(bv[i].torque,gb[i].torque)});
        need(fabs(gs[i].energy_residual)<2e-13*std::max(1.,gs[i].energy_scale)&&fabs(gs[i].torque_residual)<2e-13*std::max(1.,gs[i].torque_scale),"GPU symmetric stress identities");
        need(fabs(gb[i].work_residual)<1e-20&&fabs(gb[i].torque_residual)<1e-20,"GPU exterior work/torque identities");
    }
    need(gpu_error<2e-12,"CPU/Ada kernel agreement");
    std::cout<<std::setprecision(17)<<"{\"status\":\"COMPLETED_BASIS_ALGEBRA_AND_ADA_API_PASS_NO_PHYSICAL_RUN\",\"samples\":"<<n<<",\"synthetic_mass_min_pivot\":"<<double(minpivot)<<",\"rigid_gradient_error\":"<<double(rigid_error)<<",\"metric_point_K_difference\":"<<double(fabsl(metric_K-point_K))<<",\"metric_point_L_difference\":"<<double(fabsl(metric_L-point_L))<<",\"stress_energy_residual\":"<<stress_energy<<",\"stress_torque_residual\":"<<stress_torque<<",\"exterior_work_residual\":"<<boundary_work<<",\"exterior_torque_residual\":"<<boundary_torque<<",\"geometry_reconstruction_error\":"<<geometry_error<<",\"CPU_GPU_max_normalized_difference\":"<<gpu_error<<",\"actual_source_geometry_proven\":false,\"remap_committed\":false,\"new_physical_step\":false}\n";
    return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<'\n';return 1;}}
