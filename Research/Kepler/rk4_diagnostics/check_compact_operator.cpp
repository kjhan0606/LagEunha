// NEW finite-functional RHS fixture, not a source guard, remap or trajectory.
// Two cells/ranks, a shared MPI face, enriched nonpolynomial curved exterior,
// arbitrary projected heat, dynamic geometry modes and original resident M.
#include "material_compact_operator.h"
#include <iostream>
#include <iomanip>
#include <set>
using R=long double;namespace md=material_distributed;
static void need(bool b,const char*s){md::require(b,s);}
#include "material_shared_trace_factory.h"
#include "material_full_gpu_mass_host.h"
namespace cg=compact_geometry;namespace cv=compact_volume;namespace co=compact_operator;
const uint64_t triangle_vertex[2][3]={{3,0,1},{1,2,3}};
const double vertex[4][2]={{0,0},{1,0},{1,1},{0,1}};
const std::pair<uint64_t,uint64_t>edge_vertex[5]={{0,1},{0,3},{1,2},{1,3},{2,3}};
int edge_id(uint64_t a,uint64_t b){auto k=std::minmax(a,b);for(int j=0;j<5;j++)if(edge_vertex[j]==std::make_pair(k.first,k.second))return j;throw std::runtime_error("fixture edge");}
cg::Patch global_patch(int cell){
    cg::Patch p;p.curved=1;p.polynomial.cell=0;p.modes.component=0;p.modes.id[0]=16;p.modes.id[1]=17;p.modes.scale[0]=.02;p.modes.scale[1]=.03;
    const auto&v=triangle_vertex[cell];for(int j=0;j<3;j++)for(int k=0;k<2;k++)p.polynomial.vertex[j][k]=vertex[v[j]][k];
    compact_p3::geometry(p.polynomial);
    for(int a=0;a<10;a++){
        int index[3];gfs_material_p3_index(a,index);int which=-1,zero=-1;
        for(int k=0;k<3;k++){if(index[k]==3)which=k;if(index[k]==0)zero=k;}
        if(which>=0)p.polynomial.id[a]=v[which];
        else if(zero>=0){int i=(zero+1)%3,j=(zero+2)%3;int part=v[i]<v[j]?index[j]:index[i];p.polynomial.id[a]=4+2*edge_id(v[i],v[j])+part-1;}
        else p.polynomial.id[a]=14+cell;
    }
    STEdge e={};e.rank=2;e.eps=.047;e.A[0]=cell?1:-1;e.A[1]=cell?1:-1;e.D[0]=cell?-2:2;
    e.encode[0][0]=e.encode[1][1]=1;
    std::vector<STCoefficient>old(6);
    for(int j=0;j<4;j++){
        e.poly[j]=j;double t=j/3.,gradient[4];
        need(gfs_material_p3_radial_remainder(e.A,e.A[0]+t*e.D[0],e.A[1],e.eps,e.qnode[j],gradient),"fixture original enriched trace");
        old[j].x[0]=cell?1-t+.02*t*(1-t):t+.02*t*(1-t)*(t-.5);
        old[j].x[1]=cell?1+.03*t*(1-t)*(t-.5):-.03*t*(1-t);
    }
    e.id[0]=4;e.id[1]=5;old[4].x[0]=.04;old[4].x[1]=-.03;old[5].x[0]=.03;old[5].x[1]=.04;
    p.trace=cg::freeze_trace(e,old,cell?4:0,123,5);cg::interpolate_chart(p);return p;
}
cv::Rule gauss(int count){cv::Rule r;shared_trace::gauss(count,r.point,r.weight);return r;}
std::vector<cv::CanonicalEdge>boundary(const cg::Patch&p,int rank){
    std::vector<cv::CanonicalEdge>out;
    for(int k=0;k<3;k++){
        uint64_t a=triangle_vertex[rank][k],b=triangle_vertex[rank][(k+1)%3];bool forward=a<b;int event=edge_id(a,b);
        cv::CanonicalEdge e;e.basis=cv::from_patch(p,k,forward,event);e.side.push_back({0,forward?1:-1});e.physical_exterior=event!=3;out.push_back(e);
    }return out;
}
std::array<R,MF_EVOL_METRICS>cpu_reference(const co::Seed&s,const std::vector<STCoefficient>&c,const co::Stage&stage,std::array<md::Vec,2>&force){
    std::array<R,MF_EVOL_METRICS>m={};for(auto&f:force)f.resize(c.size());const auto&C=stage.context[0];
    for(const auto&q:s.quadrature){
        cg::Basis b;GfsMaterialPoint12 x,x0;need(cg::basis(s.patch[0],q.l,&b)&&cg::point(b,c.data(),&x)&&cg::point(b,s.generating.data(),&x0),"CPU reference completed field");
        double B[10];compact_p3::bernstein(q.l,B);double rho=compact_p3::field(s.field[0].mass,B),e0=compact_p3::field(s.field[0].internal_energy,B);
        R j=R(x.J)+C.delta,j0=R(x0.J)+C.generating_delta,jdot=R(x.J)*(x.grad_v[0]+x.grad_v[3])+C.delta_dot;
        R U=R(e0)*powl(j0/j,R(s.law.gamma)-1),P=(R(s.law.gamma)-1)*U/j,dp=P-s.law.background_pressure;
        R dm=R(q.weight)*rho,z=R(x.position[0])*x.position[0]+R(x.position[1])*x.position[1]+R(s.law.eps)*s.law.eps;
        R grav=powl(z,-1.5L),gx=-x.position[0]*grav,gy=-x.position[1]*grav;
        R value[MF_EVOL_METRICS]={q.weight,R(q.weight)*x.J,R(q.weight)*x.J*(x.grad_v[0]+x.grad_v[3]),dm,R(q.weight)*dp,R(q.weight)*U,-R(q.weight)*P*jdot,R(q.weight)*j,
            .5L*dm*(R(x.velocity[0])*x.velocity[0]+R(x.velocity[1])*x.velocity[1]),dm*(R(x.position[0])*x.velocity[1]-R(x.position[1])*x.velocity[0]),-dm/sqrtl(z),
            dm*(x.velocity[0]*gx+x.velocity[1]*gy),dm*(x.position[0]*gy-x.position[1]*gx),dm*x.position[0],dm*x.position[1],dm*x.velocity[0],dm*x.velocity[1],dm*(R(x.position[0])*x.position[0]+R(x.position[1])*x.position[1]),2*dm*(R(x.position[0])*x.velocity[0]+R(x.position[1])*x.velocity[1])};
        for(int k=0;k<MF_EVOL_METRICS;k++)m[k]+=value[k];
        double sigma[3]={-double(dp-C.mean_deltaP),0,-double(dp-C.mean_deltaP)},N[12],Nx[12],Ny[12],F[2][12]={},u=0,g=0;
        for(int a=0;a<b.n;a++){N[a]=b.N[a].v;Nx[a]=b.N[a].x;Ny[a]=b.N[a].y;}
        need(gfs_material12_stress_at_point(b.n,&x,N,Nx,Ny,rho,q.weight,sigma,double(gx),double(gy),F,&u,&g),"CPU reference stress");
        for(int a=0;a<b.n;a++)for(int d=0;d<2;d++)force[d][b.id[a]]+=F[d][a];
    }
    auto face=cv::mean_pressure_force(s.edge,stage.geometry,{R(s.law.background_pressure)+C.mean_deltaP},c.size(),s.law.background_pressure);
    for(size_t j=0;j<c.size();j++)for(int d=0;d<2;d++)force[d][j]+=face[j][d];return m;
}
R thermal_energy(const co::Seed&s,const std::vector<STCoefficient>&c){
    auto g=cv::assemble(s.edge,c,1,s.polynomial_edge_rule,s.curved_edge_rule,s.constant);R W=0,QJ=0;
    for(const auto&q:s.quadrature){cg::Basis b;GfsMaterialPoint12 x;need(cg::basis(s.patch[0],q.l,&b)&&cg::point(b,c.data(),&x),"perturbed thermal geometry");W+=q.weight;QJ+=R(q.weight)*x.J;}
    R delta=(g.cell[0].volume-QJ)/W,total=0;
    for(const auto&q:s.quadrature){cg::Basis b;GfsMaterialPoint12 x,x0;need(cg::basis(s.patch[0],q.l,&b)&&cg::point(b,c.data(),&x)&&cg::point(b,s.generating.data(),&x0),"perturbed thermal point");
        double B[10];compact_p3::bernstein(q.l,B);double e0=compact_p3::field(s.field[0].internal_energy,B);
        total+=R(q.weight)*e0*powl((R(x0.J)+s.generating_delta[0])/(R(x.J)+delta),R(s.law.gamma)-1);
    }return total;
}
int main(int argc,char**argv){MPI_Init(&argc,&argv);int rank=0,ranks=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
try{need(ranks==2,"two MPI ranks required");{
    need(cudaSetDevice(0)==cudaSuccess,"Ada device");auto p0=global_patch(0),p1=global_patch(1);std::vector<STCoefficient>global(18);std::vector<bool>seen(18);
    for(const auto&p:{p0,p1})for(int a=0;a<10;a++){
        auto x=cg::initial_polynomial_position(p,a);auto id=p.polynomial.id[a];
        if(seen[id])need(global[id].x[0]==x[0]&&global[id].x[1]==x[1],"canonical fixture shared geometry bits");
        for(int k=0;k<2;k++)global[id].x[k]=x[k];seen[id]=true;
    }
    global[16]=cg::initial_mode(p0.modes,0);global[17]=cg::initial_mode(p0.modes,1);
    auto p=rank?p1:p0;std::set<uint64_t>required;for(auto id:p.polynomial.id)required.insert(id);required.insert(16);required.insert(17);for(uint64_t id=0;id<18;id++)if(int(id%2)==rank)required.insert(id);
    std::vector<uint64_t>ids(required.begin(),required.end());std::vector<int>owners;std::map<uint64_t,size_t>local;
    for(size_t j=0;j<ids.size();j++){owners.push_back(ids[j]%2);local[ids[j]]=j;}
    md::Ownership layout(ids,owners,MPI_COMM_WORLD);for(auto&id:p.polynomial.id)id=local.at(id);for(auto&id:p.modes.id)id=local.at(id);
    co::Seed s;s.law={1.6666669845581055,1e-6,.047};s.patch={p};s.field.resize(1);s.generating.resize(ids.size());s.constant.resize(ids.size());
    for(size_t j=0;j<ids.size();j++){s.generating[j]=global[ids[j]];s.constant[j]=ids[j]<16;}
    for(int a=0;a<10;a++){s.field[0].mass[a]=1.25+.05*a;s.field[0].internal_energy[a]=1e-6*(1+.02*a+.03*rank);}
    s.edge=boundary(p,rank);s.polynomial_edge_rule=gauss(3);s.curved_edge_rule=gauss(64);
    auto rule=compact_p3::mass_rule(p.polynomial);std::vector<compact_p3::QuadraturePoint>mass_rule(rule.begin(),rule.end());
    for(const auto&q:rule)s.quadrature.push_back({{q.l[0],q.l[1],q.l[2]},q.weight});s.offset={0,s.quadrature.size()};
    auto generating_geometry=cv::assemble(s.edge,s.generating,1,s.polynomial_edge_rule,s.curved_edge_rule,s.constant);R W=0,QJ=0;
    for(const auto&q:s.quadrature){cg::Basis b;GfsMaterialPoint12 x;need(cg::basis(p,q.l,&b)&&cg::point(b,s.generating.data(),&x),"fixture generating geometry");W+=q.weight;QJ+=R(q.weight)*x.J;}
    s.generating_delta={double((generating_geometry.cell[0].volume-QJ)/W)};
    std::vector<STCoefficient>c=s.generating;
    for(size_t j=0;j<c.size();j++){
        double x=global[ids[j]].x[0],y=global[ids[j]].x[1];
        c[j].x[0]=1.01*x+.02*y+(ids[j]<16?2.:.001);c[j].x[1]=-.01*x+.99*y+(ids[j]<16?3.:-.0007);
        c[j].v[0]=ids[j]<16?.02*x-.01*y+.03*x*x:.004*(ids[j]-15);
        c[j].v[1]=ids[j]<16?.03*x+.01*y-.01*x*y:-.003*(ids[j]-15);
    }
    auto saved=c;co::Operator op(s,layout);auto stage=op.evaluate(c);
    need(memcmp(c.data(),saved.data(),c.size()*sizeof(STCoefficient))==0,"RHS does not change dynamic state");
    need(stage.geometry.external==2&&stage.geometry.partition_boundary==1,"shared MPI face is not physical exterior");
    std::array<md::Vec,2>reference_force;auto reference=cpu_reference(s,c,stage,reference_force);
    MPI_Allreduce(MPI_IN_PLACE,reference.data(),MF_EVOL_METRICS,MPI_LONG_DOUBLE,MPI_SUM,layout.comm);for(auto&f:reference_force)layout.sum(f);
    R metric_error=0,thermal_error=0,force_error=0;
    for(int k=0;k<MF_EVOL_METRICS;k++)metric_error=std::max(metric_error,fabsl(reference[k]-stage.global[k])/std::max(1.L,fabsl(reference[k])));
    for(int k:{4,5,6})thermal_error=std::max(thermal_error,fabsl(reference[k]-stage.global[k])/stage.global[5]);
    for(size_t j:layout.owned)for(int d=0;d<2;d++)force_error=std::max(force_error,fabsl(reference_force[d][j]-stage.force[d][j]));
    MPI_Allreduce(MPI_IN_PLACE,&force_error,1,MPI_LONG_DOUBLE,MPI_MAX,layout.comm);
    R heat_identity=fabsl(stage.energy_identity)/stage.global[5],torque_identity=fabsl(stage.torque_identity)/std::max(1.L,stage.torque_scale);
    need(metric_error<2e-13L&&thermal_error<2e-11L&&force_error<2e-13L&&heat_identity<1e-8L&&torque_identity<2e-13L,"completed RHS CPU/CUDA pressure/work/torque closure");
    double h=1e-5;auto plus=c,minus=c;for(size_t j=0;j<c.size();j++)for(int k=0;k<2;k++){plus[j].x[k]+=h*c[j].v[k];minus[j].x[k]-=h*c[j].v[k];}
    R fd=(thermal_energy(s,plus)-thermal_energy(s,minus))/(2*h);MPI_Allreduce(MPI_IN_PLACE,&fd,1,MPI_LONG_DOUBLE,MPI_SUM,layout.comm);
    R derivative_error=fabsl(fd-stage.global[6])/stage.global[5];need(derivative_error<2e-8L,"thermal potential finite derivative matches CUDA Udot");
    compact_mass::LocalElement e;e.component=0;e.n=12;e.mass=cg::mass(p,s.field[0],mass_rule);for(int a=0;a<10;a++)e.id[a]=p.polynomial.id[a];for(int k=0;k<2;k++)e.id[10+k]=p.modes.id[k];
    auto mass=compact_mass::prepare(layout,{e},{{0,{16,17}}});full_gpu_mass::Operator M(*mass->exact);std::array<full_gpu_mass::Solution,2>a;
    R solve_residual=0,Kdot=0;uint64_t iterations=0;
    for(int d=0;d<2;d++){
        a[d]=full_gpu_mass::solve(M,stage.force[d]);solve_residual=std::max(solve_residual,a[d].true_relative_residual);iterations+=a[d].inner_iterations;
        auto product=mass->exact->apply(a[d].x);for(size_t j:layout.owned)Kdot+=R(c[j].v[d])*product[j];
    }
    MPI_Allreduce(MPI_IN_PLACE,&Kdot,1,MPI_LONG_DOUBLE,MPI_SUM,layout.comm);R final_identity=fabsl(Kdot+stage.global[6]-stage.global[11]-stage.external_work)/stage.global[5];
    need(solve_residual<1e-16L&&final_identity<1e-8L&&M.resident,"actual original resident M solve closes continuous RHS energy");
    auto thermal=compact_thermo::evaluate(s.law,2e-6,1.2,1.1,1.1,0);
    need(thermal.valid&&thermal.energy_reference==2e-6&&fabs(thermal.pressure-s.law.background_pressure)>.1e-6,"projected thermal seed is not reset to P0");
    R face_error=stage.maximum_face_derivative_error,delta=fabsl(stage.context[0].delta);MPI_Allreduce(MPI_IN_PLACE,&delta,1,MPI_LONG_DOUBLE,MPI_MAX,layout.comm);
    if(rank==0)std::cout<<std::setprecision(17)<<"{\"status\":\"COMPACT_ENRICHED_CURVED_RHS_ORIGINAL_RESIDENT_MPI2_ADA_PASS_NO_PHYSICAL_RUN\",\"ranks\":2,\"cells\":2,\"global_owned_dofs\":"<<layout.global_owned<<",\"points_per_pass\":"<<stage.points_per_pass<<",\"nonpolynomial_geometry\":true,\"MPI_interface_treated_as_exterior\":false,\"CPU_GPU_metric_difference\":"<<double(metric_error)<<",\"CPU_GPU_thermal_difference_over_U\":"<<double(thermal_error)<<",\"CPU_GPU_force_difference\":"<<double(force_error)<<",\"RHS_energy_identity_over_U\":"<<double(heat_identity)<<",\"RHS_torque_identity\":"<<double(torque_identity)<<",\"thermal_finite_derivative_error_over_U\":"<<double(derivative_error)<<",\"face_adjoint_cofactor_difference\":"<<double(face_error)<<",\"maximum_absolute_delta\":"<<double(delta)<<",\"resident_inner_iterations_two_components\":"<<iterations<<",\"binary80_true_mass_residual\":"<<double(solve_residual)<<",\"mass_solved_energy_identity_over_U\":"<<double(final_identity)<<",\"RHS_wall_seconds\":"<<stage.wall_seconds<<",\"actual_source_geometry_proven\":false,\"conservative_remap_committed\":false,\"new_physical_step\":false}\n";
}MPI_Finalize();return 0;
}catch(const std::exception&e){std::cerr<<"rank "<<rank<<": "<<e.what()<<'\n';MPI_Abort(MPI_COMM_WORLD,1);return 1;}}
