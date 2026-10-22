#include <mpi.h>
#include <omp.h>
#include <algorithm>
#include <array>
#include <vector>
#include <map>
#include <string>
#include <set>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <stdexcept>
#include <cuda_runtime.h>
static void need(bool p,const char*s){if(!p)throw std::runtime_error(s);}
using R=long double;
#include "material_full_owned.h"
#include "material_full_initial_batch.h"
#include "material_distributed_mass.h"
namespace shared_trace_rhs {using Vec=std::vector<R>;}
namespace st=shared_trace;
#include "material_geometric_pressure_geometry.h"
namespace md=material_distributed;namespace gp=geometric_pressure;
static R value(MFDoubleDouble a){return R(a.hi)+a.lo;}
struct Device {
    void*ptr=nullptr;~Device(){mf_initial_destroy(ptr);}
};
static void config_rule(GfsEnrichConfig&cfg,double gamma,double pressure,double eps){
    cfg={};cfg.na=64;cfg.nr=28;cfg.field.ic.kind=2;cfg.field.ic.gamma=gamma;cfg.field.ic.pressure=pressure;cfg.field.ic.cx=cfg.field.ic.cy=12;cfg.field.ic.rin=2;cfg.field.ic.rout=8;cfg.field.ic.eps=eps;
    std::vector<double>x,w;shared_trace::gauss(64,x,w);for(int i=0;i<64;i++){cfg.an[i]=x[i];cfg.aw[i]=w[i];}shared_trace::gauss(28,x,w);for(int i=0;i<28;i++){cfg.rn[i]=x[i];cfg.rw[i]=w[i];}
}
static std::array<R,2> cpu_shadow(const material_full::Owned&m,size_t s,const GfsEnrichConfig&cfg,const MFInitialBlock&gpu){
    STTriangle t=m.f.triangle[s];std::array<STCoefficient,ST_MAX>c{};for(int j=0;j<t.n;j++){c[j]=m.f.coefficient[t.id[j]];t.id[j]=j;}
    R M[ST_MAX*ST_MAX]={},F[2][ST_MAX]={};double x[3],y[3];for(int j=0;j<3;j++){x[j]=t.vertex[j][0];y[j]=t.vertex[j][1];}GfsMomentPolygon poly={x,y,3,0};
    for(int j=0;j<gfs_enrich_nquad(poly,&cfg,&m.rule[s]);j++){
        STPoint q={};q.triangle=0;q.boundary=-1;need(gfs_enrich_quad(poly,&cfg,&m.rule[s],j,&q.X,&q.Y,&q.weight),"fixed full shadow quadrature");if(q.weight==0)continue;double u[4];need(q.weight>0&&gfs_enrich_reference(&cfg,q.X+12,q.Y+12,u,0,0),"physical full shadow density");q.rho=u[0];
        R ux=R(x[1])-x[0],uy=R(y[1])-y[0],vx=R(x[2])-x[0],vy=R(y[2])-y[0],dx=R(q.X)-x[0],dy=R(q.Y)-y[0],det=ux*vy-uy*vx;q.l[1]=(dx*vy-dy*vx)/det;q.l[2]=(ux*dy-uy*dx)/det;q.l[0]=1-q.l[1]-q.l[2];STJet N[ST_MAX];GfsMaterialPoint12 p;need(st_shape(t,q,N)&&st_point(t,c.data(),N,&p),"full CPU shadow field");
        R dm=R(q.weight)*q.rho;for(int a=0;a<t.n;a++)for(int b=0;b<t.n;b++)M[a*ST_MAX+b]+=dm*N[a].v*N[b].v;
        double g=pow(p.position[0]*p.position[0]+p.position[1]*p.position[1]+cfg.field.ic.eps*cfg.field.ic.eps,-1.5),sigma[3]={0,0,0};
        for(int begin=0;begin<t.n;begin+=12){int n=std::min(12,t.n-begin);double A[12],Ax[12],Ay[12],f[2][12]={},work=0,power=0;for(int a=0;a<n;a++){A[a]=N[a+begin].v;Ax[a]=N[a+begin].x;Ay[a]=N[a+begin].y;}need(gfs_material12_stress_at_point(n,&p,A,Ax,Ay,q.rho,q.weight,sigma,-p.position[0]*g,-p.position[1]*g,f,&work,&power),"original native full CPU shadow gravity");for(int a=0;a<n;a++)for(int d=0;d<2;d++)F[d][a+begin]+=f[d][a];}
    }
    R mdiff=0,mscale=0,fdiff=0,fscale=0;for(int a=0;a<t.n;a++){for(int b=0;b<t.n;b++){mdiff=std::max(mdiff,fabsl(M[a*ST_MAX+b]-value(gpu.mass[a*ST_MAX+b])));mscale=std::max(mscale,fabsl(M[a*ST_MAX+b]));}for(int d=0;d<2;d++){fdiff=std::max(fdiff,fabsl(F[d][a]-value(gpu.force[d][a])));fscale=std::max(fscale,fabsl(F[d][a]));}}
    return {mdiff/mscale,fdiff/fscale};
}
static void transform(material_full::Owned&m,size_t s,const MFInitialBlock&input,R*dest,R*force){
    R original[ST_MAX*ST_MAX];for(int j=0;j<ST_MAX*ST_MAX;j++)original[j]=value(input.mass[j]);auto&t=m.f.triangle[s];shared_trace::bubble_basis(t,original,double(value(input.metric[0])));
    R T[ST_MAX][ST_MAX]={},tmp[ST_MAX][ST_MAX]={};for(int a=0;a<t.n;a++)T[a][a]=1;
    for(int a=0;a<2;a++){int r=t.bubble_offset+a;for(int b=0;b<t.n;b++)T[r][b]=0;for(int b=0;b<2;b++){T[r][t.bubble_offset+b]=t.bubble_encode[a][b];T[r][9]-=R(t.bubble_encode[a][b])*t.bubble_shift[b];}}
    for(int a=0;a<t.n;a++)for(int b=0;b<t.n;b++)for(int k=0;k<t.n;k++)tmp[a][b]+=T[a][k]*original[k*ST_MAX+b];
    for(int a=0;a<t.n;a++)for(int b=0;b<=a;b++){R z=0;for(int k=0;k<t.n;k++)z+=tmp[a][k]*T[b][k];dest[a*ST_MAX+b]=dest[b*ST_MAX+a]=z;}
    for(int a=0;a<t.n;a++)for(int d=0;d<2;d++){R z=0;for(int k=0;k<t.n;k++)z+=T[a][k]*value(input.force[d][k]);force[d*ST_MAX+a]=z;}
    for(int d=0;d<2;d++)m.f.coefficient[t.id[9]].v[d]+=t.bubble_shift[d];for(int j=0;j<2;j++)for(int d=0;d<2;d++)m.f.coefficient[t.id[t.bubble_offset+j]].v[d]=t.bubble_decode[d][j];
}
static int initialize(const std::string&input,MPI_Comm comm,double gamma,double pressure,double eps){
    int rank,ranks;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&ranks);double begin=MPI_Wtime();int local_rank;MPI_Comm shared;MPI_Comm_split_type(comm,MPI_COMM_TYPE_SHARED,rank,MPI_INFO_NULL,&shared);MPI_Comm_rank(shared,&local_rank);MPI_Comm_free(&shared);int devices=0;need(cudaGetDeviceCount(&devices)==cudaSuccess&&devices>0,"native material visible GPUs");int device=local_rank%devices;need(cudaSetDevice(device)==cudaSuccess,"native material local-rank GPU mapping");cudaDeviceProp prop;need(cudaGetDeviceProperties(&prop,device)==cudaSuccess,"native material GPU properties");
    FILE*runtime=fopen(("full_initial_runtime_r"+std::to_string(rank)+".json").c_str(),"wx");need(runtime,"new full initial runtime metadata");fprintf(runtime,"{\"rank\":%d,\"ranks\":%d,\"OMP_threads\":%d,\"local_rank\":%d,\"visible_devices\":%d,\"selected_device\":%d,\"uuid\":\"",rank,ranks,omp_get_max_threads(),local_rank,devices,device);for(unsigned char c:prop.uuid.bytes)fprintf(runtime,"%02x",unsigned(c));fprintf(runtime,"\",\"device\":\"%s\"}\n",prop.name);material_full::durable(runtime);
    auto state=material_full::construct(input,comm);double basis_seconds=MPI_Wtime()-begin;md::Ownership layout(state.global,state.owner,comm);need(layout.global_owned==state.global_dofs,"one global owner per full scalar DOF");
    for(int component=0;component<4;component++){md::Vec v(state.global.size());for(size_t i=0;i<v.size();i++)v[i]=component<2?state.f.coefficient[i].x[component]:state.f.coefficient[i].v[component-2];auto before=v;layout.pull(v);need(v==before,"initial owned/ghost shared coefficients bitwise consistent");}
    std::string suffix=std::to_string(rank);material_full::save_basis(state,"full_unconditioned_basis_r"+suffix+".bin");
    GfsEnrichConfig cfg;config_rule(cfg,gamma,pressure,eps);Device gpu;need(mf_initial_create(state.f.triangle.data(),state.f.triangle.size(),state.f.coefficient.data(),state.f.coefficient.size(),state.rule.data(),&cfg,&gpu.ptr)==0,"persistent full GPU mass/gravity initialization");
    std::vector<size_t>shadow;for(double radius:{0.,.5,2.,8.,12.}){size_t best=0;double error=1e99;for(size_t s=0;s<state.f.triangle.size();s++){const auto&t=state.f.triangle[s];double x=0,y=0;for(int j=0;j<3;j++){x+=t.vertex[j][0]/3;y+=t.vertex[j][1]/3;}double e=fabs(hypot(x,y)-radius);if(e<error){error=e;best=s;}}shadow.push_back(best);}std::sort(shadow.begin(),shadow.end());shadow.erase(std::unique(shadow.begin(),shadow.end()),shadow.end());
    md::Vec blocks(state.f.triangle.size()*ST_MAX*ST_MAX);std::array<md::Vec,2>F={md::Vec(state.global.size()),md::Vec(state.global.size())};std::vector<std::array<R,MF_INITIAL_METRICS>>cells(state.f.cell.size());
    std::vector<MFInitialBlock>batch(MF_BATCH_TRIANGLES);std::vector<std::array<R,2*ST_MAX>>local_force(MF_BATCH_TRIANGLES);R minJ=1e99L,shadow_mass=0,shadow_force=0;uint64_t points=0,nt=state.f.triangle.size(),maxnt=0;MPI_Allreduce(&nt,&maxnt,1,MPI_UINT64_T,MPI_MAX,comm);double kernel=0,transform_seconds=0,start=MPI_Wtime();
    FILE*raw=fopen(("full_raw_initial_blocks_r"+suffix+".bin").c_str(),"wbx");FILE*progress=fopen(("full_initial_progress_r"+suffix+".tsv").c_str(),"wx");need(raw&&progress,"new fixed-mass raw blocks and progress");uint64_t header[6]={0x4d46494e49545231ULL,1,nt,sizeof(MFInitialBlock),uint64_t(rank),state.global_dofs};material_full::write(raw,header,6);fprintf(progress,"local_begin\ttriangles\tpositive_points\tkernel_seconds\tCPU_transform_seconds\n");
    for(uint64_t at=0;at<maxnt;at+=MF_BATCH_TRIANGLES){size_t n=at<nt?std::min<uint64_t>(MF_BATCH_TRIANGLES,nt-at):0;int bad=0;std::string error;
        try{if(n){double seconds=0;need(mf_initial_batch(gpu.ptr,at,n,batch.data(),&seconds)==0,"actual full GPU fixed mass and initial gravity");kernel+=seconds;material_full::write(raw,batch.data(),n);
            for(size_t j=0;j<n;j++){need(batch[j].bad_point<0,"full GPU quadrature/field failure: raw block preserves local point witness");need(batch[j].points>0&&batch[j].minJ>0,"full positive fixed mass quadrature state");points+=batch[j].points;minJ=std::min(minJ,R(batch[j].minJ));if(std::binary_search(shadow.begin(),shadow.end(),at+j)){auto compare=cpu_shadow(state,at+j,cfg,batch[j]);shadow_mass=std::max(shadow_mass,compare[0]);shadow_force=std::max(shadow_force,compare[1]);need(compare[0]<1e-11L&&compare[1]<1e-11L,"original CPU/CUDA mass/force bound on fixed full-state shadow sites");}}
            double t=MPI_Wtime();int failed=0;
            #pragma omp parallel for schedule(static) reduction(+:failed)
            for(size_t j=0;j<n;j++)try{transform(state,at+j,batch[j],blocks.data()+(at+j)*ST_MAX*ST_MAX,local_force[j].data());}catch(...){failed++;}
            need(!failed,"full fixed-mass bubble conditioning failed; raw blocks and unconditioned state retained");transform_seconds+=MPI_Wtime()-t;
            for(size_t j=0;j<n;j++){const auto&T=state.f.triangle[at+j];for(int d=0;d<2;d++)for(int a=0;a<T.n;a++)F[d][T.id[a]]+=local_force[j][d*ST_MAX+a];for(int a=0;a<MF_INITIAL_METRICS;a++)cells[T.cell][a]+=value(batch[j].metric[a]);}
            fprintf(progress,"%llu\t%zu\t%llu\t%.9g\t%.9g\n",(unsigned long long)at,n,(unsigned long long)points,kernel,transform_seconds);fflush(progress);
        }}catch(const std::exception&e){bad=1;error=e.what();}
        int any=0;MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,comm);if(any){if(bad){FILE*f=fopen(("full_initial_failure_r"+suffix+".txt").c_str(),"wx");if(f){fprintf(f,"batch_start=%llu\n%s\n",(unsigned long long)at,error.c_str());fclose(f);}fprintf(stderr,"rank%d FULL_INITIAL_REJECTED batch%llu %s\n",rank,(unsigned long long)at,error.c_str());}material_full::durable(raw);material_full::durable(progress);return 2;}
    }
    material_full::durable(raw);material_full::durable(progress);double integration_seconds=MPI_Wtime()-start;mf_initial_destroy(gpu.ptr);gpu.ptr=nullptr;
    // Preserve the fixed material operator and initialized full fields before
    // solving. This is initialization state, not an accepted evolved restart.
    material_full::save_basis(state,"full_conditioned_basis_r"+suffix+".bin");FILE*mass=fopen(("full_fixed_mass_r"+suffix+".bin").c_str(),"wbx");need(mass,"new fixed owned mass state");uint64_t mh[4]={0x4d464d4153533031ULL,1,uint64_t(blocks.size()),sizeof(R)};material_full::write(mass,mh,4);material_full::write(mass,blocks.data(),blocks.size());material_full::durable(mass);
    auto geometry=gp::geometry(state.f);std::vector<std::array<R,2>>delta(state.f.cell.size());R totalU=0,gravity_work=0,gravity_torque=0,quadratureK=0,initial_Kerror=0;
    for(size_t c=0;c<cells.size();c++){need(cells[c][1]>0&&geometry.volume[c]>0,"positive full thermodynamic source context");delta[c][0]=(geometry.volume[c]-cells[c][2])/cells[c][1];delta[c][1]=(geometry.vdot[c]-cells[c][3])/cells[c][1];totalU+=pressure*geometry.volume[c]/(gamma-1);gravity_work+=cells[c][8];gravity_torque+=cells[c][9];quadratureK+=cells[c][5];initial_Kerror+=cells[c][15];for(auto&g:geometry.G[c])for(int d=0;d<2;d++)F[d][g.first]+=R(pressure)*g.second[d];}
    R external_work=0,external_torque=0;
    for(const auto&e:geometry.edge){uint64_t poly=state.global[e.reference.poly[1]];need(poly>=state.global_vertices&&poly<state.global_vertices+2*state.global_edges,"canonical external edge DOF mapping");size_t global_edge=(poly-state.global_vertices)/2;if(!state.external_edge[global_edge])continue;need(e.side.size()==1,"one owned physical external side");
        for(int j=0;j<e.n;j++)for(int d=0;d<2;d++){R force=-R(pressure)*e.side[0].sign*e.H[j][d];int i=e.id[j];F[d][i]+=force;external_work+=state.f.coefficient[i].v[d]*force;external_torque+=(d==1?state.f.coefficient[i].x[0]:-state.f.coefficient[i].x[1])*force;}
    }
    for(int d=0;d<2;d++)layout.sum(F[d]);R aggregate[7]={totalU,gravity_work,gravity_torque,quadratureK,initial_Kerror,external_work,external_torque};MPI_Allreduce(MPI_IN_PLACE,aggregate,7,MPI_LONG_DOUBLE,MPI_SUM,comm);totalU=aggregate[0];quadratureK=aggregate[3];
    FILE*thermal=fopen(("full_generating_thermal_r"+suffix+".bin").c_str(),"wbx");need(thermal,"immutable full generating entropy context");uint64_t th[6]={0x4d46544845524d31ULL,1,uint64_t(delta.size()),sizeof(R),uint64_t(rank),1};material_full::write(thermal,th,6);double law[4]={gamma,pressure,eps,0};material_full::write(thermal,law,4);material_full::write(thermal,delta.data(),delta.size());material_full::durable(thermal);
    double solve_start=MPI_Wtime();md::Mass M(layout,state.f.triangle,std::move(blocks));std::array<md::Solution,2>solution;
    for(int d=0;d<2;d++)solution[d]=md::solve(M,F[d]);double solve_seconds=MPI_Wtime()-solve_start;
    R work=0,torque=0,K=0;
    for(int d=0;d<2;d++){md::Vec v(state.global.size()),x(state.global.size());for(size_t i=0;i<v.size();i++){v[i]=state.f.coefficient[i].v[d];x[i]=d==1?state.f.coefficient[i].x[0]:-state.f.coefficient[i].x[1];}auto Ma=M.apply(solution[d].x);work+=layout.dot(v,Ma);torque+=layout.dot(x,Ma);K+=.5L*layout.dot(v,M.apply(v));}
    // Sum Udot=-P0 Vdot and actual prescribed external work separately.
    R udot=0;for(auto v:geometry.vdot)udot-=R(pressure)*v;MPI_Allreduce(MPI_IN_PLACE,&udot,1,MPI_LONG_DOUBLE,MPI_SUM,comm);
    R energy=fabsl(work+udot-aggregate[1]-aggregate[5])/totalU,angular=fabsl(torque-aggregate[2]-aggregate[6])/(2*K),metric_error=fabsl(K-quadratureK)/totalU;
    R high[3]={shadow_mass,shadow_force,-minJ};MPI_Allreduce(MPI_IN_PLACE,high,3,MPI_LONG_DOUBLE,MPI_MAX,comm);uint64_t positive=points;MPI_Allreduce(MPI_IN_PLACE,&positive,1,MPI_UINT64_T,MPI_SUM,comm);
    bool pass=energy<1e-8L&&angular<1e-11L&&metric_error<1e-8L&&fabsl(aggregate[4])/totalU<1e-8L;
    FILE*a=fopen(("full_initial_acceleration_r"+suffix+".bin").c_str(),"wbx");need(a,"actual first full mass-equation acceleration");uint64_t ah[4]={0x4d46414343454c31ULL,1,uint64_t(state.global.size()),sizeof(R)};material_full::write(a,ah,4);for(int d=0;d<2;d++)material_full::write(a,solution[d].x.data(),solution[d].x.size());material_full::durable(a);
    if(rank==0){FILE*f=fopen("full_initial_result.json","wx");need(f,"full native initial operator result");fprintf(f,"{\"scope\":\"actual full fine_integrated shared material field, streamed GPU fixed mass and initial constant-P/Pext gravity RHS, distributed true mass solve; no accepted time step\",\"MPI_ranks\":%d,\"global_scalar_dofs\":%llu,\"global_triangles\":%llu,\"actual_positive_points\":%llu,\"min_J\":%.21Lg,\"physical_U\":%.21Lg,\"physical_K\":%.21Lg,\"initial_velocity_representation_K_error_over_U\":%.21Lg,\"mass_metric_vs_point_K_error_over_U\":%.21Lg,\"global_energy_rate_error_over_U\":%.21Lg,\"global_angular_rate_error_over_2K\":%.21Lg,\"CPU_CUDA_mass_shadow_relative\":%.21Lg,\"CPU_CUDA_force_shadow_relative\":%.21Lg,\"CG_iterations\":[%llu,%llu],\"true_equation_residual\":[%.21Lg,%.21Lg],\"basis_seconds_rank0\":%.9g,\"mass_integration_seconds_rank0\":%.9g,\"CUDA_mass_kernel_seconds_rank0\":%.9g,\"CPU_transform_seconds_rank0\":%.9g,\"factor_solve_seconds_rank0\":%.9g,\"seconds\":%.9g,\"initial_operator_algebra_pass\":%s,\"initial_curvature_checked\":false,\"accepted_steps\":0,\"long_Kepler_validated\":false}\n",ranks,(unsigned long long)layout.global_owned,(unsigned long long)state.global_triangles,(unsigned long long)positive,-high[2],totalU,K,aggregate[4]/totalU,metric_error,energy,angular,high[0],high[1],(unsigned long long)solution[0].iterations,(unsigned long long)solution[1].iterations,solution[0].true_relative_residual,solution[1].true_relative_residual,basis_seconds,integration_seconds,kernel,transform_seconds,solve_seconds,MPI_Wtime()-begin,pass?"true":"false");material_full::durable(f);}
    return pass?0:2;
}
#include "material_full_curvature_driver.h"
#include "material_full_evolution_host.h"
extern "C" int gfs_material_initialize_full(const char*input,int fcomm,double gamma,double pressure,double eps,double box){
    MPI_Comm comm=MPI_Comm_f2c(fcomm);try{need(input&&*input&&eps==.047&&box==24&&gamma>1&&pressure>0,"saved full IC initialization contract");const char*saved=getenv("GFS_MATERIAL_INITIAL_INPUT"),*evolve=getenv("GFS_MATERIAL_EVOLVE");if(evolve&&atoi(evolve)>0){need(saved&&*saved,"full evolution requires retained initialized state");return full_evolution::run(saved,comm,gamma,pressure,eps);}return saved&&saved[0]?curvature_from_saved(saved,comm,gamma,pressure,eps):initialize(input,comm,gamma,pressure,eps);}catch(const std::exception&e){int rank;MPI_Comm_rank(comm,&rank);fprintf(stderr,"rank%d FULL_INITIAL_FATAL %s\n",rank,e.what());MPI_Abort(comm,2);return 2;}
}
