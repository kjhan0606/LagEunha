// Actual smooth material evolution of the retained finite CONTROL.
// This restricted patch is not a full-disk Kepler or shock calculation.
#define GFS_GEOMETRIC_PRESSURE_LIBRARY
#include "test_geometric_pressure_v1.cpp"
#include <sstream>
#include <memory>
#include <sys/stat.h>
#include <unistd.h>
extern "C" int gp_persistent_create(const STTriangle*,size_t,size_t,const STPoint*,size_t,void**);
extern "C" int gp_persistent_evaluate(void*,const STCoefficient*,const STPhysics*,const GPInput*,GPOutput*,double*);
extern "C" void gp_persistent_destroy(void*);
static double clock_seconds(){return MPI_Wtime();}
static uint64_t hash_bytes(const void*ptr,size_t n,uint64_t h=14695981039346656037ULL){
    const auto*p=static_cast<const unsigned char*>(ptr);for(size_t i=0;i<n;i++){h^=p[i];h*=1099511628211ULL;}return h;
}
template<class T>static uint64_t hash_vector(const std::vector<T>&v){return hash_bytes(v.data(),v.size()*sizeof(T));}
static std::string read_text(const std::string&path){std::ifstream f(path,std::ios::binary);need(bool(f),"retained provenance text");return {std::istreambuf_iterator<char>(f),{}};}
struct EntropyPoint {double U0=0,jhat0=0;};
struct Fragment {uint64_t slot=0,source=0,target=0,context=0;};
struct TimeLedger {
    // Time is an offset from finite CONTROL, not from a fresh disk IC.
    double time=0,h=2e-4;uint64_t epoch=1,stage=0,fragment_epoch=1;
    R external_work=0,external_torque=0,thermal_work=0;
};
struct MaterialEvolution {
    st::Factory f;std::vector<STCoefficient> generating,stage_base;
    std::vector<STPoint> point;std::vector<EntropyPoint> entropy;
    rhs::Vec generating_delta,raw_mass;std::vector<Fragment> fragment;
    std::string overlap_provenance,target_provenance,domain_provenance;
    STPhysics physics={};TimeLedger ledger;int owner_ranks=0,rank=0;
};
static std::vector<std::string> split(const std::string&s){std::vector<std::string>v;std::istringstream in(s);std::string x;while(std::getline(in,x,','))v.push_back(x);return v;}
static void import_fragments(MaterialEvolution&m,const std::string&root){
    auto coeff=m.f.coefficient;auto basis=m.f.triangle;auto heat=m.entropy;
    m.overlap_provenance=read_text(root+"/shared_transfer_incidence_v1/all559_overlap_moments.csv");
    m.target_provenance=read_text(root+"/shared_transfer_incidence_v1/all24_completed_target_references.csv");
    const char*manifest=getenv("GFS_EVOLUTION_RETAINED_MANIFEST");need(manifest,"frozen retained-input SHA256 manifest");m.domain_provenance=read_text(manifest);
    std::istringstream in(m.overlap_provenance);std::string line;std::getline(in,line);
    std::map<std::pair<uint64_t,uint64_t>,std::set<std::string>> moments;std::set<uint64_t>sources,targets;
    while(std::getline(in,line)){auto v=split(line);need(v.size()==10&&v[8]=="True","original same-record overlap reference");
        uint64_t slot=std::stoull(v[0]),s=std::stoull(v[1]),t=std::stoull(v[2]);need(s<m.f.triangle.size(),"fragment source exists");
        auto&seen=moments[{s,t}];need(seen.insert(v[3]).second,"unique source-target moment");
        if(v[3]=="mass"){need(v[5][0]!='-'&&v[5][0]!='0',"positive retained overlap mass");m.fragment.push_back({slot,s,t,uint64_t(m.f.triangle[s].cell)});sources.insert(s);targets.insert(t);}
        if(v[3]=="physical_U")need(v[5][0]!='-'&&v[5][0]!='0',"positive retained overlap physical heat");
    }
    need(m.fragment.size()==559&&moments.size()==559&&sources.size()==384&&targets.size()==24,"complete fragment view");
    for(auto&v:moments)need(v.second.size()==12,"all twelve retained moments bound to each fragment");
    std::istringstream tin(m.target_provenance);std::getline(tin,line);std::set<uint64_t>linked;
    while(std::getline(tin,line)){auto v=split(line);need(v.size()==5&&v[2].size()==64&&targets.count(std::stoull(v[0])),"original target packet SHA binding");linked.insert(std::stoull(v[0]));}
    need(linked==targets,"all targets linked to completed packet hashes");
    need(hash_vector(coeff)==hash_vector(m.f.coefficient)&&hash_vector(basis)==hash_vector(m.f.triangle)&&hash_vector(heat)==hash_vector(m.entropy),"pure relabel leaves complete x/v/basis/heat bitwise unchanged");
}
static void restore_frozen_points(MaterialEvolution&m,const std::string&dir,int rank,int ranks){
    struct Frozen {uint64_t triangle;double X,Y,weight,rho,U;};
    std::vector<std::vector<Frozen>>by(m.f.triangle.size());
    for(int owner=0;owner<4;owner++){
        FILE*f=fopen((dir+"/frozen_physical_rank"+std::to_string(owner)+".bin").c_str(),"rb");need(f,"retained quadrature and physical U0");uint64_t h[6];
        need(fread(h,sizeof(h),1,f)==1&&h[0]==0x47465346524f5a32ULL&&h[1]==2&&h[2]==1&&h[4]==sizeof(Frozen)&&h[5]==uint64_t(owner),"original frozen point header");
        std::vector<Frozen>data(h[3]);need(fread(data.data(),sizeof(Frozen),data.size(),f)==data.size()&&fgetc(f)==EOF,"complete saved positive quadrature");fclose(f);
        for(auto&q:data){need(q.triangle<by.size()&&q.triangle%4==size_t(owner)&&q.U>0&&q.rho>0&&q.weight>0,"frozen source ownership/mass/heat");if(q.triangle%ranks==size_t(rank))by[q.triangle].push_back(q);}
    }
    for(size_t s=rank;s<by.size();s+=ranks){const auto&t=m.f.triangle[s];for(auto&q:by[s]){
        STPoint p={};p.triangle=s;p.boundary=-1;p.X=q.X;p.Y=q.Y;p.weight=q.weight;p.rho=q.rho;
        R ux=R(t.vertex[1][0])-t.vertex[0][0],uy=R(t.vertex[1][1])-t.vertex[0][1],vx=R(t.vertex[2][0])-t.vertex[0][0],vy=R(t.vertex[2][1])-t.vertex[0][1],dx=R(p.X)-t.vertex[0][0],dy=R(p.Y)-t.vertex[0][1],det=ux*vy-uy*vx;
        p.l[1]=(dx*vy-dy*vx)/det;p.l[2]=(ux*dy-uy*dx)/det;p.l[0]=1-p.l[1]-p.l[2];need(*std::min_element(p.l,p.l+3)>0,"restored reference interior point");m.point.push_back(p);m.entropy.push_back({q.U,0});
    }}
    uint64_t n=m.point.size();MPI_Allreduce(MPI_IN_PLACE,&n,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);need(n==817152,"same fixed material quadrature count without regeneration");
}
// Assemble exactly the same corrected volume functional as the CONTROL.
static State volume_context(const MaterialEvolution&m){
    State s;s.geometry=gp::geometry(m.f);s.thermal.resize(m.f.cell.size());s.input.resize(m.point.size());s.cpu.resize(m.point.size());s.gpu.resize(m.point.size());
    std::vector<GfsMaterialPoint12>field(m.point.size());int failed=0;
    #pragma omp parallel for schedule(static) reduction(+:failed)
    for(size_t i=0;i<m.point.size();i++){STJet N[ST_MAX];const auto&q=m.point[i];if(!st_shape(m.f.triangle[q.triangle],q,N)||!st_point(m.f.triangle[q.triangle],m.f.coefficient.data(),N,&field[i]))failed++;}
    MPI_Allreduce(MPI_IN_PLACE,&failed,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);need(!failed,"evolved finite positive material Jacobians");
    R sums[24][4]={};for(size_t i=0;i<m.point.size();i++){const auto&q=m.point[i];const auto&p=field[i];auto*sum=sums[m.f.triangle[q.triangle].cell];sum[0]+=q.weight;sum[1]+=R(q.weight)*p.J;sum[2]+=R(q.weight)*p.J*(R(p.grad_v[0])+p.grad_v[3]);sum[3]+=R(q.weight)*q.rho;}
    MPI_Allreduce(MPI_IN_PLACE,sums,96,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for(int c=0;c<24;c++){auto&t=s.thermal[c];t.A=sums[c][0];t.QJ=sums[c][1];t.QJdot=sums[c][2];t.mass=sums[c][3];need(t.A>0&&t.mass>0,"fixed positive source material measures");t.delta=(s.geometry.volume[c]-t.QJ)/t.A;t.delta_dot=(s.geometry.vdot[c]-t.QJdot)/t.A;}
    // Retain these already computed fields for the entropy/pressure pass.
    s.old.resize(field.size());for(size_t i=0;i<field.size();i++)s.old[i].point=field[i];return s;
}
static void bind_generating_heat(MaterialEvolution&m){
    auto s=volume_context(m);m.generating=m.f.coefficient;m.generating_delta.resize(24);
    for(int c=0;c<24;c++)m.generating_delta[c]=double(s.thermal[c].delta);
    for(size_t i=0;i<m.point.size();i++){int c=m.f.triangle[m.point[i].triangle].cell;m.entropy[i].jhat0=s.old[i].point.J+double(m.generating_delta[c]);need(m.entropy[i].jhat0>0,"actual generating J plus original geometric delta");}
}
static R boundary_pressure(const MaterialEvolution&m,const STPoint&q,const GfsMaterialPoint12&p,R delta){
    const auto&t=m.f.triangle[q.triangle];STJet N[ST_MAX];GfsMaterialPoint12 generating;
    need(st_shape(t,q,N)&&st_point(t,m.generating.data(),N,&generating),"immutable generating boundary field");
    double j0=generating.J+double(m.generating_delta[t.cell]),j=p.J+double(delta);
    need(j0>0&&j>0,"positive boundary corrected volume");
    // Initial u0 function is the old CONTROL once; later states use its
    // immutable generating geometry, never the current J as initial heat.
    double u0=m.physics.P0/(m.physics.gamma-1)*pow(generating.J,1-m.physics.gamma);
    double u=u0*pow(j/j0,1-m.physics.gamma);need(u>0&&std::isfinite(u),"positive evolved boundary thermal state");
    return (m.physics.gamma-1)*u/j;
}
#include "material_shared_evolution_assembly.h"
struct GPUContext {
    void*ptr=nullptr;size_t calls=0;double kernel=0,wall=0;
    explicit GPUContext(const MaterialEvolution&m){need(!gp_persistent_create(m.f.triangle.data(),m.f.triangle.size(),m.f.coefficient.size(),m.point.data(),m.point.size(),&ptr),"persistent native CUDA allocation");}
    ~GPUContext(){gp_persistent_destroy(ptr);}
    GPUContext(const GPUContext&)=delete;GPUContext&operator=(const GPUContext&)=delete;
};
static State prepare_evolved(const MaterialEvolution&m,GPUContext*gpu){
    auto s=volume_context(m);R projection[24]={};
    for(size_t i=0;i<m.point.size();i++){
        const auto&q=m.point[i];int c=m.f.triangle[q.triangle].cell;const auto&t=s.thermal[c];auto&in=s.input[i];
        in.deltaJ=t.delta;in.deltaJdot=t.delta_dot;double j=s.old[i].point.J+in.deltaJ;
        need(j>0&&m.entropy[i].jhat0>0,"positive evolved corrected volume");
        in.U=m.entropy[i].U0*pow(j/m.entropy[i].jhat0,1-m.physics.gamma);
        need(in.U>0&&std::isfinite(in.U),"finite smooth entropy evaluation without heat reset");
        double nu=q.weight*j,target=m.physics.P0*nu;
        double dp=(fma(m.physics.gamma-1,in.U,-target)+fma(-m.physics.P0,nu,target))/nu;
        projection[c]+=R(q.weight)*dp;
    }
    MPI_Allreduce(MPI_IN_PLACE,projection,24,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for(int c=0;c<24;c++)s.thermal[c].mean_dp=projection[c]/s.thermal[c].A;
    for(size_t i=0;i<m.point.size();i++)s.input[i].mean_deltaP=s.thermal[m.f.triangle[m.point[i].triangle].cell].mean_dp;
    #pragma omp parallel for schedule(static)
    for(size_t i=0;i<m.point.size();i++)s.cpu[i]=gp_evaluate(m.f.triangle.data(),m.f.coefficient.data(),m.point[i],m.physics,s.input[i]);
    if(gpu){double begin=clock_seconds();need(!gp_persistent_evaluate(gpu->ptr,m.f.coefficient.data(),&m.physics,s.input.data(),s.gpu.data(),&s.kernel),"actual persistent CUDA evolved stress/work");gpu->wall+=clock_seconds()-begin;gpu->kernel+=s.kernel;gpu->calls++;}
    R metrics[24][4]={};
    for(size_t i=0;i<m.point.size();i++){
        const auto&q=m.point[i];const auto&o=s.cpu[i];need(o.native.valid,"positive evolved CPU shadow");int c=m.f.triangle[q.triangle].cell;auto*t=metrics[c];t[0]+=o.nu;t[1]+=o.native.U;
        t[2]=std::max(t[2],fabsl(R(o.nu)/(R(q.weight)*(R(o.native.point.J)+s.input[i].deltaJ))-1));
        t[3]=std::max(t[3],fabsl(R(o.native.pressure)*o.nu/((R(m.physics.gamma)-1)*o.native.U)-1));
        s.minJ=std::min(s.minJ,R(o.native.point.J));s.minnu=std::min(s.minnu,R(o.nu));s.minU=std::min(s.minU,R(o.native.U));
        if(gpu){const auto&g=s.gpu[i];need(g.native.valid&&g.native.U==s.input[i].U,"actual GPU uses the same evolved heat");
            for(int d=0;d<2;d++)for(int j=0;j<m.f.triangle[q.triangle].n;j++){s.force_difference=std::max(s.force_difference,fabsl(R(o.native.force[d][j])-g.native.force[d][j]));s.force_scale=std::max(s.force_scale,fabsl(R(o.native.force[d][j])));}}
    }
    R totals[24][2],errors[24][2];for(int c=0;c<24;c++){for(int j=0;j<2;j++){totals[c][j]=metrics[c][j];errors[c][j]=metrics[c][j+2];}}
    MPI_Allreduce(MPI_IN_PLACE,totals,48,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);MPI_Allreduce(MPI_IN_PLACE,errors,48,MPI_LONG_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    for(int c=0;c<24;c++){s.thermal[c].sum_nu=totals[c][0];s.thermal[c].U=totals[c][1];s.thermal[c].max_volume_error=errors[c][0];s.thermal[c].max_pressure_error=errors[c][1];}
    R low[3]={s.minJ,s.minnu,s.minU},high[2]={s.force_difference,s.force_scale};MPI_Allreduce(MPI_IN_PLACE,low,3,MPI_LONG_DOUBLE,MPI_MIN,MPI_COMM_WORLD);MPI_Allreduce(MPI_IN_PLACE,high,2,MPI_LONG_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    s.minJ=low[0];s.minnu=low[1];s.minU=low[2];s.force_difference=high[0];s.force_scale=high[1];return s;
}
// Full versioned restart, including immutable thermal/context state and
// rollback base. FNV checks payload integrity; input SHA256 binds provenance.
struct RestartIO {
    FILE*f=nullptr;bool writing;uint64_t hash=14695981039346656037ULL;
    RestartIO(const std::string&p,bool w):writing(w){f=fopen(p.c_str(),w?"wbx":"rb");need(f,"new durable restart / real restart input");}
    ~RestartIO(){if(f)fclose(f);}
    void bytes(void*p,size_t n){if(writing)need(fwrite(p,1,n,f)==n,"restart payload write");else need(fread(p,1,n,f)==n,"restart payload read");hash=hash_bytes(p,n,hash);}
    template<class T>void scalar(T&v){bytes(&v,sizeof(v));}
    template<class T>void vec(std::vector<T>&v){uint64_t n=v.size();scalar(n);need(n<100000000,"bounded restart vector");if(!writing)v.resize(n);bytes(v.data(),n*sizeof(T));}
    void str(std::string&s){uint64_t n=s.size();scalar(n);need(n<100000000,"bounded restart provenance");if(!writing)s.resize(n);bytes(s.data(),n);}
    void finish(){if(writing){need(fwrite(&hash,sizeof(hash),1,f)==1&&fflush(f)==0&&fsync(fileno(f))==0,"durable restart checksum flush");}else{uint64_t expected;need(fread(&expected,sizeof(expected),1,f)==1&&expected==hash&&fgetc(f)==EOF,"complete verified restart payload");}need(fclose(f)==0,"restart close");f=nullptr;}
};
static void restart_payload(MaterialEvolution&m,const std::string&p,bool writing){
    RestartIO io(p,writing);uint64_t magic=0x47465345564f4c31ULL,version=1,ldsize=sizeof(R),tri=sizeof(STTriangle),point=sizeof(STPoint);
    io.scalar(magic);io.scalar(version);io.scalar(ldsize);io.scalar(tri);io.scalar(point);need(magic==0x47465345564f4c31ULL&&version==1&&ldsize==sizeof(R)&&tri==sizeof(STTriangle)&&point==sizeof(STPoint),"restart version/ABI identity");
    io.scalar(m.owner_ranks);io.scalar(m.rank);io.scalar(m.physics);io.scalar(m.ledger);
    io.vec(m.f.triangle);io.vec(m.f.cell);io.vec(m.f.dof);io.vec(m.f.coefficient);
    io.vec(m.generating);io.vec(m.generating_delta);io.vec(m.point);io.vec(m.entropy);io.vec(m.raw_mass);io.vec(m.fragment);io.vec(m.stage_base);
    io.str(m.overlap_provenance);io.str(m.target_provenance);io.str(m.domain_provenance);io.finish();
    need(m.f.coefficient.size()==3797&&m.f.triangle.size()==384&&m.f.cell.size()==24&&m.entropy.size()==m.point.size()&&m.generating.size()==m.f.coefficient.size()&&m.fragment.size()==559,"restored complete native state dimensions");
}
struct Endpoint {R mass=0,K=0,U=0,potential=0,L=0;};
static Endpoint endpoint(const MaterialEvolution&m,const State&s){
    R value[5]={};for(size_t i=0;i<m.point.size();i++){const auto&q=m.point[i];const auto&o=s.cpu[i].native;R dm=R(q.weight)*q.rho,x=o.point.position[0],y=o.point.position[1],vx=o.point.velocity[0],vy=o.point.velocity[1];value[0]+=dm;value[1]+=.5L*dm*(vx*vx+vy*vy);value[2]+=o.U;value[3]-=dm/sqrtl(x*x+y*y+R(m.physics.eps)*m.physics.eps);value[4]+=dm*(x*vy-y*vx);}
    MPI_Allreduce(MPI_IN_PLACE,value,5,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);return {value[0],value[1],value[2],value[3],value[4]};
}
struct Rate {std::vector<STCoefficient>k;R external_work=0,external_torque=0,heat=0;};
struct Timings {double prepare=0,assembly=0,solve=0,audit=0,io=0;R solve_residual=0;uint64_t rhs_calls=0;};
static void create_dir(const std::string&p,int rank){if(!rank)need(mkdir(p.c_str(),0770)==0,"fresh stage output directory");MPI_Barrier(MPI_COMM_WORLD);}
static Rate evaluate_rhs(MaterialEvolution&m,const rhs::Mass&M,GPUContext&gpu,Timings&timing,const std::string&dest,int rank,int ranks){
    double start=clock_seconds();auto s=prepare_evolved(m,&gpu);timing.prepare+=clock_seconds()-start;
    start=clock_seconds();auto load=evolved_assemble(m.f,m.point,s,m.physics,rank,ranks,m);timing.assembly+=clock_seconds()-start;
    std::array<rhs::Vec,2>a={rhs::Vec(m.f.dof.size()),rhs::Vec(m.f.dof.size())};R residual=0;
    start=clock_seconds();if(rank==0)for(int d=0;d<2;d++){if(rhs::dot(load.F[d],load.F[d])==0)continue;auto solution=rhs::solve(M,load.F[d]);a[d]=std::move(solution.x);residual=std::max(residual,solution.true_relative_residual);}
    MPI_Bcast(&residual,1,MPI_LONG_DOUBLE,0,MPI_COMM_WORLD);for(int d=0;d<2;d++)MPI_Bcast(a[d].data(),a[d].size(),MPI_LONG_DOUBLE,0,MPI_COMM_WORLD);timing.solve+=clock_seconds()-start;
    timing.solve_residual=std::max(timing.solve_residual,residual);need(residual<1e-16L,"original true mass solve residual gate");
    start=clock_seconds();auto result=evolved_audit(m.f,m.point,s,m.physics,M,load,a,false,rank,dest);timing.audit+=clock_seconds()-start;need(result.pass,"original semidiscrete energy/angular/hybrid/CUDA gates at evolved stage");
    Rate r;r.k.resize(m.f.coefficient.size());for(size_t j=0;j<r.k.size();j++)for(int d=0;d<2;d++){r.k[j].x[d]=m.f.coefficient[j].v[d];r.k[j].v[d]=double(a[d][j]);}
    R heat=0;for(auto&o:s.gpu)heat+=o.native.Udot;MPI_Allreduce(&heat,&r.heat,1,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for(size_t t=0;t<m.f.triangle.size();t++)for(int j=0;j<m.f.triangle[t].n;j++)for(int d=0;d<2;d++){int id=m.f.triangle[t].id[j];R f=load.external[(t*2+d)*ST_MAX+j];r.external_work+=m.f.coefficient[id].v[d]*f;r.external_torque+=(d==1?m.f.coefficient[id].x[0]:-m.f.coefficient[id].x[1])*f;}
    timing.rhs_calls++;return r;
}
static void displaced(MaterialEvolution&m,const std::vector<STCoefficient>&base,const Rate&r,R dt){
    for(size_t i=0;i<base.size();i++)for(int d=0;d<2;d++){m.f.coefficient[i].x[d]=double(R(base[i].x[d])+dt*r.k[i].x[d]);m.f.coefficient[i].v[d]=double(R(base[i].v[d])+dt*r.k[i].v[d]);}
}
static void advance(MaterialEvolution&m,const rhs::Mass&M,GPUContext&gpu,Timings&timing,const std::string&dest,int rank,int ranks){
    create_dir(dest,rank);m.stage_base=m.f.coefficient;auto before=m.ledger;R h=m.ledger.h;
    try{
        auto initial=prepare_evolved(m,nullptr);auto begin=endpoint(m,initial);std::array<Rate,4>k;
        for(int stage=0;stage<4;stage++){
            if(stage)displaced(m,m.stage_base,k[stage-1],stage==3?h:h/2);m.ledger.stage=stage+1;
            std::string path=dest+"/stage"+std::to_string(stage+1);create_dir(path,rank);
            // Small sparse state witness. The full immutable inputs live in
            // imported/accepted restart; no per-point stage dumps.
            FILE*f=fopen((path+"/coeff_rank"+std::to_string(rank)+".bin").c_str(),"wbx");need(f&&fwrite(m.f.coefficient.data(),sizeof(STCoefficient),m.f.coefficient.size(),f)==m.f.coefficient.size(),"stage x/v witness");fclose(f);
            k[stage]=evaluate_rhs(m,M,gpu,timing,path,rank,ranks);
        }
        const R weights[4]={1.L/6,1.L/3,1.L/3,1.L/6};
        for(size_t i=0;i<m.f.coefficient.size();i++)for(int d=0;d<2;d++){R x=m.stage_base[i].x[d],v=m.stage_base[i].v[d];for(int j=0;j<4;j++){x+=h*weights[j]*k[j].k[i].x[d];v+=h*weights[j]*k[j].k[i].v[d];}m.f.coefficient[i].x[d]=double(x);m.f.coefficient[i].v[d]=double(v);}
        R work=0,torque=0,thermal=0;for(int j=0;j<4;j++){work+=h*weights[j]*k[j].external_work;torque+=h*weights[j]*k[j].external_torque;thermal+=h*weights[j]*k[j].heat;}
        auto final=prepare_evolved(m,nullptr);auto end=endpoint(m,final);
        R energy=(end.K+end.U+end.potential)-(begin.K+begin.U+begin.potential)-work,angular=end.L-begin.L-torque,heat=end.U-begin.U-thermal;
        // Finite RK defects are separate from the RHS work identities.
        bool pass=end.mass==begin.mass&&fabsl(energy)/begin.U<1e-8L&&fabsl(angular)/(2*begin.K)<1e-11L&&fabsl(heat)/begin.U<1e-8L;
        if(rank==0){auto trace=traces(m.f,dest);pass=pass&&trace.jump==0&&trace.bubble==0;
            FILE*f=fopen((dest+"/step.json").c_str(),"w");need(f,"finite RK defect report");fprintf(f,"{\"RK4_stages\":4,\"h\":%.17g,\"CONTROL_time_before\":%.17g,\"CONTROL_time_after\":%.17g,\"mass_before\":%.21Lg,\"mass_after\":%.21Lg,\"U_before\":%.21Lg,\"U_after\":%.21Lg,\"external_work\":%.21Lg,\"external_torque\":%.21Lg,\"Udot_RK_integral\":%.21Lg,\"energy_defect_over_U\":%.21Lg,\"angular_defect_over_2K\":%.21Lg,\"thermal_increment_defect_over_U\":%.21Lg,\"shared_trace_jump\":%.17g,\"bubble_trace\":%.17g,\"finite_step_gate_pass\":%s,\"long_Kepler_validated\":false}\n",m.ledger.h,before.time,before.time+m.ledger.h,begin.mass,end.mass,begin.U,end.U,work,torque,thermal,energy/begin.U,angular/(2*begin.K),heat/begin.U,trace.jump,trace.bubble,pass?"true":"false");fclose(f);}
        MPI_Bcast(&pass,1,MPI_C_BOOL,0,MPI_COMM_WORLD);need(pass,"finite RK4 endpoint gates without energy adjustment");
        m.ledger.time+=m.ledger.h;m.ledger.epoch++;m.ledger.stage=0;m.ledger.external_work+=work;m.ledger.external_torque+=torque;m.ledger.thermal_work+=thermal;
    }catch(...){restart_payload(m,dest+"/failure_rank"+std::to_string(rank)+".bin",true);m.f.coefficient=m.stage_base;m.ledger=before;throw;}
}
template<class T>static bool same_values(const std::vector<T>&a,const std::vector<T>&b){return a.size()==b.size()&&(!a.size()||!memcmp(a.data(),b.data(),a.size()*sizeof(T)));}
// binary80 has non-value padding; mass reconstruction compares its values.
static bool same_values(const std::vector<R>&a,const std::vector<R>&b){return a==b;}
int main(int argc,char**argv){
    int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    try{
        need(argc==5&&ranks==4,"packet finite_saved root NEW_output; four original owner ranks");double start=clock_seconds();std::string root=argv[3],dest=argv[4];auto packet=curved::read(argv[1]);MaterialEvolution m;
        m.owner_ranks=ranks;m.rank=rank;m.f=load_saved(argv[2],packet);m.physics={packet.parameters[1],packet.parameters[2],packet.parameters[3],1};
        restore_frozen_points(m,root+"/geometric_pressure_v2_mpi4/finite",rank,ranks);bind_generating_heat(m);import_fragments(m,root);
        auto entropy_hash=hash_vector(m.entropy),basis_hash=hash_vector(m.f.triangle),generator_hash=hash_vector(m.generating);auto imported=prepare_evolved(m,nullptr);
        for(size_t i=0;i<m.entropy.size();i++)need(imported.input[i].U==m.entropy[i].U0,"initial smooth entropy matches every saved physical U bitwise");
        m.raw_mass=saved_raw_mass(m.f,m.point);rhs::Mass M(m.f,m.raw_mass);auto mass_hash=hash_vector(M.value);m.stage_base=m.f.coefficient;
        restart_payload(m,dest+"/imported_rank"+std::to_string(rank)+".bin",true);Timings timing;GPUContext gpu(m);
        advance(m,M,gpu,timing,dest+"/accepted_first",rank,ranks);
        std::string saved=dest+"/accepted_rank"+std::to_string(rank)+".bin";double io=clock_seconds();restart_payload(m,saved,true);timing.io+=clock_seconds()-io;
        advance(m,M,gpu,timing,dest+"/uninterrupted_second",rank,ranks);
        MaterialEvolution resumed;io=clock_seconds();restart_payload(resumed,saved,false);timing.io+=clock_seconds()-io;
        need(resumed.owner_ranks==ranks&&resumed.rank==rank&&resumed.ledger.stage==0,"restart owner/accepted epoch");
        rhs::Mass restarted_mass(resumed.f,resumed.raw_mass);need(same_values(M.value,restarted_mass.value),"restart reconstructs identical fixed material mass CSR");
        GPUContext restarted_gpu(resumed);advance(resumed,restarted_mass,restarted_gpu,timing,dest+"/restarted_second",rank,ranks);
        bool same=same_values(m.f.coefficient,resumed.f.coefficient)&&same_values(m.entropy,resumed.entropy)&&same_values(m.generating,resumed.generating)&&same_values(m.point,resumed.point)&&same_values(m.f.triangle,resumed.f.triangle)&&same_values(m.f.dof,resumed.f.dof)&&same_values(m.fragment,resumed.fragment)&&m.overlap_provenance==resumed.overlap_provenance&&m.target_provenance==resumed.target_provenance&&m.domain_provenance==resumed.domain_provenance;
        same=same&&m.ledger.time==resumed.ledger.time&&m.ledger.epoch==resumed.ledger.epoch&&m.ledger.external_work==resumed.ledger.external_work&&m.ledger.external_torque==resumed.ledger.external_torque&&m.ledger.thermal_work==resumed.ledger.thermal_work;
        need(entropy_hash==hash_vector(m.entropy)&&basis_hash==hash_vector(m.f.triangle)&&generator_hash==hash_vector(m.generating)&&mass_hash==hash_vector(M.value),"immutable entropy/basis/generator/fixed mass survived actual updates");
        int all=same;MPI_Allreduce(MPI_IN_PLACE,&all,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);need(all&&timing.rhs_calls==12&&gpu.calls+restarted_gpu.calls==12,"actual restart continuation bitmatch; exactly twelve real CUDA RHS");
        restart_payload(m,dest+"/final_rank"+std::to_string(rank)+".bin",true);
        if(rank==0)printf("{\"scope\":\"retained finite 384-source/24-cell material patch; full-field fragment view; actual RK4 and actual restart continuation\",\"MPI_ranks\":%d,\"OMP_threads\":%d,\"fragments\":559,\"fixed_material_points\":817152,\"h\":%.17g,\"actual_steps_executed\":3,\"actual_RHS_executed\":%llu,\"accepted_CONTROL_time_offset\":%.17g,\"restart_same_backend_bitmatch\":true,\"heat_reset\":false,\"spin_reservoir\":false,\"independent_target_velocity_fields\":false,\"fragment_view_epoch\":1,\"new_time_power_tessellation_claimed\":false,\"persistent_allocations_per_rank\":10,\"static_topology_point_uploads_per_context\":2,\"mass_true_relative_residual_max\":%.21Lg,\"prepare_seconds_rank0\":%.9g,\"assembly_seconds_rank0\":%.9g,\"solve_seconds_rank0\":%.9g,\"audit_seconds_rank0\":%.9g,\"restart_IO_seconds_rank0\":%.9g,\"CUDA_kernel_seconds_rank0\":%.9g,\"CUDA_transfer_and_kernel_seconds_rank0\":%.9g,\"seconds\":%.9g,\"bundle_pass\":true,\"long_Kepler_validated\":false}\n",ranks,omp_get_max_threads(),m.ledger.h,(unsigned long long)timing.rhs_calls,m.ledger.time,timing.solve_residual,timing.prepare,timing.assembly,timing.solve,timing.audit,timing.io,gpu.kernel+restarted_gpu.kernel,gpu.wall+restarted_gpu.wall,clock_seconds()-start);
    }catch(const std::exception&e){fprintf(stderr,"rank%d EVOLUTION_FAIL %s\n",rank,e.what());MPI_Abort(MPI_COMM_WORLD,2);}
    MPI_Finalize();return 0;
}
