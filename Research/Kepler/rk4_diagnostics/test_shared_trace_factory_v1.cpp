#include <mpi.h>
#include <omp.h>
#include <algorithm>
#include <array>
#include <map>
#include <set>
#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <climits>
#include <chrono>
#include <fstream>
#include <stdexcept>
static void need(bool ok,const char*message){if(!ok)throw std::runtime_error(message);}
#include "material_curved_packet.h"
#include "gfs_enriched_euler.h"
#include "material_shared_trace_factory.h"
namespace st=shared_trace;using R=long double;

static std::vector<STPoint> positive_points(const curved::Packet&p,const st::Factory&f,int rank,int ranks){
    GfsEnrichConfig cfg={};cfg.na=64;cfg.nr=28;std::vector<double>x,w;st::gauss(64,x,w);for(int i=0;i<64;i++){cfg.an[i]=x[i];cfg.aw[i]=w[i];}st::gauss(28,x,w);for(int i=0;i<28;i++){cfg.rn[i]=x[i];cfg.rw[i]=w[i];}
    cfg.field.ic.kind=2;cfg.field.ic.gamma=p.parameters[1];cfg.field.ic.pressure=p.parameters[2];cfg.field.ic.eps=p.parameters[3];cfg.field.ic.cx=cfg.field.ic.cy=12;cfg.field.ic.rin=p.parameters[4];cfg.field.ic.rout=p.parameters[5];
    std::vector<STPoint> out;
    for(size_t s=rank;s<f.triangle.size();s+=ranks){const auto&t=f.triangle[s];double px[3],py[3];for(int a=0;a<3;a++){px[a]=t.vertex[a][0];py[a]=t.vertex[a][1];}GfsMomentPolygon poly={px,py,3,0};GfsEnrichBasis basis={};basis.generator[0]=basis.generator[1]=12;GfsEnrichQuadrature quad;need(gfs_enrich_quadrature_init(poly,&cfg,&basis,&quad),"same original positive quadrature setup");
        for(int q=0;q<gfs_enrich_nquad(poly,&cfg,&quad);q++){STPoint z={};z.triangle=s;z.boundary=-1;need(gfs_enrich_quad(poly,&cfg,&quad,q,&z.X,&z.Y,&z.weight),"same original positive volume rule");if(z.weight==0)continue;need(z.weight>0,"strictly positive volume weight");double u[4];need(gfs_enrich_reference(&cfg,z.X+12,z.Y+12,u,0,0),"unchanged native reference density");z.rho=u[0];R ux=R(t.vertex[1][0])-t.vertex[0][0],uy=R(t.vertex[1][1])-t.vertex[0][1],vx=R(t.vertex[2][0])-t.vertex[0][0],vy=R(t.vertex[2][1])-t.vertex[0][1],dx=R(z.X)-t.vertex[0][0],dy=R(z.Y)-t.vertex[0][1],det=ux*vy-uy*vx;z.l[1]=(dx*vy-dy*vx)/det;z.l[2]=(ux*dy-uy*dx)/det;z.l[0]=1-z.l[1]-z.l[2];need(*std::min_element(z.l,z.l+3)>0,"interior material quadrature nodes");out.push_back(z);}
    }return out;
}
static std::vector<R> raw_mass(st::Factory&f,const std::vector<STPoint>&point,int rank,int ranks){
    std::vector<std::vector<size_t>> by(f.triangle.size());for(size_t i=0;i<point.size();i++)by[point[i].triangle].push_back(i);
    std::vector<R>M(f.triangle.size()*ST_MAX*ST_MAX),mass(f.triangle.size());
    int failed=0;
    #pragma omp parallel for schedule(static) reduction(+:failed)
    for(size_t s=0;s<f.triangle.size();s++)if(s%ranks==size_t(rank)){auto&t=f.triangle[s];R*matrix=M.data()+s*ST_MAX*ST_MAX;for(size_t ip:by[s]){STJet N[ST_MAX];if(!st_shape(t,point[ip],N)){failed++;continue;}R dm=R(point[ip].weight)*point[ip].rho;mass[s]+=dm;for(int a=0;a<t.n;a++)for(int b=0;b<t.n;b++)matrix[a*ST_MAX+b]+=dm*N[a].v*N[b].v;}}
    int global;MPI_Allreduce(&failed,&global,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);need(!global,"all canonical raw mass basis evaluations valid");
    MPI_Allreduce(MPI_IN_PLACE,M.data(),M.size(),MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);MPI_Allreduce(MPI_IN_PLACE,mass.data(),mass.size(),MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    for(size_t s=0;s<f.triangle.size();s++){st::bubble_basis(f.triangle[s],M.data()+s*ST_MAX*ST_MAX,double(mass[s]));for(int j=0;j<2;j++){double scale=0;for(int k=0;k<2;k++)scale=hypot(scale,f.triangle[s].bubble_decode[k][j]);f.min_bubble_scale=std::min(f.min_bubble_scale,scale);}}
    return M;
}
static void checkpoint(const std::string&path,const st::Factory&f,uint64_t epoch){
    FILE*out=fopen(path.c_str(),"wb");need(out,"fresh shared state checkpoint");uint64_t h[8]={0x4746535452433031ULL,1,epoch,f.dof.size(),f.triangle.size(),sizeof(STTriangle),sizeof(STCoefficient),0};need(fwrite(h,sizeof(h),1,out)==1,"shared checkpoint header");
    need(fwrite(f.triangle.data(),sizeof(STTriangle),f.triangle.size(),out)==f.triangle.size(),"checkpoint fixed canonical basis metadata");
    for(size_t i=0;i<f.dof.size();i++){const auto&d=f.dof[i];uint64_t id[4]={d.id,d.parent,uint64_t(d.kind),uint64_t(d.owner)};need(fwrite(id,sizeof(id),1,out)==1&&fwrite(&f.coefficient[i],sizeof(STCoefficient),1,out)==1,"single-owner shared coefficient checkpoint");}fclose(out);
}
static bool restart(const std::string&path,const st::Factory&expected,uint64_t epoch,std::vector<STCoefficient>&caller,std::string&witness){
    FILE*in=fopen(path.c_str(),"rb");if(!in){witness="open";return false;}auto reject=[&](const char*s){witness=s;fclose(in);return false;};uint64_t h[8];
    if(fread(h,sizeof(h),1,in)!=1||h[0]!=0x4746535452433031ULL||h[1]!=1)return reject("header");
    if(h[2]!=epoch)return reject("epoch");if(h[3]!=expected.dof.size()||h[4]!=expected.triangle.size()||h[5]!=sizeof(STTriangle)||h[6]!=sizeof(STCoefficient))return reject("dimensions");
    std::vector<STTriangle> basis(h[4]);if(fread(basis.data(),sizeof(STTriangle),basis.size(),in)!=basis.size()||std::memcmp(basis.data(),expected.triangle.data(),basis.size()*sizeof(STTriangle)))return reject("canonical_basis_metadata");
    std::vector<STCoefficient> candidate(h[3]);std::vector<bool> seen(h[3]);for(size_t i=0;i<candidate.size();i++){uint64_t id[4];STCoefficient value;if(fread(id,sizeof(id),1,in)!=1||fread(&value,sizeof(value),1,in)!=1)return reject("missing_coefficient");if(id[0]>=seen.size()||seen[id[0]])return reject("duplicate_or_invalid_global_id");auto&d=expected.dof[id[0]];if(id[1]!=d.parent||id[2]!=uint64_t(d.kind)||id[3]!=uint64_t(d.owner))return reject("owner_or_basis_identity");for(int k=0;k<2;k++)if(!std::isfinite(value.x[k])||!std::isfinite(value.v[k]))return reject("nonfinite_state");candidate[id[0]]=value;seen[id[0]]=true;}
    if(fgetc(in)!=EOF)return reject("trailing_state");fclose(in);caller.swap(candidate);witness="pass";return true;
}
static double checkpoint_tests(const std::string&dest,const st::Factory&f){
    std::string path=dest+"/shared_state.bin";checkpoint(path,f,1);std::vector<STCoefficient> caller(f.dof.size());for(auto&v:caller)for(int k=0;k<2;k++){v.x[k]=123;v.v[k]=-456;}auto before=caller;std::string witness;
    need(!restart(path,f,2,caller,witness)&&witness=="epoch"&&!std::memcmp(caller.data(),before.data(),caller.size()*sizeof(STCoefficient)),"epoch rejection preserves shared output");
    std::ifstream input(path,std::ios::binary);std::vector<char> bytes((std::istreambuf_iterator<char>(input)),std::istreambuf_iterator<char>());size_t begin=64+f.triangle.size()*sizeof(STTriangle),row=32+sizeof(STCoefficient);
    auto bad=bytes;std::memcpy(bad.data()+begin+row,bad.data()+begin,8);{std::ofstream o(dest+"/duplicate_state.bin",std::ios::binary);o.write(bad.data(),bad.size());}
    need(!restart(dest+"/duplicate_state.bin",f,1,caller,witness)&&witness=="duplicate_or_invalid_global_id"&&!std::memcmp(caller.data(),before.data(),caller.size()*sizeof(STCoefficient)),"duplicate global ID rejection preserves output");
    {std::ofstream o(dest+"/missing_state.bin",std::ios::binary);o.write(bytes.data(),bytes.size()-sizeof(STCoefficient));}
    need(!restart(dest+"/missing_state.bin",f,1,caller,witness)&&witness=="missing_coefficient"&&!std::memcmp(caller.data(),before.data(),caller.size()*sizeof(STCoefficient)),"missing coefficient rejection preserves output");
    need(restart(path,f,1,caller,witness)&&!std::memcmp(caller.data(),f.coefficient.data(),caller.size()*sizeof(STCoefficient)),"shared checkpoint roundtrip bitmatch without IC state regeneration");
    FILE*values=fopen((dest+"/shared_coefficients.bin").c_str(),"wb");need(values&&fwrite(f.coefficient.data(),sizeof(STCoefficient),f.coefficient.size(),values)==f.coefficient.size(),"canonical state value table");fclose(values);return 1;
}
struct TraceStats {double jump=0,bubble=0,reverse=0;uint64_t shared=0,queries=0;};
static TraceStats traces(const st::Factory&f,const std::string&dest){
    TraceStats out;std::map<int,std::vector<std::pair<int,int>>> edge;for(size_t s=0;s<f.triangle.size();s++)for(int a=0;a<3;a++)edge[f.triangle[s].edge[a].poly[1]].push_back({s,a});
    std::ofstream file(dest+"/canonical_trace_checks.csv");file<<"left,right,canonical_t,x_jump,v_jump,bubble_trace\n";
    for(auto&entry:edge){if(entry.second.size()!=2)continue;out.shared++;auto left=entry.second[0],right=entry.second[1];for(int j=0;j<=16;j++){double v[2][2],x[2][2],q=j/16.,b=0;for(int side=0;side<2;side++){auto index=side?right:left;const auto&t=f.triangle[index.first];const auto&e=t.edge[index.second];STPoint p={};p.triangle=index.first;p.boundary=index.second;p.X=e.A[0]+q*e.D[0];p.Y=e.A[1]+q*e.D[1];p.l[t.endpoint[index.second][0]]=1-q;p.l[t.endpoint[index.second][1]]=q;STJet N[ST_MAX];need(st_shape(t,p,N),"canonical boundary basis");b=std::max({b,fabs(N[t.bubble_offset].v),fabs(N[t.bubble_offset+1].v)});need(st_trace(e,f.coefficient.data(),q,x[side],v[side]),"single shared canonical trace evaluation");double rx[2],rv[2];need(st_trace(e,f.coefficient.data(),1-(1-q),rx,rv),"canonical orientation reversal");for(int k=0;k<2;k++)out.reverse=std::max({out.reverse,fabs(rx[k]-x[side][k]),fabs(rv[k]-v[side][k])});}double dx=0,dv=0;for(int k=0;k<2;k++){dx=std::max(dx,fabs(x[0][k]-x[1][k]));dv=std::max(dv,fabs(v[0][k]-v[1][k]));}out.jump=std::max({out.jump,dx,dv});out.bubble=std::max(out.bubble,b);out.queries++;file<<left.first<<','<<right.first<<','<<q<<','<<dx<<','<<dv<<','<<b<<'\n';}}
    return out;
}
#ifndef GFS_SHARED_TRACE_FACTORY_LIBRARY
int main(int argc,char**argv){
    int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);int code=2;
    try{need(argc==3,"original packet NEW destination");auto start=std::chrono::steady_clock::now();auto packet=curved::read(argv[1]);need(packet.record.size()==384&&packet.cells==24,"same384-source24-cell reference geometry");auto f=st::topology(packet,ranks);auto points=positive_points(packet,f,rank,ranks);auto mass=raw_mass(f,points,rank,ranks);st::initialise(f,packet.parameters[3],rank,ranks);
        STPhysics physics={packet.parameters[1],packet.parameters[2],packet.parameters[3],1};std::vector<STOutput> cpu(points.size()),gpu(points.size());
        #pragma omp parallel for schedule(static)
        for(size_t i=0;i<points.size();i++)cpu[i]=st_evaluate(f.triangle.data(),f.coefficient.data(),points[i],physics);
        double seconds=0;int cuda=st_cuda_evaluate(f.triangle.data(),f.triangle.size(),f.coefficient.data(),f.coefficient.size(),points.data(),points.size(),&physics,gpu.data(),&seconds);need(cuda==0,"native CUDA shared-state evaluator");
        R integral[24][12]={},difference[24][4]={};double stat[8]={0},minJ=DBL_MAX;uint64_t invalid=0;
        for(size_t i=0;i<points.size();i++){const auto&q=points[i];const auto&a=cpu[i];const auto&b=gpu[i];if(!a.valid||!b.valid){invalid++;continue;}int c=f.triangle[q.triangle].cell;R dm=R(q.rho)*q.weight;double X=a.point.position[0],Y=a.point.position[1],vx=a.point.velocity[0],vy=a.point.velocity[1];auto v0=st::initial_velocity({q.X,q.Y},physics.eps);
            integral[c][0]+=dm;integral[c][1]+=dm*vx;integral[c][2]+=dm*vy;integral[c][3]+=dm*(R(X)*vy-R(Y)*vx);integral[c][4]+=.5L*dm*(R(vx)*vx+R(vy)*vy);integral[c][5]+=a.U;integral[c][6]+=q.weight*a.point.J;integral[c][7]+=dm*X;integral[c][8]+=dm*Y;integral[c][9]+=dm*(R(X)*X+R(Y)*Y);integral[c][10]+=a.Udot;integral[c][11]+=a.gpower;
            difference[c][0]+=.5L*dm*((R(vx)-v0[0])*(R(vx)+v0[0])+(R(vy)-v0[1])*(R(vy)+v0[1]));difference[c][1]+=R(a.U)-R(q.weight)*physics.P0/(physics.gamma-1);difference[c][2]+=dm*((R(X)-q.X)*vy-(R(Y)-q.Y)*vx+R(q.X)*(vy-v0[1])-R(q.Y)*(vx-v0[0]));
            stat[0]=std::max({stat[0],fabs(vx-v0[0]),fabs(vy-v0[1])});stat[1]=std::max(stat[1],fabs(a.identity_energy)/(a.scale_energy+1e-300));stat[2]=std::max(stat[2],fabs(a.identity_angular)/(a.scale_angular+1e-300));minJ=std::min(minJ,a.point.J);
            for(int k=0;k<2;k++)stat[3]=std::max({stat[3],fabs(a.point.position[k]-b.point.position[k]),fabs(a.point.velocity[k]-b.point.velocity[k])});for(int k=0;k<4;k++)stat[4]=std::max(stat[4],fabs(a.point.F[k]-b.point.F[k]));for(int k=0;k<f.triangle[q.triangle].n;k++)for(int d=0;d<2;d++){stat[5]=std::max(stat[5],fabs(a.force[d][k]-b.force[d][k]));stat[6]=std::max(stat[6],fabs(a.force[d][k]));}
        }
        MPI_Allreduce(MPI_IN_PLACE,integral,24*12,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);MPI_Allreduce(MPI_IN_PLACE,difference,24*4,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);MPI_Allreduce(MPI_IN_PLACE,stat,8,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);MPI_Allreduce(MPI_IN_PLACE,&minJ,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);MPI_Allreduce(MPI_IN_PLACE,&invalid,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);uint64_t n=points.size();MPI_Allreduce(MPI_IN_PLACE,&n,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
        if(rank==0){std::string dest=argv[2];auto trace=traces(f,dest);bool restart_pass=checkpoint_tests(dest,f);R Kerr=0,Uerr=0,Lerr=0;std::ofstream cells(dest+"/initial_cell_integrals.csv");cells<<"cell,mass,U,K,L,K_error_over_U,U_relative_error,L_error_over_2K\n";for(size_t c=0;c<f.cell.size();c++){Kerr=std::max(Kerr,fabsl(difference[c][0])/integral[c][5]);Uerr=std::max(Uerr,fabsl(difference[c][1])/integral[c][5]);Lerr=std::max(Lerr,fabsl(difference[c][2])/(2*integral[c][4]));char line[600];snprintf(line,sizeof(line),"%llu,%.21Lg,%.21Lg,%.21Lg,%.21Lg,%.21Lg,%.21Lg,%.21Lg\n",(unsigned long long)f.cell[c],integral[c][0],integral[c][5],integral[c][4],integral[c][3],difference[c][0]/integral[c][5],difference[c][1]/integral[c][5],difference[c][2]/(2*integral[c][4]));cells<<line;}cells.close();
            bool pass=!invalid&&trace.jump==0&&trace.bubble==0&&trace.reverse==0&&restart_pass&&minJ>0&&Kerr<1e-8L&&Uerr<1e-10L&&Lerr<1e-11L&&stat[1]<1e-11&&stat[2]<1e-11&&stat[5]/stat[6]<1e-11;
            printf("{\"scope\":\"prospective shared-trace factory from original continuous IC and same reference mesh; initial state only, no acceleration/control/evolution\",\"MPI_ranks\":%d,\"OMP_threads\":%d,\"vertices\":%zu,\"canonical_edges\":%zu,\"rank1_edges\":%zu,\"rank2_edges\":%zu,\"P3_shared_dofs\":%zu,\"new_scalar_dofs\":%zu,\"old_broken_scalar_dofs\":%llu,\"old_nominal_conforming_scalar_dofs\":1865,\"triangles\":%zu,\"physical_cells\":%zu,\"positive_points\":%llu,\"invalid_points\":%llu,\"min_J\":%.17g,\"initial_velocity_error\":%.17g,\"K_error_over_U\":%.17Lg,\"U_relative_error\":%.17Lg,\"L_error_over_2K\":%.17Lg,\"point_energy_identity_relative\":%.17g,\"point_angular_identity_relative\":%.17g,\"CPU_CUDA_field_max_difference\":%.17g,\"CPU_CUDA_F_max_difference\":%.17g,\"CPU_CUDA_force_relative_difference\":%.17g,\"canonical_trace_jump\":%.17g,\"bubble_trace\":%.17g,\"orientation_reversal_difference\":%.17g,\"paired_edges_checked\":%llu,\"trace_queries\":%llu,\"checkpoint_and_transaction_pass\":true,\"factory_initial_gate_pass\":%s,\"new_acceleration_and_curvature_validated\":false,\"finite_control_validated\":false,\"CUDA_kernel_seconds_rank0\":%.17g,\"seconds\":%.17g}\n",ranks,omp_get_max_threads(),f.vertex.size(),f.edge.size(),f.rank1,f.rank2,f.polynomial.size(),f.dof.size(),(unsigned long long)packet.dofs,f.triangle.size(),f.cell.size(),(unsigned long long)n,(unsigned long long)invalid,minJ,stat[0],Kerr,Uerr,Lerr,stat[1],stat[2],stat[3],stat[4],stat[5]/stat[6],trace.jump,trace.bubble,trace.reverse,(unsigned long long)trace.shared,(unsigned long long)trace.queries,pass?"true":"false",seconds,std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());code=pass?0:2;
        }MPI_Bcast(&code,1,MPI_INT,0,MPI_COMM_WORLD);
    }catch(const std::exception&e){fprintf(stderr,"rank%d FAIL %s\n",rank,e.what());MPI_Abort(MPI_COMM_WORLD,2);}
    MPI_Finalize();return code;
}
#endif
