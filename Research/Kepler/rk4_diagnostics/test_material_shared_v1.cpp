#include <mpi.h>
#include <cfloat>
#include <array>
#include <algorithm>
#include <numeric>
#define GFS_ENRICH_TEST_LIBRARY
#include "test_enriched_euler_v1.cpp"
#include "gfs_material_map.h"
#include "gfs_material_triangle.h"
#include "gfs_ale_fan.h"
using Pair=std::array<double,2>;
struct Tri {GfsMaterialTriangle t;int vertex[3],cell,component;};
struct SharedMesh {
    std::vector<Pair> node;std::vector<Tri> tri;std::vector<Cell> cell;
    std::map<std::pair<int,int>,std::vector<int>> edge;int components=0;
};
static void velocity(const GfsEnrichConfig &cfg,double x,double y,double v[2],double B[4]){
    auto &ic=cfg.field.ic;if(ic.kind==1){x-=ic.cx;y-=ic.cy;v[0]=ic.vx+ic.B[0]*x+ic.B[1]*y;v[1]=ic.vy+ic.B[2]*x+ic.B[3]*y;for(int k=0;k<4;k++)B[k]=ic.B[k];}
    else{x-=ic.cx;y-=ic.cy;double z=x*x+y*y+ic.eps*ic.eps,o=pow(z,-.75),ox=-1.5*x*o/z,oy=-1.5*y*o/z;
        v[0]=-y*o;v[1]=x*o;B[0]=-y*ox;B[1]=-o-y*oy;B[2]=o+x*ox;B[3]=x*oy;}
}
static SharedMesh mesh_from(const std::map<uint64_t,Cell> &cells){
    SharedMesh m;std::map<std::pair<double,double>,int> keys;
    auto add=[&](double x,double y){auto key=std::make_pair(x,y);auto it=keys.find(key);if(it!=keys.end())return it->second;int id=m.node.size();keys[key]=id;m.node.push_back({x,y});return id;};
    for(auto &entry:cells){
        auto c=entry.second;int cid=m.cell.size();m.cell.push_back(c);GfsMomentPolygon poly={c.x.data(),c.y.data(),(int)c.x.size(),0};GfsMomentDensity g={};need(gfs_moment_geometry(poly,&g),"actual shared-mesh polygon");
        int center=add(c.gx-12+g.gx,c.gy-12+g.gy);
        for(size_t e=0;e<c.x.size();e++){
            size_t j=(e+1)%c.x.size();if(c.x[e]==c.x[j] && c.y[e]==c.y[j])continue;
            Tri t={};t.cell=cid;t.vertex[0]=center;t.vertex[1]=add(c.gx-12+c.x[e],c.gy-12+c.y[e]);t.vertex[2]=add(c.gx-12+c.x[j],c.gy-12+c.y[j]);
            for(int a=0;a<3;a++)for(int k=0;k<2;k++)t.t.x[a][k]=m.node[t.vertex[a]][k];need(gfs_material_triangle_init(&t.t),"positive shared reference triangle");m.tri.push_back(t);
        }
    }
    std::vector<int> parent(m.node.size());std::iota(parent.begin(),parent.end(),0);
    auto root=[&](int x){while(parent[x]!=x)x=parent[x];return x;};
    for(size_t i=0;i<m.tri.size();i++)for(int a=0;a<3;a++){
        int b=(a+1)%3,x=m.tri[i].vertex[a],y=m.tri[i].vertex[b];parent[root(x)]=root(y);if(x>y)std::swap(x,y);m.edge[{x,y}].push_back(i);
    }
    std::map<int,int> component;
    for(auto &t:m.tri){int r=root(t.vertex[0]);if(!component.count(r))component[r]=component.size();t.component=component[r];}
    m.components=component.size();for(auto &e:m.edge)need(e.second.size()<=2,"manifold shared edge");return m;
}
struct Point {int triangle;double x,y,w,rho,N[5],Nx[5],Ny[5];};
static void all_sum(std::vector<double> &x){MPI_Allreduce(MPI_IN_PLACE,x.data(),x.size(),MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);}
static std::vector<double> dense_solve(const std::vector<double> &M,std::vector<double> r){
    int n=r.size();std::vector<long double> L(n*n),b(r.begin(),r.end());
    for(int i=0;i<n;i++)for(int j=0;j<=i;j++){
        long double a=M[i*n+j];for(int k=0;k<j;k++)a-=L[i*n+k]*L[j*n+k];need(i!=j || a>0,"global shared material mass SPD");L[i*n+j]=i==j?sqrtl(a):a/L[j*n+j];
    }
    for(int i=0;i<n;i++){for(int j=0;j<i;j++)b[i]-=L[i*n+j]*b[j];b[i]/=L[i*n+i];}
    for(int i=n-1;i>=0;i--){for(int j=i+1;j<n;j++)b[i]-=L[j*n+i]*b[j];b[i]/=L[i*n+i];r[i]=(double)b[i];}return r;
}
static std::array<int,5> indices(const SharedMesh &mesh,const Tri &t,bool enriched){
    return {t.vertex[0],t.vertex[1],t.vertex[2],enriched?(int)mesh.node.size()+2*t.component:-1,enriched?(int)mesh.node.size()+2*t.component+1:-1};
}
static GfsMaterialMap local_map(const SharedMesh &mesh,const Tri &t,const std::vector<double> &x,const std::vector<double> &y,
    const std::vector<double> &vx,const std::vector<double> &vy,bool enriched){
    auto ids=indices(mesh,t,enriched);GfsMaterialMap m={};m.n=enriched?5:3;
    for(int a=0;a<m.n;a++){
        m.x[0][a]=x[ids[a]]-(a<3?t.t.x[0][0]:0);m.x[1][a]=y[ids[a]]-(a<3?t.t.x[0][1]:0);
        m.v[0][a]=vx[ids[a]];m.v[1][a]=vy[ids[a]];
    }return m;
}
struct CellAssembly {std::vector<double> M,fx,fy,explicit_fx,explicit_fy;double U=0,Udot=0,work=0,torque=0,external_work=0,external_torque=0,gpower=0,gtorque=0,Kdirect=0,divmax=0;};
struct SharedCaseView {
    const SharedMesh &mesh;const GfsEnrichConfig &cfg;const std::vector<Pair> &scale;const std::vector<CellAssembly> &cell;
    const std::vector<double> &x,&y,&vx,&vy,&old_ax,&old_ay;const char *name,*dest;int rank,ranks,source_mode;
};
using SharedHook=void(*)(const SharedCaseView &);static SharedHook shared_hook=nullptr;
static void run_case(SharedMesh mesh,GfsEnrichConfig cfg,const char *name,double bx,double by,int rank,int ranks,const char *dest,int source_mode=1){
    bool enriched=cfg.field.ic.kind==2;int nodes=mesh.node.size(),n=nodes+(enriched?2*mesh.components:0),nc=mesh.cell.size();
    std::vector<Point> point;std::vector<double> normal(mesh.components*3);
    for(size_t k=0;k<mesh.tri.size();k++){
        auto &t=mesh.tri[k];for(int a=0;a<3;a++){double B[4];velocity(cfg,t.t.x[a][0]+12,t.t.x[a][1]+12,t.t.v0[a],B);}
        if(t.cell%ranks!=rank)continue;
        double px[3],py[3];for(int a=0;a<3;a++){px[a]=t.t.x[a][0];py[a]=t.t.x[a][1];}
        GfsMomentPolygon poly={px,py,3,0};GfsEnrichBasis qb={};qb.generator[0]=qb.generator[1]=12;GfsEnrichQuadrature qr;
        need(gfs_enrich_quadrature_init(poly,&cfg,&qb,&qr),"shared-triangle positive IC quadrature");
        for(int j=0;j<gfs_enrich_nquad(poly,&cfg,&qr);j++){
            Point p={};p.triangle=k;need(gfs_enrich_quad(poly,&cfg,&qr,j,&p.x,&p.y,&p.w),"shared volume quadrature");if(p.w==0)continue;
            double u[4],v[2],B[4],one[2]={1,1};need(gfs_enrich_reference(&cfg,p.x+12,p.y+12,u,0,0),"same exact input density/energy");p.rho=u[0];velocity(cfg,p.x+12,p.y+12,v,B);
            gfs_material_triangle_shape(&t.t,p.x,p.y,v,B,one,enriched,p.N,p.Nx,p.Ny);
            normal[3*t.component]+=p.w*p.rho;if(enriched)for(int a=0;a<2;a++)normal[3*t.component+1+a]+=p.w*p.rho*p.N[3+a]*p.N[3+a];point.push_back(p);
        }
    }
    all_sum(normal);std::vector<Pair> scale(mesh.components,{1,1});
    if(enriched)for(int c=0;c<mesh.components;c++)for(int k=0;k<2;k++){need(normal[3*c+1+k]>0,"independent global input-velocity enrichment");scale[c][k]=sqrt(normal[3*c+1+k]/normal[3*c]);}
    std::vector<double> x(n),y(n),vx(n),vy(n);
    for(int i=0;i<nodes;i++){x[i]=mesh.node[i][0];y[i]=mesh.node[i][1];double v[2],B[4];velocity(cfg,x[i]+12,y[i]+12,v,B);vx[i]=v[0]+bx;vy[i]=v[1]+by;}
    if(enriched)for(int c=0;c<mesh.components;c++){vx[nodes+2*c]=scale[c][0];vy[nodes+2*c+1]=scale[c][1];}
    std::vector<CellAssembly> ca(nc);for(auto &a:ca){a.M.resize(n*n);a.fx.resize(n);a.fy.resize(n);a.explicit_fx.resize(n);a.explicit_fy.resize(n);}
    double velocity_error=0,point_count=point.size();
    for(auto &p:point){
        auto &t=mesh.tri[p.triangle];auto ids=indices(mesh,t,enriched);auto m=local_map(mesh,t,x,y,vx,vy,enriched);auto &a=ca[t.cell];
        if(enriched)for(int k=0;k<2;k++){p.N[3+k]/=scale[t.component][k];p.Nx[3+k]/=scale[t.component][k];p.Ny[3+k]/=scale[t.component][k];}
        GfsMaterialPoint mp;need(gfs_material_point(&m,p.N,p.Nx,p.Ny,&mp),"positive shared material map");double gx=0,gy=0;
        if(enriched && source_mode==1){double r2=p.x*p.x+p.y*p.y+.047*.047,f=pow(r2,-1.5);gx=-p.x*f;gy=-p.y*f;}
        if(enriched && source_mode==2){gx=-p.x;gy=-p.y;}
        double force[2][5]={{0}},du=0,gp=0,sigma[3]={-cfg.field.ic.pressure,0,-cfg.field.ic.pressure};
        need(gfs_material_stress_point(&m,p.N,p.Nx,p.Ny,p.rho,p.w,sigma,gx,gy,force,&du,&gp),"shared compatible stress work");
        for(int i=0;i<m.n;i++){
            a.fx[ids[i]]+=force[0][i];a.fy[ids[i]]+=force[1][i];
            a.explicit_fx[ids[i]]+=force[0][i];a.explicit_fy[ids[i]]+=force[1][i];
            for(int j=i;j<m.n;j++){double val=p.w*p.rho*(p.N[i]*p.N[j]);a.M[ids[i]*n+ids[j]]+=val;if(i!=j)a.M[ids[j]*n+ids[i]]+=val;}
        }
        double v[2],B[4];velocity(cfg,p.x+12,p.y+12,v,B);velocity_error=fmax(velocity_error,fmax(fabs(mp.velocity[0]-v[0]-bx),fabs(mp.velocity[1]-v[1]-by)));
        a.U+=p.w*cfg.field.ic.pressure/(cfg.field.ic.gamma-1);a.Udot+=du;a.gpower+=gp;a.gtorque+=p.w*p.rho*(p.x*gy-p.y*gx);
        a.Kdirect+=.5*p.w*p.rho*(pow(v[0]+bx,2)+pow(v[1]+by,2));a.divmax=fmax(a.divmax,fabs(mp.grad_v[0]+mp.grad_v[3]));
    }
    /* Complete fan traction is evaluated on each physical cell boundary,
     * including shared edges. The shared assembly cancels its two actions. */
    for(size_t k=0;k<mesh.tri.size();k++){
        auto &t=mesh.tri[k];if(t.cell%ranks!=rank)continue;auto &a=ca[t.cell];auto ids=indices(mesh,t,enriched);
        int ia=t.vertex[1],ib=t.vertex[2];double ex=mesh.node[ib][0]-mesh.node[ia][0],ey=mesh.node[ib][1]-mesh.node[ia][1],len=hypot(ex,ey),nx=ey/len,ny=-ex/len;
        auto key=std::minmax(ia,ib);auto adjacent=mesh.edge.at({key.first,key.second});bool external=adjacent.size()==1 || mesh.tri[adjacent[0]].cell==mesh.tri[adjacent.back()].cell;
        for(int q=0;q<cfg.na;q++){
            double s=cfg.an[q],w=len*cfg.aw[q],X=(1-s)*mesh.node[ia][0]+s*mesh.node[ib][0],Y=(1-s)*mesh.node[ia][1]+s*mesh.node[ib][1],v0[2],B[4],N[5],Nx[5],Ny[5],u[4];
            velocity(cfg,X+12,Y+12,v0,B);gfs_material_triangle_shape(&t.t,X,Y,v0,B,scale[t.component].data(),enriched,N,Nx,Ny);
            auto m=local_map(mesh,t,x,y,vx,vy,enriched);GfsMaterialPoint mp;need(gfs_material_point(&m,N,Nx,Ny,&mp),"shared boundary map");need(gfs_enrich_reference(&cfg,X+12,Y+12,u,0,0),"boundary density");
            GfsAlePrimitive primitive={u[0],cfg.field.ic.pressure,mp.velocity[0],mp.velocity[1]};auto fan=gfs_ale_fan_build(primitive,primitive,nx,ny,cfg.field.ic.gamma);
            auto f=gfs_ale_fan_sample(&fan,primitive.vx*nx+primitive.vy*ny);need(f.valid,"material full Euler matching trace");
            double tx=-(f.flux[1]-primitive.vx*f.flux[0]),ty=-(f.flux[2]-primitive.vy*f.flux[0]);
            for(int j=0;j<m.n;j++){a.fx[ids[j]]+=w*N[j]*tx;a.fy[ids[j]]+=w*N[j]*ty;}
            if(external)for(int j=0;j<m.n;j++){a.explicit_fx[ids[j]]+=w*N[j]*tx;a.explicit_fy[ids[j]]+=w*N[j]*ty;}
            a.work+=w*(mp.velocity[0]*tx+mp.velocity[1]*ty);a.torque+=w*(X*ty-Y*tx);
            if(external){a.external_work+=w*(mp.velocity[0]*tx+mp.velocity[1]*ty);a.external_torque+=w*(X*ty-Y*tx);}
        }
    }
    for(auto &a:ca){all_sum(a.M);all_sum(a.fx);all_sum(a.fy);all_sum(a.explicit_fx);all_sum(a.explicit_fy);std::vector<double> data={a.U,a.Udot,a.work,a.torque,a.gpower,a.gtorque,a.Kdirect,a.external_work,a.external_torque};all_sum(data);
        a.U=data[0];a.Udot=data[1];a.work=data[2];a.torque=data[3];a.gpower=data[4];a.gtorque=data[5];a.Kdirect=data[6];a.external_work=data[7];a.external_torque=data[8];MPI_Allreduce(MPI_IN_PLACE,&a.divmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);}
    MPI_Allreduce(MPI_IN_PLACE,&velocity_error,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);MPI_Allreduce(MPI_IN_PLACE,&point_count,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    std::vector<double> M(n*n),fx(n),fy(n);for(auto &a:ca){for(int i=0;i<n*n;i++)M[i]+=a.M[i];for(int i=0;i<n;i++){fx[i]+=a.fx[i];fy[i]+=a.fy[i];}}
    auto ax=dense_solve(M,fx),ay=dense_solve(M,fy);double solve_error=0,force_scale=0;
    if(shared_hook){SharedCaseView view={mesh,cfg,scale,ca,x,y,vx,vy,ax,ay,name,dest,rank,ranks,source_mode};shared_hook(view);return;}
    for(int i=0;i<n;i++){
        double a=0,b=0;for(int j=0;j<n;j++){a+=M[i*n+j]*ax[j];b+=M[i*n+j]*ay[j];}solve_error=fmax(solve_error,fmax(fabs(a-fx[i]),fabs(b-fy[i])));force_scale=fmax(force_scale,fmax(fabs(fx[i]),fabs(fy[i])));
    }
    double trace_v=0,trace_a=0,trace_x=0;
    for(auto &edge:mesh.edge)if(edge.second.size()==2)for(double s:{.13,.5,.87}){
        Pair X={(1-s)*mesh.node[edge.first.first][0]+s*mesh.node[edge.first.second][0],(1-s)*mesh.node[edge.first.first][1]+s*mesh.node[edge.first.second][1]};GfsMaterialPoint out[2],acc[2];
        for(int side=0;side<2;side++){
            auto &t=mesh.tri[edge.second[side]];double v[2],B[4],N[5],Nx[5],Ny[5];velocity(cfg,X[0]+12,X[1]+12,v,B);gfs_material_triangle_shape(&t.t,X[0],X[1],v,B,scale[t.component].data(),enriched,N,Nx,Ny);
            auto m=local_map(mesh,t,x,y,vx,vy,enriched);need(gfs_material_point(&m,N,Nx,Ny,&out[side]),"matching shared velocity");
            auto am=local_map(mesh,t,x,y,ax,ay,enriched);need(gfs_material_point(&am,N,Nx,Ny,&acc[side]),"matching shared acceleration");
            for(int k=0;k<2;k++)out[side].position[k]+=t.t.x[0][k];
        }
        for(int k=0;k<2;k++){trace_v=fmax(trace_v,fabs(out[0].velocity[k]-out[1].velocity[k]));trace_a=fmax(trace_a,fabs(acc[0].velocity[k]-acc[1].velocity[k]));trace_x=fmax(trace_x,fabs(out[0].position[k]-out[1].position[k]));}
    }
    if(rank==0){
        double global_energy=0,global_angular=0,totalK=0,totalU=0,max_e_over_U=0,max_L_over_K=0;
        printf("{\"case\":\"%s\",\"ranks\":%d,\"nodes\":%d,\"components\":%d,\"scalar_DOFs\":%d,\"points\":%.0f,\"velocity_error\":%.17g,\"mass_solve_relative\":%.17g,\"shared_position_jump\":%.17g,\"shared_velocity_jump\":%.17g,\"shared_acceleration_jump\":%.17g,\"cells\":[",name,ranks,nodes,mesh.components,n,point_count,velocity_error,solve_error/(force_scale+1.e-300),trace_x,trace_v,trace_a);
        for(int c=0;c<nc;c++){
            auto &a=ca[c];double K=0,kdot=0,Ldot=0;for(int i=0;i<n;i++)for(int j=0;j<n;j++){
                K+=.5*a.M[i*n+j]*(vx[i]*vx[j]+vy[i]*vy[j]);kdot+=a.M[i*n+j]*(vx[i]*ax[j]+vy[i]*ay[j]);
                Ldot+=a.M[i*n+j]*(x[i]*ay[j]-y[i]*ax[j]+vx[i]*vy[j]-vy[i]*vx[j]);
            }
            double er=kdot+a.Udot-a.work-a.gpower,lr=Ldot-a.torque-a.gtorque,defect_work=0,defect_torque=0;
            for(int i=0;i<n;i++){
                double rx=-a.fx[i],ry=-a.fy[i];for(int j=0;j<n;j++){rx+=a.M[i*n+j]*ax[j];ry+=a.M[i*n+j]*ay[j];}
                defect_work+=vx[i]*rx+vy[i]*ry;defect_torque+=x[i]*ry-y[i]*rx;
            }
            global_energy+=er;global_angular+=lr;totalK+=K;totalU+=a.U;max_e_over_U=fmax(max_e_over_U,fabs(er)/a.U);max_L_over_K=fmax(max_L_over_K,fabs(lr)/(2*K));
            if(c)printf(",");printf("{\"id\":%llu,\"K\":%.17g,\"U\":%.17g,\"K_initial_error_over_U\":%.17g,\"Euler_cell_energy_residual_over_U\":%.17g,\"Euler_cell_angular_residual_over_2K\":%.17g,\"Pdot_over_P0_max\":%.17g,\"mass_load_defect_work_over_U\":%.17g,\"mass_load_defect_torque_over_2K\":%.17g,\"energy_defect_identity_error_over_U\":%.17g,\"angular_defect_identity_error_over_2K\":%.17g}",
                (unsigned long long)mesh.cell[c].id,K,a.U,(K-a.Kdirect)/a.U,er/a.U,lr/(2*K),cfg.field.ic.gamma*a.divmax,defect_work/a.U,defect_torque/(2*K),(er-defect_work)/a.U,(lr-defect_torque)/(2*K));
        }
        printf("],\"global_energy_residual_over_U\":%.17g,\"global_angular_residual_over_2K\":%.17g,\"max_cell_energy_residual_over_U\":%.17g,\"max_cell_angular_residual_over_2K\":%.17g}",global_energy/totalU,global_angular/(2*totalK),max_e_over_U,max_L_over_K);
        std::string path=std::string(dest)+"/"+name+"_state.bin";FILE *f=fopen(path.c_str(),"wb");need(f,"new shared-state audit");uint64_t header[2]={0x47465348414c3031ULL,(uint64_t)n};fwrite(header,sizeof(header),1,f);
        for(auto *a:{&x,&y,&vx,&vy,&ax,&ay})need(fwrite(a->data(),8,n,f)==(size_t)n,"shared state export");fclose(f);fflush(stdout);
    }
}
#ifndef GFS_MATERIAL_SHARED_LIBRARY
int main(int argc,char **argv){
    MPI_Init(&argc,&argv);int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);need(argc==3 || argc==4,"native polygon directory, new output directory, optional source_control");
    if(argc==4)need(std::string(argv[3])=="source_control","known source-control option");
    std::set<uint64_t> ids={17894,20502,39878,44945,67927,67934,67935,67942,67943,67992,67993,68000,68001,68002,68525,68526,68527,68534,68535,68584,68585,68592,68593,68600};auto mesh=mesh_from(load(argv[1],ids));
    if(rank==0)printf("{\"scope\":\"globally shared continuous material kinematics and consistent mass solve on actual24-cell patch; initial smooth RHS only; local physical Euler traction/work ledgers are measured separately from global conservation\",\"cases\":[");
    int cases=argc==4?3:3;
    for(int test=0;test<cases;test++){
        auto cfg=config(32,20);const char *name=test==0?"shear":test==1?"Kepler":"Kepler_perturb";int source=1;
        if(argc==4){name=test==0?"Kepler_profile_gravity_off":test==1?"Kepler_profile_linear_source":"Kepler_mass_defect_check";source=test==0?0:test==1?2:1;}
        else if(test==0){cfg.field.ic.kind=1;cfg.field.ic.rho=1;cfg.field.ic.vx=.2;cfg.field.ic.vy=.1;cfg.field.ic.B[1]=1.2;}
        if(rank==0 && test)printf(",\n");run_case(mesh,cfg,name,argc==3 && test==2?.013:0,argc==3 && test==2?-.007:0,rank,ranks,argv[2],source);
    }
    if(rank==0)printf("]}\n");MPI_Finalize();return 0;
}
#endif
