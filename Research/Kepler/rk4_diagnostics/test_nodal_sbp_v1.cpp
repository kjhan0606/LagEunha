#include <cfloat>
#include <array>
#include <algorithm>
#define GFS_ENRICH_TEST_LIBRARY
#include "test_enriched_euler_v1.cpp"
#include "gfs_sbp_flux.h"
#include "gfs_sbp_polygon.h"
#include "gfs_sbp_upwind.h"

struct Node {double x,y,m,ex,ey;};
struct Operator {
    std::vector<Node> p;std::vector<double> qx,qy;GfsSbpMatrixAudit audit;double area;
};
using State=std::array<double,4>;
static bool advected_wave=false;
static bool central_control=false;
static bool patch_only=false;
static bool upwind_candidate=false;
static Operator make_operator(const Cell &c){
    GfsMomentPolygon poly={c.x.data(),c.y.data(),(int)c.x.size(),0};GfsMomentDensity geom={};
    need(gfs_moment_geometry(poly,&geom),"actual polygon geometry");Operator o;o.area=geom.area;
    const double an[3]={.5*(1-sqrt(3./5)),.5,.5*(1+sqrt(3./5))},aw[3]={5./18,4./9,5./18};
    const double rn[3]={.5*(1-1/sqrt(5.)),.5*(1+1/sqrt(5.)),1},rw[3]={(5-sqrt(5.))/24,(5+sqrt(5.))/24,1./12};
    for(int e=0;e<poly.n;e++){
        int j=(e+1)%poly.n;double ex=poly.x[j]-poly.x[e],ey=poly.y[j]-poly.y[e];
        if(ex==0 && ey==0)continue;
        double area2=(poly.x[e]-geom.gx)*(poly.y[j]-geom.gy)-(poly.y[e]-geom.gy)*(poly.x[j]-geom.gx);
        need(area2>0,"positive centroid fan triangle");
        for(int ir=0;ir<3;ir++)for(int ia=0;ia<3;ia++){
            double edge_x=(1-an[ia])*poly.x[e]+an[ia]*poly.x[j],edge_y=(1-an[ia])*poly.y[e]+an[ia]*poly.y[j];
            Node p={geom.gx+rn[ir]*(edge_x-geom.gx),geom.gy+rn[ir]*(edge_y-geom.gy),area2*rw[ir]*aw[ia],0,0};
            if(ir==2){p.ex=aw[ia]*ey;p.ey=-aw[ia]*ex;}o.p.push_back(p);
        }
    }
    int n=o.p.size();std::vector<double> m(n),ex(n),ey(n),v(6*n),vx(6*n),vy(6*n),scratch(12*n);
    GfsDgBasis b={};b.gx=geom.gx;b.gy=geom.gy;b.l00=geom.l00;b.l10=geom.l10;b.l11=geom.l11;
    for(int i=0;i<n;i++){
        m[i]=o.p[i].m;ex[i]=o.p[i].ex;ey[i]=o.p[i].ey;
        gfs_dg_monomial(&b,2,o.p[i].x,o.p[i].y,v.data()+6*i,vx.data()+6*i,vy.data()+6*i);
    }
    o.qx.resize(n*n);o.qy.resize(n*n);
    need(gfs_sbp_polygon_p2(n,m.data(),ex.data(),ey.data(),v.data(),vx.data(),vy.data(),o.qx.data(),o.qy.data(),scratch.data(),&o.audit),"native positive-diagonal P2 SBP construction");
    return o;
}
static GfsSbpPrimitive primitive(const State &u,double gamma){
    return {u[0],gfs_dg_pressure(u.data(),gamma),u[1]/u[0],u[2]/u[0]};
}
static State conservative(GfsSbpPrimitive a,double gamma){
    return {a.rho,a.rho*a.vx,a.rho*a.vy,a.p/(gamma-1)+.5*a.rho*(a.vx*a.vx+a.vy*a.vy)};
}
static int diagnostic_volume_flux(GfsSbpPrimitive l,GfsSbpPrimitive r,double nx,double ny,double wn,double gamma,double f[4]){
    if(!central_control)return gfs_sbp_volume_flux(l,r,nx,ny,wn,gamma,f);
    auto ul=conservative(l,gamma),ur=conservative(r,gamma);double a=l.vx*nx+l.vy*ny,b=r.vx*nx+r.vy*ny;
    double fl[4]={l.rho*(a-wn),l.rho*l.vx*(a-wn)+l.p*nx,l.rho*l.vy*(a-wn)+l.p*ny,ul[3]*(a-wn)+l.p*a};
    double fr[4]={r.rho*(b-wn),r.rho*r.vx*(b-wn)+r.p*nx,r.rho*r.vy*(b-wn)+r.p*ny,ur[3]*(b-wn)+r.p*b};
    for(int c=0;c<4;c++)f[c]=.5*(fl[c]+fr[c]);return 1;
}
static std::vector<std::array<double,8>> primitive_gradients(const Operator &o,const std::vector<State> &u,double gamma){
    int n=o.p.size();std::vector<std::array<double,8>> gradient(n);
    for(int i=0;i<n;i++)for(int j=0;j<n;j++){
        auto p=primitive(u[j],gamma);double q[4]={p.rho,p.p,p.vx,p.vy};
        for(int c=0;c<4;c++){gradient[i][c]+=o.qx[i*n+j]*q[c]/o.p[i].m;gradient[i][4+c]+=o.qy[i*n+j]*q[c]/o.p[i].m;}
    }return gradient;
}
static GfsSbpPrimitive midpoint_trace(GfsSbpPrimitive q,const std::array<double,8> &gradient,double dx,double dy){
    q.rho+=gradient[0]*dx+gradient[4]*dy;q.p+=gradient[1]*dx+gradient[5]*dy;
    q.vx+=gradient[2]*dx+gradient[6]*dy;q.vy+=gradient[3]*dx+gradient[7]*dy;return q;
}
static void volume_rhs(const Operator &o,const std::vector<State> &u,const std::vector<std::array<double,2>> &w,double gamma,std::vector<State> &r,std::vector<State> *correction=nullptr){
    int n=o.p.size();r.assign(n,State{});if(correction)correction->assign(n,State{});
    std::vector<std::array<double,8>> gradient;if(upwind_candidate)gradient=primitive_gradients(o,u,gamma);
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
        double sx=o.qx[i*n+j],sy=o.qy[i*n+j],wn=.5*((w[i][0]+w[j][0])*sx+(w[i][1]+w[j][1])*sy),f[4];
        need(diagnostic_volume_flux(primitive(u[i],gamma),primitive(u[j],gamma),sx,sy,wn,gamma,f),"positive nodal volume flux");
        if(upwind_candidate){
            double dx=.5*(o.p[j].x-o.p[i].x),dy=.5*(o.p[j].y-o.p[i].y),df[4];
            auto l=midpoint_trace(primitive(u[i],gamma),gradient[i],dx,dy),rr=midpoint_trace(primitive(u[j],gamma),gradient[j],-dx,-dy);
            if(!gfs_sbp_upwind_correction(l,rr,sx,sy,wn,gamma,df)){
                fprintf(stderr,"inadmissible reconstructed trace or HLL state: nodes %d,%d rho %.17g %.17g P %.17g %.17g\n",i,j,l.rho,rr.rho,l.p,rr.p);need(false,"volume upwind admissibility");
            }
            for(int c=0;c<4;c++){f[c]+=df[c];if(correction){(*correction)[i][c]-=2*df[c];(*correction)[j][c]+=2*df[c];}}
        }
        for(int c=0;c<4;c++){r[i][c]-=2*f[c];r[j][c]+=2*f[c];}
    }
}
static double pressure_dot(const State &u,const State &dot,double gamma){
    double vx=u[1]/u[0],vy=u[2]/u[0];return (gamma-1)*(dot[3]-vx*dot[1]-vy*dot[2]+.5*(vx*vx+vy*vy)*dot[0]);
}
static void self_boundary(const Operator &o,const std::vector<State> &u,const std::vector<std::array<double,2>> &w,double gamma,std::vector<State> &r){
    for(size_t i=0;i<o.p.size();i++){
        auto &p=o.p[i];double a=hypot(p.ex,p.ey);if(a==0)continue;auto v=primitive(u[i],gamma);
        GfsAlePrimitive q={v.rho,v.p,v.vx,v.vy};
        auto f=gfs_ale_flux(q,q,p.ex/a,p.ey/a,w[i][0],w[i][1],gamma);need(f.valid,"same-state full Euler HLLC boundary SAT");
        for(int c=0;c<4;c++)r[i][c]-=a*f.flux[c];
    }
}
static void algebra(){
    double sym=0,consistent=0,kep=0,boost_error=0,radau=0;uint64_t seed=173;
    auto random=[&](){seed=seed*6364136223846793005ULL+1442695040888963407ULL;return (seed>>11)*0x1.0p-53;};
    for(int k=0;k<12000;k++){
        GfsSbpPrimitive l={exp(8*random()-4),exp(10*random()-8),6*random()-3,6*random()-3};
        GfsSbpPrimitive r={exp(8*random()-4),exp(10*random()-8),6*random()-3,6*random()-3};
        double nx=2*random()-1,ny=2*random()-1,w=2*random()-1,gamma=k%2?1.4:5./3,f[4],g[4],ff[4];
        need(gfs_sbp_volume_flux(l,r,nx,ny,w,gamma,f),"volume flux sample");
        need(gfs_sbp_volume_flux(r,l,nx,ny,w,gamma,g),"symmetric volume sample");
        for(int c=0;c<4;c++)sym=fmax(sym,fabs(f[c]-g[c])/(1+fabs(f[c])));
        double dvx=r.vx-l.vx,dvy=r.vy-l.vy;
        double kd=dvx*f[1]+dvy*f[2]-.5*(r.vx*r.vx+r.vy*r.vy-l.vx*l.vx-l.vy*l.vy)*f[0];
        double expected=.5*(l.p+r.p)*(dvx*nx+dvy*ny);
        kep=fmax(kep,fabs(kd-expected)/(1+fabs(kd)+fabs(f[3])));
        need(gfs_sbp_volume_flux(l,l,nx,ny,w,gamma,ff),"Euler ALE consistency sample");
        auto u=conservative(l,gamma);double vn=l.vx*nx+l.vy*ny;
        double exact[4]={l.rho*(vn-w),l.rho*l.vx*(vn-w)+l.p*nx,l.rho*l.vy*(vn-w)+l.p*ny,(u[3]+l.p)*vn-w*u[3]};
        for(int c=0;c<4;c++)consistent=fmax(consistent,fabs(ff[c]-exact[c])/(1+fabs(exact[c])));
        double bx=.013,by=-.007;l.vx+=bx;l.vy+=by;r.vx+=bx;r.vy+=by;
        need(gfs_sbp_volume_flux(l,r,nx,ny,w+bx*nx+by*ny,gamma,g),"boosted volume flux");
        double transformed[4]={f[0],f[1]+bx*f[0],f[2]+by*f[0],f[3]+bx*f[1]+by*f[2]+.5*(bx*bx+by*by)*f[0]};
        for(int c=0;c<4;c++)boost_error=fmax(boost_error,fabs(g[c]-transformed[c])/(1+fabs(transformed[c])));
    }
    double rn[3]={.5*(1-1/sqrt(5.)),.5*(1+1/sqrt(5.)),1},rw[3]={(5-sqrt(5.))/24,(5+sqrt(5.))/24,1./12};
    for(int k=0;k<=4;k++){double s=0;for(int j=0;j<3;j++)s+=rw[j]*pow(rn[j],k);radau=fmax(radau,fabs(s-1./(k+2)));}
    printf("\"algebra\":{\"samples\":12000,\"symmetry\":%.17g,\"Euler_ALE_consistency\":%.17g,\"kinetic_identity\":%.17g,\"boost\":%.17g,\"positive_radial_rule_moment_error\":%.17g},\n",sym,consistent,kep,boost_error,radau);
    need(sym<1.e-13 && consistent<1.e-12 && kep<1.e-12 && boost_error<1.e-12 && radau<1.e-14,"two-point and cubature algebra gates");
}
static std::set<uint64_t> neighbors(const char *path){
    std::set<uint64_t> core={67943,68000,68527,68584,67935,68001,68526,68592},ids=core;
    for(auto id:{17894ULL,20502ULL,44945ULL,39878ULL})ids.insert(id);
    FILE *f=fopen((std::string(path)+"/enriched_pairs_r0.bin").c_str(),"rb");need(f,"canonical native face list");
    uint64_t h[4];need(fread(h,sizeof(h),1,f)==1 && h[0]==0x4746534543503031ULL && h[3]==144,"canonical face packet schema");
    std::vector<char> row(144);for(uint64_t k=0;k<h[2];k++){
        need(fread(row.data(),144,1,f)==1,"canonical face row");uint64_t a,b;memcpy(&a,row.data(),8);memcpy(&b,row.data()+8,8);
        if(core.count(a))ids.insert(b);if(core.count(b))ids.insert(a);
    }fclose(f);return ids;
}
static std::map<uint64_t,State> reference_integrals(const char *path,const std::set<uint64_t> &ids){
    FILE *f=fopen((std::string(path)+"/enriched_project_r0.bin").c_str(),"rb");need(f,"native integrated IC packet");
    uint64_t h[8];need(fread(h,sizeof(h),1,f)==1 && h[5]==1528,"native integrated IC schema");
    std::map<uint64_t,State> out;std::vector<char> row(1528);
    for(uint64_t k=0;k<h[2];k++){
        need(fread(row.data(),1528,1,f)==1,"native initial integral row");uint64_t id;memcpy(&id,row.data(),8);
        if(ids.count(id)){State a;memcpy(a.data(),row.data()+1360,32);out[id]=a;}
    }fclose(f);need(out.size()==ids.size(),"all reference integrals");return out;
}
static void patch(const Cell &cell,const State &reference){
    auto o=make_operator(cell);auto cfg=config(16,16);int n=o.p.size();double gamma=cfg.field.ic.gamma,P0=cfg.field.ic.pressure;
    printf("{\"id\":%llu,\"nodes\":%d,\"minimum_mass_fraction\":%.17g,\"SBP\":%.17g,\"row_sum\":%.17g,\"P2_exactness\":%.17g,\"cubature_compatibility\":%.17g,",
        (unsigned long long)cell.id,n,o.audit.minimum_mass/o.area,o.audit.skew_error,o.audit.rowsum_error,o.audit.polynomial_error,o.audit.compatibility);
    need(o.audit.minimum_mass>0 && o.audit.skew_error<1.e-13 && o.audit.rowsum_error<1.e-11 && o.audit.polynomial_error<1.e-11 && o.audit.compatibility<1.e-11,"actual polygon positive SBP gates");
    std::vector<State> u(n),r;std::vector<std::array<double,2>> w(n);
    State integral={};double initial_p=0;
    for(int i=0;i<n;i++){
        double q[4];need(gfs_enrich_reference(&cfg,cell.gx+o.p[i].x,cell.gy+o.p[i].y,q,0,0),"same original continuous Kepler IC");
        for(int c=0;c<4;c++){u[i][c]=q[c];integral[c]+=o.p[i].m*q[c];}
        initial_p=fmax(initial_p,fabs(gfs_dg_pressure(q,gamma)/P0-1));
    }
    printf("\"initial_P_error\":%.17g,\"initial_mass_quadrature_relative\":%.17g,\"initial_energy_quadrature_over_thermal\":%.17g,\"states\":[",initial_p,
        (integral[0]-reference[0])/reference[0],(integral[3]-reference[3])/(o.area*P0/(gamma-1)));
    for(int test=0;test<3;test++){
        auto states=u;double bx=test?.013:0,by=test?-.007:0;
        double trace=test==2?.04:0;
        for(int i=0;i<n;i++){
            transform(states[i].data(),bx,by);w[i]={0,0};
            if(test==2)w[i]={.02+.07*o.p[i].x-.03*o.p[i].y,-.01+.02*o.p[i].x-.03*o.p[i].y};
        }
        std::vector<State> correction;volume_rhs(o,states,w,gamma,r,&correction);self_boundary(o,states,w,gamma,r);
        double chain=0,pdotmax=0,divmax=0,gravity=0,mean=0,l1=0,gcl=0;
        double angular=0,angular_boundary=0,angular_scale=0,angular_correction=0,correction_pressure=0;
        double correction_total[4]={},correction_norm[4]={};
        for(int i=0;i<n;i++){
            auto &p=o.p[i];double X=cell.gx+p.x-12,Y=cell.gy+p.y-12;auto v=primitive(states[i],gamma);
            angular+=X*r[i][2]-Y*r[i][1]+p.m*(w[i][0]*states[i][2]-w[i][1]*states[i][1]);
            angular_correction+=X*correction[i][2]-Y*correction[i][1];
            correction_pressure=fmax(correction_pressure,fabs(pressure_dot(states[i],correction[i],gamma))/(P0*p.m));
            for(int c=0;c<4;c++){correction_total[c]+=correction[i][c];correction_norm[c]+=fabs(correction[i][c]);}
            angular_scale+=p.m*v.rho*(v.vx*v.vx+v.vy*v.vy);
            double area=hypot(p.ex,p.ey);if(area==0)continue;
            GfsAlePrimitive q={v.rho,v.p,v.vx,v.vy};auto f=gfs_ale_flux(q,q,p.ex/area,p.ey/area,w[i][0],w[i][1],gamma);
            need(f.valid,"actual angular face flux");angular_boundary+=area*(X*f.flux[2]-Y*f.flux[1]);
        }
        for(int i=0;i<n;i++){
            double div=0,divw=0;for(int j=0;j<n;j++){
                auto v=primitive(states[j],gamma);div+=o.qx[i*n+j]*v.vx+o.qy[i*n+j]*v.vy;
                divw+=o.qx[i*n+j]*w[j][0]+o.qy[i*n+j]*w[j][1];
            }div/=o.p[i].m;divw/=o.p[i].m;gcl=fmax(gcl,fabs(divw-trace));
            State dot={};for(int c=0;c<4;c++)dot[c]=r[i][c]/o.p[i].m-trace*states[i][c];
            double X=cell.gx+o.p[i].x-12,Y=cell.gy+o.p[i].y-12,a=pow(X*X+Y*Y+.047*.047,-1.5);
            double fx=-states[i][0]*X*a,fy=-states[i][0]*Y*a,vx=states[i][1]/states[i][0],vy=states[i][2]/states[i][0];
            State source={0,fx,fy,vx*fx+vy*fy};gravity=fmax(gravity,fabs(pressure_dot(states[i],source,gamma))/P0);
            for(int c=0;c<4;c++)dot[c]+=source[c];
            double pd=pressure_dot(states[i],dot,gamma)/P0,predicted=-gamma*div;
            chain=fmax(chain,fabs(pd-predicted));pdotmax=fmax(pdotmax,fabs(pd));divmax=fmax(divmax,fabs(div));
            mean+=o.p[i].m*pd;l1+=o.p[i].m*fabs(pd);
        }
        if(test)printf(",");
        double correction_conservation=0;for(int c=0;c<4;c++)if(correction_norm[c]>0)correction_conservation=fmax(correction_conservation,fabs(correction_total[c])/correction_norm[c]);
        printf("{\"perturb\":[%.17g,%.17g],\"affine_mesh\":%s,\"chain_rule_error_over_P0\":%.17g,\"max_Pdot_over_P0\":%.17g,\"mean_Pdot_over_P0\":%.17g,\"L1_Pdot_over_P0\":%.17g,\"max_divergence_error\":%.17g,\"gravity_thermal_residual_over_P0\":%.17g,\"affine_metric_GCL\":%.17g,\"angular_volume_plus_boundary_residual\":%.17g,\"angular_residual_over_kinetic_scale\":%.17g,\"upwind_angular_correction\":%.17g,\"upwind_angular_correction_over_kinetic_scale\":%.17g,\"upwind_max_abs_pressure_rate_over_P0\":%.17g,\"upwind_conservative_sum_relative\":%.17g}",
            bx,by,test==2?"true":"false",chain,pdotmax,mean/o.area,l1/o.area,divmax,gravity,gcl,angular+angular_boundary,(angular+angular_boundary)/angular_scale,angular_correction,angular_correction/angular_scale,correction_pressure,correction_conservation);
    }printf("]}");
}
static void upwind_affine_contract(const std::map<uint64_t,Cell> &cells){
    double trace_error=0,correction=0,antisymmetry=0,boost=0;int pairs=0;
    for(auto &entry:cells){
        auto o=make_operator(entry.second);int n=o.p.size();double gamma=5./3;std::vector<State> u(n);
        auto affine=[](double x,double y){return GfsSbpPrimitive{1+.1*x-.2*y,1.e-6*(1+.05*x-.04*y),.2+.3*x-.7*y,.1+.4*x+.02*y};};
        for(int i=0;i<n;i++)u[i]=conservative(affine(o.p[i].x,o.p[i].y),gamma);
        auto gradient=primitive_gradients(o,u,gamma);
        for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
            double dx=.5*(o.p[j].x-o.p[i].x),dy=.5*(o.p[j].y-o.p[i].y),sx=o.qx[i*n+j],sy=o.qy[i*n+j],norm=hypot(sx,sy);if(norm==0)continue;
            auto l=midpoint_trace(primitive(u[i],gamma),gradient[i],dx,dy),r=midpoint_trace(primitive(u[j],gamma),gradient[j],-dx,-dy);
            auto q=affine(o.p[i].x+dx,o.p[i].y+dy);double wn=.023*sx-.017*sy,f[4];
            double error=fmax(fabs(l.rho-q.rho),fmax(fabs(l.p-q.p)/1.e-6,fmax(fabs(l.vx-q.vx),fabs(l.vy-q.vy))));
            error=fmax(error,fmax(fabs(r.rho-q.rho),fmax(fabs(r.p-q.p)/1.e-6,fmax(fabs(r.vx-q.vx),fabs(r.vy-q.vy)))));
            trace_error=fmax(trace_error,error);need(gfs_sbp_upwind_correction(l,r,sx,sy,wn,gamma,f),"affine full Euler HLL correction");
            for(int c=0;c<4;c++)correction=fmax(correction,fabs(f[c])/norm);
            /* Nonmatching admissible traces, same wave-speed rule. */
            r.rho*=1.03;r.p*=1.04;r.vx+=.007;r.vy-=.004;double a[4],b[4],br[4];
            need(gfs_sbp_upwind_correction(l,r,sx,sy,wn,gamma,a),"nonmatching HLL trace");
            need(gfs_sbp_upwind_correction(r,l,-sx,-sy,-wn,gamma,b),"opposite HLL trace");
            double bx=.013,by=-.007;l.vx+=bx;l.vy+=by;r.vx+=bx;r.vy+=by;
            need(gfs_sbp_upwind_correction(l,r,sx,sy,wn+bx*sx+by*sy,gamma,br),"boosted reconstructed HLL correction");
            double expected[4]={a[0],a[1]+bx*a[0],a[2]+by*a[0],a[3]+bx*a[1]+by*a[2]+.5*(bx*bx+by*by)*a[0]};
            for(int c=0;c<4;c++){antisymmetry=fmax(antisymmetry,fabs(a[c]+b[c])/norm);boost=fmax(boost,fabs(br[c]-expected[c])/norm);}pairs++;
        }
    }
    printf("\"upwind_affine_contract\":{\"pairs\":%d,\"affine_primitive_trace_error\":%.17g,\"matching_trace_correction_per_metric\":%.17g,\"opposite_pair_flux_error_per_metric\":%.17g,\"boost_correction_error_per_metric\":%.17g},\n",pairs,trace_error,correction,antisymmetry,boost);
    need(trace_error<1.e-9 && correction<1.e-10 && antisymmetry<1.e-11 && boost<1.e-11,"upwind affine/pair/boost algebra");
}
struct Mesh {std::vector<Operator> cells;std::vector<int> offsets;std::vector<std::pair<int,int>> faces;int nodes=0;};
static Mesh density_mesh(){
    Mesh mesh;std::map<std::pair<long long,long long>,std::vector<int>> edges;
    for(int j=0;j<2;j++)for(int i=0;i<2;i++){
        Cell c={};c.x={i*.5,(i+1)*.5,(i+1)*.5,i*.5};c.y={j*.5,j*.5,(j+1)*.5,(j+1)*.5};auto o=make_operator(c);
        mesh.offsets.push_back(mesh.nodes);
        for(size_t k=0;k<o.p.size();k++){
            auto &p=o.p[k];if(p.ex==0 && p.ey==0)continue;
            long long x=llround(p.x*1.e12)%1000000000000LL,y=llround(p.y*1.e12)%1000000000000LL;
            edges[{x,y}].push_back(mesh.nodes+k);
        }mesh.nodes+=o.p.size();mesh.cells.push_back(o);
    }
    for(auto &e:edges){need(e.second.size()==2,"periodic canonical nodal face pair");mesh.faces.push_back({e.second[0],e.second[1]});}
    return mesh;
}
static void density_rhs(const Mesh &mesh,const std::vector<State> &u,std::vector<State> &r){
    const double gamma=1.4,wx=advected_wave?0:.8,wy=advected_wave?0:-.3;r.assign(mesh.nodes,State{});
    for(size_t c=0;c<mesh.cells.size();c++){
        auto &o=mesh.cells[c];int off=mesh.offsets[c],n=o.p.size();
        std::vector<State> local(u.begin()+off,u.begin()+off+n),local_r;std::vector<std::array<double,2>> w(n,{wx,wy});
        volume_rhs(o,local,w,gamma,local_r);for(int i=0;i<n;i++)r[off+i]=local_r[i];
    }
    std::vector<Node> p;for(auto &c:mesh.cells)p.insert(p.end(),c.p.begin(),c.p.end());
    for(auto pair:mesh.faces){
        int i=pair.first,j=pair.second;double a=hypot(p[i].ex,p[i].ey);
        need(hypot(p[i].ex+p[j].ex,p[i].ey+p[j].ey)<1.e-14,"canonical opposite face measures");
        auto l=primitive(u[i],gamma),rr=primitive(u[j],gamma);
        GfsAlePrimitive left={l.rho,l.p,l.vx,l.vy},right={rr.rho,rr.p,rr.vx,rr.vy};
        auto f=gfs_ale_flux(left,right,p[i].ex/a,p[i].ey/a,wx,wy,gamma);need(f.valid,"full HLLC/HLL boundary flux");
        for(int c=0;c<4;c++){r[i][c]-=a*f.flux[c];r[j][c]+=a*f.flux[c];}
    }
    for(int i=0;i<mesh.nodes;i++)for(int c=0;c<4;c++)r[i][c]/=p[i].m;
}
static GfsSbpPrimitive primitive_increment(const State &u,const State &du,double gamma){
    auto v=primitive(u,gamma);return {du[0],pressure_dot(u,du,gamma),(du[1]-v.vx*du[0])/v.rho,(du[2]-v.vy*du[0])/v.rho};
}
static State conservative_increment(GfsSbpPrimitive v,GfsSbpPrimitive d,double gamma){
    return {d.rho,v.vx*d.rho+v.rho*d.vx,v.vy*d.rho+v.rho*d.vy,
        d.p/(gamma-1)+.5*(v.vx*v.vx+v.vy*v.vy)*d.rho+v.rho*(v.vx*d.vx+v.vy*d.vy)};
}
static State euler_A_jvp(GfsSbpPrimitive v,GfsSbpPrimitive d,double nx,double ny,double wn,double gamma){
    double vn=v.vx*nx+v.vy*ny,a=vn-wn,dvn=d.vx*nx+d.vy*ny;
    auto du=conservative_increment(v,d,gamma);auto u=conservative(v,gamma);
    double fm=a*d.rho+v.rho*dvn;
    return {fm,v.vx*fm+v.rho*a*d.vx+nx*d.p,v.vy*fm+v.rho*a*d.vy+ny*d.p,
        a*du[3]+(u[3]+v.p)*dvn+vn*d.p};
}
static State hllc_matching_jvp(GfsSbpPrimitive v,GfsSbpPrimitive dl,GfsSbpPrimitive dr,
    double nx,double ny,double wn,double gamma){
    auto left=euler_A_jvp(v,dl,nx,ny,wn,gamma),right=euler_A_jvp(v,dr,nx,ny,wn,gamma);
    GfsSbpPrimitive d={dr.rho-dl.rho,dr.p-dl.p,dr.vx-dl.vx,dr.vy-dl.vy};
    double vn=v.vx*nx+v.vy*ny,a=vn-wn,cs=sqrt(gamma*v.p/v.rho),dvn=d.vx*nx+d.vy*ny;
    double ap=.5*(d.p+v.rho*cs*dvn),am=.5*(d.p-v.rho*cs*dvn);
    GfsSbpPrimitive entropy={d.rho-d.p/(cs*cs),0,0,0},tangent={0,0,d.vx-nx*dvn,d.vy-ny*dvn};
    GfsSbpPrimitive plus={ap/(cs*cs),ap,nx*ap/(v.rho*cs),ny*ap/(v.rho*cs)};
    GfsSbpPrimitive minus={am/(cs*cs),am,-nx*am/(v.rho*cs),-ny*am/(v.rho*cs)};
    auto se=conservative_increment(v,entropy,gamma),st=conservative_increment(v,tangent,gamma);
    auto sp=conservative_increment(v,plus,gamma),sm=conservative_increment(v,minus,gamma);State result;
    for(int c=0;c<4;c++)result[c]=.5*(left[c]+right[c])- .5*(fabs(a)*(se[c]+st[c])+fabs(a+cs)*sp[c]+fabs(a-cs)*sm[c]);
    return result;
}
static void general_density_jvp(const Mesh &mesh,const std::vector<State> &u,const std::vector<State> &du,std::vector<State> &r){
    double gamma=1.4,wx=advected_wave?0:.8,wy=advected_wave?0:-.3;r.assign(mesh.nodes,State{});std::vector<Node> points;
    for(size_t k=0;k<mesh.cells.size();k++){
        auto &o=mesh.cells[k];int off=mesh.offsets[k],n=o.p.size();points.insert(points.end(),o.p.begin(),o.p.end());
        for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
            int a=off+i,b=off+j;double sx=o.qx[i*n+j],sy=o.qy[i*n+j],f[4];
            if(central_control){
                auto l=euler_A_jvp(primitive(u[a],gamma),primitive_increment(u[a],du[a],gamma),sx,sy,wx*sx+wy*sy,gamma);
                auto rr=euler_A_jvp(primitive(u[b],gamma),primitive_increment(u[b],du[b],gamma),sx,sy,wx*sx+wy*sy,gamma);
                for(int c=0;c<4;c++)f[c]=.5*(l[c]+rr[c]);
            }else need(gfs_sbp_volume_jvp(primitive(u[a],gamma),primitive(u[b],gamma),primitive_increment(u[a],du[a],gamma),
                primitive_increment(u[b],du[b],gamma),sx,sy,wx*sx+wy*sy,gamma,f),"exact native volume JVP");
            for(int c=0;c<4;c++){r[a][c]-=2*f[c];r[b][c]+=2*f[c];}
        }
    }
    for(auto pair:mesh.faces){
        int a=pair.first,b=pair.second;double len=hypot(points[a].ex,points[a].ey),nx=points[a].ex/len,ny=points[a].ey/len;
        for(int c=0;c<4;c++)need(fabs(u[a][c]-u[b][c])<1.e-13,"identical physical traces for exact HLLC variation");
        auto f=hllc_matching_jvp(primitive(u[a],gamma),primitive_increment(u[a],du[a],gamma),primitive_increment(u[b],du[b],gamma),nx,ny,wx*nx+wy*ny,gamma);
        for(int c=0;c<4;c++){r[a][c]-=len*f[c];r[b][c]+=len*f[c];}
    }
    for(int i=0;i<mesh.nodes;i++)for(int c=0;c<4;c++)r[i][c]/=points[i].m;
}
/* Exact first variation of the SAME split-volume/HLLC operator about the
 * stationary translating density wave. Primitive order is rho,vx,vy,P.
 * At identical traces and relative normal speed zero, the full HLLC fan's
 * derivative is the exact linear acoustic Riemann solver. Entropy/tangential
 * characteristics have speed zero. This is verified against recorded full
 * conservative-state finite differences below; it is not a replacement flux.
 */
static void density_primitive_jvp(const Mesh &mesh,const std::vector<State> &base,
    const std::vector<State> &du,std::vector<State> &r){
    const double gamma=1.4,P=1;r.assign(mesh.nodes,State{});std::vector<Node> points;
    for(size_t k=0;k<mesh.cells.size();k++){
        auto &o=mesh.cells[k];int n=o.p.size(),off=mesh.offsets[k];points.insert(points.end(),o.p.begin(),o.p.end());
        for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
            int a=off+i,b=off+j;double nx=o.qx[i*n+j],ny=o.qy[i*n+j];
            double dvn=.5*((du[a][1]+du[b][1])*nx+(du[a][2]+du[b][2])*ny),dp=.5*(du[a][3]+du[b][3]);
            double fm=.5*(base[a][0]+base[b][0])*dvn,fp=gamma*P*dvn;
            r[a][0]-=2*fm;r[b][0]+=2*fm;
            r[a][1]-=2*nx*dp/base[a][0];r[b][1]+=2*nx*dp/base[b][0];
            r[a][2]-=2*ny*dp/base[a][0];r[b][2]+=2*ny*dp/base[b][0];
            r[a][3]-=2*fp;r[b][3]+=2*fp;
        }
    }
    for(auto pair:mesh.faces){
        int a=pair.first,b=pair.second;auto &p=points[a];double len=hypot(p.ex,p.ey),nx=p.ex/len,ny=p.ey/len;
        need(fabs(base[a][0]-base[b][0])<1.e-13,"matching continuous density-wave traces for exact HLLC derivative");
        double rho=.5*(base[a][0]+base[b][0]),cs=sqrt(gamma*P/rho);
        double va=du[a][1]*nx+du[a][2]*ny,vb=du[b][1]*nx+du[b][2]*ny;
        double fm=rho*.5*(va+vb)-.5*(du[b][3]-du[a][3])/cs;
        double ps=.5*(du[a][3]+du[b][3])-.5*rho*cs*(vb-va);
        r[a][0]-=len*fm;r[b][0]+=len*fm;
        r[a][1]-=p.ex*ps/base[a][0];r[b][1]+=p.ex*ps/base[b][0];
        r[a][2]-=p.ey*ps/base[a][0];r[b][2]+=p.ey*ps/base[b][0];
        r[a][3]-=len*gamma*P/rho*fm;r[b][3]+=len*gamma*P/rho*fm;
    }
    for(int i=0;i<mesh.nodes;i++)for(int c=0;c<4;c++)r[i][c]/=points[i].m;
}
static void matrix_file(const std::string &path,int n,const std::vector<double> &J){
    FILE *f=fopen(path.c_str(),"wb");need(f,"new Jacobian file");uint64_t h[2]={0x4746535342504a31ULL,(uint64_t)n};
    need(fwrite(h,sizeof(h),1,f)==1 && fwrite(J.data(),8,J.size(),f)==J.size(),"Jacobian packet write");fclose(f);
}
static void density_wave(const char *dest){
    auto mesh=density_mesh();std::vector<State> u(mesh.nodes),base;
    for(size_t k=0;k<mesh.cells.size();k++)for(size_t i=0;i<mesh.cells[k].p.size();i++){
        auto &p=mesh.cells[k].p[i];GfsSbpPrimitive v={1+.6*sin(2*acos(-1.)*p.x),1.,.8,-.3};u[mesh.offsets[k]+i]=conservative(v,1.4);
    }
    density_rhs(mesh,u,base);double maxbase=0,mass=0,cons[4]={};
    for(size_t k=0;k<mesh.cells.size();k++)for(size_t i=0;i<mesh.cells[k].p.size();i++)for(int c=0;c<4;c++){
        double rate=base[mesh.offsets[k]+i][c];maxbase=fmax(maxbase,fabs(rate));cons[c]+=mesh.cells[k].p[i].m*rate;
    }
    for(double z:cons)mass=fmax(mass,fabs(z));
    printf("\"density_wave\":{\"scope\":\"actual nodal P2 polygon SBP plus full HLLC boundaries; periodic 2x2 squares; fixed spatial metrics; no time evolution or weight solve\",\"mesh_velocity\":[%.17g,%.17g],\"nodes\":%d,\"variables\":%d,\"canonical_face_nodes\":%zu,\"base_RHS_max\":%.17g,\"conservative_rate_residual\":%.17g,\"jacobians\":[",advected_wave?0:.8,advected_wave?0:-.3,mesh.nodes,4*mesh.nodes,mesh.faces.size(),maxbase,mass);
    int dim=4*mesh.nodes;std::vector<double> exact(dim*dim),exact_conservative(dim*dim);
    for(int j=0;j<dim;j++){
        std::vector<State> du(mesh.nodes),dot;du[j/4][j%4]=1;density_primitive_jvp(mesh,u,du,dot);
        for(int i=0;i<dim;i++)exact[i*dim+j]=dot[i/4][i%4];
        du.assign(mesh.nodes,State{});int b=j/4,k=j%4;double rho=u[b][0],vx=.8,vy=-.3;
        du[b][0]=k==0?1:0;du[b][1]=((k==1?1:0)-vx*du[b][0])/rho;du[b][2]=((k==2?1:0)-vy*du[b][0])/rho;
        du[b][3]=.4*((k==3?1:0)-vx*(k==1?1:0)-vy*(k==2?1:0)+.5*(vx*vx+vy*vy)*du[b][0]);
        density_primitive_jvp(mesh,u,du,dot);
        for(int i=0;i<mesh.nodes;i++){
            double r=u[i][0];State dc={dot[i][0],vx*dot[i][0]+r*dot[i][1],vy*dot[i][0]+r*dot[i][2],
                dot[i][3]/.4+.5*(vx*vx+vy*vy)*dot[i][0]+r*(vx*dot[i][1]+vy*dot[i][2])};
            for(int c=0;c<4;c++)exact_conservative[(4*i+c)*dim+j]=dc[c];
        }
        /* Independent generic native flux JVP also checks the stationary
         * primitive derivation and supplies the advected state's Jacobian. */
        du.assign(mesh.nodes,State{});du[j/4][j%4]=1;general_density_jvp(mesh,u,du,dot);
        for(int i=0;i<dim;i++){
            if(!advected_wave)need(fabs(exact_conservative[i*dim+j]-dot[i/4][i%4])<1.e-10,"two independent stationary linearizations");
            exact_conservative[i*dim+j]=dot[i/4][i%4];
        }
    }
    if(!advected_wave)matrix_file(std::string(dest)+"/density_jacobian_exact_primitive.bin",dim,exact);
    matrix_file(std::string(dest)+"/density_jacobian_exact_conservative.bin",dim,exact_conservative);
    FILE *metrics=fopen((std::string(dest)+"/density_acoustic_metric.bin").c_str(),"wb");need(metrics,"new acoustic energy metric");
    for(size_t k=0;k<mesh.cells.size();k++)for(size_t i=0;i<mesh.cells[k].p.size();i++){
        double m=mesh.cells[k].p[i].m,rho=u[mesh.offsets[k]+i][0],h[3]={m*rho,m*rho,m/1.4};
        need(fwrite(h,8,3,metrics)==3,"acoustic metric write");
    }fclose(metrics);
    for(int t=0;t<3;t++){
        double step=pow(10.,-4-t);std::vector<double> J(dim*dim);std::vector<State> plus,minus;
        for(int j=0;j<dim;j++){
            auto pert=u;int node=j/4,component=j%4;double h=step*fmax(1.,fabs(u[node][component]));
            pert[node][component]+=h;density_rhs(mesh,pert,plus);pert[node][component]=u[node][component]-h;density_rhs(mesh,pert,minus);
            for(int i=0;i<dim;i++)J[i*dim+j]=(plus[i/4][i%4]-minus[i/4][i%4])/(2*h);
        }
        matrix_file(std::string(dest)+"/density_jacobian_"+std::to_string(t)+".bin",dim,J);
        double error=0;for(size_t i=0;i<J.size();i++)error=fmax(error,fabs(J[i]-exact_conservative[i]));
        if(t)printf(",");printf("{\"relative_FD_step\":%.17g,\"file\":\"density_jacobian_%d.bin\",\"maximum_error_vs_analytic\":%.17g}",step,t,error);
    }printf("]}\n");
}
int main(int argc,char **argv){
    need(argc==4 || argc==5,"original polygon packet, native coupled IC, new output directory, optional advected");
    if(argc==5){std::string choice(argv[4]);need(choice=="advected" || choice=="advectedcentral" || choice=="patch" || choice=="patchcentral" || choice=="patchupwind","known spatial audit variation");
        upwind_candidate=choice=="patchupwind";patch_only=choice=="patch" || choice=="patchcentral" || upwind_candidate;advected_wave=!patch_only;central_control=choice=="advectedcentral" || choice=="patchcentral";}
    printf("{\"scope\":\"native C/C++ actual-polygon nodal SBP spatial/algebra audit only; no evolved state, no arbitrary heat or viscosity\",\"central_Euler_diagnostic_control\":%s,\"reconstructed_HLL_volume_candidate\":%s,\n",central_control?"true":"false",upwind_candidate?"true":"false");
    algebra();auto ids=neighbors(argv[2]);auto cells=load(argv[1],ids);auto integrals=reference_integrals(argv[2],ids);
    if(upwind_candidate)upwind_affine_contract(cells);
    printf("\"actual_cells\":[");bool first=true;for(auto &c:cells){if(!first)printf(",\n");first=false;patch(c.second,integrals[c.first]);fflush(stdout);}printf("],\n");
    if(patch_only)printf("\"density_wave\":{\"scope\":\"not repeated; refer to frozen v3 split-flux and v4 central-control matrices\",\"jacobians\":[]}\n");
    else density_wave(argv[3]);printf("}\n");return 0;
}
