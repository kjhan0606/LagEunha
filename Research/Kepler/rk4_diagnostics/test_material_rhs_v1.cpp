#include <cfloat>
#include <array>
#include <algorithm>
#define GFS_ENRICH_TEST_LIBRARY
#include "test_enriched_euler_v1.cpp"
#include "gfs_ale_fan.h"
#include "gfs_material_map.h"
struct MaterialBasis {GfsMomentDensity g;double gx,gy,L[5][5]={{0}},velocity_scale;int n;};
struct MaterialSample {double x,y,w,rho,u[4],raw[5],rx[5],ry[5];};
static double power_remainder(double z){
    const double a=-.75;if(fabs(z)>.125)return expm1(a*log1p(z))-a*z;
    double term=.5*a*(a-1)*z*z,sum=term;
    for(int k=3;k<60;k++){term*=((a-k+1)/k)*z;double old=sum;sum+=term;if(sum==old)break;}return sum;
}
static void velocity_data(const GfsEnrichConfig &cfg,double X,double Y,double v[2],double B[4]){
    auto &ic=cfg.field.ic;X-=ic.cx;Y-=ic.cy;
    if(ic.kind==1){
        v[0]=ic.vx+ic.B[0]*X+ic.B[1]*Y;v[1]=ic.vy+ic.B[2]*X+ic.B[3]*Y;
        for(int i=0;i<4;i++)B[i]=ic.B[i];return;
    }
    double z=X*X+Y*Y+ic.eps*ic.eps,om=pow(z,-.75),ox=-1.5*X*om/z,oy=-1.5*Y*om/z;
    v[0]=-Y*om;v[1]=X*om;B[0]=-Y*ox;B[1]=-om-Y*oy;B[2]=om+X*ox;B[3]=X*oy;
}
static void raw_shape(const MaterialBasis &b,const GfsEnrichConfig &cfg,double x,double y,double f[5],double fx[5],double fy[5]){
    double u=(x-b.g.gx)/b.g.l00,v=(y-b.g.gy-b.g.l10*u)/b.g.l11;
    f[0]=1;f[1]=u;f[2]=v;fx[0]=fy[0]=0;fx[1]=1/b.g.l00;fy[1]=0;fx[2]=-b.g.l10/(b.g.l00*b.g.l11);fy[2]=1/b.g.l11;
    if(b.n==3)return;
    double Cx=b.gx+b.g.gx-cfg.field.ic.cx,Cy=b.gy+b.g.gy-cfg.field.ic.cy;
    double rx=x-b.g.gx,ry=y-b.g.gy,z0=Cx*Cx+Cy*Cy+cfg.field.ic.eps*cfg.field.ic.eps;
    double z=(2*(Cx*rx+Cy*ry)+rx*rx+ry*ry)/z0,om=pow(z0,-.75),delta=om*expm1(-.75*log1p(z));
    double nonlinear=om*(power_remainder(z)-.75*(rx*rx+ry*ry)/z0);
    f[3]=(-ry*delta-Cy*nonlinear)/b.velocity_scale;f[4]=(rx*delta+Cx*nonlinear)/b.velocity_scale;
    double vc[2],Bc[4],vp[2],Bp[4];velocity_data(cfg,b.gx+b.g.gx,b.gy+b.g.gy,vc,Bc);velocity_data(cfg,b.gx+x,b.gy+y,vp,Bp);
    fx[3]=(Bp[0]-Bc[0])/b.velocity_scale;fy[3]=(Bp[1]-Bc[1])/b.velocity_scale;
    fx[4]=(Bp[2]-Bc[2])/b.velocity_scale;fy[4]=(Bp[3]-Bc[3])/b.velocity_scale;
}
static void shape(const MaterialBasis &b,const GfsEnrichConfig &cfg,double x,double y,double N[5],double Nx[5],double Ny[5]){
    double f[5],fx[5],fy[5];raw_shape(b,cfg,x,y,f,fx,fy);
    for(int i=0;i<b.n;i++){
        N[i]=f[i];Nx[i]=fx[i];Ny[i]=fy[i];
        for(int j=0;j<i;j++){N[i]-=b.L[i][j]*N[j];Nx[i]-=b.L[i][j]*Nx[j];Ny[i]-=b.L[i][j]*Ny[j];}
        N[i]/=b.L[i][i];Nx[i]/=b.L[i][i];Ny[i]/=b.L[i][i];
    }
}
static void solve_small(long double A[5][5],long double rhs[5],int n){
    long double L[5][5]={{0}};
    for(int i=0;i<n;i++)for(int j=0;j<=i;j++){
        long double z=A[i][j];for(int k=0;k<j;k++)z-=L[i][k]*L[j][k];need(i!=j || z>0,"positive actual material mass/Gram pivot");L[i][j]=i==j?sqrtl(z):z/L[j][j];
    }
    for(int i=0;i<n;i++){for(int j=0;j<i;j++)rhs[i]-=L[i][j]*rhs[j];rhs[i]/=L[i][i];}
    for(int i=n-1;i>=0;i--){for(int j=i+1;j<n;j++)rhs[i]-=L[j][i]*rhs[j];rhs[i]/=L[i][i];}
}
static MaterialBasis build_basis(const Cell &cell,const GfsEnrichConfig &cfg,std::vector<MaterialSample> &points){
    MaterialBasis b;b.gx=cell.gx;b.gy=cell.gy;b.n=cfg.field.ic.kind==1?3:5;
    GfsMomentPolygon polygon={cell.x.data(),cell.y.data(),(int)cell.x.size(),0};need(gfs_moment_geometry(polygon,&b.g),"material reference geometry");
    b.velocity_scale=pow(pow(cell.gx+b.g.gx-12,2)+pow(cell.gy+b.g.gy-12,2)+.047*.047,-.75)*fmax(b.g.l00,b.g.l11);
    GfsEnrichBasis qb={};qb.generator[0]=cell.gx;qb.generator[1]=cell.gy;GfsEnrichQuadrature qr;
    need(gfs_enrich_quadrature_init(polygon,&cfg,&qb,&qr),"same-cusp positive reference quadrature");
    long double mass=0;
    for(int k=0;k<gfs_enrich_nquad(polygon,&cfg,&qr);k++){
        MaterialSample p={};need(gfs_enrich_quad(polygon,&cfg,&qr,k,&p.x,&p.y,&p.w),"positive material point");if(p.w==0)continue;
        need(gfs_enrich_reference(&cfg,cell.gx+p.x,cell.gy+p.y,p.u,0,0),"unchanged exact IC reference");p.rho=p.u[0];
        raw_shape(b,cfg,p.x,p.y,p.raw,p.rx,p.ry);mass+=(long double)p.w*p.rho;points.push_back(p);
    }
    b.L[0][0]=1;
    for(int k=1;k<b.n;k++){
        for(int pass=0;pass<2;pass++){
            long double G[5][5]={{0}},r[5]={};
            for(auto &p:points){
                double q[5];for(int i=0;i<k;i++){q[i]=p.raw[i];for(int j=0;j<i;j++)q[i]-=b.L[i][j]*q[j];q[i]/=b.L[i][i];}
                double residual=p.raw[k];for(int i=0;i<k;i++)residual-=b.L[k][i]*q[i];long double w=(long double)p.w*p.rho/mass;
                for(int i=0;i<k;i++){r[i]+=w*q[i]*residual;for(int j=0;j<k;j++)G[i][j]+=w*q[i]*q[j];}
            }
            solve_small(G,r,k);for(int i=0;i<k;i++)b.L[k][i]+=(double)r[i];
        }
        long double norm=0;
        for(auto &p:points){
            double q[5];for(int i=0;i<k;i++){q[i]=p.raw[i];for(int j=0;j<i;j++)q[i]-=b.L[i][j]*q[j];q[i]/=b.L[i][i];}
            double r=p.raw[k];for(int i=0;i<k;i++)r-=b.L[k][i]*q[i];norm+=(long double)p.w*p.rho*r*r/mass;
        }need(norm>0,"independent kinematic enrichment");b.L[k][k]=sqrt((double)norm);
    }return b;
}
static GfsMaterialMap initial_map(const MaterialBasis &b,const GfsEnrichConfig &cfg,double bx,double by){
    GfsMaterialMap m={};m.n=b.n;double c[2][5]={{b.g.gx,b.g.l00,0,0,0},{b.g.gy,b.g.l10,b.g.l11,0,0}};
    double vc[2],B[4];velocity_data(cfg,b.gx+b.g.gx,b.gy+b.g.gy,vc,B);
    double v[2][5]={{vc[0]+bx,B[0]*b.g.l00+B[1]*b.g.l10,B[1]*b.g.l11,b.n==5?b.velocity_scale:0,0},
        {vc[1]+by,B[2]*b.g.l00+B[3]*b.g.l10,B[3]*b.g.l11,0,b.n==5?b.velocity_scale:0}};
    for(int k=0;k<2;k++)for(int i=0;i<m.n;i++)for(int j=0;j<=i;j++){m.x[k][j]+=c[k][i]*b.L[i][j];m.v[k][j]+=v[k][i]*b.L[i][j];}return m;
}
static void assess_material(const Cell &cell,GfsEnrichConfig cfg,const char *kind,double bx,double by,const std::array<double,4> *native=nullptr){
    std::vector<MaterialSample> points;auto b=build_basis(cell,cfg,points);auto m=initial_map(b,cfg,bx,by);
    double M[5][5]={{0}},force[2][5]={{0}},Udot=0,Gpower=0,Gtorque=0,Kdirect=0,U=0,area=0,Vdot=0,Vfd=0,minJ=1;
    double field=0,gradient=0,pressure=0,coords=0,mass=0,power_abs=0;const double h=1.e-5;
    for(auto &p:points){
        double N[5],Nx[5],Ny[5];shape(b,cfg,p.x,p.y,N,Nx,Ny);GfsMaterialPoint mp;need(gfs_material_point(&m,N,Nx,Ny,&mp),"positive current material map");
        gfs_material_mass_point(m.n,N,p.w*p.rho,M);double gx=0,gy=0;
        if(cfg.field.ic.kind==2){double X=b.gx+p.x-12,Y=b.gy+p.y-12,a=pow(X*X+Y*Y+.047*.047,-1.5);gx=-X*a;gy=-Y*a;}
        double sigma[3]={-cfg.field.ic.pressure,0,-cfg.field.ic.pressure};
        need(gfs_material_stress_point(&m,N,Nx,Ny,p.rho,p.w,sigma,gx,gy,force,&Udot,&Gpower),"compatible material stress/work point");
        double X=b.gx+mp.position[0]-12,Y=b.gy+mp.position[1]-12;
        Gtorque+=p.w*p.rho*(X*gy-Y*gx);mass+=p.w*p.rho;
        Kdirect+=.5*p.w*p.rho*(pow(p.u[1]/p.rho+bx,2)+pow(p.u[2]/p.rho+by,2));U+=p.w*cfg.field.ic.pressure/(cfg.field.ic.gamma-1);area+=p.w;
        field=fmax(field,fmax(fabs(mp.velocity[0]-p.u[1]/p.rho-bx),fabs(mp.velocity[1]-p.u[2]/p.rho-by)));
        coords=fmax(coords,fmax(fabs(mp.position[0]-p.x),fabs(mp.position[1]-p.y)));
        double div=mp.grad_v[0]+mp.grad_v[3];Vdot+=p.w*mp.J*div;
        if(cfg.field.ic.kind==2)gradient=fmax(gradient,fabs(div));
        double rho=p.rho/mp.J,e=cfg.field.ic.pressure/((cfg.field.ic.gamma-1)*p.rho),edot=-cfg.field.ic.pressure*div/rho;
        double pdot=(cfg.field.ic.gamma-1)*(rho*edot-e*rho*div);pressure=fmax(pressure,fabs(pdot)/cfg.field.ic.pressure);
        GfsMaterialMap plus=m,minus=m;for(int k=0;k<2;k++)for(int a=0;a<m.n;a++){plus.x[k][a]+=h*m.v[k][a];minus.x[k][a]-=h*m.v[k][a];}
        GfsMaterialPoint pp,pm;need(gfs_material_point(&plus,N,Nx,Ny,&pp) && gfs_material_point(&minus,N,Nx,Ny,&pm),"finite mapped volume control");
        minJ=fmin(minJ,fmin(pp.J,pm.J));Vfd+=p.w*(pp.J-pm.J)/(2*h);
        power_abs+=p.w*p.rho*(fabs(mp.velocity[0]*gx)+fabs(mp.velocity[1]*gy))+p.w*fabs(cfg.field.ic.pressure*div);
        for(int a=0;a<m.n;a++)power_abs+=p.w*cfg.field.ic.pressure*(fabs(m.v[0][a]*Nx[a])+fabs(m.v[1][a]*Ny[a]));
    }
    double boundary_power=0,boundary_torque=0,boundary_volume=0,fan_mass=0,fan_work=0;
    faces(cell,cfg,[&](double x,double y,double nx,double ny,double w){
        double N[5],Nx[5],Ny[5],u[4];shape(b,cfg,x,y,N,Nx,Ny);GfsMaterialPoint p;need(gfs_material_point(&m,N,Nx,Ny,&p),"material boundary point");
        need(gfs_enrich_reference(&cfg,b.gx+x,b.gy+y,u,0,0),"material boundary initial density");
        GfsAlePrimitive q={u[0],cfg.field.ic.pressure,p.velocity[0],p.velocity[1]};double wn=q.vx*nx+q.vy*ny;
        auto fan=gfs_ale_fan_build(q,q,nx,ny,cfg.field.ic.gamma);auto f=gfs_ale_fan_sample(&fan,wn);need(f.valid,"complete Euler fan on matching material subface");
        double traction[2]={-(f.flux[1]-q.vx*f.flux[0]),-(f.flux[2]-q.vy*f.flux[0])};
        double work=p.velocity[0]*traction[0]+p.velocity[1]*traction[1];
        for(int a=0;a<m.n;a++){force[0][a]+=w*N[a]*traction[0];force[1][a]+=w*N[a]*traction[1];}
        boundary_power+=w*work;boundary_torque+=w*((b.gx+p.position[0]-12)*traction[1]-(b.gy+p.position[1]-12)*traction[0]);
        boundary_volume+=w*wn;fan_mass=fmax(fan_mass,fabs(f.flux[0]));
        double E=q.p/(cfg.field.ic.gamma-1)+.5*q.rho*(q.vx*q.vx+q.vy*q.vy);
        fan_work=fmax(fan_work,fabs(work+f.flux[3]-E/q.rho*f.flux[0]));
    });
    double Kdot=0,Ldot=0;for(int a=0;a<m.n;a++){
        Kdot+=m.v[0][a]*force[0][a]+m.v[1][a]*force[1][a];
        double x=m.x[0][a]+(a==0?b.gx-12:0),y=m.x[1][a]+(a==0?b.gy-12:0);
        Ldot+=x*force[1][a]-y*force[0][a];
        for(int j=0;j<m.n;j++)Ldot+=M[a][j]*(m.v[0][a]*m.v[1][j]-m.v[1][a]*m.v[0][j]);
    }
    double K=gfs_material_kinetic(&m,M),energy_error=Kdot+Udot-boundary_power-Gpower;
    double native_mass=0,native_energy=0;
    if(native){auto u=*native;transform(u.data(),bx,by);native_mass=(mass-u[0])/u[0];native_energy=(K+U-u[3])/U;}
    double acceleration[2][5]={{0}},mass_residual=0,force_norm=0,Kdot_acceleration=0;
    for(int k=0;k<2;k++){
        long double gram[5][5],r[5]={};for(int a=0;a<m.n;a++){r[a]=force[k][a];for(int j=0;j<m.n;j++)gram[a][j]=M[a][j];}solve_small(gram,r,m.n);
        for(int a=0;a<m.n;a++)acceleration[k][a]=(double)r[a];
        for(int a=0;a<m.n;a++){
            double f=0;for(int j=0;j<m.n;j++)f+=M[a][j]*acceleration[k][j];mass_residual=fmax(mass_residual,fabs(f-force[k][a]));force_norm=fmax(force_norm,fabs(force[k][a]));Kdot_acceleration+=m.v[k][a]*f;
        }
    }
    GfsMaterialMoments moments;need(gfs_material_moments(&m,M,&moments),"actual represented intrinsic angular/kinetic moments");
    double intrinsic_rate=Ldot-(b.gx-12+moments.center[0])*force[1][0]+(b.gy-12+moments.center[1])*force[0][0];
    printf("{\"id\":%llu,\"case\":\"%s\",\"perturb\":[%.17g,%.17g],\"basis\":%d,\"points\":%zu,\"initial_velocity_error\":%.17g,\"coordinate_error\":%.17g,\"K_minus_reference_over_U\":%.17g,\"energy_work_residual\":%.17g,\"energy_work_relative\":%.17g,\"angular_residual\":%.17g,\"angular_residual_over_2K\":%.17g,\"max_Kepler_divergence_error\":%.17g,\"max_abs_Pdot_over_P0\":%.17g,\"volume_vs_boundary_rate_over_area\":%.17g,\"finite_map_volume_rate_error_over_area\":%.17g,\"minimum_finite_J\":%.17g,\"matching_fan_mass_residual\":%.17g,\"matching_fan_work_residual\":%.17g,\"mass\":%.17g,\"K\":%.17g,\"U\":%.17g,\"native_initial_reference\":%s,\"native_initial_mass_relative_error\":%.17g,\"native_initial_energy_error_over_U\":%.17g,\"mass_matrix_solve_relative\":%.17g,\"mass_solve_kinetic_power_residual\":%.17g,\"intrinsic_angular\":%.17g,\"intrinsic_angular_rate\":%.17g,\"internal_kinetic\":%.17g,\"rigid_rotational_kinetic\":%.17g,\"resolved_strain_remainder\":%.17g}",
        (unsigned long long)cell.id,kind,bx,by,m.n,points.size(),field,coords,(K-Kdirect)/U,energy_error,energy_error/(power_abs+fabs(boundary_power)+1.e-300),Ldot-boundary_torque-Gtorque,(Ldot-boundary_torque-Gtorque)/(2*K),gradient,pressure,(Vdot-boundary_volume)/area,(Vdot-Vfd)/area,minJ,fan_mass,fan_work,mass,K,U,native?"true":"false",native_mass,native_energy,mass_residual/(force_norm+1.e-300),Kdot_acceleration-Kdot,moments.intrinsic_angular,intrinsic_rate,moments.internal_kinetic,moments.rigid_rotation_kinetic,moments.resolved_remainder);
}
static std::map<uint64_t,std::array<double,4>> native_integrals(const char *directory,const std::set<uint64_t> &ids){
    std::string path=std::string(directory)+"/enriched_project_r0.bin";FILE *f=fopen(path.c_str(),"rb");need(f,"native initial integral reference");
    uint64_t h[8];need(fread(h,sizeof(h),1,f)==1 && h[5]==1528,"native initial packet schema");std::vector<char> row(h[5]);std::map<uint64_t,std::array<double,4>> out;
    for(uint64_t k=0;k<h[2];k++){need(fread(row.data(),row.size(),1,f)==1,"native initial packet row");uint64_t id;memcpy(&id,row.data(),8);
        if(ids.count(id)){std::array<double,4> u;memcpy(u.data(),row.data()+1360,32);out[id]=u;}}
    fclose(f);need(out.size()==ids.size(),"all native initial references");return out;
}
int main(int argc,char **argv){
    need(argc==2 || argc==3,"original native polygon packet directory and optional native initial integrals for finer Kepler-only check");
    std::set<uint64_t> ids={17894,20502,39878,44945,67927,67934,67935,67942,67943,67992,67993,68000,68001,68002,68525,68526,68527,68534,68535,68584,68585,68592,68593,68600};
    auto cells=load(argv[1],ids);std::map<uint64_t,std::array<double,4>> native;if(argc==3)native=native_integrals(argv[2],ids);
    printf("{\"scope\":\"cell-local compatible smooth material RHS and finite-map volume controls; no global H1 assembly, shock stress, actual remap or timestep certified\",\"cells\":[");bool first=true;
    for(auto &entry:cells){
        for(int test=argc==3?3:0;test<5;test++){
            auto cfg=argc==3?config(64,28):config(32,20);const char *kind=test==0?"rigid":test==1?"simple_shear":test==2?"compression":test==3?"Kepler":"Kepler_perturb";
            if(test<3){auto &ic=cfg.field.ic;ic.kind=1;ic.rho=1;ic.vx=.2;ic.vy=.1;
                if(test==0){ic.B[1]=-.7;ic.B[2]=.7;}if(test==1)ic.B[1]=1.2;if(test==2){ic.B[0]=-.2;ic.B[3]=-.1;}}
            if(!first)printf(",\n");first=false;assess_material(entry.second,cfg,kind,test==4?.013:0,test==4?-.007:0,argc==3?&native[entry.first]:nullptr);fflush(stdout);
        }
    }printf("]}\n");return 0;
}
