#define GFS_ENRICH_TEST_LIBRARY
#include "test_enriched_euler_v1.cpp"
struct NativeCell {GfsEnrichBasis b;GfsEnrichState s;double rhs[4][10],geo[4][10];};
static void load_native(const char *path,const std::set<uint64_t> &ids,
    std::map<uint64_t,NativeCell> &cells,double perturb[2]){
    std::string root(path);FILE *f=fopen((root+"/enriched_project_r0.bin").c_str(),"rb");need(f,"native projection file");
    uint64_t h[8];need(fread(h,sizeof(h),1,f)==1 && h[5]==1528,"projection schema");
    std::vector<char> row(h[5]);
    for(uint64_t k=0;k<h[2];k++){
        need(fread(row.data(),row.size(),1,f)==1,"projection row");uint64_t id;memcpy(&id,row.data(),8);if(!ids.count(id))continue;
        NativeCell c={};memcpy(&c.b.geometry.area,row.data()+24,6*8);memcpy(c.b.generator,row.data()+72,2*8);
        memcpy(c.b.scale,row.data()+88,4*8);memcpy(c.b.L,row.data()+120,100*8);memcpy(c.b.mean,row.data()+920,10*8);
        memcpy(&c.b.minimum_norm,row.data()+1000,4*8);memcpy(&c.b.rank,row.data()+1032,2*4);memcpy(c.s.c,row.data()+1040,40*8);cells[id]=c;
    }fclose(f);
    f=fopen((root+"/enriched_coupled.bin").c_str(),"rb");need(f,"native coupled file");
    need(fread(h,sizeof(h),1,f)==1 && h[5]==688,"coupled schema");need(fread(perturb,8,2,f)==2,"perturbation metadata");row.resize(688);
    for(uint64_t k=0;k<h[2];k++){
        need(fread(row.data(),row.size(),1,f)==1,"coupled row");uint64_t id;memcpy(&id,row.data(),8);if(!ids.count(id))continue;
        need(cells.count(id),"projection contains selected cell");memcpy(cells[id].rhs,row.data()+48,40*8);memcpy(cells[id].geo,row.data()+368,40*8);
        cells[id].s=boost(cells[id].s,perturb[0],perturb[1]);
    }fclose(f);need(cells.size()==ids.size(),"all native reference cells");
}
static void projection_solve(const GfsEnrichBasis &b,double G[10][10],const double rhs[4][10],double dot[4][10]){
    int active[10],n=0;for(int k=0;k<10;k++)if(b.L[k][k]>0)active[n++]=k;
    double L[10][10]={};
    for(int i=0;i<n;i++)for(int j=0;j<=i;j++){
        double a=G[active[i]][active[j]];for(int k=0;k<j;k++)a=fma(-L[i][k],L[j][k],a);
        if(i==j){need(a>0,"independent active Gram positivity");L[i][i]=sqrt(a);}else L[i][j]=a/L[j][j];
    }
    for(int c=0;c<4;c++){
        for(int k=0;k<10;k++)dot[c][k]=0;double z[10];
        for(int i=0;i<n;i++){z[i]=rhs[c][active[i]];for(int j=0;j<i;j++)z[i]-=L[i][j]*z[j];z[i]/=L[i][i];}
        for(int i=n-1;i>=0;i--){for(int j=i+1;j<n;j++)z[i]-=L[j][i]*z[j];z[i]/=L[i][i];dot[c][active[i]]=z[i];}
    }
}
struct ThermalDerivative {
    double max=0,min=DBL_MAX,maxvalue=-DBL_MAX,true_max=0,minp=DBL_MAX,maxperror=0;
    double weight=0,L1=0,signed_rate=0,volume_max=0,volume_min=DBL_MAX,worstx=0,worsty=0;int worstkind=0;
};
static void thermal_point(const NativeCell &a,const GfsEnrichConfig &cfg,const double dot[4][10],
    double bx,double by,double x,double y,ThermalDerivative &out,double weight=0,int kind=0){
    double q[10],u[4],d[4]={};need(gfs_enrich_q(&a.b,&cfg,10,x,y,q,0,0),"thermal derivative basis");
    need(gfs_enrich_evaluate(&a.b,&cfg,&a.s,x,y,u),"thermal derivative state");
    for(int c=0;c<4;c++)for(int k=0;k<10;k++)d[c]+=dot[c][k]*(k?q[k]-a.b.mean[k]:1);
    double vx=u[1]/u[0],vy=u[2]/u[0],factor=(cfg.field.ic.gamma-1)/cfg.field.ic.pressure;
    double pdot=factor*(d[3]-vx*d[1]-vy*d[2]+.5*(vx*vx+vy*vy)*d[0]);
    if(fabs(pdot)>out.max){out.max=fabs(pdot);out.worstx=x;out.worsty=y;out.worstkind=kind;}
    out.min=fmin(out.min,pdot);out.maxvalue=fmax(out.maxvalue,pdot);
    if(weight>0){out.weight+=weight;out.L1+=weight*fabs(pdot);out.signed_rate+=weight*pdot;
        out.volume_max=fmax(out.volume_max,fabs(pdot));out.volume_min=fmin(out.volume_min,pdot);}
    double P=gfs_dg_pressure(u,cfg.field.ic.gamma);out.minp=fmin(out.minp,P);out.maxperror=fmax(out.maxperror,fabs(P/cfg.field.ic.pressure-1));
    double X=a.b.generator[0]+x-cfg.field.ic.cx,Y=a.b.generator[1]+y-cfg.field.ic.cy;
    if(X==0 && Y==0)return;double inverse=pow(X*X+Y*Y+cfg.field.ic.eps*cfg.field.ic.eps,-1.5),exact[4];
    analytic_rhs(cfg,a.b.generator[0]+x,a.b.generator[1]+y,bx,by,-X*inverse,-Y*inverse,exact);
    out.true_max=fmax(out.true_max,fabs(factor*(exact[3]-vx*exact[1]-vy*exact[2]+.5*(vx*vx+vy*vy)*exact[0])));
}
int main(int argc,char **argv){
    need(argc==3,"coupled case and original native polygon directories");
    std::set<uint64_t> ids={67943,68000,68527,68584,67935,68001,68526,68592,17894,20502,44945,39878};
    auto polygons=load(argv[2],ids);std::map<uint64_t,NativeCell> data;double perturb[2];load_native(argv[1],ids,data,perturb);
    printf("{\"scope\":\"independent selected-cell strong Euler projection versus actual native moving-cell moment RHS minus geometric transport; no time evolution\",\"perturb\":[%.17g,%.17g],\"cells\":[\n",perturb[0],perturb[1]);bool first=true;
    for(auto &entry:data){
        uint64_t id=entry.first;auto &cell=polygons[id];auto &native=entry.second;
        GfsMomentPolygon polygon={cell.x.data(),cell.y.data(),(int)cell.x.size(),0};
        for(auto rule:std::vector<std::pair<int,int>>{{32,16},{64,28}}){
            auto cfg=config(rule.first,rule.second);GfsEnrichQuadrature qr;need(gfs_enrich_quadrature_init(polygon,&cfg,&native.b,&qr),"reference quadrature");
            double G[10][10]={},strong[4][10]={},numerical[4][10]={};
            for(int n=0;n<gfs_enrich_nquad(polygon,&cfg,&qr);n++){
                double x,y,w;need(gfs_enrich_quad(polygon,&cfg,&qr,n,&x,&y,&w),"reference volume point");if(w==0)continue;
                double q[10];need(gfs_enrich_q(&native.b,&cfg,10,x,y,q,0,0),"reference basis");
                double X=cell.gx+x-12,Y=cell.gy+y-12,inverse=pow(X*X+Y*Y+.047*.047,-1.5),u[4];
                analytic_rhs(cfg,cell.gx+x,cell.gy+y,perturb[0],perturb[1],-X*inverse,-Y*inverse,u);
                for(int k=0;k<10;k++){
                    double a=k?q[k]-native.b.mean[k]:1;
                    for(int j=0;j<=k;j++)G[k][j]+=w*a*(j?q[j]-native.b.mean[j]:1);
                    for(int c=0;c<4;c++)strong[c][k]+=w*a*u[c];
                }
            }
            for(int k=0;k<10;k++)for(int j=0;j<k;j++)G[j][k]=G[k][j];
            double error=0,rate=sqrt(2*native.b.scale[3]/native.b.scale[0])/fmin(native.b.geometry.l00,native.b.geometry.l11);
            for(int c=0;c<4;c++)for(int k=0;k<10;k++){
                numerical[c][k]=native.rhs[c][k]-native.geo[c][k];
                error=fmax(error,fabs(numerical[c][k]-strong[c][k])/(native.b.geometry.area*native.b.scale[c]*rate));
            }
            double exact_dot[4][10],native_dot[4][10];projection_solve(native.b,G,strong,exact_dot);projection_solve(native.b,G,numerical,native_dot);
            ThermalDerivative pe,pn;auto check=config(48,24);GfsEnrichQuadrature checkrule;
            need(gfs_enrich_quadrature_init(polygon,&check,&native.b,&checkrule),"thermal assessment rule");
            for(int n=0;n<gfs_enrich_nquad(polygon,&check,&checkrule);n++){
                double x,y,w;need(gfs_enrich_quad(polygon,&check,&checkrule,n,&x,&y,&w),"thermal assessment volume");if(w==0)continue;
                thermal_point(native,check,exact_dot,perturb[0],perturb[1],x,y,pe,w,0);thermal_point(native,check,native_dot,perturb[0],perturb[1],x,y,pn,w,0);
            }
            for(size_t k=0;k<cell.x.size();k++){
                thermal_point(native,check,exact_dot,perturb[0],perturb[1],cell.x[k],cell.y[k],pe,0,1);
                thermal_point(native,check,native_dot,perturb[0],perturb[1],cell.x[k],cell.y[k],pn,0,1);
            }
            faces(cell,check,[&](double x,double y,double nx,double ny,double w){
                thermal_point(native,check,exact_dot,perturb[0],perturb[1],x,y,pe,0,2);thermal_point(native,check,native_dot,perturb[0],perturb[1],x,y,pn,0,2);
            });
            if(!first)printf(",\n");first=false;
            printf("{\"id\":%llu,\"angular_n\":%d,\"radial_n\":%d,\"normalized_moment_RHS_error\":%.17g,\"initial_max_relative_P_error\":%.17g,\"independent_projected_max_abs_Pdot_over_P0\":%.17g,\"native_max_abs_Pdot_over_P0\":%.17g,\"native_min_Pdot_over_P0\":%.17g,\"native_max_Pdot_over_P0\":%.17g,\"unprojected_exact_max_abs_Pdot_over_P0\":%.17g,\"native_volume_L1_Pdot_over_P0\":%.17g,\"native_volume_mean_Pdot_over_P0\":%.17g,\"native_volume_max_abs_Pdot_over_P0\":%.17g,\"native_volume_min_Pdot_over_P0\":%.17g,\"independent_projected_volume_L1_Pdot_over_P0\":%.17g,\"worst_relative_coordinates\":[%.17g,%.17g],\"worst_kind\":%d}\n",
                (unsigned long long)id,cfg.na,cfg.nr,error,pn.maxperror,pe.max,pn.max,pn.min,pn.maxvalue,pe.true_max,
                pn.L1/pn.weight,pn.signed_rate/pn.weight,pn.volume_max,pn.volume_min,pe.L1/pe.weight,pn.worstx,pn.worsty,pn.worstkind);fflush(stdout);
        }
    }printf("]}\n");return 0;
}
