/* Native Cartesian cold birth RHS; no PDE time update. */
#define main archived_fv_evolution_main
#include "fv_native_evolution.c"
#undef main
int main(int argc,char **argv){
    if(argc!=2)return 2;int side=atoi(argv[1]);if(side!=64&&side!=128)return 2;
    GfsFvMesh m={0};m.dim=2;m.n=side*side;m.gamma=(double)(float)(5./3);m.cell=calloc(m.n,sizeof(*m.cell));
    for(int a=0;a<3;a++){m.L[a]=1;m.bc[a]=GFS_FV_PERIODIC;}
    double projection;initialize(&m,side,COLD_SHEAR,&projection);if(gfs_fv_geometry(&m)||valid(&m))return 3;
    GfsFvWorkspace *w=gfs_fv_workspace_create();GfsFvRate *rates=calloc(m.n,sizeof(*rates));GfsFvRhsStats stats;
    int code=gfs_fv_native_rhs(w,&m,rates,&stats);if(code)return 4;
    double thermal_rate=0,pressure_rate_L1=0;
    for(int i=0;i<m.n;i++){
        GfsFvPrimitive p;gfs_fv_primitive(2,m.cell[i].Q,m.geometry[i].V,m.gamma,&p);double internal=gfs_fv_internal(2,m.cell[i].Q),udot=rates[i].dQ[4],v2=0;
        for(int a=0;a<2;a++){udot-=p.u[a]*rates[i].dQ[a+1];v2+=p.u[a]*p.u[a];}udot+=.5*v2*rates[i].dQ[0];thermal_rate+=udot;
        pressure_rate_L1+=(m.gamma-1)*fabs(udot-internal*rates[i].Vdot/m.geometry[i].V)/1e-6;
    }
    double pair=fabs(stats.kinetic_rate-stats.kinetic_pressure_work);
    int pass=fabs(thermal_rate)<1e-12&&pressure_rate_L1<1e-6&&pair<1e-12&&stats.stress_thermal_rate<1e-12&&fabs(stats.entropy_balance_defect)<1e-10;
    printf("{\"scope\":\"native KEP cold birth RHS, no PDE time update\",\"side\":%d,\"thermal_rate\":%.17g,\"relative_pressure_rate_L1\":%.17g,\"kinetic_pair_defect\":%.17g,\"compression_thermal_rate\":%.17g,\"true_weight_residual\":%.17g,\"PCG\":%d,\"actual_geometry_entropy_balance\":%.17g,\"pass\":%s}\n",side,thermal_rate,pressure_rate_L1,pair,stats.stress_thermal_rate,stats.root_residual,stats.cg,stats.entropy_balance_defect,pass?"true":"false");
    gfs_fv_workspace_free(w);gfs_fv_geometry_free(&m);free(m.cell);free(rates);return pass?0:1;
}
