/* Common mean-approach stress: frozen states and exact affine h scaling. */
#define main archived_native_main
#include "fv_native_evolution.c"
#undef main
#include "../gfs_fv_native.c"
static int saved_probe(const char *file,int wall){
    FILE *in=fopen(file,"rb");FvCheckpoint h;
    if(!in||fread(&h,sizeof(h),1,in)!=1||h.magic!=magic||h.version!=8||h.dim!=1||h.cell_bytes!=sizeof(GfsFvCell)||h.count<2||h.count>4096)return 2;
    GfsFvMesh m={0};m.dim=1;m.n=h.count;m.gamma=h.gamma;memcpy(m.bc,h.bc,sizeof(m.bc));memcpy(m.lo,h.lo,sizeof(m.lo));memcpy(m.L,h.L,sizeof(m.L));
    m.cell=calloc(m.n,sizeof(*m.cell));if(!m.cell||fread(m.cell,sizeof(*m.cell),m.n,in)!=(size_t)m.n||fgetc(in)!=EOF||fclose(in)||gfs_fv_geometry(&m))return 2;
    GfsFvCell *saved=malloc(m.n*sizeof(*saved));memcpy(saved,m.cell,m.n*sizeof(*saved));
    GfsFvWorkspace *w=gfs_fv_workspace_create();GfsFvRate *rate=calloc(m.n,sizeof(*rate));GfsFvRhsStats s;
    if(gfs_fv_native_rhs(w,&m,rate,&s))return 3;
    int chosen=-1;double gap=HUGE_VAL;
    for(int k=0;k<w->ne;k++){
        FvEdge *e=w->edge+k;
        if(wall){if(e->i==0&&e->j<0)chosen=k;}
        else if(e->j>=0){double d=fabs(m.cell[e->i].q[0]-m.cell[e->j].q[0]);if(d<gap){gap=d;chosen=k;}}
    }
    if(chosen<0)return 3;FvEdge *e=w->edge+chosen;FvFanPoint *p=w->point+e->offset;
    double heat=-(e->j<0?.5:1)*p->weight*p->compression*p->delta;
    double kinetic=fabs(s.kinetic_rate-s.kinetic_pressure_work+s.kinetic_boundary_flux);
    int unchanged=!memcmp(saved,m.cell,m.n*sizeof(*saved));
    int pass=unchanged&&heat>=0&&kinetic<1e-10*(1+fabs(s.kinetic_rate))&&fabs(s.entropy_balance_defect)<1e-10*(1+fabs(s.entropy_rate));
    if(wall){
        pass&=e->flux[0]==0&&e->flux[4]==0&&p->delta<0&&p->fraction==1;
#ifdef GFS_NATIVE_MEAN_APPROACH_STRESS
        pass&=p->compression>0&&heat>0;
#else
        pass&=p->compression==0;
#endif
    }
    double expected_viscous=HUGE_VAL;
    for(int i=0;i<m.n;i++){
        double damping=0;for(int k=w->offset[i];k<w->offset[i+1];k++)damping+=w->edge[w->incident[k]].viscous;
        if(damping>0)expected_viscous=fmin(expected_viscous,.4*m.cell[i].Q[0]/damping);
    }
    double dt_error=fabs(expected_viscous-s.viscous_dt)/(1+fabs(expected_viscous));pass&=dt_error<1e-14;
    printf("{\"case\":\"%s\",\"source_version\":8,\"operator_version\":%llu,\"pair\":[%d,%d],\"mean_jump\":%.17g,\"chi\":%.17g,\"Pi\":%.17g,\"physical_stress_heat\":%.17g,\"boundary_mass_flux\":%.17g,\"boundary_energy_flux\":%.17g,\"kinetic_identity_defect\":%.17g,\"entropy_identity_defect\":%.17g,\"viscous_dt_coefficient_error\":%.17g,\"state_and_birth_parameters_unchanged\":%s,\"accepted_steps\":0,\"pass\":%s}\n",wall?"frozen_noh_wall":"frozen_doublefan_closest_pair",(unsigned long long)checkpoint_version,e->i,e->j,p->delta,p->fraction,p->compression,heat,e->flux[0],e->flux[4],kinetic,s.entropy_balance_defect,dt_error,unchanged?"true":"false",pass?"true":"false");
    free(saved);free(rate);gfs_fv_workspace_free(w);gfs_fv_geometry_free(&m);free(m.cell);return pass?0:1;
}
static int compression_scaling(int d){
    double measured[2],linear[2],quadratic[2];int failed=0;
    for(int trial=0;trial<2;trial++){
        GfsFvMesh m={0};m.dim=d;m.n=1;m.gamma=1.4;double L=trial?.5:1;
        for(int a=0;a<3;a++){m.L[a]=L;if(a<d)m.n*=8;}
        m.cell=calloc(m.n,sizeof(*m.cell));int center=0,stride=1;
        for(int a=0;a<d;a++){center+=3*stride;stride*=8;}
        for(int i=0;i<m.n;i++){int id=i;for(int a=0;a<d;a++){m.cell[i].q[a]=(id%8+.5)*L/8;id/=8;}m.cell[i].id=i;}
        if(gfs_fv_geometry(&m))return 2;
        for(int i=0;i<m.n;i++){
            GfsFvCell *c=m.cell+i;double V=m.geometry[i].V,u=-(c->q[0]+m.geometry[i].centroid[0]);
            c->m0=c->Q[0]=V;c->mu=V/pow(V,2./d);c->Q[1]=V*u;c->Q[4]=V/(m.gamma-1)+.5*V*u*u;
        }
        GfsFvWorkspace *w=gfs_fv_workspace_create();GfsFvRate *rate=calloc(m.n,sizeof(*rate));GfsFvRhsStats s;
        if(gfs_fv_native_rhs(w,&m,rate,&s))return 3;
        int chosen=-1;for(int k=0;k<w->ne;k++)if(w->edge[k].i==center&&w->edge[k].j==center+1)chosen=k;
        if(chosen<0)return 3;FvFanPoint *p=w->point+w->edge[chosen].offset;
        double h=L/8;measured[trial]=p->compression;linear[trial]=sqrt(m.gamma)*h;quadratic[trial]=measured[trial]-linear[trial];
        failed|=fabs(measured[trial]-(sqrt(m.gamma)*h+2*h*h))>1e-12||fabs(p->fraction-1)>1e-12;
        free(rate);gfs_fv_workspace_free(w);gfs_fv_geometry_free(&m);free(m.cell);
    }
    failed|=fabs(linear[1]/linear[0]-.5)>1e-12||fabs(quadratic[1]/quadratic[0]-.25)>1e-11;
    printf("{\"case\":\"affine_compression_local_h_scaling\",\"dim\":%d,\"Pi_h\":%.17g,\"Pi_h_over2\":%.17g,\"linear_ratio\":%.17g,\"quadratic_ratio\":%.17g,\"local_operator_not_PDE_convergence\":true,\"pass\":%s}\n",d,measured[0],measured[1],linear[1]/linear[0],quadratic[1]/quadratic[0],failed?"false":"true");return failed;
}
int main(int argc,char **argv){
    if(argc!=3)return 2;int failed=saved_probe(argv[1],1)|saved_probe(argv[2],0);
#ifdef GFS_NATIVE_MEAN_APPROACH_STRESS
    for(int d=1;d<=3;d++)failed|=compression_scaling(d);
#endif
    return failed;
}
