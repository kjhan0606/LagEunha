/* Read-only use of the ACTUAL native implementation and its diagnostic
 * arrays. No surrogate reconstruction and no late-state continuation. */
#define main archived_native_main
#include "fv_native_evolution.c"
#undef main
#include "../gfs_fv_native.c"
#ifndef GFS_AUDIT_SOURCE_VERSION
#define GFS_AUDIT_SOURCE_VERSION 2
#endif
static int inspect(GfsFvMesh *m,const char *name,const char *csv,int mode){
    GfsFvWorkspace *w=gfs_fv_workspace_create();GfsFvRate *r=calloc(m->n,sizeof(*r));GfsFvRhsStats stats;
    if(gfs_fv_native_rhs(w,m,r,&stats))return 2;
    double *old=calloc(m->n,sizeof(double)),*updated=calloc(m->n,sizeof(double));
    double H0=0,H1=0,predicted=0,Hbulk=0,Hresidual=0,orthogonal=0,fixture_error=0;int checked=0;
    for(int k=0;k<w->ne;k++){
        FvEdge *e=w->edge+k;int winner=e->j>=0&&w->stress_fraction[e->j]>w->stress_fraction[e->i]?e->j:e->i;
        double den=w->parallel_energy[winner]+w->perpendicular_energy[winner];
        double bulk=den>0?w->bulk_fraction[winner]*w->parallel_energy[winner]/den:0,residual=den>0?w->perpendicular_energy[winner]/den:0;
        for(int q=0;q<e->count;q++){
            FvFanPoint *p=w->point+e->offset+q;double factor=(e->j>=0?1:.5)*p->weight*(-p->delta);
            double before=p->original_compression*factor,after=p->compression*factor;
            H0+=before;H1+=after;predicted+=p->fraction*before;Hbulk+=bulk*before;Hresidual+=residual*before;
            old[e->i]+=(e->j>=0?.5:1)*before;updated[e->i]+=(e->j>=0?.5:1)*after;
            if(e->j>=0){old[e->j]+=.5*before;updated[e->j]+=.5*after;}
        }
    }
    FILE *out=csv?fopen(csv,"wx"):NULL;if(csv&&!out)return 3;
    if(out)fprintf(out,"id,Eparallel,Eperpendicular,Edata,compression_fraction,chi,old_heat_share,new_heat_share\n");
    for(int i=0;i<m->n;i++){
        orthogonal=fmax(orthogonal,fabs(w->data_energy[i]-w->parallel_energy[i]-w->perpendicular_energy[i]-w->cross_energy[i])/(1+w->data_energy[i]));
        if(out)fprintf(out,"%llu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",(unsigned long long)m->cell[i].id,w->parallel_energy[i],w->perpendicular_energy[i],w->data_energy[i],w->bulk_fraction[i],w->stress_fraction[i],old[i],updated[i]);
        int interior=1;for(int a=0;a<m->dim;a++)if(m->cell[i].q[a]<.25||m->cell[i].q[a]>.75)interior=0;
        if(mode>=0&&(mode==5||interior)){
            double expected=mode>=4?1:0;fixture_error=fmax(fixture_error,fabs(w->stress_fraction[i]-expected));checked++;
        }
    }
    if(out&&fclose(out))return 3;
    double kinetic=fabs(stats.kinetic_rate-stats.kinetic_pressure_work+stats.kinetic_boundary_flux);
    int pass=orthogonal<1e-11&&fabs(H1-predicted)<1e-12*(1+H0)&&fabs(H1-Hbulk-Hresidual)<1e-12*(1+H0)&&fabs(H1-stats.stress_thermal_rate)<1e-12*(1+H0)&&kinetic<1e-10&&fabs(stats.entropy_balance_defect)<1e-10*(1+fabs(stats.entropy_rate));
    if(mode>=0)pass&=checked>0&&fixture_error<1e-12;
    if(mode==5)pass&=H0>0&&fabs(H1/H0-1)<1e-12;
    printf("{\"case\":\"%s\",\"dim\":%d,\"cells\":%d,\"checked_cells\":%d,\"chi_fixture_error\":%.17g,\"LSQ_data_identity_with_cross_term_error\":%.17g,\"original_heat_rate\":%.17g,\"filtered_heat_rate\":%.17g,\"same_state_predicted_heat_rate\":%.17g,\"compression_component_heat_rate\":%.17g,\"residual_component_heat_rate\":%.17g,\"kinetic_pair_defect\":%.17g,\"actual_geometry_entropy_balance\":%.17g,\"pass\":%s}\n",name,m->dim,m->n,checked,fixture_error,orthogonal,H0,H1,predicted,Hbulk,Hresidual,kinetic,stats.entropy_balance_defect,pass?"true":"false");
    free(old);free(updated);free(r);gfs_fv_workspace_free(w);return pass?0:1;
}
#ifndef GFS_STRESS_AUDIT_ENTRY
#define GFS_STRESS_AUDIT_ENTRY main
#endif
int GFS_STRESS_AUDIT_ENTRY(int argc,char **argv){
    if(argc==2&&!strcmp(argv[1],"fixtures")){
        const char *names[]={"translation","affine_shear","rotation","traceless_strain","planar_compression","alternating"};int failed=0;
        for(int d=1;d<=3;d++)for(int mode=0;mode<6;mode++){
            if(d==1&&mode>=1&&mode<=3)continue;
            GfsFvMesh m={0};m.dim=d;m.n=1;m.gamma=1.4;for(int a=0;a<3;a++)m.L[a]=1;
            for(int a=0;a<d;a++)m.n*=8;m.cell=calloc(m.n,sizeof(*m.cell));
            for(int i=0;i<m.n;i++){int index=i;for(int a=0;a<d;a++){m.cell[i].q[a]=(index%8+.5)/8;index/=8;}m.cell[i].id=i;}
            if(gfs_fv_geometry(&m))return 2;
            for(int i=0;i<m.n;i++){
                GfsFvCell *c=m.cell+i;double x[3],u[3]={0};for(int a=0;a<3;a++)x[a]=c->q[a]+m.geometry[i].centroid[a];
                if(mode==0)for(int a=0;a<d;a++)u[a]=.5;
                if(mode==1)u[0]=x[1];if(mode==2){u[0]=-x[1];u[1]=x[0];}
                if(mode==3){u[0]=x[0];u[1]=-x[1];}if(mode==4)u[0]=-x[0];if(mode==5)u[0]=i%2?-.5:.5;
                c->Q[0]=m.geometry[i].V;c->m0=c->Q[0];c->mu=c->m0/pow(m.geometry[i].V,2./d);c->Q[4]=m.geometry[i].V/(m.gamma-1);
                for(int a=0;a<d;a++){c->Q[a+1]=c->Q[0]*u[a];c->Q[4]+=.5*c->Q[0]*u[a]*u[a];}
            }
            failed|=inspect(&m,names[mode],NULL,mode);gfs_fv_geometry_free(&m);free(m.cell);
        }return failed;
    }
    if(argc!=3)return 2;FILE *in=fopen(argv[1],"rb");FvCheckpoint h;
    if(!in||fread(&h,sizeof(h),1,in)!=1||h.magic!=magic||h.version!=GFS_AUDIT_SOURCE_VERSION||h.cell_bytes!=sizeof(GfsFvCell)||h.dim!=2||h.count>1000000||(h.problem!=COLD_SHEAR&&h.problem!=GRESHO))return 2;
    GfsFvMesh m={0};m.dim=h.dim;m.n=h.count;m.gamma=h.gamma;memcpy(m.bc,h.bc,sizeof(m.bc));memcpy(m.lo,h.lo,sizeof(m.lo));memcpy(m.L,h.L,sizeof(m.L));m.cell=calloc(m.n,sizeof(*m.cell));
    if(!m.cell||fread(m.cell,sizeof(*m.cell),m.n,in)!=(size_t)m.n||fgetc(in)!=EOF||fclose(in)||gfs_fv_geometry(&m))return 2;
    int code=inspect(&m,h.problem==COLD_SHEAR?"frozen_cold64":"frozen_gresho64",argv[2],-1);gfs_fv_geometry_free(&m);free(m.cell);return code;
}
