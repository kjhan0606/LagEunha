/* Read-only thermal chain-rule and analytic physical-flux comparison.
 * This is not a modification of the failed FV candidate. */
#define main archived_fv_evolution_main
#include "fv_native_evolution.c"
#undef main
int main(int argc,char **argv){
    if(argc!=3)return 2;FILE *in=fopen(argv[1],"rb");FvCheckpoint h;
    if(!in||fread(&h,sizeof(h),1,in)!=1||fclose(in)||h.magic!=magic||h.dim!=2||h.count>1000000||h.problem!=COLD_SHEAR)return 2;
    GfsFvMesh m={0};m.dim=2;m.n=h.count;m.gamma=h.gamma;memcpy(m.bc,h.bc,sizeof(m.bc));memcpy(m.lo,h.lo,sizeof(m.lo));memcpy(m.L,h.L,sizeof(m.L));m.cell=calloc(m.n,sizeof(*m.cell));
    double time,boundary[5],initial[5];uint64_t steps;
    if(!m.cell||load_state(argv[1],&m,h.side,h.problem,&time,&steps,boundary,initial)||gfs_fv_geometry(&m))return 3;
    GfsFvRate *rate=calloc(m.n,sizeof(*rate));GfsFvRhsStats stats;GfsFvWorkspace *work=gfs_fv_workspace_create();
    int code=gfs_fv_native_rhs(work,&m,rate,&stats);if(code){fprintf(stderr,"Readonly native RHS failed: %d\n",code);return 4;}
    FILE *out=fopen(argv[2],"wx");if(!out)return 5;
    fprintf(out,"id,P,rho_dot,Vdot,native_Udot,analytic_flux_proxy_Udot,native_relative_Pdot,analytic_flux_relative_Pdot\n");
    double Udot=0,exact_Udot=0,Pdot_L1=0,exact_Pdot_L1=0,rho_dot_L1=0,volume_rate_L1=0,energy_rate=0,difference_L1=0;
    for(int i=0;i<m.n;i++){
        GfsFvGeometry *g=m.geometry+i;GfsFvPrimitive p;gfs_fv_primitive(2,m.cell[i].Q,g->V,m.gamma,&p);double exact[5]={0};
        for(int k=0;k<g->faces;k++){
            GfsFvNativeFace *f=g->face+k;int j=f->j;
            if(f->mirror)return 6;
            for(int q=0;q<f->count;q++){
                GfsFvPoint *point=g->point+f->point+q;double x[3],wn=(rate[i].z-rate[j].z)/(2*f->distance),un=0,v2=0;
                for(int a=0;a<3;a++)x[a]=m.cell[i].q[a]+point->x[a]-floor(m.cell[i].q[a]+point->x[a]);GfsFvPrimitive physical;analytic(COLD_SHEAR,x,&physical);
                for(int a=0;a<2;a++){
                    wn+=.5*f->normal[a]*(rate[i].qdot[a]+rate[j].qdot[a])-(rate[j].qdot[a]-rate[i].qdot[a])*(point->x[a]-.5*f->r[a])/f->distance;
                    un+=physical.u[a]*f->normal[a];v2+=physical.u[a]*physical.u[a];
                }
                double energy=physical.p/(m.gamma-1)+.5*physical.rho*v2,massflux=physical.rho*(un-wn),A=point->weight;
                exact[0]-=A*massflux;
                for(int a=0;a<2;a++)exact[a+1]-=A*(massflux*physical.u[a]+physical.p*f->normal[a]);
                exact[4]-=A*((energy+physical.p)*un-wn*energy);
            }
        }
        double v2=0,thermal=rate[i].dQ[4],physical=exact[4];
        for(int a=0;a<2;a++){v2+=p.u[a]*p.u[a];thermal-=p.u[a]*rate[i].dQ[a+1];physical-=p.u[a]*exact[a+1];}
        thermal+=.5*v2*rate[i].dQ[0];physical+=.5*v2*exact[0];
        double internal=gfs_fv_internal(2,m.cell[i].Q),pdot=(m.gamma-1)*(thermal-internal*rate[i].Vdot/g->V)/g->V;
        double epdot=(m.gamma-1)*(physical-internal*rate[i].Vdot/g->V)/g->V,rhodot=(rate[i].dQ[0]-p.rho*rate[i].Vdot)/g->V;
        fprintf(out,"%llu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",(unsigned long long)m.cell[i].id,p.p,rhodot,rate[i].Vdot,thermal,physical,pdot/p.p,epdot/p.p);
        Udot+=thermal;exact_Udot+=physical;Pdot_L1+=g->V*fabs(pdot)/1e-6;exact_Pdot_L1+=g->V*fabs(epdot)/1e-6;rho_dot_L1+=g->V*fabs(rhodot);volume_rate_L1+=fabs(rate[i].Vdot);energy_rate+=rate[i].dQ[4];difference_L1+=fabs(thermal-physical);
    }
    int bad=fclose(out);gfs_fv_workspace_free(work);gfs_fv_geometry_free(&m);free(m.cell);free(rate);
    printf("{\"scope\":\"readonly native thermal chain rule and analytic physical face flux on same saved geometry and mesh velocity; analytic flux proxy is not an evolved state\",\"time\":%.17g,\"side\":%d,\"thermal_rate\":%.17g,\"analytic_flux_proxy_thermal_rate\":%.17g,\"relative_pressure_rate_L1_over_initial_P\":%.17g,\"analytic_flux_proxy_relative_pressure_rate_L1_over_initial_P\":%.17g,\"thermal_rate_difference_L1\":%.17g,\"density_rate_L1\":%.17g,\"volume_rate_L1\":%.17g,\"total_energy_rate\":%.17g}\n",time,h.side,Udot,exact_Udot,Pdot_L1,exact_Pdot_L1,difference_L1,rho_dot_L1,volume_rate_L1,energy_rate);return bad?1:0;
}
