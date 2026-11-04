/* Same accepted state, same RHS, different time integrators and step sizes. */
#define main archived_suite_main
#include "../../Hydro1DExam/laguerre_sod.c"
#undef main
#include <stdint.h>
typedef struct {uint64_t magic;int n,pid;double time,work;} SavedState;
#define PROBE_MAGIC UINT64_C(0x574f433144303031)
static void prepare(int pid) {
    setup_problem(pid);set_method_gfs();unified1d=1;pair_1d=0;mode=VORONOI;
}
static double gap(int *idx) {
    double d=1e30;*idx=-1;
    for(int i=0;i<N-1;i++)if(P[i+1].x-P[i].x<d){d=P[i+1].x-P[i].x;*idx=i;}
    return d;
}
static double energy(void) {
    double s=0;for(int i=0;i<N;i++)s+=P[i].m*(P[i].e+.5*P[i].v*P[i].v);return s;
}
int main(int argc,char **argv) {
    if(argc<2)return 2;
    const int save=!strcmp(argv[1],"save"),noh=!strcmp(argv[1],"noh");
    int rk4=argc>2&&!strcmp(argv[2],"rk4");
    double factor=argc>3?atof(argv[3]):1,t=0;
    if(!(factor>0))return 2;
    prepare(noh?PROB_NOH:PROB_SOD);
    if(save||noh)init(200);
    else {
        if(argc<5)return 2;
        FILE *f=fopen(argv[4],"rb");SavedState h;
        if(!f||fread(&h,sizeof(h),1,f)!=1||h.magic!=PROBE_MAGIC||h.n<2||h.n>NMAX)return 2;
        prepare(h.pid);N=h.n;meand=(x1_dom-x0_dom)/N;t=h.time;bc_eng_net=h.work;
        if(fread(P,sizeof(Pt),N,f)!=(size_t)N)return 2;
        fclose(f);
    }
    if(!gfs_unified1d_rhs())return 2;
    double e0=energy(),w0=bc_eng_net;
    FILE *trace=fopen("trajectory.csv","wx");
    if(!trace)return 2;
    fprintf(trace,"t,dt,min_gap,i,closing,relative_acc_pressure,relative_acc_dissipation,volume_i,volume_j,U_i,U_j,P_i,P_j,velocity_i,velocity_j\n");
    const double goal=noh?.6:.025;
    int collapsed=0,steps=0;
    while(t<goal-1e-14){
        if(!gfs_unified1d_rhs())return 2;
        int i;double d=gap(&i),dt=factor*get_dt();
        if(i<1||i>N-3){if(!noh&&!save)return 2;}
        double ai=i>0?.5*(P[i-1].p-P[i+1].p)/P[i].m:0;
        double aj=i+2<N?.5*(P[i].p-P[i+2].p)/P[i+1].m:0;
        fprintf(trace,"%.17g,%.17g,%.17g,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
            t,dt,d,i,P[i].v-P[i+1].v,aj-ai,P[i+1].av-P[i].av-(aj-ai),
            P[i].vol,P[i+1].vol,P[i].m*P[i].e,P[i+1].m*P[i+1].e,
            P[i].p,P[i+1].p,P[i].v,P[i+1].v);
        if(save&&d<.1*meand){
            FILE *f=fopen("accepted_sod_checkpoint.bin","wbx");
            SavedState h={PROBE_MAGIC,N,problem_id,t,bc_eng_net};
            if(!f||fwrite(&h,sizeof(h),1,f)!=1||fwrite(P,sizeof(Pt),N,f)!=(size_t)N)return 2;
            fclose(f);break;
        }
        if(!(d>1e-11*meand)||dt<1e-12||steps>=200000){collapsed=1;break;}
        if(t+dt>goal)dt=goal-t;
        int tries=0;
        while((rk4?rk4_step(dt):gfs_unified1d_step(dt))<0){
            dt*=.5;if(++tries>24||dt<1e-12){collapsed=1;break;}
        }
        if(collapsed)break;
        t+=dt;steps++;
    }
    fclose(trace);
    int i;double d=gap(&i);
    printf("{\"integrator\":\"%s\",\"factor\":%.17g,\"problem\":\"%s\",\"time\":%.17g,\"steps\":%d,\"collapsed\":%s,\"checkpoint_saved\":%s,\"min_gap\":%.17g,\"relative_energy_residual\":%.17g,\"thermal_floors\":%ld}\n",
        rk4?"rk4":"kdk",factor,noh?"noh":"sod",t,steps,collapsed?"true":"false",save?"true":"false",d,(energy()-e0-(bc_eng_net-w0))/e0,n_efloor_1d);
    return collapsed?1:0;
}
