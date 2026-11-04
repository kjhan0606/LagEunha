/* Native FV--ALE evolution harness using the common C/OpenMP operator and
 * original project geometry. This is a new experimental state format, not
 * a reinterpretation of accepted U64/variational checkpoints. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <float.h>
#include <omp.h>
#include "../gfs_fv_native.h"
enum {CONTACT,MOVING_CONTACT,SOD,GRESHO,COLD_SHEAR,ALTERNATING,BLAST,COLLISION,SHU_OSHER,NOH,LAX,DOUBLE_FAN,CONTACT_OPEN,MOVING_CONTACT_OPEN};
#include "fv_original_cases.h"
static int fixture_side;
typedef struct {
    uint64_t magic,version,cell_bytes,count,steps;int dim,side,problem,bc[3];
    double lo[3],L[3],gamma,time,boundary[5],initial[5];
} FvCheckpoint;
static const uint64_t magic=UINT64_C(0x4656414c45303031);
#if defined(GFS_NATIVE_GLOBAL_VELOCITY)
static const uint64_t checkpoint_version=4;
static const char *core_name="cell_mean_KEP_EC_ALE_LSQ_stress_global_physical_velocity";
#elif defined(GFS_NATIVE_KEP_FILTER)
static const uint64_t checkpoint_version=3;
static const char *core_name="cell_mean_KEP_EC_ALE_LSQ_compression_residual_stress";
#elif defined(GFS_NATIVE_KEP)
static const uint64_t checkpoint_version=2;
static const char *core_name="cell_mean_KEP_EC_ALE_compression_stress";
#else
static const uint64_t checkpoint_version=1;
static const char *core_name="reconstructed_HLLC_HLL_ALE";
#endif
static const uint64_t diagnostic_trial_flag=UINT64_C(0x8000000000000000);
static int write_state_kind(const char *path,GfsFvMesh *m,int side,int problem,double time,
        uint64_t steps,const double boundary[5],const double initial[5],int trial){
    FvCheckpoint h={0};h.magic=magic;h.version=checkpoint_version|(trial?diagnostic_trial_flag:0);h.cell_bytes=sizeof(GfsFvCell);h.count=m->n;h.steps=steps;h.dim=m->dim;h.side=side;h.problem=problem;
    memcpy(h.bc,m->bc,sizeof(h.bc));memcpy(h.lo,m->lo,sizeof(h.lo));memcpy(h.L,m->L,sizeof(h.L));h.gamma=m->gamma;h.time=time;memcpy(h.boundary,boundary,sizeof(h.boundary));memcpy(h.initial,initial,sizeof(h.initial));
    FILE *out=fopen(path,"wbx");if(!out)return 1;
    int bad=fwrite(&h,sizeof(h),1,out)!=1||fwrite(m->cell,sizeof(*m->cell),m->n,out)!=(size_t)m->n;
    if(fclose(out))bad=1;return bad;
}
static int write_state(const char *path,GfsFvMesh *m,int side,int problem,double time,
        uint64_t steps,const double boundary[5],const double initial[5]){
    return write_state_kind(path,m,side,problem,time,steps,boundary,initial,0);
}
static void failure_receipt(const char *stage,double time,double dt,int code,const GfsFvRhsStats *s){
    FILE *out=fopen("failure_stage.json","wx");if(!out)return;char number[64];
    if(isfinite(s->root_residual))snprintf(number,sizeof(number),"%.17g",s->root_residual);else strcpy(number,"null");
    fprintf(out,"{\"stage\":\"%s\",\"accepted_time\":%.17g,\"trial_dt\":%.17g,\"rhs_status\":%d,\"actual_PCG_iterations\":%d,\"true_scaled_residual\":%s,\"PCG_cap\":64,\"diagnostic_trial_is_restartable\":false}\n",stage,time,dt,code,s->cg,number);fclose(out);
}
typedef struct {
    double seconds;long long calls,points,hops,changed,switches;int maximum_hops,maximum_changed,owner_switch_retries;
} QueryLedger;
static void add_queries(QueryLedger *q,const GfsFvRhsStats *s,int n){
#ifdef GFS_NATIVE_GLOBAL_VELOCITY
    q->seconds+=s->velocity_query_seconds;q->calls++;q->points+=n;q->hops+=s->total_query_hops;q->changed+=s->velocity_owner_changed;q->switches+=s->velocity_owner_switches;
    if(s->maximum_query_hops>q->maximum_hops)q->maximum_hops=s->maximum_query_hops;
    if(s->velocity_owner_changed>q->maximum_changed)q->maximum_changed=s->velocity_owner_changed;
#else
    (void)q;(void)s;(void)n;
#endif
}
static int load_state(const char *path,GfsFvMesh *m,int side,int problem,double *time,
        uint64_t *steps,double boundary[5],double initial[5]){
    FILE *in=fopen(path,"rb");FvCheckpoint h;
    if(!in||fread(&h,sizeof(h),1,in)!=1||h.magic!=magic||h.version!=checkpoint_version||h.cell_bytes!=sizeof(GfsFvCell)||h.count!=(uint64_t)m->n||h.dim!=m->dim||h.side!=side||h.problem!=problem||h.gamma!=m->gamma||memcmp(h.bc,m->bc,sizeof(h.bc))||memcmp(h.lo,m->lo,sizeof(h.lo))||memcmp(h.L,m->L,sizeof(h.L)))return 1;
    if(fread(m->cell,sizeof(*m->cell),m->n,in)!=(size_t)m->n||fgetc(in)!=EOF||fclose(in))return 1;
    for(int i=0;i<m->n;i++)if(m->cell[i].id!=(uint64_t)i||!(m->cell[i].mu>0)||!(m->cell[i].m0>0))return 1;
    *time=h.time;*steps=h.steps;memcpy(boundary,h.boundary,sizeof(h.boundary));memcpy(initial,h.initial,sizeof(h.initial));return 0;
}
static void analytic(int problem,const double x[3],GfsFvPrimitive *p){
    if(problem>=BLAST){FvOriginalCase spec;if(original_case(problem,&spec))original_point(&spec,x[0],p);return;}
    *p=(GfsFvPrimitive){1,1,{0,0,0}};
    if(problem==CONTACT||problem==MOVING_CONTACT){p->rho=x[0]<.5?1:.125;p->u[0]=problem==MOVING_CONTACT?1:0;}
    if(problem==SOD){p->rho=x[0]<.5?1:.125;p->p=x[0]<.5?1:.1;}
    if(problem==COLD_SHEAR){p->p=1e-6;p->u[0]=sin(2*M_PI*x[1]);}
    if(problem==ALTERNATING)p->u[0]=((int)floor(fixture_side*x[0]))%2?-1:1;
    if(problem==GRESHO){
        double dx=x[0]-.5,dy=x[1]-.5,r=hypot(dx,dy),speed=r<=.2?5*r:r<=.4?2-5*r:0;
        p->p=r<=.2?5+12.5*r*r:r<=.4?9+12.5*r*r-20*r+4*log(5*r):3+4*log(2.);
        if(r>0){p->u[0]=-speed*dy/r;p->u[1]=speed*dx/r;}
    }
}
static void initialize(GfsFvMesh *m,int side,int problem,double *projection){
    const double xg[8]={-.9602898564975363,-.7966664774136267,-.5255324099163290,-.1834346424956498,.1834346424956498,.5255324099163290,.7966664774136267,.9602898564975363};
    const double wg[8]={.1012285362903763,.2223810344533745,.3137066458778873,.3626837833783620,.3626837833783620,.3137066458778873,.2223810344533745,.1012285362903763};
    double V=pow(1./side,m->dim),projected=0;
    for(int i=0;i<m->n;i++){
        GfsFvCell *c=m->cell+i;int cell=i;
        for(int a=0;a<m->dim;a++){c->q[a]=(cell%side+.5)/side;cell/=side;}c->id=i;
        double pressure=0,kinetic=0;
        for(int k=m->dim==3?0:3;k<(m->dim==3?8:4);k++)for(int j=m->dim>=2?0:3;j<(m->dim>=2?8:4);j++)for(int q=0;q<8;q++){
            int index[3]={q,j,k};double x[3]={0},weight=V;
            for(int a=0;a<m->dim;a++){x[a]=c->q[a]+.5*xg[index[a]]/side;weight*=.5*wg[index[a]];}
            GfsFvPrimitive p;analytic(problem,x,&p);c->Q[0]+=weight*p.rho;pressure+=weight*p.p/(m->gamma-1);
            double v2=0;for(int a=0;a<m->dim;a++){c->Q[1+a]+=weight*p.rho*p.u[a];v2+=p.u[a]*p.u[a];}kinetic+=.5*weight*p.rho*v2;
        }
        double mean_kinetic=0;for(int a=0;a<m->dim;a++)mean_kinetic+=.5*c->Q[a+1]*c->Q[a+1]/c->Q[0];
        c->Q[4]=pressure+mean_kinetic;c->m0=c->Q[0];c->mu=c->m0/pow(V,2./m->dim);projected+=kinetic-mean_kinetic;
    }*projection=projected;
}
static int valid(GfsFvMesh *m){
    for(int i=0;i<m->n;i++){GfsFvPrimitive p;if(!gfs_fv_primitive(m->dim,m->cell[i].Q,m->geometry[i].V,m->gamma,&p)||!isfinite(m->cell[i].W))return 1;}return 0;
}
static void totals(const GfsFvMesh *m,double q[5]){
    for(int c=0;c<5;c++){long double sum=0;for(int i=0;i<m->n;i++)sum+=m->cell[i].Q[c];q[c]=(double)sum;}
}
static int write_table(const char *path,const GfsFvMesh *m){
    FILE *out=fopen(path,"wx");if(!out)return 1;
    fprintf(out,"# id x y z W mass px py pz E V rho P vx_mean vy_mean vz_mean cx_offset cy_offset cz_offset mu m0 residual_integral\n");
    for(int i=0;i<m->n;i++){
        const GfsFvCell *c=m->cell+i;GfsFvPrimitive p;if(!gfs_fv_primitive(m->dim,c->Q,m->geometry[i].V,m->gamma,&p))return 1;
        fprintf(out,"%llu",(unsigned long long)c->id);for(int a=0;a<3;a++)fprintf(out," %.17g",c->q[a]);fprintf(out," %.17g",c->W);
        for(int a=0;a<5;a++)fprintf(out," %.17g",c->Q[a]);fprintf(out," %.17g %.17g %.17g",m->geometry[i].V,p.rho,p.p);
        for(int a=0;a<3;a++)fprintf(out," %.17g",p.u[a]);for(int a=0;a<3;a++)fprintf(out," %.17g",m->geometry[i].centroid[a]);fprintf(out," %.17g %.17g %.17g\n",c->mu,c->m0,c->residual_integral);
    }return fclose(out)!=0;
}
static void reference_point(int problem,const double q[3],double x[3]){
    for(int a=0;a<3;a++)x[a]=problem>=BLAST?q[a]:q[a]-floor(q[a]);
}
int main(int argc,char **argv){
    if(argc<5||argc>7){fprintf(stderr,"usage: fv_native dim side contact|moving_contact|sod|gresho|cold|alternating|blast|collision|shu|noh|lax|doublefan|contact_open|moving_contact_open target [dt_factor [restart]]\n");return 2;}
    GfsFvMesh m={0};m.dim=atoi(argv[1]);int side=atoi(argv[2]),problem=-1;double target=atof(argv[4]),factor=argc>5?atof(argv[5]):1;
    const char *names[]={"contact","moving_contact","sod","gresho","cold","alternating","blast","collision","shu","noh","lax","doublefan","contact_open","moving_contact_open"};for(int k=0;k<14;k++)if(!strcmp(argv[3],names[k]))problem=k;
    if(m.dim<1||m.dim>3||side<4||side>512||problem<0||!(target>0)||!(factor>0)||factor>1||((problem==GRESHO||problem==COLD_SHEAR)&&m.dim!=2)||(problem>=BLAST&&m.dim!=1))return 2;
    m.n=1;for(int a=0;a<m.dim;a++)m.n*=side;if(m.n>1000000)return 2;
    fixture_side=side;
    m.gamma=problem==GRESHO||problem==COLD_SHEAR?(double)(float)(5./3):1.4;
    for(int a=0;a<3;a++){m.L[a]=1;m.bc[a]=GFS_FV_PERIODIC;}if(problem==SOD)m.bc[0]=GFS_FV_REFLECT;
    FvOriginalCase original;int padding[2]={0},guard_multiplier=1;
    if(problem>=BLAST){
        const char *value=getenv("GFS_TEST_GUARD_MULTIPLIER");if(value){char *end;long number=strtol(value,&end,10);if(*end||number<1||number>2)return 2;guard_multiplier=number;}
        if(!original_case(problem,&original)||original_mesh(&m,side,&original,guard_multiplier,padding))return 2;
    }
    m.cell=calloc(m.n,sizeof(*m.cell));if(!m.cell)return 2;
    double time=0,projection=0,boundary[5]={0},initial[5]={0};uint64_t steps=0,initial_steps=0;
    if(argc==7){if(load_state(argv[6],&m,side,problem,&time,&steps,boundary,initial))return 3;initial_steps=steps;}
    else{if(problem>=BLAST){if(initialize_original(&m,side,&original,padding,&projection))return 4;}else initialize(&m,side,problem,&projection);totals(&m,initial);}
    if(problem>=BLAST&&original_receipt(&original,side,m.n,padding,guard_multiplier))return 5;
    if(target<=time||gfs_fv_geometry(&m)||valid(&m))return 4;
    if(write_state("initial.fv",&m,side,problem,time,steps,boundary,initial)||write_table("initial.dat",&m))return 5;
    if(getenv("GFS_TEST_INIT_ONLY")){
        FILE *check=fopen("initialization_check.json","wx");if(!check)return 5;
        fprintf(check,"{\"scope\":\"initialization only; no PDE evolution\",\"problem\":\"%s\",\"cells\":%d,\"gamma\":%.17g,\"physical_resolution\":%d,\"time\":%.17g,\"bc_flags\":[%d,%d,%d],\"initial_computational_interval\":[%.17g,%.17g],\"evolution_performed\":false}\n",names[problem],m.n,m.gamma,side,time,m.bc[0],m.bc[1],m.bc[2],m.lo[0],m.lo[0]+m.L[0]);
        if(fclose(check))return 5;gfs_fv_geometry_free(&m);free(m.cell);return 0;
    }
    FILE *history=fopen("history.csv","wx");if(!history)return 5;
    fprintf(history,"time,dt,steps,retries,mass,momentum_x,momentum_y,momentum_z,energy,mass_ledger_error,energy_ledger_error,weight_mass_invariant,integrated_residual_invariant,pressure_relative_L1_from_initial_profile,mean_velocity_vs_initial_node_diagnostic_L1,minP,minV,minMass,GCL_L1,root_residual,Newton,PCG,physical_velocity_node_L1_from_initial_profile,kinetic_RK_defect,kinetic_pair_defect,stress_thermal_rate\n");
    GfsFvCell *saved=malloc(m.n*sizeof(*saved));GfsFvRate *r0=calloc(m.n,sizeof(*r0)),*r1=calloc(m.n,sizeof(*r1));double *V0=malloc(m.n*sizeof(double));
    double (*node_velocity)[3]=malloc(m.n*sizeof(*node_velocity));
    if(!saved||!r0||!r1||!V0||!node_velocity)return 12;GfsFvWorkspace *work=gfs_fv_workspace_create();
    double start=omp_get_wtime(),gcl_cumulative=0,budget_cumulative=0,max_ledger=0,max_link=0,last_output=-1,recon_time=0,fan_time=0,solve_time=0,flux_time=0;
    int failed=0,total_retries=0,max_newton=0,max_cg=0;
    double kinetic_time_defect_sum=0,maximum_kinetic_pair_defect=0,integrated_stress_heat=0,maximum_entropy_balance_defect=0;
    QueryLedger query={0};
    while(time<target){
        GfsFvRhsStats s0,s1;int code=gfs_fv_native_rhs(work,&m,r0,&s0);add_queries(&query,&s0,m.n);if(code){fprintf(stderr,"FV RHS0 failed code=%d t=%.17g\n",code,time);failure_receipt("RHS0 accepted state",time,0,code,&s0);failed=1;break;}
        recon_time+=s0.reconstruction_seconds;fan_time+=s0.fan_seconds;solve_time+=s0.solve_seconds;flux_time+=s0.flux_seconds;
        double dt=fmin(target-time,factor*fmin(s0.wave_dt,s0.admissible_dt));
        memcpy(saved,m.cell,m.n*sizeof(*saved));double kinetic0=0;
        for(int i=0;i<m.n;i++){V0[i]=m.geometry[i].V;for(int a=0;a<m.dim;a++)kinetic0+=.5*saved[i].Q[a+1]*saved[i].Q[a+1]/saved[i].Q[0];}
        int accepted=0,retry=0;
        for(;retry<=8;retry++){
            if(!(dt>1e-14)||!isfinite(dt))break;
            for(int i=0;i<m.n;i++){
                m.cell[i]=saved[i];for(int a=0;a<m.dim;a++)m.cell[i].q[a]+=dt*r0[i].qdot[a];m.cell[i].W+=dt*r0[i].z;
                for(int c=0;c<5;c++)m.cell[i].Q[c]+=dt*r0[i].dQ[c];
            }
            if(gfs_fv_geometry(&m)||valid(&m)){dt*=.5;continue;}
            code=gfs_fv_native_rhs(work,&m,r1,&s1);
            add_queries(&query,&s1,m.n);
            recon_time+=s1.reconstruction_seconds;fan_time+=s1.fan_seconds;solve_time+=s1.solve_seconds;flux_time+=s1.flux_seconds;
            if(code){fprintf(stderr,"FV RHS1 failed code=%d t=%.17g dt=%.17g\n",code,time,dt);failure_receipt("RHS1 unaccepted Euler trial",time,dt,code,&s1);if(write_state_kind("failure_trial.fv",&m,side,problem,time+dt,steps,boundary,initial,1))return 5;failed=1;break;}
            if(dt>factor*fmin(s1.wave_dt,s1.admissible_dt)*(1+32*DBL_EPSILON)){if(s1.velocity_owner_switches)query.owner_switch_retries++;dt*=.5;continue;}
            int cone_bad=0;
            for(int i=0;i<m.n;i++){
                double fe[5];for(int c=0;c<5;c++)fe[c]=m.cell[i].Q[c]+dt*r1[i].dQ[c];if(!(gfs_fv_internal(m.dim,fe)>0))cone_bad=1;
                for(int a=0;a<m.dim;a++)m.cell[i].q[a]=.5*(saved[i].q[a]+m.cell[i].q[a]+dt*r1[i].qdot[a]);
                m.cell[i].W=.5*(saved[i].W+m.cell[i].W+dt*r1[i].z);
                for(int c=0;c<5;c++)m.cell[i].Q[c]=.5*(saved[i].Q[c]+fe[c]);
                m.cell[i].residual_integral=saved[i].residual_integral+.5*dt*(r0[i].residual+r1[i].residual);
            }
            if(cone_bad||gfs_fv_geometry(&m)||valid(&m)){if(s1.velocity_owner_switches)query.owner_switch_retries++;dt*=.5;continue;}
            accepted=1;break;
        }
        if(!accepted){memcpy(m.cell,saved,m.n*sizeof(*saved));if(gfs_fv_geometry(&m))return 8;failed=1;fprintf(stderr,"FV step rejected t=%.17g dt=%.17g retry=%d\n",time,dt,retry);break;}
        total_retries+=retry;steps++;time+=dt;if(target-time<16*DBL_EPSILON*target)time=target;
        for(int c=0;c<5;c++)boundary[c]+=.5*dt*(s0.boundary_flux[c]+s1.boundary_flux[c]);
        double sum[5];totals(&m,sum);double pnorm=0,vnorm=0,minP=HUGE_VAL,minV=HUGE_VAL,minM=HUGE_VAL,link=0,link_res=0,gcl=0,budget=0,kinetic1=0;
        for(int i=0;i<m.n;i++){
            GfsFvCell *c=m.cell+i;GfsFvGeometry *g=m.geometry+i;GfsFvPrimitive p,ref;gfs_fv_primitive(m.dim,c->Q,g->V,m.gamma,&p);
            double x[3];reference_point(problem,c->q,x);analytic(problem,x,&ref);
            pnorm+=g->V*fabs(p.p-ref.p)/ref.p;double dv=0;for(int a=0;a<m.dim;a++)dv+=pow(p.u[a]-ref.u[a],2);vnorm+=g->V*sqrt(dv);
            minP=fmin(minP,p.p);minV=fmin(minV,g->V);minM=fmin(minM,c->Q[0]);
            double invariant=c->Q[0]+c->mu*(c->W-c->W0)-c->m0;link=fmax(link,fabs(invariant)/c->m0);link_res=fmax(link_res,fabs(invariant+c->residual_integral)/c->m0);
            double defect=g->V-V0[i]-.5*dt*(r0[i].Vdot+r1[i].Vdot);gcl+=fabs(defect);budget+=fabs(c->Q[4]*defect/g->V);
            for(int a=0;a<m.dim;a++)kinetic1+=.5*c->Q[a+1]*c->Q[a+1]/c->Q[0];
        }
        double kinetic_defect=kinetic1-kinetic0-.5*dt*(s0.kinetic_rate+s1.kinetic_rate);
        double pair_defect=fmax(fabs(s0.kinetic_rate-s0.kinetic_pressure_work+s0.kinetic_boundary_flux),fabs(s1.kinetic_rate-s1.kinetic_pressure_work+s1.kinetic_boundary_flux));
        kinetic_time_defect_sum+=fabs(kinetic_defect);maximum_kinetic_pair_defect=fmax(maximum_kinetic_pair_defect,pair_defect);
        integrated_stress_heat+=.5*dt*(s0.stress_thermal_rate+s1.stress_thermal_rate);
        maximum_entropy_balance_defect=fmax(maximum_entropy_balance_defect,fmax(fabs(s0.entropy_balance_defect),fabs(s1.entropy_balance_defect)));
#ifdef GFS_NATIVE_KEP
        if(pair_defect>1e-10*(1+fabs(s0.kinetic_rate)+fabs(s1.kinetic_rate))){failed=1;fprintf(stderr,"KEP pair identity failed t=%.17g defect=%.17g\n",time,pair_defect);}
#endif
        gcl_cumulative+=gcl;budget_cumulative+=budget;max_link=fmax(max_link,link_res);
        double ledger=0;for(int c=0;c<5;c++)ledger=fmax(ledger,fabs(sum[c]+boundary[c]-initial[c])/(1+fabs(initial[c])));max_ledger=fmax(max_ledger,ledger);
        if(ledger>1e-10||link_res>1e-9){failed=1;fprintf(stderr,"FV conservation failure t=%.17g ledger=%.17g invariant=%.17g\n",time,ledger,link_res);}
        max_newton=max_newton>s0.newton?max_newton:s0.newton;max_newton=max_newton>s1.newton?max_newton:s1.newton;
        max_cg=max_cg>s0.cg?max_cg:s0.cg;max_cg=max_cg>s1.cg?max_cg:s1.cg;
        if(time-last_output>=.01||time>=target||failed){
            if(gfs_fv_native_velocity(work,&m,node_velocity)){failed=1;break;}
            double node_error=0;
            for(int i=0;i<m.n;i++){
                double x[3],error=0;reference_point(problem,m.cell[i].q,x);GfsFvPrimitive ref;analytic(problem,x,&ref);
                for(int a=0;a<m.dim;a++)error+=pow(node_velocity[i][a]-ref.u[a],2);node_error+=m.geometry[i].V*sqrt(error);
            }
            fprintf(history,"%.17g,%.17g,%llu,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%d,%.17g,%.17g,%.17g,%.17g\n",time,dt,(unsigned long long)steps,total_retries,sum[0],sum[1],sum[2],sum[3],sum[4],sum[0]+boundary[0]-initial[0],sum[4]+boundary[4]-initial[4],link,link_res,pnorm,vnorm,minP,minV,minM,gcl,fmax(s0.root_residual,s1.root_residual),s0.newton+s1.newton,s0.cg+s1.cg,node_error,kinetic_defect,pair_defect,.5*(s0.stress_thermal_rate+s1.stress_thermal_rate));fflush(history);last_output=time;
            fprintf(stderr,"FV progress dim=%d N=%d t=%.9g/%.9g step=%llu dt=%.3g root=%.3g initial_profile_Perror=%.6g\n",m.dim,m.n,time,target,(unsigned long long)steps,dt,fmax(s0.root_residual,s1.root_residual),pnorm);
        }
        if(failed||steps-initial_steps>100000){failed=1;break;}
    }
    double seconds=omp_get_wtime()-start;
    if(write_state("final.fv",&m,side,problem,time,steps,boundary,initial)||write_table("final.dat",&m))return 5;
    if(!gfs_fv_native_velocity(work,&m,node_velocity)){
        FILE *node=fopen("final_velocity.dat","wx");if(!node)return 5;
        fprintf(node,"# id reconstructed_physical_vx_at_q vy_at_q vz_at_q\n");
        for(int i=0;i<m.n;i++)fprintf(node,"%llu %.17g %.17g %.17g\n",(unsigned long long)m.cell[i].id,node_velocity[i][0],node_velocity[i][1],node_velocity[i][2]);
        if(fclose(node))return 5;
    }
    FILE *out=fopen("result.json","wx");if(!out)return 5;
    fprintf(out,"{\"scope\":\"native C/OpenMP FV--ALE with actual original Voro geometry; experimental format\",\"dim\":%d,\"side\":%d,\"cells\":%d,\"problem\":\"%s\",\"time\":%.17g,\"target\":%.17g,\"steps\":%llu,\"new_steps\":%llu,\"retries\":%d,\"failed\":%s,\"wall_seconds\":%.17g,\"geometry_seconds\":%.17g,\"reconstruction_seconds\":%.17g,\"fan_seconds\":%.17g,\"solve_seconds\":%.17g,\"flux_seconds\":%.17g,\"maximum_ledger_error\":%.17g,\"maximum_integrated_residual_invariant_error\":%.17g,\"cumulative_GCL_L1\":%.17g,\"cumulative_GCL_energy_budget\":%.17g,\"initial_continuum_kinetic_projection_difference\":%.17g,\"maximum_Newton\":%d,\"maximum_CG_total_per_rhs\":%d,\"maximum_geometry_candidates\":%d,\"maximum_volume_closure\":%.17g,\"floor_or_repair_count\":0}\n",m.dim,side,m.n,names[problem],time,target,(unsigned long long)steps,(unsigned long long)(steps-initial_steps),total_retries,failed?"true":"false",seconds,m.geometry_seconds,recon_time,fan_time,solve_time,flux_time,max_ledger,max_link,gcl_cumulative,budget_cumulative,projection,max_newton,max_cg,m.maximum_candidates,m.maximum_volume_closure);fclose(out);fclose(history);
    FILE *kinetic=fopen("kinetic_ledger.json","wx");if(!kinetic)return 5;
    fprintf(kinetic,"{\"core\":\"%s\",\"checkpoint_version\":%llu,\"new_segment_only\":true,\"has_pair_pressure_identity\":%s,\"kinetic_RK_defect_abs_sum\":%.17g,\"maximum_kinetic_pair_defect\":%.17g,\"integrated_compression_thermal_production\":%.17g,\"maximum_semidiscrete_entropy_balance_defect\":%.17g}\n",core_name,(unsigned long long)checkpoint_version,checkpoint_version>=2?"true":"false",kinetic_time_defect_sum,maximum_kinetic_pair_defect,integrated_stress_heat,maximum_entropy_balance_defect);
    if(fclose(kinetic))return 5;
    FILE *queries=fopen("velocity_query.json","wx");if(!queries)return 5;
    double minimum=HUGE_VAL,minimum_ratio=HUGE_VAL;
    for(int i=0;i<m.n;i++)for(int k=0;k<m.geometry[i].faces;k++){
        GfsFvNativeFace *f=m.geometry[i].face+k;if(f->mirror)continue;
        minimum=fmin(minimum,f->distance);minimum_ratio=fmin(minimum_ratio,f->distance/sqrt(m.cell[i].m0/m.cell[i].mu));
    }
    fprintf(queries,"{\"core\":\"%s\",\"scope\":\"RHS calls including rejected trials; query seconds are a subset of reconstruction time\",\"query_seconds\":%.17g,\"RHS_calls\":%lld,\"points\":%lld,\"owner_different_from_own_count\":%lld,\"owner_different_fraction\":%.17g,\"owner_transition_count\":%lld,\"maximum_different_owners_per_call\":%d,\"average_hops\":%.17g,\"maximum_hops\":%d,\"rollback_attempts_with_owner_transition\":%d,\"minimum_final_face_neighbor_generator_distance\":%.17g,\"minimum_final_face_neighbor_distance_over_birth_h\":%.17g}\n",core_name,query.seconds,query.calls,query.points,query.changed,query.points?(double)query.changed/query.points:0,query.switches,query.maximum_changed,query.points?(double)query.hops/query.points:0,query.maximum_hops,query.owner_switch_retries,minimum,minimum_ratio);
    if(fclose(queries))return 5;
    gfs_fv_workspace_free(work);gfs_fv_geometry_free(&m);free(m.cell);free(saved);free(r0);free(r1);free(V0);free(node_velocity);return failed||time!=target?1:0;
}
