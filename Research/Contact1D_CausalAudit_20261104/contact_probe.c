/* Audit wrapper. The copied production solver is included unchanged. */
#define main archived_hydro_suite_main
#include "source/Hydro1DExam/laguerre_sod.c"
#undef main

static void measure(const char *path, double *out) {
    double l1=0, li=0, le=0, pmax=0, vmax=0, amax=0, esum=0, msum=0;
    FILE *f=fopen(path,"w");
    if(!f) { perror(path); exit(2); }
    fprintf(f,"# x rho v p e w vol mass face_left face_right acceleration\n");
    for(int i=0;i<N;i++) {
        double a=face[i], b=face[i+1], v=b-a;
        double overlap=fmax(0.0,fmin(b,0.5)-fmax(a,0.0));
        double rr=0.125+0.875*overlap/v;
        double sample=P[i].x<0.5?1.0:0.125;
        double left=fmax(0.0,fmin(b,0.5)-a);
        double right=fmax(0.0,b-fmax(a,0.5));
        l1+=fabs(P[i].rho-sample)*v;
        li+=fabs(P[i].rho-rr)*v;
        le+=fabs(P[i].rho-1.0)*left+fabs(P[i].rho-0.125)*right;
        pmax=fmax(pmax,fabs(P[i].p-1));
        vmax=fmax(vmax,fabs(P[i].v)); amax=fmax(amax,fabs(P[i].av));
        esum+=P[i].m*(P[i].e+0.5*P[i].v*P[i].v); msum+=P[i].m;
        fprintf(f,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
            P[i].x,P[i].rho,P[i].v,P[i].p,P[i].e,P[i].w,P[i].vol,P[i].m,a,b,P[i].av);
    }
    fclose(f);
    out[0]=l1; out[1]=li; out[2]=le; out[3]=pmax; out[4]=vmax; out[5]=amax;
    out[6]=esum; out[7]=msum; out[8]=face[0]; out[9]=face[N];
}

int main(int argc,char **argv) {
    if(argc!=6) { fprintf(stderr,"usage: probe method scenario N end_time dt_factor\n"); return 2; }
    const char *method=argv[1], *scenario=argv[2];
    int n=atoi(argv[3]); double end=atof(argv[4]), factor=atof(argv[5]);
    if(n<10||n>NMAX-2||!(factor>0)||end<0) return 2;
    setup_problem(PROB_CONTACT);
    if(!strcmp(method,"gfs")) { set_method_gfs(); mode=VORONOI; }
    else if(!strcmp(method,"hllc")) { set_method_hllc_cd10(); mode=VORONOI; }
    else if(!strcmp(method,"ns")) { set_method_ns_av(); mode=VORONOI; }
    else if(!strcmp(method,"sph")) { set_method_sph(); mode=SPH_MODE; }
    else return 2;
    ic_type=!strcmp(scenario,"equalvolume")?0:1;
    init(n);
    if(!strcmp(scenario,"pressure_only")) {
        for(int i=0;i<N;i++) P[i].e=1.0/(GM1*P[i].rho);
        update_derived();
    } else if(!strcmp(scenario,"fv_projected")) {
        if(mode==SPH_MODE) return 2;
        for(int i=0;i<N;i++) {
            double a=face[i], b=face[i+1];
            double left=fmax(0.0,fmin(b,0.5)-fmax(a,0.0));
            double right=fmax(0.0,fmin(b,1.0)-fmax(a,0.5));
            P[i].m=left+0.125*right;
            P[i].e=(b-a)/(GM1*P[i].m);
        }
        update_derived();
    } else if(!strcmp(scenario,"laguerre_aligned")) {
        if(n%9||mode==SPH_MODE) return 2;
        direct_hll=1; rk_face=1; xf_ready=1;
        xf_face[0]=0;
        for(int i=0;i<N;i++) {
            double rho=i<NL_actual?1.0:0.125;
            xf_face[i+1]=xf_face[i]+P[i].m/rho;
            P[i].e=1.0/(GM1*rho);
        }
        xf_face[N]=1.0; update_derived();
    } else if(strcmp(scenario,"equalvolume")&&strcmp(scenario,"original")) return 2;
    for(int i=0;i<N;i++) P[i].K=P[i].p/pow(P[i].rho,GAM);
    compute_rhs();
    double initial[10],final[10];
    measure("initial.dat",initial);
    double xinit[NMAX],massinit[NMAX];
    for(int i=0;i<N;i++) {xinit[i]=P[i].x;massinit[i]=P[i].m;}
    n_pair_1d=n_efloor_1d=n_retry_1d=n_wlim_1d=0; bc_eng_net=0;
    double t=0,dtmin=1e30; int steps=0,failed=0;
    while(t<end-1e-14) {
        double dt=factor*get_dt(); if(t+dt>end) dt=end-t;
        int tries=0;
        while(rk4_step(dt)<0) {
            dt*=0.5; n_retry_1d++;
            if(++tries>24||dt<1e-12) {failed=1;break;}
        }
        if(failed||dt<1e-12||steps>=200000) {failed=1;break;}
        ++steps;
        dtmin=fmin(dtmin,dt); t+=dt;
    }
    compute_rhs(); measure("final.dat",final);
    double xchange=0,mchange=0;
    for(int i=0;i<N;i++) {
        xchange=fmax(xchange,fabs(P[i].x-xinit[i]));
        mchange=fmax(mchange,fabs(P[i].m-massinit[i]));
    }
    printf("{\"method\":\"%s\",\"scenario\":\"%s\",\"N\":%d,\"NL\":%d,\"NR\":%d,\"t\":%.17g,\"dt_factor\":%.17g,\"steps\":%d,\"failed\":%d,\"initial\":[",
        method,scenario,N,NL_actual,NR_actual,t,factor,steps,failed);
    for(int k=0;k<10;k++) printf("%s%.17g",k?",":"",initial[k]);
    printf("],\"final\":[");
    for(int k=0;k<10;k++) printf("%s%.17g",k?",":"",final[k]);
    printf("],\"relative_energy_residual\":%.17g,\"boundary_energy\":%.17g,\"max_position_change\":%.17g,\"max_mass_change\":%.17g,\"pair_calls\":%ld,\"floor_hits\":%ld,\"retries\":%ld}\n",
        (final[6]-initial[6]-bc_eng_net)/initial[6],bc_eng_net,xchange,mchange,
        n_pair_1d,n_efloor_1d,n_retry_1d);
    return 0;
}
