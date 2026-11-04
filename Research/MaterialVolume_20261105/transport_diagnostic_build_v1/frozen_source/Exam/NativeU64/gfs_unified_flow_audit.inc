/* Read-only test observables. No force, timestep, or state correction. */
static void gfs_unified_flow_audit(SimParameters *simpar,double time,double gamma)
{
    static int enabled=-1,gresho,kh;static double next;static FILE *out;
    if(enabled<0){
        const char *flag=getenv("GFS_UNIFIED_FLOW_AUDIT"),*ic=getenv("EUNHA_IC");
        enabled=flag&&!strcmp(flag,"1");gresho=ic&&!strcmp(ic,"gresho");kh=ic&&!strcmp(ic,"kh");
        if(enabled&&!gresho&&!kh)MPI_Abort(MPI_COMM(simpar),99);
    }
    if(!enabled)return;
    const char *stop=getenv("HYDRO_TSTOP");
    if(time<next&&(!stop||time<atof(stop)))return;
    if(!out){
        out=fopen("unified_flow.csv","wx");if(!out)MPI_Abort(MPI_COMM(simpar),99);
        fprintf(out,"time,mass,Ekin,U,momentum_x,momentum_y,Lz_raw,density_L1_from_initial_profile,pressure_L1_from_initial_profile,velocity_L1_from_initial_profile,radial_velocity_L1,mode_sin,mode_cos,mode_amplitude,Ekin_y,material_cohort_overlap_width\n");
    }
    const int bins=128;static double *bin_mass,*bin_colour;static int threads;
    if(kh&&!bin_mass){
        threads=1;
#ifdef _OPENMP
        threads=omp_get_max_threads();
#endif
        bin_mass=calloc((size_t)threads*bins,sizeof(double));bin_colour=calloc((size_t)threads*bins,sizeof(double));
        if(!bin_mass||!bin_colour)MPI_Abort(MPI_COMM(simpar),99);
    }
    if(kh){memset(bin_mass,0,(size_t)threads*bins*sizeof(double));memset(bin_colour,0,(size_t)threads*bins*sizeof(double));}
    double lx=SIMBOX(simpar).x.max-SIMBOX(simpar).x.min,ly=SIMBOX(simpar).y.max-SIMBOX(simpar).y.min;
    double birth_volume=lx*ly/NX(simpar)/NY(simpar),xmid=.5*(SIMBOX(simpar).x.max+SIMBOX(simpar).x.min),ymid=.5*(SIMBOX(simpar).y.max+SIMBOX(simpar).y.min);
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(simpar);int n=VORO_NP(simpar);
    long double M=0,K=0,U=0,px=0,py=0,J=0,drho=0,dP=0,dv=0,vr=0,sine=0,cosine=0,Ky=0;
#pragma omp parallel for schedule(static) reduction(+:M,K,U,px,py,J,drho,dP,dv,vr,sine,cosine,Ky)
    for(int i=0;i<n;i++){
        double x=p[i].x,y=p[i].y,m=p[i].mass,V=p[i].weight_thermal.V,u=p[i].weight_thermal.U,vx=p[i].vx,vy=p[i].vy;
        double dx=x-xmid,dy=y-ymid,r=hypot(dx,dy),ref_rho=1,ref_P=10,ref_vx=0,ref_vy=0;
        if(gresho){
            double speed=r<=.2?5*r:r<=.4?2-5*r:0;
            ref_P=r<=.2?5+12.5*r*r:r<=.4?9+12.5*r*r-20*r+4*log(5*r):3+4*log(2.);
            if(r>0){ref_vx=-speed*dy/r;ref_vy=speed*dx/r;vr+=V*fabs((dx*vx+dy*vy)/r);}
        }else{
            double profile=.5*(tanh((y-.5)/.05)-tanh((y-1.5)/.05));ref_rho=1+profile;ref_vx=-1+2*profile;
            ref_vy=.01*sin(2*M_PI*x)*(exp(-pow((y-.5)/.2,2))+exp(-pow((y-1.5)/.2,2)));
            sine+=m*vy*sin(2*M_PI*x/lx);cosine+=m*vy*cos(2*M_PI*x/lx);
            /* A passive material label from the actual FV birth mass.
             * The overlap width measures spatial interleaving, not molecular diffusion. */
            double colour=fmin(1,fmax(0,m/birth_volume-1));int bin=(int)floor(y/ly*bins),tid=0;
            if(bin<0)bin=0;if(bin>=bins)bin=bins-1;
#ifdef _OPENMP
            tid=omp_get_thread_num();
#endif
            bin_mass[tid*bins+bin]+=m;bin_colour[tid*bins+bin]+=m*colour;
        }
        M+=m;K+=.5L*m*(vx*vx+vy*vy);Ky+=.5L*m*vy*vy;U+=u;px+=m*vx;py+=m*vy;J+=m*(dx*vy-dy*vx);
        drho+=V*fabs(m/V-ref_rho);dP+=V*fabs((gamma-1)*u/V-ref_P);dv+=V*hypot(vx-ref_vx,vy-ref_vy);
    }
    double width=0;
    if(kh)for(int b=0;b<bins;b++){
        double m=0,c=0;for(int t=0;t<threads;t++){m+=bin_mass[t*bins+b];c+=bin_colour[t*bins+b];}
        if(m>0){double mean=c/m;width+=ly/bins*4*mean*(1-mean);}
    }
    double A=2*(double)sine/(double)M,B=2*(double)cosine/(double)M;
    fprintf(out,"%.17g,%.17Lg,%.17Lg,%.17Lg,%.17Lg,%.17Lg,%.17Lg,%.17Lg,%.17Lg,%.17Lg,%.17Lg,%.17g,%.17g,%.17g,%.17Lg,%.17g\n",time,M,K,U,px,py,J,drho/(lx*ly),dP/(lx*ly),dv/(lx*ly),vr/(lx*ly),A,B,hypot(A,B),Ky,width);
    if(fflush(out))MPI_Abort(MPI_COMM(simpar),99);
    next=time+.05;
}
