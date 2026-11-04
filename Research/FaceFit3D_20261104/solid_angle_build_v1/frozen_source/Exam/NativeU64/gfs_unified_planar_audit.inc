/* Opt-in observation of a square, material-column planar IC. Never changes
 * the integration state. IDs encode row*N+column only in this diagnostic. */
static void gfs_unified_planar_audit(SimParameters *simpar,const char *stage,
        double time,double gamma)
{
    const char *e=getenv("GFS_UNIFIED_PLANAR_AUDIT_N");if(!e)return;
    int side=atoi(e),n=VORO_NP(simpar);
    if(side<2||(long long)side*side!=n||!gfs_unified_periodic_mode()){
        fprintf(stderr,"[UNIFIED_PLANAR] incompatible diagnostic IC\n");MPI_Abort(MPI_COMM(simpar),99);
    }
    enum {FIELDS=7};static double *minimum,*maximum;static int capacity;
    if(side>capacity){
        free(minimum);free(maximum);
        minimum=malloc((size_t)side*FIELDS*sizeof(double));maximum=malloc((size_t)side*FIELDS*sizeof(double));
        if(!minimum||!maximum)MPI_Abort(MPI_COMM(simpar),99);capacity=side;
    }
    for(int k=0;k<side*FIELDS;k++){minimum[k]=HUGE_VAL;maximum[k]=-HUGE_VAL;}
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(simpar);
    double ly=SIMBOX(simpar).y.max-SIMBOX(simpar).y.min,ylo=SIMBOX(simpar).y.min;
    double maxvy=0,maxay=0;unsigned long long idvy=0,iday=0;
    static FILE *trace;const char *steps_env=getenv("GFS_UNIFIED_PLANAR_TRACE_STEPS");
    int tracing=steps_env&&gfs_rk4_context_step<=atoi(steps_env);
    if(tracing&&!trace){
        trace=fopen("unified_planar_trace.csv","wx");
        if(!trace)MPI_Abort(MPI_COMM(simpar),99);
        fprintf(trace,"step,time,stage,id,x,y,m,vx,vy,U,V,ax,ay\n");
    }
    for(int i=0;i<n;i++){
        unsigned long long id=(unsigned long long)PINDX(p+i);
        if(id>=(unsigned long long)n)MPI_Abort(MPI_COMM(simpar),99);
        int col=id%side,row=id/side;
        double V=p[i].weight_thermal.V,U=p[i].weight_thermal.U;
        double q[FIELDS]={p[i].x,remainder(p[i].y-ylo-(row+.5)*ly/side,ly),
            p[i].vx,p[i].vy,U,V,(gamma-1)*U/V};
        for(int a=0;a<FIELDS;a++){
            int k=col*FIELDS+a;minimum[k]=fmin(minimum[k],q[a]);maximum[k]=fmax(maximum[k],q[a]);
        }
        if(fabs(p[i].vy)>maxvy){maxvy=fabs(p[i].vy);idvy=id;}
        if(fabs(p[i].ay)>maxay){maxay=fabs(p[i].ay);iday=id;}
        if(tracing)fprintf(trace,"%d,%.17g,%s,%llu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
            gfs_rk4_context_step,time,stage,id,p[i].x,p[i].y,p[i].mass,p[i].vx,p[i].vy,U,V,p[i].ax,p[i].ay);
    }
    double span[FIELDS]={0};
    for(int col=0;col<side;col++)for(int a=0;a<FIELDS;a++){
        int k=col*FIELDS+a;span[a]=fmax(span[a],maximum[k]-minimum[k]);
    }
    if(tracing&&fflush(trace))MPI_Abort(MPI_COMM(simpar),99);
    fprintf(stderr,"[UNIFIED_PLANAR] step=%d t=%.17g stage=%s column_span=(%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e) max_vy=%.9e id_vy=%llu max_ay=%.9e id_ay=%llu\n",
        gfs_rk4_context_step,time,stage,span[0],span[1],span[2],span[3],span[4],span[5],span[6],maxvy,idvy,maxay,iday);
}
