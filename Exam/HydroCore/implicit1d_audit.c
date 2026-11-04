/* Direct conservation and time-refinement checks of the actual 1D driver. */
#define main archived_hydro1d_main
#include "../../Hydro1DExam/laguerre_sod.c"
#undef main
static long double conserved_energy(void)
{
    long double E=0;
    for(int i=0;i<N;i++)E+=(long double)P[i].m*(P[i].e+.5L*P[i].v*P[i].v);
    return E;
}
static int algebra(int walls,double boost,double h)
{
    static double v0[NMAX],e0[NMAX],m0[NMAX],x0[NMAX];
    setup_problem(PROB_CONTACT);set_method_gfs();mode=VORONOI;
    unified1d=unified_implicit=1;pair_1d=0;
    bc_left=bc_right=walls?BC_REFLECT:BC_TRANSMIT;
    init(64);
    long double p0=0,p1=0,dE=0,expected_wall_impulse=0;
    for(int i=0;i<N;i++){
        P[i].v=boost+.3*sin(1.7*i)+.2*cos(.63*i);
        v0[i]=P[i].v;e0[i]=P[i].e;m0[i]=P[i].m;x0[i]=P[i].x;
        p0+=(long double)m0[i]*v0[i];
    }
    if(!gfs_u1_implicit_viscosity(h))return 2;
    double max_residual=0,min_heat=1e30,max_heat_error=0,scale=0;
    for(int i=0;i<N;i++){
        double v=P[i].v,dv=v-v0[i],mv=m0[i]*dv;
        double dissip=0,expected_heat=.5*m0[i]*dv*dv;
        if(i&&v<P[i-1].v){
            double g=v-P[i-1].v;
            dissip+=gfs_u1_D[i-1]*g;
            expected_heat+=.5*h*gfs_u1_D[i-1]*g*g;
        }
        if(i<N-1&&v>P[i+1].v){
            double g=v-P[i+1].v;
            dissip+=gfs_u1_D[i]*g;
            expected_heat+=.5*h*gfs_u1_D[i]*g*g;
        }
        if(i==0&&v<0){dissip+=gfs_u1_Dwall[0]*v;expected_heat+=h*gfs_u1_Dwall[0]*v*v;expected_wall_impulse-=h*gfs_u1_Dwall[0]*v;}
        if(i==N-1&&v>0){dissip+=gfs_u1_Dwall[1]*v;expected_heat+=h*gfs_u1_Dwall[1]*v*v;expected_wall_impulse-=h*gfs_u1_Dwall[1]*v;}
        double heat=m0[i]*(P[i].e-e0[i]);
        min_heat=fmin(min_heat,heat);max_heat_error=fmax(max_heat_error,fabs(heat-expected_heat));
        max_residual=fmax(max_residual,fabs(mv+h*dissip));
        scale=fmax(scale,fmax(fabs(mv),fabs(h*dissip)));
        dE+=(long double)m0[i]*(.5L*(v+v0[i])*(v-v0[i])+(P[i].e-e0[i]));
        p1+=(long double)m0[i]*v;
        if(P[i].m!=m0[i]||P[i].x!=x0[i])return 2;
    }
    long double pdefect=p1-p0-expected_wall_impulse;
    double tolerance=4096*DBL_EPSILON*fmax(1.,fabs(boost));
    int ok=fabs(dE)<tolerance&&fabs(pdefect)<tolerance&&min_heat>=-DBL_EPSILON&&max_residual<tolerance&&max_heat_error<tolerance;
    printf("{\"scope\":\"actual fixed-position 1D viscosity half-step\",\"walls\":%d,\"boost\":%.17g,\"h\":%.17g,\"energy_defect\":%.17Lg,\"momentum_defect_with_wall\":%.17Lg,\"max_equation_residual\":%.17g,\"equation_scale\":%.17g,\"min_heat\":%.17g,\"max_local_heat_error\":%.17g,\"pass\":%s}\n",walls,boost,h,dE,pdefect,max_residual,scale,min_heat,max_heat_error,ok?"true":"false");
    return !ok;
}
static int evolve(int pid,int np,double factor,int implicit)
{
    setup_problem(pid);set_method_gfs();mode=VORONOI;
    unified1d=1;unified_implicit=implicit;pair_1d=0;bc_eng_net=0;
    init(np);long double E0=conserved_energy();
    double t=0;int steps=0,retries=0;
    while(t<tend-1e-14){
        double dt=fmin(factor*get_dt(),tend-t);int tries=0;
        while((implicit?gfs_unified1d_implicit_step(dt):gfs_unified1d_step(dt))<0){
            dt*=.5;tries++;retries++;
            if(tries>24||dt<1e-12)return 2;
        }
        if(!(dt>0)||!isfinite(dt)||steps>=1000000)return 2;
        t+=dt;steps++;
    }
    FILE *out=fopen("final.csv","wx");if(!out)return 2;
    fprintf(out,"x,V,rho,v,p,U\n");
    for(int i=0;i<N;i++)fprintf(out,"%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",P[i].x,P[i].vol,P[i].rho,P[i].v,P[i].p,P[i].m*P[i].e);
    if(fclose(out))return 2;
    long double dE=conserved_energy()-E0-bc_eng_net;
    printf("{\"scope\":\"full 1D evolution\",\"problem_id\":%d,\"initial_N\":%d,\"final_N\":%d,\"dt_factor\":%.17g,\"implicit\":%d,\"time\":%.17g,\"target\":%.17g,\"steps\":%d,\"retries\":%d,\"relative_energy_error\":%.17Lg}\n",pid,np,N,factor,implicit,t,tend,steps,retries,fabsl(dE)/fabsl(E0));
    return 0;
}
int main(int argc,char **argv)
{
    if(argc==2&&!strcmp(argv[1],"algebra")){
        int failed=0;
        const double hs[]={1e-4,.01,1.};
        for(int walls=0;walls<2;walls++)for(int b=0;b<(walls?1:2);b++)for(int k=0;k<3;k++)failed|=algebra(walls,b?1000:0,hs[k]);
        return failed;
    }
    if(argc!=6||strcmp(argv[1],"evolve"))return 2;
    return evolve(atoi(argv[2]),atoi(argv[3]),atof(argv[4]),atoi(argv[5]));
}
