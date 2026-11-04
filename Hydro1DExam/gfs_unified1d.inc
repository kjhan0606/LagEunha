/*
 * First dimensional integration of the shared interior operator.
 * Fixed masses, no steering. Initial conserved quantities are integrated
 * over the actual uniform Voronoi cells, not sampled then reinterpreted.
 * KDK kicks use midpoint RK2 for pressure plus candidate dissipative force.
 * The drift uses the exact adiabatic volume map. Not yet shock-qualified.
 */
#include "../Exam/gfs_velocity_fit.h"

static double gfs_u1_heat[NMAX],gfs_u1_slope[NMAX],gfs_u1_rate[NMAX];
static double gfs_u1_D[NMAX],gfs_u1_Dwall[2],gfs_u1_pressure_acc[NMAX];
static int gfs_u1_rhs_cached=0;
static double gfs_unified1d_visc_dt(void)
{
    double maxrate=0;
    for(int i=0;i<N;i++)maxrate=fmax(maxrate,gfs_u1_rate[i]);
    return maxrate>0?CFL/maxrate:1e30;
}
static int gfs_unified1d_geometry(void)
{
    if(N<2)return 0;
    for(int i=1;i<N;i++){
        if(!(P[i].x>P[i-1].x))return 0;
        face[i]=.5*(P[i-1].x+P[i].x);
    }
    face[0]=bc_left==BC_REFLECT?x0_dom:1.5*P[0].x-.5*P[1].x;
    face[N]=bc_right==BC_REFLECT?x1_dom:1.5*P[N-1].x-.5*P[N-2].x;
    for(int i=0;i<N;i++){
        double V=face[i+1]-face[i];
        if(!(V>0)||!(P[i].m>0)||!(P[i].e>0)||!isfinite(V)||!isfinite(P[i].e))return 0;
        P[i].vol=V;P[i].rho=P[i].m/V;P[i].p=GM1*P[i].rho*P[i].e;
        P[i].c=sqrt(GAM*GM1*P[i].e);P[i].w=0;
        if(!isfinite(P[i].p)||!isfinite(P[i].c))return 0;
    }
    return 1;
}

static int gfs_unified1d_rhs(void)
{
    gfs_u1_rhs_cached=0;
    if(!gfs_unified1d_geometry())return 0;
    gfs_u1_Dwall[0]=gfs_u1_Dwall[1]=0;
    for(int i=0;i<N;i++){
        GfsVelocityFit fit;double G[3][3];
        if(!gfs_velocity_fit_init(&fit,1))return 0;
        if(i>0){
            double r[3]={P[i-1].x-P[i].x,0,0},dv[3]={P[i-1].v-P[i].v,0,0};
            if(!gfs_velocity_fit_add_face(&fit,r,dv,1))return 0;
        }else if(bc_left==BC_REFLECT){
            double r[3]={2*(x0_dom-P[i].x),0,0},dv[3]={-2*P[i].v,0,0};
            if(!gfs_velocity_fit_add_face(&fit,r,dv,1))return 0;
        }
        if(i<N-1){
            double r[3]={P[i+1].x-P[i].x,0,0},dv[3]={P[i+1].v-P[i].v,0,0};
            if(!gfs_velocity_fit_add_face(&fit,r,dv,1))return 0;
        }else if(bc_right==BC_REFLECT){
            double r[3]={2*(x1_dom-P[i].x),0,0},dv[3]={-2*P[i].v,0,0};
            if(!gfs_velocity_fit_add_face(&fit,r,dv,1))return 0;
        }
        if(!gfs_velocity_fit_solve(&fit,G,NULL))return 0;
        gfs_u1_slope[i]=G[0][0];
        P[i].av=0;gfs_u1_heat[i]=0;gfs_u1_rate[i]=0;
        gfs_u1_pressure_acc[i]=0;gfs_u1_D[i]=0;
    }
    vf[0]=bc_left==BC_REFLECT?0:1.5*P[0].v-.5*P[1].v;
    vf[N]=bc_right==BC_REFLECT?0:1.5*P[N-1].v-.5*P[N-2].v;
    pf[0]=P[0].p;pf[N]=P[N-1].p;
    for(int i=0;i<N-1;i++){
        double r[3]={P[i+1].x-P[i].x,0,0},c[3]={.5*r[0],0,0};
        double vi[3]={P[i].v,0,0},vj[3]={P[i+1].v,0,0};
        double jump[3]={P[i+1].v-.5*r[0]*gfs_u1_slope[i+1]
                      -P[i].v-.5*r[0]*gfs_u1_slope[i],0,0};
        GfsVariationalFace b;double force[3];
        if(!gfs_variational_face(1,r,c,1,&b))return 0;
        gfs_variational_pressure(&b,P[i].p,P[i+1].p,force);
        double q=gfs_variational_pair_pressure(&b,vi,vj,jump,
                P[i].rho,P[i].c,P[i+1].rho,P[i+1].c,
                gfs_variational_compression_fraction(gfs_u1_slope[i],0,gfs_u1_slope[i+1],0),
                P[i].vol,P[i+1].vol);
        P[i].av+=(force[0]-q)/P[i].m;
        P[i+1].av+=(force[0]+q)/P[i+1].m;
        gfs_u1_pressure_acc[i]+=force[0]/P[i].m;
        gfs_u1_pressure_acc[i+1]+=force[0]/P[i+1].m;
        double heat=.5*q*(P[i].v-P[i+1].v);
        gfs_u1_heat[i]+=heat/P[i].m;gfs_u1_heat[i+1]+=heat/P[i+1].m;
        if(q>0&&P[i].v>P[i+1].v){
            double damping=2*q/(P[i].v-P[i+1].v);
            gfs_u1_D[i]=.5*damping;
            gfs_u1_rate[i]+=damping/P[i].m;gfs_u1_rate[i+1]+=damping/P[i+1].m;
        }
        vf[i+1]=.5*(P[i].v+P[i+1].v);
    }
    for(int side=0;side<2;side++){
        int i=side?N-1:0,reflect=side?bc_right==BC_REFLECT:bc_left==BC_REFLECT;
        if(reflect){
            double wall=side?x1_dom:x0_dom;
            double r[3]={2*(wall-P[i].x),0,0},c[3]={wall-P[i].x,0,0};
            double vi[3]={P[i].v,0,0},vj[3]={-P[i].v,0,0};
            double at_face=P[i].v+gfs_u1_slope[i]*(wall-P[i].x);
            double jump[3]={-2*at_face,0,0};GfsVariationalFace b;
            if(!gfs_variational_face(1,r,c,1,&b))return 0;
            double q=gfs_variational_pair_pressure(&b,vi,vj,jump,
                    P[i].rho,P[i].c,P[i].rho,P[i].c,
                    gfs_variational_compression_fraction(gfs_u1_slope[i],0,gfs_u1_slope[i],0),
                    P[i].vol,P[i].vol);
            P[i].av-=q*b.normal[0]/P[i].m;
            gfs_u1_heat[i]+=q*P[i].v*b.normal[0]/P[i].m;
            if(q>0&&P[i].v*b.normal[0]>0)
            {
                gfs_u1_Dwall[side]=q/(P[i].v*b.normal[0]);
                gfs_u1_rate[i]+=2*gfs_u1_Dwall[side]/P[i].m;
            }
            pf[side?N:0]+=q;
        }
    }
    for(int i=0;i<N;i++){
        double rate=vf[i+1]-vf[i];
        P[i].ae=-P[i].p*rate/P[i].m+gfs_u1_heat[i];
        divv[i]=rate/P[i].vol;
    }
    gfs_u1_rhs_cached=1;
    return 1;
}

static void gfs_u1_integrate(double a,double b,double rho,double v,double p,
        int sine,double *mass,double *mom,double *energy)
{
    if(!(b>a))return;
    double m=rho*(b-a);
    if(sine)m+=.04*(cos(5*a)-cos(5*b));
    *mass+=m;*mom+=m*v;*energy+=p*(b-a)/GM1+.5*m*v*v;
}
static void gfs_unified1d_init(int n)
{
    gfs_u1_rhs_cached=0;
    if(n<4||n>NMAX-2){fprintf(stderr,"unified 1D invalid N\n");exit(2);}
    N=n;meand=(x1_dom-x0_dom)/n;NL_actual=NR_actual=0;
    for(int i=0;i<n;i++){
        memset(P+i,0,sizeof(Pt));double a=x0_dom+i*meand,b=a+meand,m=0,mom=0,E=0;
        P[i].x=.5*(a+b);
        gfs_u1_integrate(a,fmin(b,xmid1),rhoL_ic,vL_ic,pL_ic,0,&m,&mom,&E);
        if(n_regions==3){
            gfs_u1_integrate(fmax(a,xmid1),fmin(b,xmid2),rhoM_ic,vM_ic,pM_ic,0,&m,&mom,&E);
            gfs_u1_integrate(fmax(a,xmid2),b,rhoR_ic,vR_ic,pR_ic,0,&m,&mom,&E);
        }else{
            gfs_u1_integrate(fmax(a,xmid1),b,rhoR_ic,vR_ic,pR_ic,
                problem_id==PROB_SHUOSHER,&m,&mom,&E);
        }
        P[i].m=m;P[i].v=mom/m;P[i].e=(E-.5*m*P[i].v*P[i].v)/m;
        if(P[i].x<xmid1)NL_actual++;else NR_actual++;
    }
    memset(divv,0,sizeof(divv));memset(divv_old,0,sizeof(divv_old));
    memset(alpha_cd,0,sizeof(alpha_cd));memset(vsig_max,0,sizeof(vsig_max));
    if(!gfs_unified1d_geometry()){fprintf(stderr,"unified 1D invalid initial state\n");exit(2);}
}

static int gfs_u1_kick(double h)
{
    static double v0[NMAX],e0[NMAX];
    if(!gfs_u1_rhs_cached&&!gfs_unified1d_rhs())return 0;
    for(int i=0;i<N;i++){
        v0[i]=P[i].v;e0[i]=P[i].e;
        P[i].v=v0[i]+.5*h*P[i].av;P[i].e=e0[i]+.5*h*gfs_u1_heat[i];
    }
    gfs_u1_rhs_cached=0;
    if(!gfs_unified1d_rhs())return 0;
    for(int i=0;i<N;i++){
        P[i].v=v0[i]+h*P[i].av;P[i].e=e0[i]+h*gfs_u1_heat[i];
    }
    gfs_u1_rhs_cached=0;
    return gfs_unified1d_geometry();
}
static double gfs_u1_boundary_power(void)
{
    const double vl=bc_left==BC_REFLECT?0:1.5*P[0].v-.5*P[1].v;
    const double vr=bc_right==BC_REFLECT?0:1.5*P[N-1].v-.5*P[N-2].v;
    return P[0].p*vl-P[N-1].p*vr;
}
static int gfs_unified1d_step(double dt)
{
    static Pt saved[NMAX];static double oldV[NMAX],oldU[NMAX];
    const int n0=N;const double work0=bc_eng_net;
    memcpy(saved,P,N*sizeof(Pt));
    if(!gfs_u1_kick(.5*dt))goto reject;
    double power0=gfs_u1_boundary_power();
    for(int i=0;i<N;i++){oldV[i]=P[i].vol;oldU[i]=P[i].m*P[i].e;P[i].x+=dt*P[i].v;}
    gfs_u1_rhs_cached=0;
    if(!gfs_unified1d_geometry())goto reject;
    for(int i=0;i<N;i++){
        double U;
        if(!gfs_variational_adiabatic(oldU[i],oldV[i],P[i].vol,GAM,&U))goto reject;
        P[i].e=U/P[i].m;
    }
    if(!gfs_unified1d_geometry())goto reject;
    bc_eng_net+=.5*dt*(power0+gfs_u1_boundary_power());
    if(!gfs_u1_kick(.5*dt))goto reject;
    inflow_particles();
    gfs_u1_rhs_cached=0;
    if(!gfs_unified1d_geometry())goto reject;
    return 0;
reject:
    gfs_u1_rhs_cached=0;
    N=n0;memcpy(P,saved,N*sizeof(Pt));bc_eng_net=work0;
    (void)gfs_unified1d_geometry();
    return -1;
}
