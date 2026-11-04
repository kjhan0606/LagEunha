/* Coefficients are lagged once at the start of each viscosity half-step.
 * Active signs are solved at the endpoint. This viscosity time update is
 * first order, L-stable, and compression-only. Inviscid KDK remains second
 * order. Heat includes the positive BE time-discretization dissipation.
 */
static int gfs_u1_implicit_viscosity(double h)
{
    static double v0[NMAX],dv[NMAX],diag[NMAX],lower[NMAX],rhs[NMAX],heat[NMAX];
    static unsigned char active[NMAX],next[NMAX];
    int wall_active[2],new_wall[2];
    if(!gfs_unified1d_rhs())return 0;
    for(int i=0;i<N;i++){v0[i]=P[i].v;dv[i]=0;heat[i]=0;}
    for(int i=0;i<N-1;i++)active[i]=gfs_u1_D[i]>0;
    wall_active[0]=gfs_u1_Dwall[0]>0;wall_active[1]=gfs_u1_Dwall[1]>0;
    int converged=0;
    for(int iteration=0;iteration<32;iteration++){
        for(int i=0;i<N;i++){
            double dl=i&&active[i-1]?gfs_u1_D[i-1]:0;
            double dr=i<N-1&&active[i]?gfs_u1_D[i]:0;
            double wall=(i==0&&wall_active[0]?gfs_u1_Dwall[0]:0)+
                        (i==N-1&&wall_active[1]?gfs_u1_Dwall[1]:0);
            diag[i]=P[i].m+h*(dl+dr+wall);lower[i]=-h*dl;
            rhs[i]=-h*(dl*(v0[i]-(i?v0[i-1]:0))+
                         dr*(v0[i]-(i<N-1?v0[i+1]:0))+wall*v0[i]);
        }
        for(int i=1;i<N;i++){
            if(!(diag[i-1]>0))return 0;
            double ratio=lower[i]/diag[i-1];
            diag[i]-=ratio*lower[i];rhs[i]-=ratio*rhs[i-1];
        }
        if(!(diag[N-1]>0))return 0;
        dv[N-1]=rhs[N-1]/diag[N-1];
        for(int i=N-2;i>=0;i--){
            if(!(diag[i]>0))return 0;
            double upper=-h*(active[i]?gfs_u1_D[i]:0);
            dv[i]=(rhs[i]-upper*dv[i+1])/diag[i];
        }
        int changed=0;
        for(int i=0;i<N-1;i++){
            next[i]=gfs_u1_D[i]>0 && (v0[i]-v0[i+1])+(dv[i]-dv[i+1])>0;
            if(next[i]!=active[i])changed=1;
        }
        new_wall[0]=gfs_u1_Dwall[0]>0&&v0[0]+dv[0]<0;
        new_wall[1]=gfs_u1_Dwall[1]>0&&v0[N-1]+dv[N-1]>0;
        if(new_wall[0]!=wall_active[0]||new_wall[1]!=wall_active[1])changed=1;
        if(!changed){converged=1;break;}
        memcpy(active,next,N-1);wall_active[0]=new_wall[0];wall_active[1]=new_wall[1];
    }
    if(!converged)return 0;
    for(int i=0;i<N;i++)heat[i]=.5*P[i].m*dv[i]*dv[i];
    for(int i=0;i<N-1;i++)if(active[i]){
        double delta=(v0[i]+dv[i])-(v0[i+1]+dv[i+1]);
        double value=.5*h*gfs_u1_D[i]*delta*delta;
        heat[i]+=value;heat[i+1]+=value;
    }
    if(wall_active[0])heat[0]+=h*gfs_u1_Dwall[0]*pow(v0[0]+dv[0],2);
    if(wall_active[1])heat[N-1]+=h*gfs_u1_Dwall[1]*pow(v0[N-1]+dv[N-1],2);
    for(int i=0;i<N;i++){
        if(!isfinite(dv[i])||!(heat[i]>=0)||!isfinite(heat[i]))return 0;
        P[i].v=v0[i]+dv[i];P[i].e+=heat[i]/P[i].m;
    }
    gfs_u1_rhs_cached=0;return gfs_unified1d_geometry();
}
static int gfs_unified1d_implicit_step(double dt)
{
    static Pt saved[NMAX];static double V0[NMAX],U0[NMAX];
    int n0=N;double work0=bc_eng_net;
    memcpy(saved,P,N*sizeof(Pt));
    if(!gfs_u1_rhs_cached&&!gfs_unified1d_rhs())goto reject;
    for(int i=0;i<N;i++)P[i].v+=.5*dt*gfs_u1_pressure_acc[i];
    gfs_u1_rhs_cached=0;
    if(!gfs_u1_implicit_viscosity(.5*dt))goto reject;
    double power0=gfs_u1_boundary_power();
    for(int i=0;i<N;i++){V0[i]=P[i].vol;U0[i]=P[i].m*P[i].e;P[i].x+=dt*P[i].v;}
    gfs_u1_rhs_cached=0;
    if(!gfs_unified1d_geometry())goto reject;
    for(int i=0;i<N;i++){
        double U;
        if(!gfs_variational_adiabatic(U0[i],V0[i],P[i].vol,GAM,&U))goto reject;
        P[i].e=U/P[i].m;
    }
    if(!gfs_unified1d_geometry())goto reject;
    bc_eng_net+=.5*dt*(power0+gfs_u1_boundary_power());
    if(!gfs_u1_implicit_viscosity(.5*dt)||!gfs_unified1d_rhs())goto reject;
    for(int i=0;i<N;i++)P[i].v+=.5*dt*gfs_u1_pressure_acc[i];
    gfs_u1_rhs_cached=0;
    inflow_particles();
    if(!gfs_unified1d_geometry())goto reject;
    return 0;
reject:
    N=n0;memcpy(P,saved,N*sizeof(Pt));bc_eng_net=work0;
    gfs_u1_rhs_cached=0;(void)gfs_unified1d_geometry();return -1;
}
