#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include "../gfs_viscosity_implicit.h"
#define main archived_hydro1d_main
#include "../../Hydro1DExam/laguerre_sod.c"
#undef main
static int one_dimensional_match(int walls,double h)
{
    setup_problem(PROB_CONTACT);set_method_gfs();mode=VORONOI;
    unified1d=unified_implicit=1;pair_1d=0;
    bc_left=bc_right=walls?BC_REFLECT:BC_TRANSMIT;init(64);
    for(int i=0;i<N;i++)P[i].v=.3*sin(1.7*i)+.2*cos(.63*i);
    if(!gfs_unified1d_rhs())return 2;
    double mass[64],x[64],v[64],U[64];GfsViscosityEdge edges[66];int ne=0;
    for(int i=0;i<N;i++){
        mass[i]=P[i].m;x[i]=P[i].x;v[i]=P[i].v;U[i]=P[i].m*P[i].e;
        if(i<N-1)edges[ne++]=(GfsViscosityEdge){i,i+1,{1,0,0},gfs_u1_D[i]};
    }
    if(walls){
        edges[ne++]=(GfsViscosityEdge){0,-1,{-1,0,0},gfs_u1_Dwall[0]};
        edges[ne++]=(GfsViscosityEdge){N-1,-1,{1,0,0},gfs_u1_Dwall[1]};
    }
    GfsViscosityWorkspace *w=gfs_viscosity_workspace_create();GfsViscosityStats s;
    int rc=gfs_viscosity_be(w,1,N,ne,edges,mass,x,v,U,h,&s);
    gfs_viscosity_workspace_free(w);
    if(rc||!gfs_u1_implicit_viscosity(h))return 2;
    double dv=0,dU=0;
    for(int i=0;i<N;i++){dv=fmax(dv,fabs(v[i]-P[i].v));dU=fmax(dU,fabs(U[i]-P[i].m*P[i].e));}
    int pass=dv<1e-11&&dU<1e-12&&fabs(s.energy_defect)<1e-12;
    printf("{\"scope\":\"common graph versus actual 1D tridiagonal viscosity\",\"dim\":1,\"walls\":%d,\"h\":%.17g,\"iterations\":%d,\"active_iterations\":%d,\"max_velocity_difference\":%.17g,\"max_U_difference\":%.17g,\"energy_defect\":%.17g,\"pass\":%s}\n",walls,h,s.iterations,s.active_iterations,dv,dU,s.energy_defect,pass?"true":"false");
    return !pass;
}
static int multidimensional(int dim,int walls,double boost,double h)
{
    enum {C=64,E=4096};double m[C],U[C],x[3*C],v[3*C];GfsViscosityEdge edges[E];
    int n=dim==2?16:64,ne=0;
    for(int i=0;i<n;i++){
        int ix=i%4,iy=(i/4)%4,iz=i/16;
        m[i]=.2+.01*i;U[i]=.01;
        x[dim*i]=ix*.3;x[dim*i+1]=iy*.3;if(dim==3)x[dim*i+2]=iz*.3;
        for(int a=0;a<dim;a++)v[dim*i+a]=boost+.2*sin(1.7*i+.9*a);
    }
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
        double r[3]={0},d2=0;
        for(int a=0;a<dim;a++){r[a]=x[dim*j+a]-x[dim*i+a];d2+=r[a]*r[a];}
        if(d2>.19)continue;
        double norm=sqrt(d2);
        for(int a=0;a<dim;a++)r[a]/=norm;
        edges[ne]=(GfsViscosityEdge){i,j,{r[0],r[1],r[2]},.2+.01*(ne%13)};ne++;
    }
    if(walls)for(int i=0;i<n;i++)for(int a=0;a<dim;a++)for(int side=0;side<2;side++){
        double coord=x[dim*i+a];if((side&&coord<.89)||(!side&&coord>.01))continue;
        edges[ne]=(GfsViscosityEdge){i,-1,{0,0,0},.5};edges[ne].normal[a]=side?1:-1;ne++;
    }
    GfsViscosityWorkspace *w=gfs_viscosity_workspace_create();GfsViscosityStats s;
    int rc=gfs_viscosity_be(w,dim,n,ne,edges,m,x,v,U,h,&s);
    gfs_viscosity_workspace_free(w);if(rc)return 2;
    double dp=0,dJ=0;
    for(int a=0;a<3;a++){dp=fmax(dp,fabs(s.momentum_defect[a]));dJ=fmax(dJ,fabs(s.angular_defect[a]));}
    double tol=4096*DBL_EPSILON*fmax(1,fabs(boost));
    int pass=fabs(s.energy_defect)<tol&&dp<tol&&dJ<tol&&s.min_heat>=0;
    printf("{\"scope\":\"common graph central viscosity conservation\",\"dim\":%d,\"walls\":%d,\"boost\":%.17g,\"h\":%.17g,\"iterations\":%d,\"active_iterations\":%d,\"energy_defect\":%.17g,\"max_momentum_defect\":%.17g,\"max_angular_defect\":%.17g,\"min_heat\":%.17g,\"equation_residual\":%.17g,\"pass\":%s}\n",dim,walls,boost,h,s.iterations,s.active_iterations,s.energy_defect,dp,dJ,s.min_heat,s.relative_equation_residual,pass?"true":"false");
    return !pass;
}
int main(void)
{
    int failed=0;double hs[]={1e-4,.01,1};
    for(int walls=0;walls<2;walls++)for(int k=0;k<3;k++)failed|=one_dimensional_match(walls,hs[k]);
    for(int dim=2;dim<=3;dim++)for(int walls=0;walls<2;walls++)
        for(int b=0;b<(walls?1:2);b++)for(int k=0;k<3;k++)failed|=multidimensional(dim,walls,b?1000:0,hs[k]);
    return failed;
}
