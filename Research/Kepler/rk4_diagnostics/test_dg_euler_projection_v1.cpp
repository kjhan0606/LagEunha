#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "gfs_dg_euler.h"
static void need(bool ok,const char *s){if(!ok){fprintf(stderr,"FAIL %s\n",s);exit(1);}}
static void gauss(int n,double *x,double *w){
    for(int i=0;i<(n+1)/2;i++){
        long double z=cosl(acosl(-1.L)*(i+.75L)/(n+.5L)),d=0;
        for(int it=0;it<80;it++){
            long double a=1,b=z;for(int k=2;k<=n;k++){long double c=((2*k-1)*z*b-(k-1)*a)/k;a=b;b=c;}
            d=n*(z*b-a)/(z*z-1);long double t=z-b/d;
            if(fabsl(t-z)<4*LDBL_EPSILON){z=t;break;}z=t;
        }
        long double a=1,b=z;for(int k=2;k<=n;k++){long double c=((2*k-1)*z*b-(k-1)*a)/k;a=b;b=c;}
        d=n*(z*b-a)/(z*z-1);double ww=1/((1-z*z)*d*d);
        x[i]=(1-z)*.5L;x[n-1-i]=(1+z)*.5L;w[i]=w[n-1-i]=ww;
    }
}
int main(){
    double x[6]={-.8,.45,1,.7,-.4,-.9},y[6]={-.4,-.7,-.1,.7,.8,.2};
    GfsMomentPolygon polygon={x,y,6,0};
    double max_pressure=0,max_mean=0,max_gcl=0,max_gauss=0;
    for(int order=2;order<=6;order+=2)for(int variant=0;variant<2;variant++){
        GfsDgConfig cfg={};cfg.order=order;cfg.ng[0]=order+2;cfg.ng[1]=order+4;cfg.ng[2]=order+3;
        for(int r=0;r<3;r++){
            gauss(cfg.ng[r],cfg.node[r],cfg.weight[r]);
            for(int power=0;power<2*cfg.ng[r];power++){
                long double v=0;for(int k=0;k<cfg.ng[r];k++)v+=(long double)cfg.weight[r][k]*powl(cfg.node[r][k],power);
                max_gauss=fmax(max_gauss,fabsl(v-1.L/(power+1)));
            }
        }
        cfg.ic.kind=1;cfg.ic.rho=1.2;cfg.ic.pressure=1;cfg.ic.gamma=5./3.;cfg.ic.vx=.2;cfg.ic.vy=-.3;
        if(variant==0){cfg.ic.B[0]=.1;cfg.ic.B[1]=-.7;cfg.ic.B[2]=.9;cfg.ic.B[3]=-.03;}
        else {cfg.density_slope[0]=.2;cfg.density_slope[1]=-.1;}
        GfsDgProjectionAudit audit={};GfsDgPolynomial coarse={};
        need(gfs_dg_project(polygon,&cfg,0,0,&audit,&coarse),"Euler polynomial projection");
        double integrated[4]={},qsum=0;
        for(int k=0;k<polygon.n*cfg.ng[2]*cfg.ng[2];k++){
            double xx,yy,w,u[4];gfs_dg_quad(polygon,&audit.basis,&cfg,2,k,&xx,&yy,&w);
            gfs_dg_evaluate(&audit.basis,&audit.polynomial,order,xx,yy,u);
            max_pressure=fmax(max_pressure,fabs(gfs_dg_pressure(u,cfg.ic.gamma)-1));qsum+=w;
            for(int c=0;c<4;c++)integrated[c]+=w*u[c];
        }
        for(int c=0;c<4;c++)max_mean=fmax(max_mean,fabs(integrated[c]-audit.integral[c]));
        need(fabs(qsum-audit.basis.area)<1.e-13,"independent polygon quadrature");
        /* Uniform Euler state, affine mesh velocity: weak form agrees with
         * Reynolds transport for every modal moment, without evolving it. */
        double rhs[4][GFS_DG_MAX_BASIS]={},reference[4][GFS_DG_MAX_BASIS]={};
        double state[4]={1.2,.24,-.36,1.5+.5*1.2*(.2*.2+.3*.3)};
        for(int k=0;k<polygon.n*cfg.ng[2]*cfg.ng[2];k++){
            double xx,yy,w,phi[GFS_DG_MAX_BASIS],dx[GFS_DG_MAX_BASIS],dy[GFS_DG_MAX_BASIS];
            gfs_dg_quad(polygon,&audit.basis,&cfg,2,k,&xx,&yy,&w);
            need(gfs_dg_volume_point(&audit.basis,order,state,cfg.ic.gamma,xx,yy,0,0,w,rhs),"Euler volume point");
            gfs_dg_monomial(&audit.basis,order,xx,yy,phi,dx,dy);
            double wx=.1+.03*xx-.2*yy,wy=-.2+.4*xx-.07*yy;
            for(int j=0;j<gfs_dg_nb(order);j++)for(int c=0;c<4;c++)
                reference[c][j]+=w*state[c]*(wx*dx[j]+wy*dy[j]-.04*(j?phi[j]-audit.basis.mean[j]:1));
        }
        for(int edge=0;edge<polygon.n;edge++){
            int next=(edge+1)%polygon.n;double sx=y[next]-y[edge],sy=x[edge]-x[next];
            double length=hypot(sx,sy),nx=sx/length,ny=sy/length;
            for(int q=0;q<cfg.ng[2];q++){
                double t=cfg.node[2][q],xx=(1-t)*x[edge]+t*x[next],yy=(1-t)*y[edge]+t*y[next];
                double wx=.1+.03*xx-.2*yy,wy=-.2+.4*xx-.07*yy,fx[4],fy[4],f[4];
                need(gfs_dg_euler_flux(state,cfg.ic.gamma,fx,fy),"Euler face state");
                for(int c=0;c<4;c++)f[c]=fx[c]*nx+fy[c]*ny-state[c]*(wx*nx+wy*ny);
                gfs_dg_face_point(&audit.basis,order,xx,yy,f,length*cfg.weight[2][q],rhs);
            }
        }
        for(int c=0;c<4;c++)for(int j=0;j<gfs_dg_nb(order);j++)max_gcl=fmax(max_gcl,fabs(rhs[c][j]-reference[c][j]));
    }
    need(max_gauss<1.e-15,"Gauss polynomial moments");need(max_pressure<5.e-11,"constant-pressure recovery");
    need(max_mean<5.e-13,"conserved cell integrals");need(max_gcl<2.e-10,"affine mesh Reynolds transport");
    printf("{\"status\":\"pass\",\"max_order\":6,\"gauss_moment_error\":%.17g,\"pressure_error\":%.17g,\"conserved_integral_error\":%.17g,\"fixed_basis_ALE_weak_form_error\":%.17g,\"scope\":\"polynomial projection and instantaneous Euler weak form; no moving-basis timestep or nonlinear stability claim\"}\n",max_gauss,max_pressure,max_mean,max_gcl);
}
