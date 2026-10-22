#include <cstdio>
#include <cfloat>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include "gfs_ale_fan.h"
static void need(bool b,const char *s){if(!b){fprintf(stderr,"FAIL %s\n",s);exit(1);}}
static double scale4(const double *u,const double *v){double a=1;for(int c=0;c<4;c++)a=fmax(a,fmax(fabs(u[c]),fabs(v[c])));return a;}
static std::vector<GfsWeightPoint> points;
static const int N=8;
static void evaluate(const double z[N],const double mu[N],double R[N],double J[N][N]){
    for(int i=0;i<N;i++){R[i]=mu[i]*z[i];for(int j=0;j<N;j++)J[i][j]=i==j?mu[i]:0;}
    for(auto &p:points){
        auto f=gfs_weight_flux(&p,z[p.i],z[p.j]);need(f.valid,"graph fan sample");
        double v=p.weight*f.flux[0],a=gfs_weight_jacobian(&p,z[p.i],z[p.j]);need(a>0,"positive face derivative");
        R[p.i]-=v;R[p.j]+=v;J[p.i][p.i]+=a;J[p.j][p.j]+=a;J[p.i][p.j]-=a;J[p.j][p.i]-=a;
    }
}
static double norm(const double *x){double a=0;for(int i=0;i<N;i++)a+=x[i]*x[i];return sqrt(a);}
static void solve(double A[N][N],const double rhs[N],double x[N]){
    for(int i=0;i<N;i++)for(int j=0;j<=i;j++){
        double v=A[i][j];for(int k=0;k<j;k++)v=fma(-A[i][k],A[j][k],v);
        if(i==j){need(v>0,"SPD Cholesky pivot");A[i][i]=sqrt(v);}else A[i][j]=v/A[j][j];
    }
    for(int i=0;i<N;i++){x[i]=rhs[i];for(int j=0;j<i;j++)x[i]-=A[i][j]*x[j];x[i]/=A[i][i];}
    for(int i=N-1;i>=0;i--){for(int j=i+1;j<N;j++)x[i]-=A[j][i]*x[j];x[i]/=A[i][i];}
}
int main(){
    const double gamma=5./3.;
    const GfsAlePrimitive states[][2]={
        {{1,1,.2,-.3},{1,1,.2,-.3}},{{1,1,0,0},{.125,.1,0,0}},
        {{.001,1.e-6,.8,-.3},{.0011,1.1e-6,.79,-.31}},
        {{1,.4,-10,0},{1,.4,10,0}},{{.2,.1,2,-3},{2,3,-1,1}},
        {{1,1.e-6,3.2,-1.9},{.8,2.e-6,3.21,-1.91}}};
    double max_flux=0,max_deriv=0,max_continuity=0,max_boost=0,max_monotone=0;unsigned modes=0;int samples=0,fallback=0;
    for(auto &lr:states)for(int angle=0;angle<5;angle++){
        double theta=.37*angle,nx=cos(theta),ny=sin(theta);auto fan=gfs_ale_fan_build(lr[0],lr[1],nx,ny,gamma);
        need(fan.valid,"full lab fan");if(fan.hll)fallback++;
        double span=fan.speed[2]-fan.speed[0],lo=fan.speed[0]-span,hi=fan.speed[2]+span,prev=DBL_MAX;
        for(int k=0;k<=500;k++){
            double w=lo+(hi-lo)*k/500;auto a=gfs_ale_fan_sample(&fan,w);need(a.valid,"full sample");modes|=1u<<a.mode;samples++;
            max_monotone=fmax(max_monotone,(a.flux[0]-prev)/(1+fabs(a.flux[0])));prev=a.flux[0];
            auto ref=gfs_ale_flux(lr[0],lr[1],nx,ny,w*nx,w*ny,gamma);need(ref.valid,"independent translated fan flux");
            if(a.mode==ref.mode)for(int c=0;c<4;c++)max_flux=fmax(max_flux,fabs(a.flux[c]-ref.flux[c])/(1+fabs(ref.flux[c])));
            double h=1.e-5*span;auto am=gfs_ale_fan_sample(&fan,w-h),ap=gfs_ale_fan_sample(&fan,w+h);
            if(am.region==a.region && ap.region==a.region)
                max_deriv=fmax(max_deriv,fabs((ap.flux[0]-am.flux[0])/(2*h)+a.density)/(1+a.density));
            GfsAlePrimitive bl=lr[0],br=lr[1];double bx=2.1,by=-1.3;
            bl.vx+=bx;bl.vy+=by;br.vx+=bx;br.vy+=by;
            auto bf=gfs_ale_fan_build(bl,br,nx,ny,gamma);need(bf.valid,"boosted fan");
            auto b=gfs_ale_fan_sample(&bf,w+bx*nx+by*ny);double t[4]={a.flux[0],a.flux[1]+bx*a.flux[0],a.flux[2]+by*a.flux[0],
                a.flux[3]+bx*a.flux[1]+by*a.flux[2]+.5*(bx*bx+by*by)*a.flux[0]};
            for(int c=0;c<4;c++)max_boost=fmax(max_boost,fabs(b.flux[c]-t[c])/(1+fabs(t[c])));
        }
        for(int k=0;k<3;k++){
            if(fan.hll && k==1)continue;double w=fan.speed[k];
            int left=k==0?0:(fan.hll?1:k),right=k==2?3:(fan.hll?1:k+1);
            for(int c=0;c<4;c++){
                double a=fma(-w,fan.u[left][c],fan.f[left][c]),b=fma(-w,fan.u[right][c],fan.f[right][c]);
                max_continuity=fmax(max_continuity,fabs(a-b)/(1+fabs(a)+fabs(b)));
            }
        }
    }
    double mu[N],z[N],R[N],J[N][N];for(int i=0;i<N;i++){mu[i]=.5+.1*i;z[i]=7*sin(1.7*i);}
    for(int i=0;i<N;i++)for(int j=i+1;j<N;j++)if(j==i+1 || (i+j)%3==0){
        int which=(i+j)%6;double theta=.17*(i+2*j);GfsWeightPoint p={};p.i=i;p.j=j;p.weight=.2+.03*j;p.inv2d=.6+.04*i;p.w0=.3*(i-j);
        p.fan=gfs_ale_fan_build(states[which][0],states[which][1],cos(theta),sin(theta),gamma);need(p.fan.valid,"graph precompute");points.push_back(p);
    }
    evaluate(z,mu,R,J);double max_fd=0,max_symmetry=0,min_rayleigh=DBL_MAX,gauge_identity=0;
    for(int j=0;j<N;j++){
        double zp[N],zm[N],rp[N],rm[N],tmp[N][N];for(int i=0;i<N;i++)zp[i]=zm[i]=z[i];double h=1.e-6;zp[j]+=h;zm[j]-=h;
        evaluate(zp,mu,rp,tmp);evaluate(zm,mu,rm,tmp);
        for(int i=0;i<N;i++){max_fd=fmax(max_fd,fabs((rp[i]-rm[i])/(2*h)-J[i][j])/(1+fabs(J[i][j])));max_symmetry=fmax(max_symmetry,fabs(J[i][j]-J[j][i]));}
    }
    for(int trial=0;trial<100;trial++){
        double x[N],lhs=0,bound=0;for(int i=0;i<N;i++)x[i]=sin(.31*(i+1)*(trial+1));
        for(int i=0;i<N;i++){bound+=mu[i]*x[i]*x[i];for(int j=0;j<N;j++)lhs+=x[i]*J[i][j]*x[j];}
        min_rayleigh=fmin(min_rayleigh,lhs/bound);need(lhs>=bound*(1-1.e-13),"positive graph Rayleigh lower bound");
    }
    int iterations=0,backtracks=0;
    for(;iterations<40;iterations++){
        evaluate(z,mu,R,J);double rn=norm(R);if(rn<1.e-12)break;
        double direction[N],minus[N];for(int i=0;i<N;i++)minus[i]=-R[i];solve(J,minus,direction);
        double factor=1;bool accepted=false;
        for(int retry=0;retry<40;retry++){
            double trial[N],res[N],mat[N][N];for(int i=0;i<N;i++)trial[i]=z[i]+factor*direction[i];evaluate(trial,mu,res,mat);
            if(norm(res)<rn*(1-1.e-4*factor)){for(int i=0;i<N;i++)z[i]=trial[i];accepted=true;break;}
            factor*=.5;backtracks++;
        }need(accepted,"nonlinear line search");
    }
    evaluate(z,mu,R,J);double sumR=0,sumMuZ=0;for(int i=0;i<N;i++){sumR+=R[i];sumMuZ+=mu[i]*z[i];}
    gauge_identity=fabs(sumR-sumMuZ);
    need(norm(R)<1.e-12 && fabs(sumMuZ)<1.e-12,"coupled mass/weight root and weighted gauge");
    need(max_flux<2.e-8 && max_deriv<1.e-8 && max_boost<2.e-8 && max_continuity<1.e-12 && max_monotone<1.e-13,"fan physical contracts");
    need(max_fd<1.e-7 && max_symmetry==0 && (modes&(1u<<5)),"Jacobian including HLL fallback");
    printf("{\"scope\":\"native C++/CUDA shared lab fan and bounded nonlinear graph contract; no evolution\",\"samples\":%d,\"mode_bitmask\":%u,\"HLL_fans\":%d,\"max_flux_discrepancy\":%.17g,\"max_mass_derivative_error\":%.17g,\"max_fan_continuity_error\":%.17g,\"max_boost_error\":%.17g,\"max_monotonicity_violation\":%.17g,\"graph_Jacobian_FD_error\":%.17g,\"graph_symmetry_error\":%.17g,\"min_Rayleigh_over_mu_bound\":%.17g,\"Newton_iterations\":%d,\"backtracks\":%d,\"final_residual\":%.17g,\"weighted_gauge_rate\":%.17g,\"canonical_sum_identity_error\":%.17g}\n",samples,modes,fallback,max_flux,max_deriv,max_continuity,max_boost,max_monotone,max_fd,max_symmetry,min_rayleigh,iterations,backtracks,norm(R),sumMuZ,gauge_identity);
}
