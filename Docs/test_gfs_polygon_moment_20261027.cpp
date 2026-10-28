#include "gfs_polygon_moment_20261027.h"
#include <cmath>
#include <cstdio>

static int require(bool ok,const char *message) {
    if(!ok) std::fprintf(stderr,"FAIL: %s\n",message);
    return ok ? 0 : 1;
}

static gfs_polygon_moment polygon(const double *x,const double *y,int n,double mass,bool reverse=false) {
    double x1[8],y1[8],x2[8],y2[8];
    for(int k=0;k<n;++k) {
        int a=reverse ? n-1-k : k;
        int b=reverse ? (n-2-k+n)%n : (k+1)%n;
        x1[k]=x[a];y1[k]=y[a];x2[k]=x[b];y2[k]=y[b];
    }
    return gfs_polygon_moments(x1,y1,x2,y2,0,n,mass);
}

int main() {
    int failed=0;
    const double sx[4]={-.5,.5,.5,-.5},sy[4]={-.5,-.5,.5,.5};
    auto s=polygon(sx,sy,4,2.);
    failed+=require(s.status==0 && std::fabs(s.area-1.)<1e-14 &&
                    std::fabs(s.cx)<1e-14 && std::fabs(s.cy)<1e-14 &&
                    std::fabs(s.inertia-1./3.)<1e-14,
                    "unit square area, centroid and inertia");
    const double rx[4]={-.7,1.3,1.3,-.7},ry[4]={-.7,-.7,.3,.3};
    auto r=polygon(rx,ry,4,3.);
    auto reverse=polygon(rx,ry,4,3.,true);
    failed+=require(r.status==0 && reverse.status==0 && std::fabs(r.area-2.)<1e-14 &&
                    std::fabs(r.cx-.3)<1e-14 && std::fabs(r.cy+.2)<1e-14 &&
                    std::fabs(r.inertia-1.25)<1e-14 &&
                    std::fabs(reverse.inertia-r.inertia)<1e-14,
                    "offset rectangle and reversed edge orientation");
    const double tx[3]={0.,2.,0.},ty[3]={0.,0.,1.};
    auto t=polygon(tx,ty,3,3.);
    failed+=require(t.status==0 && std::fabs(t.area-1.)<1e-14 &&
                    std::fabs(t.cx-2./3.)<1e-14 && std::fabs(t.cy-1./3.)<1e-14 &&
                    std::fabs(t.inertia-5./6.)<1e-14,
                    "triangle analytical centroid and inertia");
    const double farx[4]={99.5,100.5,100.5,99.5};
    const double fary[4]={-70.5,-70.5,-69.5,-69.5};
    auto far=polygon(farx,fary,4,2.);
    failed+=require(far.status==0 && std::fabs(far.cx-100.)<1e-10 &&
                    std::fabs(far.cy+70.)<1e-10 && std::fabs(far.inertia-1./3.)<1e-10,
                    "inertia remains centered when the generator is far from the cell");
    const double broken_x1[3]={0.,1.,0.},broken_y1[3]={0.,0.,1.};
    const double broken_x2[3]={1.,0.,.1},broken_y2[3]={0.,1.,0.};
    auto broken=gfs_polygon_moments(broken_x1,broken_y1,broken_x2,broken_y2,0,3,1.);
    failed+=require(broken.status==-2,"open polygon is rejected");
    auto negative=polygon(sx,sy,4,-1.);
    failed+=require(negative.status!=0,"nonpositive mass is rejected");
    if(failed) return 1;
    std::printf("PASS: polygon area, mass centroid, inertia, orientation and invalid geometry\n");
    return 0;
}
