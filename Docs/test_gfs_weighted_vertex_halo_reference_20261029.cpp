#include "gfs_weighted_vertex_halo_reference_20261029.h"
#include <algorithm>
#include <cassert>
#include <cstdio>
#include <limits>
#include <random>

static bool cutter(long double px,long double py,long double ux,long double uy,
                   long double wp,long double qx,long double qy,long double wq,
                   bool clamped) {
    long double rx=qx-px,ry=qy-py,d2=rx*rx+ry*ry;
    assert(d2>0);
    long double rhs=.5L*(d2+wp-wq);
    if (clamped) rhs=std::max((long double).05*d2,
                              std::min((long double).95*d2,rhs));
    return ux*rx+uy*ry>rhs;
}

int main() {
    auto safe=gfs_weighted_vertex_halo_disk_inside_reference;
    // delta=9 and |u|^2=16 give an exactly tangent radius of five.
    assert(!safe(0,0,4,0,0,9,-1,-5,9,5));
    assert(safe(0,0,4,0,0,9,-1.001,-5.001,9.001,5.001));
    assert(safe(0,0,4,0,-4,5,-1.001,-5.001,9.001,5.001));
    assert(!safe(0,0,0,0,1,0,-1,-1,1,1));
    assert(!safe(0,0,0,0,NAN,0,-1,-1,1,1));
    assert(!safe(0,0,0,0,0,INFINITY,-1,-1,1,1));
    assert(!safe(0,0,0,0,-1e308,1e308,-1,-1,1,1));
    assert(!safe(0,0,1e200,0,0,1,-1e201,-1e201,1e201,1e201));
    assert(!safe(0,0,0,0,0,1,1,-1,-1,1));
    assert(!safe(NAN,0,0,0,0,1,-1,-1,1,1));
    assert(safe(0,0,1e-200,0,0,1e-320,-1e-150,-1e-150,1e-150,1e-150));
    assert(safe(1e15,-1e15,4,0,0,9,1e15-2,-1e15-6,1e15+10,-1e15+6));
    assert(!safe(1e15,-1e15,4,0,0,9,1e15-1,-1e15-5,1e15+9,-1e15+5));
    assert(!safe(0,0,4,0,0,9,-1,-5,9,std::nextafter(5.,0.)));

    // Explicit unclamped, lower/upper-clamped and transition-adjacent
    // planes. q=(2,0) lies outside every certified rectangle below.
    unsigned clamp_cases=0;
    for(double wq:{0.,1.,3.8,-4.,3.6,std::nextafter(3.6,0.),
                   std::nextafter(3.6,INFINITY),-3.6,
                   std::nextafter(-3.6,0.),std::nextafter(-3.6,-INFINITY)}) {
        double wm=std::max(0.,wq),pad=sqrt(.01+wm)+.001;
        assert(pad<2.);
        assert(safe(0,0,0,.1,0,wm,-pad,.1-pad,pad,.1+pad));
        assert(!cutter(0,0,0,.1,0,2,0,wq,false));
        assert(!cutter(0,0,0,.1,0,2,0,wq,true));
        clamp_cases++;
    }

    // A stale global upper bound can pass although a missing site cuts.
    // Updating Wmax to the omitted site's current weight rejects the halo.
    assert(safe(0,0,.1,0,0,0,-.5,-.5,.5,.5));
    assert(cutter(0,0,.1,0,0,.6,0,.36,false));
    assert(cutter(0,0,.1,0,0,.6,0,.36,true));
    assert(!safe(0,0,.1,0,0,.36,-.5,-.5,.5,.5));
    // A whole finite triangle, not just one isolated vertex, exposes
    // the same missing-current-global-bound error without a clamp.
    for(double vy:{-.5,.5})
        assert(safe(0,0,0,vy,0,0,-3,-3,2.1,3));
    assert(safe(0,0,1,0,0,0,-3,-3,2.1,3));
    assert(cutter(0,0,1,0,0,3,0,4,false));
    assert(cutter(0,0,1,0,0,3,0,4,true));
    assert(!safe(0,0,1,0,0,4,-3,-3,2.1,3));
    // An extremely low omitted weight activates the upper clamp. It is
    // still harmless when the current uniform-weight bound certifies D.
    assert(!cutter(0,0,.1,0,0,.6,0,-100,true));

    std::mt19937_64 rng(2026102903);
    std::uniform_real_distribution<double> u(-1.,1.),z(0.,1.);
    unsigned zero_reductions=0,accepted=0,rejected=0,outside_planes=0;
    for (unsigned k=0;k<100000;k++) {
        double px=1e6*u(rng),py=1e6*u(rng),dx=10*u(rng),dy=10*u(rng);
        double pad=hypot(dx,dy)*(1+1e-7*u(rng));
        double xmin=(px+dx)-pad,xmax=(px+dx)+pad;
        double ymin=(py+dy)-pad,ymax=(py+dy)+pad;
        int ordinary=gfs_vertex_halo_disk_inside(px,py,dx,dy,xmin,ymin,xmax,ymax);
        assert(safe(px,py,dx,dy,0,0,xmin,ymin,xmax,ymax)==ordinary);
        assert(safe(px,py,dx,dy,-4,-4,xmin,ymin,xmax,ymax)==ordinary);
        zero_reductions++;
    }
    for (unsigned k=0;k<20000;k++) {
        double px=1000*u(rng),py=1000*u(rng),dx=2*u(rng),dy=2*u(rng);
        double wp=10*u(rng),wmax=wp+4*z(rng);
        double r=sqrt(dx*dx+dy*dy+(wmax-wp));
        double pad=r*(1+1e-7*u(rng));
        double xmin=(px+dx)-pad,xmax=(px+dx)+pad;
        double ymin=(py+dy)-pad,ymax=(py+dy)+pad;
        if (!safe(px,py,dx,dy,wp,wmax,xmin,ymin,xmax,ymax)) {rejected++;continue;}
        accepted++;
        long double vx=(long double)px+dx,vy=(long double)py+dy;
        long double r2=(long double)dx*dx+(long double)dy*dy+(long double)wmax-wp;
        for(long double gap:{vx-xmin,(long double)xmax-vx,vy-ymin,(long double)ymax-vy})
            assert(gap>0 && gap*gap>r2);
        // Independent oracle: evaluate actual outside-site halfplanes,
        // both standard and clamped, rather than reusing disk acceptance.
        for(unsigned n=0;n<24;n++) {
            double qx=xmin+(xmax-xmin)*z(rng),qy=ymin+(ymax-ymin)*z(rng);
            double gap=(.001+z(rng))*std::max(1.,r);
            switch(n%4) {case 0:qx=xmin-gap;break;case 1:qx=xmax+gap;break;
                         case 2:qy=ymin-gap;break;default:qy=ymax+gap;break;}
            double wq=n%3==0?wmax:(n%3==1?wp:wmax-1000*z(rng));
            assert(!cutter(px,py,dx,dy,wp,qx,qy,wq,false));
            assert(!cutter(px,py,dx,dy,wp,qx,qy,wq,true));
            outside_planes++;
        }
    }
    assert(accepted>8000 && rejected>8000);
    printf("{\"zero_weight_reduction_cases\":%u,\"weighted_boundary_cases\":20000,"
           "\"accepted_checked_in_long_double\":%u,\"rejected\":%u,"
           "\"outside_sites_checked_standard_and_clamped\":%u,"
           "\"explicit_clamp_transition_cases\":%u,"
           "\"tangent_invalid_overflow_stale_bound_cases\":\"pass\","
           "\"native_solver_changed\":false}\n",
           zero_reductions,accepted,rejected,outside_planes,clamp_cases);
}
