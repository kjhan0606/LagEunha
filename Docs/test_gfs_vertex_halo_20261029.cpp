#include "gfs_vertex_halo.h"
#include <cassert>
#include <cstdio>
#include <random>

int main() {
    for(double x: {-.1,.1}) for(double y: {-1.,1.})
        assert(gfs_vertex_halo_disk_inside(0,0,x,y,-1.5,-2.1,1.5,2.1));
    assert(!gfs_vertex_halo_disk_inside(0,0,.1,1.,-.5,-2.1,.5,2.1));
    assert(!gfs_vertex_halo_disk_inside(0,0,3,4,-2,-1,8,9)); // exact tangent
    assert(gfs_vertex_halo_disk_inside(0,0,3,4,-2.001,-1.001,8.001,9.001));
    assert(!gfs_vertex_halo_disk_inside(NAN,0,0,0,-1,-1,1,1));
    assert(!gfs_vertex_halo_disk_inside(0,0,0,0,-1,-1,INFINITY,1));
    assert(!gfs_vertex_halo_disk_inside(0,0,0,0,1,-1,-1,1));
    std::mt19937_64 random(20261029);
    std::uniform_real_distribution<double> u(-1.,1.);
    unsigned passed=0, rejected=0;
    for(unsigned k=0;k<100000;k++) {
        double px=1e6*u(random), py=1e6*u(random);
        double dx=10*u(random),dy=10*u(random), r=hypot(dx,dy);
        double pad=r*(1.+1e-7*u(random));
        double xmin=(px+dx)-pad,xmax=(px+dx)+pad;
        double ymin=(py+dy)-pad,ymax=(py+dy)+pad;
        int safe=gfs_vertex_halo_disk_inside(px,py,dx,dy,xmin,ymin,xmax,ymax);
        if(safe) {
            long double vx=(long double)px+dx,vy=(long double)py+dy;
            long double r2=(long double)dx*dx+(long double)dy*dy;
            for(long double gap: {vx-xmin,(long double)xmax-vx,vy-ymin,(long double)ymax-vy})
                assert(gap>0 && gap*gap>r2);
            passed++;
        } else rejected++;
    }
    assert(passed>40000 && rejected>40000);
    printf("PASS: exact tangent/invalid inputs rejected; %u accepted disks checked with long-double inequalities, %u rejected\n",passed,rejected);
}
