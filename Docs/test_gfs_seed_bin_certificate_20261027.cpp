#include "gfs_seed_bin_certificate_20261027.h"
#include <cmath>
#include <cstdio>
#include <random>

int main() {
    std::mt19937_64 rng(0x20261027);
    std::uniform_real_distribution<double> unit(0.,1.);
    unsigned long long skipped=0;
    for (int trial=0; trial<500000; ++trial) {
        const int mx=10,my=10;
        const int cx=(int)(unit(rng)*mx),cy=(int)(unit(rng)*my);
        const double cell=.1, xmin=-.5,ymin=-.5;
        const double x=xmin+unit(rng),y=ymin+unit(rng);
        double xx=xmin+(cx+unit(rng))*cell, yy=ymin+(cy+unit(rng))*cell;
        if (cx==0 && unit(rng)<.5) xx=xmin-unit(rng);
        if (cx==mx-1 && unit(rng)<.5) xx=xmin+mx*cell+unit(rng);
        if (cy==0 && unit(rng)<.5) yy=ymin-unit(rng);
        if (cy==my-1 && unit(rng)<.5) yy=ymin+my*cell+unit(rng);
        const double wi=unit(rng)*.1, wj=unit(rng)*.1;
        const double maxw=std::max(wi,wj)+unit(rng)*.1;
        const double radius2=unit(rng)*.2;
        if (!gfs_seed_bin_cannot_cut(x,y,cx,cy,mx,my,xmin,ymin,cell,wi,maxw,radius2)) continue;
        ++skipped;
        const double dx=xx-x,dy=yy-y,d2=dx*dx+dy*dy;
        double f=.5+.5*(wi-wj)/d2;
        f=std::max(.05,std::min(.95,f));
        if (!(f*f*d2>radius2)) {
            std::fprintf(stderr,"unsafe certificate trial=%d\n",trial);
            return 1;
        }
    }
    if (skipped<1000 || gfs_seed_bin_cannot_cut(0,0,5,5,10,10,-.5,-.5,.1,0,0,NAN)) {
        std::fprintf(stderr,"insufficient skips or nonfinite fallback failed\n");
        return 1;
    }
    std::printf("PASS 500000 bounded/edge and weighted trials; certified skips=%llu\n",skipped);
}
