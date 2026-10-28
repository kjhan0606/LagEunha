#include "gfs_spin_affine_face_exchange_20261027.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>

static bool near(double x,double y){return std::fabs(x-y)<1e-12*(1.+std::fabs(y));}
static double angular(const gfs_spin_affine_cell &a){
    const auto &s=a.mechanical;return s.x*s.py-s.y*s.px+s.spin;
}

int main(){
    gfs_spin_affine_cell a{{-1.,0.,1.,.1,.2,2.,.1,1.2},.25};
    gfs_spin_affine_cell b{{1.,0.,1.3,-.1,.05,2.3,-.2,1.6},.3};
    const auto x=a,y=b;
    gfs_spin_affine_transfer f{{.1,.4,.03,.05,.03,.02,1.2,1.6},-.02,.02};
    if(gfs_spin_affine_face_exchange(&a,&b,f))return 1;
    if(!near(a.mechanical.m+b.mechanical.m,x.mechanical.m+y.mechanical.m) ||
       !near(a.mechanical.px+b.mechanical.px,x.mechanical.px+y.mechanical.px) ||
       !near(a.mechanical.py+b.mechanical.py,x.mechanical.py+y.mechanical.py) ||
       !near(a.mechanical.energy+b.mechanical.energy,x.mechanical.energy+y.mechanical.energy) ||
       !near(a.mechanical.m*a.mechanical.x+b.mechanical.m*b.mechanical.x,
             x.mechanical.m*x.mechanical.x+y.mechanical.m*y.mechanical.x) ||
       !near(a.mechanical.m*a.mechanical.y+b.mechanical.m*b.mechanical.y,
             x.mechanical.m*x.mechanical.y+y.mechanical.m*y.mechanical.y) ||
       !near(angular(a)+angular(b),angular(x)+angular(y)) ||
       !near(a.nonrigid_energy+b.nonrigid_energy,x.nonrigid_energy+y.nonrigid_energy))return 2;
    if(!(gfs_spin_affine_thermal(a)>0.) || !(gfs_spin_affine_thermal(b)>0.))return 3;
    const auto beforea=a,beforeb=b;
    f.nonrigid_delta_i=100.;
    if(gfs_spin_affine_face_exchange(&a,&b,f)==0 || std::memcmp(&a,&beforea,sizeof(a)) ||
       std::memcmp(&b,&beforeb,sizeof(b)))return 4;
    std::mt19937_64 rng(20261027);
    std::uniform_real_distribution<double> unit(-1.,1.);
    for(int k=0;k<50000;k++) {
        gfs_spin_affine_cell c{{unit(rng),unit(rng),1.5+.2*unit(rng),
                                 .1*unit(rng),.1*unit(rng),6.,.1*unit(rng),
                                 1.+.2*unit(rng)},.3+.05*unit(rng)};
        gfs_spin_affine_cell d{{unit(rng),unit(rng),1.5+.2*unit(rng),
                                 .1*unit(rng),.1*unit(rng),6.,.1*unit(rng),
                                 1.+.2*unit(rng)},.3+.05*unit(rng)};
        const auto c0=c,d0=d;
        const double nr=.01*unit(rng);
        gfs_spin_affine_transfer t{{unit(rng),unit(rng),.03*unit(rng),
                                    .03*unit(rng),.03*unit(rng),.03*unit(rng),
                                    c.mechanical.inertia,d.mechanical.inertia},nr,-nr};
        const double ux=.7,uy=-.4,u2=ux*ux+uy*uy;
        auto boost=[&](gfs_spin_affine_cell &z) {
            auto &s=z.mechanical;
            s.energy+=ux*s.px+uy*s.py+.5*s.m*u2;
            s.px+=s.m*ux;s.py+=s.m*uy;
        };
        auto cb=c,db=d;boost(cb);boost(db);
        auto tb=t;
        tb.mechanical.denergy+=ux*t.mechanical.dpx+uy*t.mechanical.dpy+.5*t.mechanical.dm*u2;
        tb.mechanical.dpx+=t.mechanical.dm*ux;tb.mechanical.dpy+=t.mechanical.dm*uy;
        if(gfs_spin_affine_face_exchange(&c,&d,t))return 5;
        if(gfs_spin_affine_face_exchange(&cb,&db,tb))return 7;
        if(!near(c.mechanical.m+d.mechanical.m,c0.mechanical.m+d0.mechanical.m) ||
           !near(c.mechanical.px+d.mechanical.px,c0.mechanical.px+d0.mechanical.px) ||
           !near(c.mechanical.py+d.mechanical.py,c0.mechanical.py+d0.mechanical.py) ||
           !near(c.mechanical.energy+d.mechanical.energy,c0.mechanical.energy+d0.mechanical.energy) ||
           !near(c.mechanical.m*c.mechanical.x+d.mechanical.m*d.mechanical.x,
                 c0.mechanical.m*c0.mechanical.x+d0.mechanical.m*d0.mechanical.x) ||
           !near(c.mechanical.m*c.mechanical.y+d.mechanical.m*d.mechanical.y,
                 c0.mechanical.m*c0.mechanical.y+d0.mechanical.m*d0.mechanical.y) ||
           !near(angular(c)+angular(d),angular(c0)+angular(d0)) ||
           !near(c.nonrigid_energy+d.nonrigid_energy,c0.nonrigid_energy+d0.nonrigid_energy))return 6;
        if(!near(cb.mechanical.x,c.mechanical.x) || !near(cb.mechanical.y,c.mechanical.y) ||
           !near(db.mechanical.x,d.mechanical.x) || !near(db.mechanical.y,d.mechanical.y) ||
           !near(cb.mechanical.spin,c.mechanical.spin) || !near(db.mechanical.spin,d.mechanical.spin) ||
           !near(cb.mechanical.px,c.mechanical.px+c.mechanical.m*ux) ||
           !near(cb.mechanical.py,c.mechanical.py+c.mechanical.m*uy) ||
           !near(db.mechanical.px,d.mechanical.px+d.mechanical.m*ux) ||
           !near(db.mechanical.py,d.mechanical.py+d.mechanical.m*uy) ||
           !near(gfs_spin_affine_thermal(cb),gfs_spin_affine_thermal(c)) ||
           !near(gfs_spin_affine_thermal(db),gfs_spin_affine_thermal(d)))return 8;
    }
    std::puts("PASS: 50000 pair invariants, Galilean boosts and positivity rollback");
}
