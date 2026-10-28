#include "gfs_weight_motion_remap_20261027.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>

static bool near(double x,double y){return std::fabs(x-y)<1e-12*(1.+std::fabs(y));}
static double angular(const gfs_spin_affine_cell &a){
    const auto &s=a.mechanical;return s.x*s.py-s.y*s.px+s.spin;
}

int main(){
    const double vel=.4,u=2.,etot=u+.5*vel*vel;
    gfs_spin_affine_cell a{{-.5,0.,1.,vel,0.,etot,0.,1.},0.};
    gfs_spin_affine_cell b{{.5,0.,1.,vel,0.,etot,0.,1.},0.};
    const auto olda=a,oldb=b;
    gfs_weight_motion_input in{0.,0.,1.,1.,0.,0.,.2,0.,0.,1.,1.,1.,
                                {1.,vel,0.,etot,0.},{1.,vel,0.,etot,0.}};
    gfs_spin_affine_transfer f;
    if(!gfs_weight_motion_transfer(in,&f) || !near(f.mechanical.dm,-.1))return 1;
    if(gfs_spin_affine_face_exchange(&a,&b,f))return 2;
    if(!near(a.mechanical.m,1.1) || !near(b.mechanical.m,.9) ||
       !near(a.mechanical.m/1.1,1.) || !near(b.mechanical.m/.9,1.) ||
       !near(gfs_spin_affine_thermal(a)/a.mechanical.m,u) ||
       !near(gfs_spin_affine_thermal(b)/b.mechanical.m,u) ||
       !near(a.mechanical.energy+b.mechanical.energy,olda.mechanical.energy+oldb.mechanical.energy) ||
       !near(angular(a)+angular(b),angular(olda)+angular(oldb)))return 3;
    in.wdot_i=0.;in.wdot_j=0.;
    if(!gfs_weight_motion_transfer(in,&f) || f.mechanical.dm!=0. ||
       f.mechanical.dpx!=0. || f.mechanical.dpy!=0. || f.mechanical.denergy!=0.)return 4;
    in.wdot_i=.3;in.wdot_j=.3;
    if(!gfs_weight_motion_transfer(in,&f) || f.mechanical.dm!=0.)return 5;
    in.wdot_i=0.;in.wdot_j=.2;in.state_i.density=2.;
    if(!gfs_weight_motion_transfer(in,&f) || !near(f.mechanical.dm,.2))return 6;
    const auto oldf=f;
    in.site_distance=0.;
    if(gfs_weight_motion_transfer(in,&f) || std::memcmp(&f,&oldf,sizeof(f)))return 7;
    std::mt19937_64 rng(20261027);
    std::uniform_real_distribution<double> v(-1.,1.);
    for(int k=0;k<50000;k++){
        gfs_weight_motion_input q{v(rng),v(rng),.5+.2*v(rng),1.+.2*v(rng),
            v(rng),v(rng),v(rng),v(rng),v(rng),.1,1.,1.,
            {1.+.2*v(rng),v(rng),v(rng),3.,.1},
                                      {1.+.2*v(rng),v(rng),v(rng),3.,.1}};
        gfs_spin_affine_transfer x,y,z;
        if(!gfs_weight_motion_transfer(q,&x))return 8;
        auto swap=q;swap.wdot_i=q.wdot_j;swap.wdot_j=q.wdot_i;
        swap.weight_i=q.weight_j;swap.weight_j=q.weight_i;
        swap.state_i=q.state_j;swap.state_j=q.state_i;
        if(!gfs_weight_motion_transfer(swap,&y) ||
           !near(x.mechanical.dm,-y.mechanical.dm) ||
           !near(x.mechanical.dpx,-y.mechanical.dpx) ||
           !near(x.mechanical.dpy,-y.mechanical.dpy) ||
           !near(x.mechanical.denergy,-y.mechanical.denergy) ||
           !near(x.nonrigid_delta_i,y.nonrigid_delta_j))return 9;
        const double ux=.7,uy=-.4,u2=ux*ux+uy*uy;
        auto boost=q;
        for(auto *s:{&boost.state_i,&boost.state_j}){
            s->specific_energy+=ux*s->vx+uy*s->vy+.5*u2;
            s->vx+=ux;s->vy+=uy;
        }
        if(!gfs_weight_motion_transfer(boost,&z) ||
           !near(z.mechanical.dm,x.mechanical.dm) ||
           !near(z.mechanical.dpx,x.mechanical.dpx+ux*x.mechanical.dm) ||
           !near(z.mechanical.dpy,x.mechanical.dpy+uy*x.mechanical.dm) ||
           !near(z.mechanical.denergy,x.mechanical.denergy+ux*x.mechanical.dpx+
                                     uy*x.mechanical.dpy+.5*u2*x.mechanical.dm))return 10;
    }
    gfs_weight_motion_input sep{0.,0.,1.,1.,.2,0.,0.,0.,.4,.1,1.,1.,
                                {2.,0.,0.,3.,0.},{1.,0.,0.,3.,0.}};
    if(!gfs_weight_motion_transfer(sep,&f) || !near(f.mechanical.dm,.008))return 11;
    sep.weight_i=sep.weight_j=0.;
    if(!gfs_weight_motion_transfer(sep,&f) || f.mechanical.dm!=0.)return 12;
    std::puts("PASS: W0 reduction, both power-plane speed terms, donor sign, 50000 face swaps/boosts and rollback");
}
