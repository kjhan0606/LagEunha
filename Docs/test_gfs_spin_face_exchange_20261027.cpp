#include "gfs_spin_face_exchange_20261027.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>

static int require(bool ok,const char *message) {
    if(!ok) std::fprintf(stderr,"FAIL: %s\n",message);
    return ok ? 0 : 1;
}

static double angular(const gfs_spin_cell &a,const gfs_spin_cell &b) {
    return a.x*a.py-a.y*a.px+a.spin+b.x*b.py-b.y*b.px+b.spin;
}

static bool conserved(const gfs_spin_cell &a0,const gfs_spin_cell &b0,
                      const gfs_spin_cell &a,const gfs_spin_cell &b) {
    const double tol=1e-12;
    return std::fabs((a.m+b.m)-(a0.m+b0.m))<tol &&
           std::fabs((a.m*a.x+b.m*b.x)-(a0.m*a0.x+b0.m*b0.x))<tol &&
           std::fabs((a.m*a.y+b.m*b.y)-(a0.m*a0.y+b0.m*b0.y))<tol &&
           std::fabs((a.px+b.px)-(a0.px+b0.px))<tol &&
           std::fabs((a.py+b.py)-(a0.py+b0.py))<tol &&
           std::fabs((a.energy+b.energy)-(a0.energy+b0.energy))<tol &&
           std::fabs(angular(a,b)-angular(a0,b0))<tol;
}

int main() {
    int failed=0;
    /* A non-centered set of three masses undergoing rigid rotation gives
     * K=|p|²/(2m)+S²/(2I) only about its actual mass centroid. */
    const double pm[3]={1.,2.,1.},px[3]={0.,2.,0.},py[3]={0.,0.,2.};
    const double center_x=1.,center_y=.5,bulk_x=.3,bulk_y=-.2,omega=.7;
    double moment=0.,kinetic=0.;
    for(int k=0;k<3;++k) {
        double dx=px[k]-center_x,dy=py[k]-center_y;
        double vx=bulk_x-omega*dy,vy=bulk_y+omega*dx;
        moment+=pm[k]*(dx*dx+dy*dy);
        kinetic+=.5*pm[k]*(vx*vx+vy*vy);
    }
    const gfs_spin_cell rigid={center_x,center_y,4.,4.*bulk_x,4.*bulk_y,
                               kinetic+3.,moment*omega,moment};
    failed+=require(std::fabs(gfs_spin_thermal_energy(rigid)-3.)<1e-12,
                    "rigid-spin energy closure uses the actual mass centroid");
    const gfs_spin_cell a0={0.,0.,2.,.4,.2,10.,.01,.2};
    const gfs_spin_cell b0={1.,0.,1.5,-.1,.3,8.,-.02,.15};
    gfs_spin_cell a=a0,b=b0;
    const gfs_spin_face_transfer pressure={.45,.2,0.,.02,.03,.01,.2,.15};
    int status=gfs_spin_face_exchange(&a,&b,pressure);
    failed+=require(status==0 && conserved(a0,b0,a,b),
                    "off-center pressure exchange conserves mass, momentum, energy and angular momentum");
    failed+=require(a.m==a0.m && b.m==b0.m && a.spin!=a0.spin && b.spin!=b0.spin,
                    "zero mass flux preserves each mass while off-center momentum changes spin");
    failed+=require(gfs_spin_thermal_energy(a)>0. && gfs_spin_thermal_energy(b)>0.,
                    "thermal energy remains positive after rotational energy is included");

    gfs_spin_cell translated_a=a0,translated_b=b0;
    translated_a.x+=5.;translated_a.y-=3.;translated_b.x+=5.;translated_b.y-=3.;
    gfs_spin_face_transfer translated=pressure;translated.cx+=5.;translated.cy-=3.;
    status=gfs_spin_face_exchange(&translated_a,&translated_b,translated);
    failed+=require(status==0 && std::fabs(translated_a.spin-a.spin)<1e-14 &&
                    std::fabs(translated_b.spin-b.spin)<1e-14,
                    "spin moments are invariant under coordinate translation");

    a=a0;b=b0;
    const double dm=.02;
    const gfs_spin_face_transfer advect={.45,.2,dm,.008,-.004,.03,
                                         a0.inertia*(a0.m-dm)/a0.m,
                                         b0.inertia*(b0.m+dm)/b0.m};
    status=gfs_spin_face_exchange(&a,&b,advect);
    failed+=require(status==0 && conserved(a0,b0,a,b) && a.m<a0.m && b.m>b0.m &&
                    (a.x!=a0.x || a.y!=a0.y) && (b.x!=b0.x || b.y!=b0.y),
                    "mass, centroid, momentum, spin and energy transport conserves pair budgets");
    translated_a=a0;translated_b=b0;
    translated_a.x+=5.;translated_a.y-=3.;translated_b.x+=5.;translated_b.y-=3.;
    translated=advect;translated.cx+=5.;translated.cy-=3.;
    status=gfs_spin_face_exchange(&translated_a,&translated_b,translated);
    failed+=require(status==0 && std::fabs(translated_a.spin-a.spin)<1e-14 &&
                    std::fabs(translated_b.spin-b.spin)<1e-14 &&
                    std::fabs((translated_a.x-5.)-a.x)<1e-14 &&
                    std::fabs((translated_b.y+3.)-b.y)<1e-14,
                    "mass-transfer centroid and spin are translation invariant");

    a=a0;b=b0;
    gfs_spin_face_transfer invalid=advect;invalid.dm=3.;
    status=gfs_spin_face_exchange(&a,&b,invalid);
    failed+=require(status==-2 && std::memcmp(&a,&a0,sizeof a)==0 &&
                    std::memcmp(&b,&b0,sizeof b)==0,
                    "nonpositive mass rejects atomically");
    invalid=pressure;invalid.denergy=100.;
    status=gfs_spin_face_exchange(&a,&b,invalid);
    failed+=require(status==-3 && std::memcmp(&a,&a0,sizeof a)==0 &&
                    std::memcmp(&b,&b0,sizeof b)==0,
                    "nonpositive thermal energy rejects atomically");
    invalid=pressure;invalid.dpx=std::numeric_limits<double>::quiet_NaN();
    status=gfs_spin_face_exchange(&a,&b,invalid);
    failed+=require(status==-1 && std::memcmp(&a,&a0,sizeof a)==0 &&
                    std::memcmp(&b,&b0,sizeof b)==0,
                    "nonfinite transfer rejects atomically");
    if(failed) return 1;
    std::printf("PASS: pair mass, momentum, energy, spin, angular momentum and rejection gates\n");
    return 0;
}
