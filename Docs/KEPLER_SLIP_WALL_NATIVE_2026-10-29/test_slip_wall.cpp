#include "Exam/kepler_slip_wall.h"
#include <cassert>
#include <cstdio>
#include <cstring>
#include <limits>

static bool near(double a, double b, double scale=1.) {
    return fabs(a-b) <= 256.*DBL_EPSILON*fmax(scale,fmax(fabs(a),fabs(b)));
}

int main() {
    GfsSlipImpulse hit;
    /* Plane: oblique incidence preserves tangential velocity and kinetic E. */
    assert(gfs_slip_plane_impact(2.,.3,0.,0.,1.,0.,2.,3.,4.,5.,&hit)==GFS_SLIP_REFLECTED);
    assert(hit.vx_after==-4. && hit.vy_after==5.);
    assert(hit.wall_jx==24. && hit.wall_jy==0.);
    assert(near(hit.wall_torque_impulse,-7.2));
    assert(near(hit.fluid_kinetic_roundoff,0.));
    const double l_before=3.*(2.*5.-.3*4.);
    const double l_after=3.*(2.*hit.vy_after-.3*hit.vx_after);
    assert(near(l_after-l_before+hit.wall_torque_impulse,0.));

    /* Incoming and tangent states receive no impulse; no velocity clamp. */
    assert(gfs_slip_plane_impact(2.,0.,0.,0.,1.,0.,2.,1.,-4.,5.,&hit)==GFS_SLIP_NO_IMPULSE);
    assert(hit.vx_after==-4. && hit.vy_after==5. && hit.wall_jx==0.);
    assert(gfs_slip_circle_impact(2.,0.,0.,0.,2.,1.,0.,7.,&hit)==GFS_SLIP_NO_IMPULSE);
    assert(hit.vx_after==0. && hit.vy_after==7.);

    /* Sweep non-axis-aligned circular impacts. Center angular momentum and
     * kinetic energy are invariant, and fluid+wall linear impulse closes. */
    for (int k=0;k<1024;++k) {
        const double theta=(k+.37)*6.2831853071795864769/1024.;
        const double nx=cos(theta),ny=sin(theta),mass=.017;
        const double rx=2.75*nx,ry=2.75*ny;
        const double vx=.43*nx-.71*ny,vy=.43*ny+.71*nx;
        assert(gfs_slip_circle_impact(3.+rx,3.+ry,3.,3.,2.75,mass,vx,vy,&hit)==GFS_SLIP_REFLECTED);
        assert(near(vx*vx+vy*vy,hit.vx_after*hit.vx_after+hit.vy_after*hit.vy_after));
        assert(near(mass*(hit.vx_after-vx)+hit.wall_jx,0.));
        assert(near(mass*(hit.vy_after-vy)+hit.wall_jy,0.));
        assert(near(rx*vy-ry*vx,rx*hit.vy_after-ry*hit.vx_after));
        assert(near(hit.wall_torque_impulse,0.));
        assert(near(hit.fluid_kinetic_roundoff,0.));
    }

    /* A polygon facet generally exerts torque: account for it, rather than
     * incorrectly claiming circle angular momentum conservation. */
    const double theta=.3,nx=cos(theta),ny=sin(theta),a=2.75*cos(3.141592653589793/1024.);
    const double rx=a*nx-.001*ny,ry=a*ny+.001*nx;
    const double vx=.4*nx-.7*ny,vy=.4*ny+.7*nx;
    assert(gfs_slip_plane_impact(3.+rx,3.+ry,3.,3.,nx,ny,a,.2,vx,vy,&hit)==GFS_SLIP_REFLECTED);
    assert(fabs(hit.wall_torque_impulse)>1.e-5);
    assert(near(.2*(rx*(hit.vy_after-vy)-ry*(hit.vx_after-vx))+hit.wall_torque_impulse,0.));
    assert(near(hit.fluid_kinetic_roundoff,0.));

    /* Invalid/off-boundary events fail without modifying the output. */
    const GfsSlipImpulse before=hit;
    assert(gfs_slip_circle_impact(2.01,0.,0.,0.,2.,1.,1.,0.,&hit)==GFS_SLIP_INVALID);
    assert(memcmp(&before,&hit,sizeof(hit))==0);
    assert(gfs_slip_plane_impact(2.,0.,0.,0.,2.,0.,2.,1.,1.,0.,&hit)==GFS_SLIP_INVALID);
    assert(gfs_slip_circle_impact(2.,0.,0.,0.,2.,0.,1.,0.,&hit)==GFS_SLIP_INVALID);
    assert(gfs_slip_circle_impact(2.,0.,0.,0.,2.,1.,std::numeric_limits<double>::infinity(),0.,&hit)==GFS_SLIP_INVALID);
    puts("PASS: native slip impulse, 1024 circular oblique impacts, polygon torque ledger, tangent/incoming and invalid events");
    return 0;
}
