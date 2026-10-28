/* Compile the actual shared header as C, matching native Exam/exam.c. */
#include "Exam/kepler_slip_wall.h"
#include <assert.h>
#include <stdio.h>

int main(void) {
    GfsSlipImpulse hit;
    assert(gfs_slip_plane_impact(2.,.3,0.,0.,1.,0.,2.,3.,4.,5.,&hit)==GFS_SLIP_REFLECTED);
    assert(hit.vx_after==-4. && hit.vy_after==5.);
    assert(hit.wall_jx==24. && hit.wall_jy==0.);
    assert(fabs(hit.wall_torque_impulse+7.2)<1.e-14);
    assert(fabs(hit.fluid_kinetic_roundoff)<1.e-14);
    puts("PASS: C native header compatibility and oblique impulse/torque");
    return 0;
}
