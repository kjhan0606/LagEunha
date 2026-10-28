#include "Exam/gfs_wall_dense.h"
#include <assert.h>
#include <stdio.h>

int main(void) {
    double root=0.;
    assert(gfs_wall_first_outward((GfsWallCubic){-.25,1.,0.,0.},&root)==1);
    assert(fabs(root-.25)<1.e-14);
    /* Accelerated out-and-back: final position inside, first event still found. */
    assert(gfs_wall_first_outward((GfsWallCubic){-.01,.1,-.1,0.},&root)==1);
    assert(fabs(root-(1.-sqrt(.6))*.5)<1.e-14);
    assert(gfs_wall_first_outward((GfsWallCubic){-.08,.66,-1.5,1.},&root)==1);
    assert(fabs(root-.2)<1.e-13);
    /* Euler predictor crosses, actual constant-acceleration trajectory turns. */
    GfsWallCubic p=gfs_wall_cubic(-.1,.3,-.15,-.15,-.6);
    assert(p.a0+.5*p.a1>0.);
    assert(gfs_wall_first_outward(p,&root)==0);
    assert(gfs_wall_first_outward((GfsWallCubic){-.25,1.,-1.,0.},&root)==0);
    /* Inward motion from a tiny positive contact offset is not another hit. */
    assert(gfs_wall_first_outward((GfsWallCubic){1.e-16,-1.,0.,0.},&root)==0);
    assert(gfs_wall_first_outward((GfsWallCubic){NAN,1.,0.,0.},&root)==-1);
    puts("PASS: native RK4 dense first event, out-and-back, three roots, turn-back and grazing");
    return 0;
}
