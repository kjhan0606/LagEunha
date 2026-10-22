#ifndef MATERIAL_SHARED_WEIGHTED_ANGULAR_H
#define MATERIAL_SHARED_WEIGHTED_ANGULAR_H
namespace shared_target {
// SAME t=(q_lo+eta*(q_hi-q_lo))/s. Collect affine s dependence and
// weights at the eta center BEFORE bounding. No graph/field DOF changes.
static S weighted_angular(const sd::AffineAngular&lo,const sd::AffineAngular&hi,S s,S eta){
 S center(eta.n),left=S(1)-center;
 S slope=left*lo.slope+center*hi.slope;
 S intercept=left*lo.intercept+center*hi.intercept;
 S correction=left*lo.correction+center*hi.correction;
 S width_slope=hi.slope-lo.slope;
 S width_intercept=hi.intercept-lo.intercept;
 S width_correction=hi.correction-lo.correction;
 return slope+(intercept+correction)/s+S(0,eta.e)*(width_slope+(width_intercept+width_correction)/s);
}
}
#endif
