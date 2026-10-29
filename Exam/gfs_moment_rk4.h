#ifndef GFS_MOMENT_RK4_H
#define GFS_MOMENT_RK4_H
#include "gfs_moment.h"

/* Per-cell scratch independent of the RK4 storage used by the tagged moment
 * payload. Four RHS evaluations must use mutually consistent mesh/face times.
 * An external orchestrator must collectively reject on any rank/cell failure
 * before writing any candidate into the committed particle/checkpoint. */
typedef struct {
    GfsMomentState base,stage,sum;
    double dt;
    int phase,fixed_mass;
} GfsMomentRK4;

GFS_MOM_HD static inline int gfs_moment_rk4_begin(
    GfsMomentRK4 *r,const GfsMomentState *accepted,double dt,int fixed_mass){
    GfsMomentPrimitive primitive;
    if(!(dt>0) || !isfinite(dt) || (fixed_mass!=0 && fixed_mass!=1) ||
       !gfs_moment_primitive(accepted,&primitive))return 0;
    r->base=*accepted;r->stage=*accepted;gfs_moment_zero(&r->sum);
    r->dt=dt;r->phase=0;r->fixed_mass=fixed_mass;return 1;
}

/* Advances one classical RK4 stage. Failure changes neither r nor accepted.
 * Each stage checks positive-definite mass moments and positive thermal
 * residual. A finished result is still tentative until collective acceptance. */
GFS_MOM_HD static inline int gfs_moment_rk4_add_rhs(GfsMomentRK4 *r,const GfsMomentState *rhs){
    if(r->phase<0 || r->phase>=4 || !(r->dt>0) || !isfinite(r->dt))return 0;
    if(r->fixed_mass && rhs->a[GFM_M]!=0)return 0;
    for(int k=0;k<GFM_NSTATE;k++)if(!isfinite(rhs->a[k]))return 0;
    GfsMomentState next=r->base,sum=r->sum;
    for(int k=0;k<GFM_NSTATE;k++){
        double d=rhs->a[k];
        if(r->phase==0){sum.a[k]=d;next.a[k]+=.5*r->dt*d;}
        else if(r->phase==1){sum.a[k]+=2*d;next.a[k]+=.5*r->dt*d;}
        else if(r->phase==2){sum.a[k]+=2*d;next.a[k]+=r->dt*d;}
        else next.a[k]+=(r->dt/6.)*(sum.a[k]+d);
    }
    GfsMomentPrimitive primitive;
    if(!gfs_moment_primitive(&next,&primitive))return 0;
    r->sum=sum;r->stage=next;r->phase++;return 1;
}

GFS_MOM_HD static inline int gfs_moment_rk4_candidate(const GfsMomentRK4 *r,GfsMomentState *candidate){
    if(r->phase!=4)return 0;
    *candidate=r->stage;return 1;
}
#endif
