#include "gfs_circle_observer.h"
#include <string.h>

static int finite_state(GfsCircleReplayState s){
    return isfinite(s.x)&&isfinite(s.y)&&isfinite(s.vx)&&isfinite(s.vy);
}
static double gap(GfsCircleReplayState s,GfsCircleObserverOptions o){
    return hypot(s.x-o.cx,s.y-o.cy)-o.radius;
}
static double speed(GfsCircleReplayState s,GfsCircleObserverOptions o){
    double x=s.x-o.cx,y=s.y-o.cy,r=hypot(x,y);
    return r>0.?(x*s.vx+y*s.vy)/r:0.;
}
double gfs_circle_observer_tolerance(GfsCircleReplayState s,GfsCircleObserverOptions o){
    return 32.*DBL_EPSILON*fmax(o.radius,fabs(s.x-o.cx)+fabs(s.y-o.cy));
}
static int stop(GfsCircleObserver *s,int reason){s->action=GFS_CIRCLE_OBSERVER_STOP;s->reason=reason;return s->action;}
static int retry(GfsCircleObserver *s,double h,int reason){
    if(!(h>0.)||!isfinite(h))return stop(s,GFS_CIRCLE_OBSERVER_TIME_LIMIT);
    if(s->trials>=s->options.max_trials)return stop(s,GFS_CIRCLE_OBSERVER_LIMIT);
    s->ceiling=h;s->reason=reason;s->action=GFS_CIRCLE_OBSERVER_RETRY;return s->action;
}
static int smaller(GfsCircleObserver *s,double h,int reason){
    s->refining=0;s->stalled=0;return retry(s,.5*h,reason);
}
void gfs_circle_observer_init(GfsCircleObserver *s,GfsCircleObserverOptions o,
    size_t count,uint64_t epoch,double ceiling){
    memset(s,0,sizeof(*s));s->options=o;s->count=count;s->epoch=epoch;
    s->selected_id=UINT64_MAX;s->action=GFS_CIRCLE_OBSERVER_RETRY;s->ceiling=ceiling;
    if(!count||!isfinite(ceiling)||!(ceiling>0.)||!isfinite(o.cx)||!isfinite(o.cy)||
       !isfinite(o.radius)||!(o.radius>0.)||!isfinite(o.minimum_normal_speed)||
       !(o.minimum_normal_speed>=0.)||o.max_trials<1||o.max_trials>128||
       o.max_dense_nodes<1||o.max_dense_nodes>4096||o.max_dense_depth<0||o.max_dense_depth>48)
        stop(s,GFS_CIRCLE_OBSERVER_INVALID);
}
int gfs_circle_observer_initial(GfsCircleObserver *s,GfsCircleObservationView view,void *ctx){
    if(s->action!=GFS_CIRCLE_OBSERVER_RETRY||s->initial_checked||!view)
        return stop(s,GFS_CIRCLE_OBSERVER_INVALID);
    s->zero_contacts=0;
    for(size_t i=0;i<s->count;i++){
        GfsCircleObservation p;
        if(view(ctx,i,&p)||!finite_state(p.start)||p.id==UINT64_MAX)
            return stop(s,GFS_CIRCLE_OBSERVER_INVALID);
        double g=gap(p.start,s->options),tol=gfs_circle_observer_tolerance(p.start,s->options);
        if(!isfinite(g)||g>tol)return stop(s,GFS_CIRCLE_OBSERVER_INVALID);
        if(fabs(g)<=tol && speed(p.start,s->options)>s->options.minimum_normal_speed){
            if(!s->zero_contacts)s->selected_id=p.id;
            s->zero_contacts++;
        }
    }
    s->initial_checked=1;
    if(s->zero_contacts){s->action=GFS_CIRCLE_OBSERVER_ZERO;s->reason=GFS_CIRCLE_OBSERVER_INITIAL;}
    return s->action;
}
static int valid_trial(GfsCircleObservation p){
    if(!finite_state(p.start)||!finite_state(p.trial.end))return 0;
    for(int k=0;k<4;k++)if(!isfinite(p.trial.kx[k])||!isfinite(p.trial.ky[k]))return 0;
    double dx=(p.trial.kx[0]+2.*p.trial.kx[1]+2.*p.trial.kx[2]+p.trial.kx[3])/6.;
    double dy=(p.trial.ky[0]+2.*p.trial.ky[1]+2.*p.trial.ky[2]+p.trial.ky[3])/6.;
    return isfinite(dx)&&isfinite(dy)&&
        fabs(p.trial.end.x-(p.start.x+dx))<=32.*DBL_EPSILON*fmax(1.,fabs(p.start.x)+fabs(dx))&&
        fabs(p.trial.end.y-(p.start.y+dy))<=32.*DBL_EPSILON*fmax(1.,fabs(p.start.y)+fabs(dy));
}
static GfsCircleDenseResult candidate(GfsCircleObservation p,double h,GfsCircleObserverOptions o){
    double x=p.start.x-o.cx,y=p.start.y-o.cy;
    double bound=(5./24.)*hypot(p.trial.kx[0],p.trial.ky[0])+
        (hypot(p.trial.kx[1],p.trial.ky[1])+hypot(p.trial.kx[2],p.trial.ky[2]))/3.+
        hypot(p.trial.kx[3],p.trial.ky[3])/6.;
    double guard=64.*DBL_EPSILON*fmax(o.radius,fabs(x)+fabs(y)+bound);
    if(isfinite(bound)&&hypot(x,y)+bound+guard<o.radius){
        GfsCircleDenseResult no;
        memset(&no,0,sizeof(no));no.status=GFS_CIRCLE_DENSE_NONE;no.theta_hi=1.;return no;
    }
    GfsWallCubic a=gfs_wall_cubic(x,p.trial.kx[0],p.trial.kx[1],p.trial.kx[2],p.trial.kx[3]);
    GfsWallCubic b=gfs_wall_cubic(y,p.trial.ky[0],p.trial.ky[1],p.trial.ky[2],p.trial.ky[3]);
    GfsCircleDenseOptions options={o.radius,h,gfs_circle_observer_tolerance(p.start,o),
        o.minimum_normal_speed,o.max_dense_nodes,o.max_dense_depth};
    return gfs_circle_dense_first_candidate(a,b,options);
}
static int shorten_before(GfsCircleObserver *s,double h,GfsCircleDenseResult d,int reason){
    s->refining=0;s->stalled=0;
    double theta=d.theta_lo>0.?d.theta_lo:d.theta_hi;
    if(!(theta>0.))return smaller(s,h,GFS_CIRCLE_OBSERVER_UNCERTAIN);
    return retry(s,.5*theta*h,reason);
}
static int probe_first_exit(GfsCircleObserver *s,double h,GfsCircleDenseResult d,int reason){
    /* A small real-RK probe beyond the certified first dense root can obtain
     * an actual outside bracket without a Zeno sequence of accepted steps
     * approaching contact from below. This is only a proposed time, never
     * a projected position or an event certificate. Every participant is
     * rescanned on the resulting actual native trial.
     */
    if(d.dense_sign_bracket_certified && d.theta_hi>0. && d.theta_hi<.75){
        s->refining=0;s->stalled=0;
        return retry(s,fmin(.875,1.125*d.theta_hi)*h,reason);
    }
    return shorten_before(s,h,d,reason);
}
static int next_bracket(GfsCircleObserver *s,double seed){
    double h=seed;
    if(s->stalled>=2){h=s->lo+.5*(s->hi-s->lo);s->stalled=0;s->bisections++;}
    if(!(h>s->lo && h<s->hi)||!isfinite(h))
        h=s->lo-s->glo*(s->hi-s->lo)/(s->ghi-s->glo);
    if(!(h>s->lo && h<s->hi)||!isfinite(h)){
        h=s->lo+.5*(s->hi-s->lo);s->bisections++;
    }
    if(!(h>s->lo && h<s->hi))return stop(s,GFS_CIRCLE_OBSERVER_TIME_LIMIT);
    return retry(s,h,GFS_CIRCLE_OBSERVER_BRACKET);
}
int gfs_circle_observer_observe(GfsCircleObserver *s,uint64_t epoch,double h,
    GfsCircleObservationView view,void *ctx){
    if(epoch!=s->epoch)return stop(s,GFS_CIRCLE_OBSERVER_BASE_CHANGED);
    if(s->action!=GFS_CIRCLE_OBSERVER_RETRY||!s->initial_checked||!view||
       !isfinite(h)||!(h>0.)||h>s->ceiling*(1.+64.*DBL_EPSILON))
        return stop(s,GFS_CIRCLE_OBSERVER_INVALID);
    if(s->trials>=s->options.max_trials)return stop(s,GFS_CIRCLE_OBSERVER_LIMIT);
    s->trials++;
    int have=0,active_found=0;
    GfsCircleObservation best,active;
    GfsCircleDenseResult best_d,active_d;
    memset(&best,0,sizeof(best));memset(&active,0,sizeof(active));
    memset(&best_d,0,sizeof(best_d));memset(&active_d,0,sizeof(active_d));
    double second_lo=INFINITY;
    for(size_t i=0;i<s->count;i++){
        GfsCircleObservation p;
        if(view(ctx,i,&p)||!valid_trial(p)||p.id==UINT64_MAX)
            return stop(s,GFS_CIRCLE_OBSERVER_INVALID);
        double initial_gap=gap(p.start,s->options),tol=gfs_circle_observer_tolerance(p.start,s->options);
        if(initial_gap>tol)return stop(s,GFS_CIRCLE_OBSERVER_BASE_CHANGED);
        GfsCircleDenseResult d=candidate(p,h,s->options);
        if(d.status==GFS_CIRCLE_DENSE_UNRESOLVED)
            return smaller(s,h,GFS_CIRCLE_OBSERVER_UNCERTAIN);
        if(d.status==GFS_CIRCLE_DENSE_NONE && gap(p.trial.end,s->options)>
            gfs_circle_observer_tolerance(p.trial.end,s->options))
            return stop(s,GFS_CIRCLE_OBSERVER_INVALID);
        if(s->refining && p.id==s->selected_id){
            if(memcmp(&p.start,&s->selected_start,sizeof(p.start)))
                return stop(s,GFS_CIRCLE_OBSERVER_BASE_CHANGED);
            active=p;active_d=d;active_found=1;
        }
        if(d.status==GFS_CIRCLE_DENSE_CANDIDATE){
            if(d.reason==GFS_CIRCLE_DENSE_START_CONTACT)
                return stop(s,GFS_CIRCLE_OBSERVER_BASE_CHANGED); /* must initial-scan first */
            /* A tangent/weak window cannot seed an actual impact bracket.
             * A weak endpoint during refinement is allowed only after a
             * first-root sign witness plus an actual outside endpoint were
             * already retained for this selected particle. Otherwise stop
             * explicitly rather than turning repeated approach into a false
             * initial zero-time impulse at a merely tangent trajectory.
             */
            if(!d.dense_sign_bracket_certified && !(s->refining && p.id==s->selected_id))
                return stop(s,GFS_CIRCLE_OBSERVER_UNCERTAIN);
            if(!have || d.theta_hi<best_d.theta_hi){
                if(have)second_lo=fmin(second_lo,best_d.theta_lo);
                best=p;best_d=d;have=1;
            }else second_lo=fmin(second_lo,d.theta_lo);
        }
    }
    if(s->refining && !active_found)return stop(s,GFS_CIRCLE_OBSERVER_BASE_CHANGED);
    if(have && (!s->refining||best.id!=s->selected_id)){
        if(s->refining)s->switches++;
        s->refining=0;
        double g=gap(best.trial.end,s->options);
        if(!(g>0.))return probe_first_exit(s,h,best_d,GFS_CIRCLE_OBSERVER_OUT_AND_BACK);
        double start_gap=gap(best.start,s->options);
        if(!(start_gap<0.))
            return shorten_before(s,h,best_d,GFS_CIRCLE_OBSERVER_INWARD_START);
        s->selected_id=best.id;s->selected_start=best.start;s->refining=1;
        s->lo=0.;s->glo=start_gap;s->hi=h;s->ghi=g;s->stalled=0;
        return next_bracket(s,.5*(best_d.theta_lo+best_d.theta_hi)*h);
    }
    if(!s->refining){
        if(have)return stop(s,GFS_CIRCLE_OBSERVER_INVALID);
        s->action=GFS_CIRCLE_OBSERVER_QUIET;s->reason=GFS_CIRCLE_OBSERVER_NO_CANDIDATE;return s->action;
    }
    double g=gap(active.trial.end,s->options);
    double tol=gfs_circle_observer_tolerance(active.trial.end,s->options);
    double vn=speed(active.trial.end,s->options);
    double time_band=2.*tol/fmax(fabs(vn),fmax(s->options.minimum_normal_speed,DBL_MIN));
    if(fabs(g)<=tol){
        if(have && second_lo<=best_d.theta_hi)
            return stop(s,GFS_CIRCLE_OBSERVER_OVERLAP);
        if(have && best_d.theta_hi*h<h-time_band)
            return shorten_before(s,h,best_d,GFS_CIRCLE_OBSERVER_OUT_AND_BACK);
        if(!(vn>s->options.minimum_normal_speed))
            return smaller(s,h,GFS_CIRCLE_OBSERVER_GRAZING);
        s->action=GFS_CIRCLE_OBSERVER_CONTACT;s->reason=GFS_CIRCLE_OBSERVER_REFINED;return s->action;
    }
    if(g<0. && active_d.status==GFS_CIRCLE_DENSE_CANDIDATE)
        return probe_first_exit(s,h,active_d,GFS_CIRCLE_OBSERVER_OUT_AND_BACK);
    if(!have && g<0. && (h<s->ceiling*(1.-64.*DBL_EPSILON)||h<=s->lo)){
        s->action=GFS_CIRCLE_OBSERVER_QUIET;s->reason=GFS_CIRCLE_OBSERVER_NO_CANDIDATE;return s->action;
    }
    double width=s->hi-s->lo;
    if(g<0.){if(!(h>s->lo))return smaller(s,h,GFS_CIRCLE_OBSERVER_UNCERTAIN);s->lo=h;s->glo=g;}
    else {if(!(h<s->hi))return smaller(s,h,GFS_CIRCLE_OBSERVER_UNCERTAIN);s->hi=h;s->ghi=g;}
    if(s->hi-s->lo>.75*width)s->stalled++;else s->stalled=0;
    return next_bracket(s,NAN);
}
int gfs_circle_observer_reject_stage(GfsCircleObserver *s,uint64_t epoch,double h){
    if(epoch!=s->epoch)return stop(s,GFS_CIRCLE_OBSERVER_BASE_CHANGED);
    if(s->action!=GFS_CIRCLE_OBSERVER_RETRY||!s->initial_checked||!isfinite(h)||!(h>0.)||
       h>s->ceiling*(1.+64.*DBL_EPSILON))
        return stop(s,GFS_CIRCLE_OBSERVER_INVALID);
    s->trials++;
    return smaller(s,h,GFS_CIRCLE_OBSERVER_FAILED_STAGE);
}
