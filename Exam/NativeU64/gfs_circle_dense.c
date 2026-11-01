#include "gfs_circle_dense.h"
#include <string.h>
#if defined(__x86_64__) || defined(__i386__)
#include <xmmintrin.h>
#endif

/* Compile this small translation unit with strict IEEE arithmetic, separately
 * from the original host solver. No fast-math, reassociation or contraction.
 * Outward nextafter bounds include each rounded primitive operation.
 */
typedef struct { double lo,hi; } Interval;
typedef struct {
    Interval gap[7],derivative[6];
    double lo,hi;
    int depth;
} Node;
static const unsigned char choose[7][7]={
    {1},{1,1},{1,2,1},{1,3,3,1},{1,4,6,4,1},
    {1,5,10,10,5,1},{1,6,15,20,15,6,1}};

static Interval point(double x){Interval z={x,x};return z;}
static int is_zero(Interval a){return a.lo==0. && a.hi==0.;}
static int finite_interval(Interval a){return isfinite(a.lo)&&isfinite(a.hi)&&a.lo<=a.hi;}
static Interval add(Interval a,Interval b){
    if(!finite_interval(a)||!finite_interval(b)){Interval bad={NAN,NAN};return bad;}
    if(is_zero(a))return b;
    if(is_zero(b))return a;
    if(a.lo==a.hi && b.lo==b.hi && a.lo==-b.lo)return point(0.);
    Interval c={nextafter(a.lo+b.lo,-INFINITY),nextafter(a.hi+b.hi,INFINITY)};return c;
}
static Interval multiply(Interval a,Interval b){
    if(!finite_interval(a)||!finite_interval(b)){Interval bad={NAN,NAN};return bad;}
    if(is_zero(a)||is_zero(b))return point(0.);
    if(a.lo==1. && a.hi==1.)return b;
    if(b.lo==1. && b.hi==1.)return a;
    double v[4]={a.lo*b.lo,a.lo*b.hi,a.hi*b.lo,a.hi*b.hi};
    double lo=v[0],hi=v[0];
    for(int k=1;k<4;k++){lo=fmin(lo,v[k]);hi=fmax(hi,v[k]);}
    Interval c={nextafter(lo,-INFINITY),nextafter(hi,INFINITY)};return c;
}
static Interval divide_positive(Interval a,double b){
    if(is_zero(a))return a;
    if(b==1.)return a;
    if(a.lo==a.hi && a.lo==b)return point(1.);
    if(a.lo==a.hi && a.lo==-b)return point(-1.);
    Interval c={nextafter(a.lo/b,-INFINITY),nextafter(a.hi/b,INFINITY)};return c;
}
static Interval scale(Interval a,double b){return multiply(a,point(b));}
static Interval evaluate(const Interval *p,int degree,double theta){
    Interval value=p[degree],t=point(theta);
    for(int k=degree-1;k>=0;k--)value=add(multiply(value,t),p[k]);
    return value;
}
static Interval circle_gap(const Interval *x,const Interval *y,double theta){
    /* Preserve the shared-coordinate dependency when refining a point.
     * Expanding x^2+y^2 first needlessly widens cancellation near a root.
     */
    Interval a=evaluate(x,3,theta),b=evaluate(y,3,theta);
    return add(add(multiply(a,a),multiply(b,b)),point(-1.));
}
static void bernstein(const Interval *power,int degree,Interval *out){
    for(int k=0;k<=degree;k++){
        out[k]=point(0.);
        for(int i=0;i<=k;i++){
            Interval weight=divide_positive(point(choose[k][i]),choose[degree][i]);
            out[k]=add(out[k],multiply(weight,power[i]));
        }
    }
}
static void split(const Interval *p,int degree,Interval *left,Interval *right){
    Interval t[7];for(int k=0;k<=degree;k++)t[k]=p[k];
    left[0]=t[0];right[degree]=t[degree];
    for(int level=1;level<=degree;level++){
        for(int k=0;k<=degree-level;k++)t[k]=scale(add(t[k],t[k+1]),.5);
        left[level]=t[0];right[degree-level]=t[degree-level];
    }
}
static double minimum_lower(const Interval *p,int n){
    double x=p[0].lo;for(int k=1;k<n;k++)x=fmin(x,p[k].lo);return x;
}
static double maximum_upper(const Interval *p,int n){
    double x=p[0].hi;for(int k=1;k<n;k++)x=fmax(x,p[k].hi);return x;
}
static GfsCircleDenseResult unresolved(GfsCircleDenseResult r,int reason){
    r.status=GFS_CIRCLE_DENSE_UNRESOLVED;r.reason=reason;return r;
}
static GfsCircleDenseResult contact_window(GfsCircleDenseResult r,double theta,
    Interval gap,double derivative_lower,double bracket_lo,double bracket_hi,int reason){
    double error=fmax(fabs(gap.lo),fabs(gap.hi));
    double dt=nextafter(error/derivative_lower,INFINITY);
    r.theta_lo=fmax(bracket_lo,nextafter(theta-dt,-INFINITY));
    r.theta_hi=fmin(bracket_hi,nextafter(theta+dt,INFINITY));
    r.status=GFS_CIRCLE_DENSE_CANDIDATE;r.reason=reason;return r;
}

/* Keep the strict predicate behind a call boundary: the native Intel driver
 * enables FTZ/DAZ at startup. Moving arithmetic across the scoped MXCSR
 * change would invalidate subnormal interval endpoints. No solver-wide
 * floating-point option or caller exception flag is changed.
 */
#if defined(__GNUC__) || defined(__clang__)
__attribute__((noinline))
#endif
static GfsCircleDenseResult first_candidate_strict(
    GfsWallCubic x,GfsWallCubic y,GfsCircleDenseOptions o){
    GfsCircleDenseResult result;
    memset(&result,0,sizeof(result));
    result.theta_hi=1.;
    if(!(o.radius>0.) || !(o.step>0.) || !(o.radial_tolerance>0.) ||
       !(o.radial_tolerance<o.radius) || !(o.minimum_normal_speed>=0.) ||
       !isfinite(o.radius)||!isfinite(o.step)||!isfinite(o.radial_tolerance)||
       !isfinite(o.minimum_normal_speed)||o.max_nodes<1||o.max_nodes>4096||
       o.max_depth<0||o.max_depth>48)
        return unresolved(result,GFS_CIRCLE_DENSE_INVALID_INPUT);
    volatile double tiny=DBL_MIN;
    volatile double half=tiny*.5;
    if(FLT_RADIX!=2 || DBL_MANT_DIG!=53 || !(half>0.))
        return unresolved(result,GFS_CIRCLE_DENSE_FP_ENVIRONMENT);

    double xr[4]={x.a0,x.a1,x.a2,x.a3},yr[4]={y.a0,y.a1,y.a2,y.a3};
    Interval a[4],b[4],power[7],derivative[6];
    for(int k=0;k<4;k++){
        if(!isfinite(xr[k])||!isfinite(yr[k]))return unresolved(result,GFS_CIRCLE_DENSE_INVALID_INPUT);
        a[k]=divide_positive(point(xr[k]),o.radius);b[k]=divide_positive(point(yr[k]),o.radius);
    }
    for(int k=0;k<7;k++)power[k]=point(0.);
    for(int i=0;i<4;i++)for(int j=0;j<4;j++){
        power[i+j]=add(power[i+j],add(multiply(a[i],a[j]),multiply(b[i],b[j])));
    }
    power[0]=add(power[0],point(-1.));
    for(int k=1;k<7;k++)derivative[k-1]=scale(power[k],k);
    Interval u=divide_positive(point(o.radial_tolerance),o.radius);
    Interval low=add(scale(u,-2.),multiply(u,u)),high=add(scale(u,2.),multiply(u,u));
    Interval floor=scale(multiply(divide_positive(multiply(point(o.step),point(o.minimum_normal_speed)),o.radius),add(point(1.),u)),2.);
    for(int k=0;k<7;k++)if(!finite_interval(power[k]))return unresolved(result,GFS_CIRCLE_DENSE_INVALID_INPUT);
    if(!finite_interval(low)||!finite_interval(high)||!finite_interval(floor))return unresolved(result,GFS_CIRCLE_DENSE_INVALID_INPUT);
    if(power[0].lo>high.hi)return unresolved(result,GFS_CIRCLE_DENSE_OUTSIDE_START);
    if(power[0].hi>high.lo)return unresolved(result,GFS_CIRCLE_DENSE_UNCERTAIN_START);
    if(power[0].lo>=low.hi && power[0].hi<=high.lo && derivative[0].lo>floor.hi){
        result.status=GFS_CIRCLE_DENSE_CANDIDATE;result.reason=GFS_CIRCLE_DENSE_START_CONTACT;
        result.theta_lo=result.theta_hi=0.;return result;
    }

    Node stack[50];int count=1;
    bernstein(power,6,stack[0].gap);bernstein(derivative,5,stack[0].derivative);
    stack[0].lo=0.;stack[0].hi=1.;stack[0].depth=0;
    while(count){
        if(result.nodes>=o.max_nodes){
            result.theta_lo=stack[count-1].lo;result.theta_hi=stack[count-1].hi;
            return unresolved(result,GFS_CIRCLE_DENSE_NODE_LIMIT);
        }
        Node node=stack[--count];result.nodes++;
        result.max_depth=node.depth>result.max_depth?node.depth:result.max_depth;
        result.theta_lo=node.lo;result.theta_hi=node.hi;
        for(int k=0;k<7;k++)if(!finite_interval(node.gap[k]))
            return unresolved(result,GFS_CIRCLE_DENSE_ROUNDING_LIMIT);
        for(int k=0;k<6;k++)if(!finite_interval(node.derivative[k]))
            return unresolved(result,GFS_CIRCLE_DENSE_ROUNDING_LIMIT);
        double gmin=minimum_lower(node.gap,7),gmax=maximum_upper(node.gap,7);
        double qmin=minimum_lower(node.derivative,6),qmax=maximum_upper(node.derivative,6);
        if(!isfinite(gmin)||!isfinite(gmax)||!isfinite(qmin)||!isfinite(qmax))
            return unresolved(result,GFS_CIRCLE_DENSE_ROUNDING_LIMIT);
        /* Monotonically inward intervals cannot contain an outward event,
         * including the tolerated initial contact immediately after reflection.
         */
        if(qmax<=0.)continue;
        int end_contact=0;Interval end=point(0.);
        if(node.gap[6].lo<=high.hi && node.gap[6].hi>=low.lo){
            end=circle_gap(a,b,node.hi);
            end_contact=end.lo>=low.hi && end.hi<=high.lo && node.derivative[5].lo>floor.hi;
        }
        if(qmin>floor.hi && node.gap[0].hi<0. && end_contact)
            return contact_window(result,node.hi,end,qmin,node.lo,node.hi,GFS_CIRCLE_DENSE_CONTACT_WINDOW);
        if(gmax<=0. && !end_contact)continue;
        if(qmin>floor.hi && node.gap[0].hi<0. && node.gap[6].lo>0.){
            /* Certified monotone dense sign bracket. Interval Horner
             * refinement avoids building a deep subdivision tree at a simple hit.
             */
            double lo=node.lo,hi=node.hi;
            for(int k=0;k<64;k++){
                double mid=lo+.5*(hi-lo);
                if(mid==lo||mid==hi)return unresolved(result,GFS_CIRCLE_DENSE_ROUNDING_LIMIT);
                Interval value=circle_gap(a,b,mid);result.refinement_evaluations++;
                if(value.lo>=low.hi && value.hi<=high.lo){
                    result.dense_sign_bracket_certified=1;
                    return contact_window(result,mid,value,qmin,node.lo,node.hi,GFS_CIRCLE_DENSE_SIGN_BRACKET);
                }
                if(value.hi<0.)lo=mid;
                else if(value.lo>0.)hi=mid;
                else return unresolved(result,GFS_CIRCLE_DENSE_ROUNDING_LIMIT);
            }
            return unresolved(result,GFS_CIRCLE_DENSE_ROUNDING_LIMIT);
        }
        if(gmin>=low.hi && gmax<=high.lo && qmin>floor.hi){
            result.status=GFS_CIRCLE_DENSE_CANDIDATE;
            result.reason=GFS_CIRCLE_DENSE_CONTACT_WINDOW;return result;
        }
        if(result.nodes>=o.max_nodes)return unresolved(result,GFS_CIRCLE_DENSE_NODE_LIMIT);
        if(node.depth>=o.max_depth)return unresolved(result,GFS_CIRCLE_DENSE_DEPTH_LIMIT);
        double mid=node.lo+.5*(node.hi-node.lo);
        if(mid==node.lo||mid==node.hi||count+2>50)return unresolved(result,GFS_CIRCLE_DENSE_ROUNDING_LIMIT);
        Node left,right;
        split(node.gap,6,left.gap,right.gap);
        split(node.derivative,5,left.derivative,right.derivative);
        left.lo=node.lo;left.hi=mid;left.depth=node.depth+1;
        right.lo=mid;right.hi=node.hi;right.depth=node.depth+1;
        stack[count++]=right;stack[count++]=left; /* left interval is processed first */
    }
    result.status=GFS_CIRCLE_DENSE_NONE;result.reason=GFS_CIRCLE_DENSE_EXCLUDED;
    result.theta_lo=0.;result.theta_hi=1.;return result;
}

GfsCircleDenseResult gfs_circle_dense_first_candidate(
    GfsWallCubic x,GfsWallCubic y,GfsCircleDenseOptions o){
#if defined(__x86_64__) || defined(__i386__)
    const unsigned int saved=_mm_getcsr();
    _mm_setcsr(saved & ~0x8040u); /* Gradual underflow inside this call only. */
    GfsCircleDenseResult result=first_candidate_strict(x,y,o);
    _mm_setcsr(saved);
    return result;
#else
    /* Unsupported arithmetic remains an explicit unresolved result. */
    return first_candidate_strict(x,y,o);
#endif
}
