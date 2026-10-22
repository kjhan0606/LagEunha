#ifndef MATERIAL_CURVED_FOLD_CHART_H
#define MATERIAL_CURVED_FOLD_CHART_H
#include "material_curved_branch_bounds.h"
#include "material_curved_topology.h"
namespace curved_fold_chart {
using R=long double;using I=curved_bounds::I;using Box=curved_events::Box;using Plane=curved_events::Plane;using W=curved_roots::Window;
struct Family {
    curved_compiled::Map map;Plane plane;Box critical;I radius,fr,ftt;
    int i=-1,j=-1,chart=0,reason=0;R eps=0,contraction=0,position_error=0,gradient_error=0,validation_width[2]={};
    bool unique=false,inside=false,outside=false,valid=false;
};
static Family certify(const curved::Record &r,const curved_branches::FixedMap &fixed,R eps,const std::vector<curved_partition::Site> &site,const std::vector<R> &weight,int i,int j,const curved_bounds::Enclosure &bound){
    Family out;out.i=i;out.j=j;out.eps=eps;out.position_error=bound.native_position_envelope;out.gradient_error=bound.native_gradient_envelope;
    curved_events::Context nominal(r,eps,site,weight,0,0),family(r,eps,site,weight,bound.native_position_envelope,bound.native_gradient_envelope);
    for(int k=0;k<2;k++){for(int a=0;a<10;a++)nominal.map.position[k].a[a]=I(fixed.position[k][a]);for(int d=0;d<2;d++)nominal.map.position_B[k][d]=I(fixed.B[k][d]);}
    auto scalar=curved_events::plane(site,weight,i,j);scalar.n[0]=I(site[j].x-site[i].x);scalar.n[1]=I(site[j].y-site[i].y);scalar.c=I((site[j].x*site[j].x+site[j].y*site[j].y-site[i].x*site[i].x-site[i].y*site[i].y+weight[i]-weight[j])/2);
    auto original=curved_events::plane(site,weight,i,j);curved_events::System ns={scalar,{},true},fs={original,{},true};R nx=scalar.n[0].lo,ny=scalar.n[1].lo,c=scalar.c.lo,n2=nx*nx+ny*ny;
    if(!(n2>0)){out.reason=1;return out;}R x[2]={c*nx/n2,c*ny/n2};
    for(int iteration=0;iteration<12;iteration++){
        Box p;p.x[0]=I(x[0]);p.x[1]=I(x[1]);I f[2],J[2][2];curved_events::system_value(nominal,ns,p,f,J);R a=curved_events::mid(J[0][0]),b=curved_events::mid(J[0][1]),v=curved_events::mid(J[1][0]),d=curved_events::mid(J[1][1]),det=a*d-b*v;
        if(det==0 || !std::isfinite(det)){out.reason=2;return out;}R f0=curved_events::mid(f[0]),f1=curved_events::mid(f[1]),dx=(d*f0-b*f1)/det,dy=(-v*f0+a*f1)/det;x[0]-=dx;x[1]-=dy;
        if(!std::isfinite(x[0]) || !std::isfinite(x[1])){out.reason=3;return out;}if(fabsl(dx)+fabsl(dy)<64*LDBL_EPSILON*(1+hypotl(x[0],x[1])))break;
    }
    Box box,center;for(int k=0;k<2;k++)center.x[k]=I(x[k]);I f[2],J[2][2],fc[2],unused[2][2];curved_events::system_value(family,fs,center,fc,unused);R a0=curved_events::mid(unused[0][0]),b0=curved_events::mid(unused[0][1]),v0=curved_events::mid(unused[1][0]),d0=curved_events::mid(unused[1][1]),det0=a0*d0-b0*v0;if(det0==0 || !std::isfinite(det0)){out.reason=4;return out;}R C0[2][2]={{d0/det0,-b0/det0},{-v0/det0,a0/det0}};
    /* One a-priori box from the actual residual uncertainty. No repeated
     * box-width search: strict Krawczyk inclusion remains mandatory. */
    for(int k=0;k<2;k++){R error=0;for(int z=0;z<2;z++)error+=fabsl(C0[k][z])*curved_bounds::absmax(fc[z]);out.validation_width[k]=std::max(1e-12L,8*error+512*LDBL_EPSILON*(1+fabsl(x[k])));if(!std::isfinite(out.validation_width[k])){out.reason=4;return out;}box.x[k]=I(x[k]-out.validation_width[k],x[k]+out.validation_width[k]);}curved_events::system_value(family,fs,box,f,J);R a=curved_events::mid(J[0][0]),b=curved_events::mid(J[0][1]),v=curved_events::mid(J[1][0]),d=curved_events::mid(J[1][1]),det=a*d-b*v;
    if(det==0 || !std::isfinite(det)){out.reason=4;return out;}R C[2][2]={{d/det,-b/det},{-v/det,a/det}};bool interior=true;Box K;
    for(int k=0;k<2;k++){K.x[k]=center.x[k];R row=0;for(int z=0;z<2;z++)K.x[k]=K.x[k]-I(C[k][z])*fc[z];for(int l=0;l<2;l++){I E(k==l);for(int z=0;z<2;z++)E=E-I(C[k][z])*J[z][l];row+=curved_bounds::absmax(E);K.x[k]=K.x[k]+E*(box.x[l]-center.x[l]);}out.contraction=std::max(out.contraction,row);interior=interior && K.x[k].lo>box.x[k].lo && K.x[k].hi<box.x[k].hi;}
    out.unique=interior && out.contraction<1;if(!out.unique){out.reason=5;return out;}out.critical=K;out.radius=curved_events::radial(K);out.map=family.map;out.plane=original;
    if(out.radius.lo<=0){out.reason=7;return out;}auto jet=curved_events::projected(curved_events::evaluate(family.map,K,eps),original,0);I X=K.x[0],Y=K.x[1];out.fr=(X*jet.d[0]+Y*jet.d[1])/out.radius;out.ftt=curved_bounds::square(Y)*jet.dd[0][0]-I(2)*X*Y*jet.dd[0][1]+curved_bounds::square(X)*jet.dd[1][1]-X*jet.d[0]-Y*jet.d[1];
    out.inside=true;out.outside=false;for(int k=0;k<3;k++){I predicate=fixed.source_exact_n[k][0]*X+fixed.source_exact_n[k][1]*Y-fixed.source_exact_c[k];out.inside=out.inside && predicate.lo>0;out.outside=out.outside || predicate.hi<0;}
    out.chart=X.lo>0?1:(X.hi<0?-1:((Y.lo>0 || Y.hi<0)?1:0));out.valid=out.chart!=0 && curved_events::excludes(out.fr) && curved_events::excludes(out.ftt) && (out.fr*out.ftt).hi<0;out.reason=out.valid?0:6;return out;
}
static Box polar_box(I radius,W u,int chart){I U(u.lo,u.hi),D=I(1)+curved_bounds::square(U);Box b;b.x[0]=I(chart)*radius*(I(1)-curved_bounds::square(U))/D;b.x[1]=I(2*chart)*radius*U/D;for(int k=0;k<2;k++){b.x[k].lo=std::max(b.x[k].lo,-radius.hi);b.x[k].hi=std::min(b.x[k].hi,radius.hi);}return b;}
struct Control {bool complete=false;int reason=0;W cone={-1,1};I critical;R min_radial=0,min_curvature=0,C=0,numerical=0;};
static Control control(const Family &f,const curved::Record &r,const curved_branches::FixedMap &map,const std::vector<curved_partition::Site> &site,const std::vector<R> &weight,R maximum){
    Control out;out.critical=f.radius;if(!f.valid){out.reason=1;return out;}if(maximum<=f.radius.hi){out.reason=6;return out;}W cone;if(!curved_topology::cone(r,f.chart,cone)){out.reason=2;return out;}
    I critical_u=I(f.chart)*f.critical.x[1]/(f.radius+I(f.chart)*f.critical.x[0]);R lo=std::min(cone.lo,critical_u.lo),hi=std::max(cone.hi,critical_u.hi),pad=(hi-lo)/16+1e-6L;
    /* This certificate neighborhood may cross u=+/-1. Quadrature charts
     * still partition [-1,1]; a chart-edge fold is never assumed to have an
     * interior angular critical point in that half-circle. */
    W neighborhood={lo-pad,hi+pad};I full_radius(f.radius.lo,maximum);R keep_lo=neighborhood.hi,keep_hi=neighborhood.lo;bool retained=false;
    for(int k=0;k<8;k++){R a=k==0?neighborhood.lo:neighborhood.lo+(neighborhood.hi-neighborhood.lo)*k/8,b=k==7?neighborhood.hi:neighborhood.lo+(neighborhood.hi-neighborhood.lo)*(k+1)/8;auto box=polar_box(full_radius,{a,b},f.chart);auto jet=curved_events::projected(curved_events::evaluate(f.map,box,f.eps),f.plane,f.position_error);I ft=-box.x[1]*jet.d[0]+box.x[0]*jet.d[1];ft=curved_bounds::widen(ft,f.gradient_error*(curved_bounds::absmax(f.plane.n[0])+curved_bounds::absmax(f.plane.n[1]))*(curved_bounds::absmax(box.x[0])+curved_bounds::absmax(box.x[1])));if(curved_events::excludes(jet.v) && curved_events::excludes(ft))continue;retained=true;keep_lo=std::min(keep_lo,a);keep_hi=std::max(keep_hi,b);}
    if(!retained){out.reason=7;return out;}out.cone={keep_lo,keep_hi};int orientation=f.ftt.lo>0?1:-1;R radial_min=1e99L,curvature_min=1e99L;
    /* Fixed inclusion grid, not a convergence/order/cap scan. It encloses
     * the entire neighborhood and radial path from the critical point. */
    for(int k=0;k<8;k++)for(int z=0;z<4;z++){R a=k==0?keep_lo:keep_lo+(keep_hi-keep_lo)*k/8,b=k==7?keep_hi:keep_lo+(keep_hi-keep_lo)*(k+1)/8,rl=z==0?f.radius.lo:f.radius.lo+(maximum-f.radius.lo)*z/4,rh=z==3?maximum:f.radius.lo+(maximum-f.radius.lo)*(z+1)/4;I radius(rl,rh);auto box=polar_box(radius,{a,b},f.chart);auto jet=curved_events::projected(curved_events::evaluate(f.map,box,f.eps),f.plane,0);I X=box.x[0],Y=box.x[1],fr=I(-orientation)*(X*jet.d[0]+Y*jet.d[1])/radius,ftt=I(orientation)*(curved_bounds::square(Y)*jet.dd[0][0]-I(2)*X*Y*jet.dd[0][1]+curved_bounds::square(X)*jet.dd[1][1]-X*jet.d[0]-Y*jet.d[1]);radial_min=std::min(radial_min,fr.lo);curvature_min=std::min(curvature_min,ftt.lo);}
    if(curvature_min<=0){out.reason=3;return out;}if(radial_min<=0){out.reason=4;return out;}
    for(int z=0;z<4;z++){I radius(z==0?f.radius.lo:f.radius.lo+(maximum-f.radius.lo)*z/4,z==3?maximum:f.radius.lo+(maximum-f.radius.lo)*(z+1)/4);auto angular=[&](R U){auto box=polar_box(radius,{U,U},f.chart);auto jet=curved_events::projected(curved_events::evaluate(f.map,box,f.eps),f.plane,0);return I(orientation)*(-box.x[1]*jet.d[0]+box.x[0]*jet.d[1]);};if(!(angular(keep_lo).hi<0 && angular(keep_hi).lo>0)){out.reason=5;return out;}}
    out.min_radial=radial_min;out.min_curvature=curvature_min;out.C=curved_bounds::root(I(2)*I(radial_min)*I(curvature_min)).lo;curved_branches::Curve curve={-1,-1,f.i,f.j};auto ray=curved_branch_bounds::representative_ray(r,map,full_radius,f.chart,f.eps);auto p=curved_branches::polynomial(ray,curve,site,weight);auto dp=curved_branches::derivative(p);R root_error=curved_bounds::absmax(curved_branches::hull(dp,out.cone))*128*LDBL_EPSILON*(1+std::max(fabsl(out.cone.lo),fabsl(out.cone.hi))),u_factor=powl(std::max({1.L,fabsl(out.cone.lo),fabsl(out.cone.hi)}),6.L);out.numerical=curved_bounds::up(curved_branch_bounds::arithmetic(map,curve,maximum,f.eps,site,weight)*u_factor+root_error);R expand=curved_bounds::up(2*out.numerical/out.min_radial);out.critical=I(std::max(0.L,curved_bounds::down(f.radius.lo-expand)),curved_bounds::up(f.radius.hi+expand));out.complete=out.C>0;return out;
}
/* Strong convexity and a nonzero radial derivative imply
 * |f_theta| >= C sqrt(r-H) for each root, with H an UPPER BOUND on all
 * relevant critical radii, including arithmetic uncertainty. Its reciprocal
 * is integrated analytically; no denominator floor or infinite node bound. */
static R reciprocal_measure(const Control &c,R lo,R hi){need(c.complete && lo>=c.critical.hi,"fold coarea starts beyond the full protected event box");R H=c.critical.hi,a=sqrtl(lo-H),b=sqrtl(hi-H),width=(hi-lo)/(a+b);return curved_bounds::up(2*width*(H+(a*a+a*b+b*b)/3)/c.C);}
}
#endif
