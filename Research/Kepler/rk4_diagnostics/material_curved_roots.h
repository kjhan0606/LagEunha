#ifndef MATERIAL_CURVED_ROOTS_H
#define MATERIAL_CURVED_ROOTS_H
#include "material_curved_bounds.h"
#include "material_curved_compiled.h"
namespace curved_roots {
using R=long double;using I=curved_bounds::I;using curved_bounds::absmax;
struct P {int degree;std::array<I,7> c;explicit P(int n=0):degree(n),c{} {}};
static P constant(I x){P p;p.c[0]=x;return p;}
static P add(const P &a,const P &b){P p(std::max(a.degree,b.degree));for(int i=0;i<=p.degree;i++)p.c[i]=a.c[i]+b.c[i];return p;}
static P scale(const P &a,I x){P p(a.degree);for(int i=0;i<=a.degree;i++)p.c[i]=a.c[i]*x;return p;}
static P sub(const P &a,const P &b){return add(a,scale(b,I(-1)));}
static P mul(const P &a,const P &b){need(a.degree+b.degree<=6,"interval polynomial degree, no truncation");P p(a.degree+b.degree);for(int i=0;i<=a.degree;i++)for(int j=0;j<=b.degree;j++)p.c[i+j]=p.c[i+j]+a.c[i]*b.c[j];return p;}
static I value(const P &a,I x){I y=a.c[a.degree];for(int i=a.degree-1;i>=0;i--)y=y*x+a.c[i];return y;}
static int choose(int n,int k){if(k<0 || k>n)return 0;int a=1;for(int j=1;j<=k;j++)a=a*(n-j+1)/j;return a;}
static I power(I a,int n){I p(1);for(int i=0;i<n;i++)p=p*a;return p;}
using Bernstein=std::array<I,7>;
static Bernstein bernstein(const P &p,R left,R right){
    I a(left),w=I(right)-a;std::array<I,7> translated={};Bernstein b={};
    for(int j=0;j<=p.degree;j++)for(int i=j;i<=p.degree;i++)translated[j]=translated[j]+p.c[i]*I(choose(i,j))*power(a,i-j)*power(w,j);
    for(int k=0;k<=p.degree;k++)for(int j=0;j<=k;j++)b[k]=b[k]+translated[j]*I(choose(k,j))/I(choose(p.degree,j));return b;
}
static int sign(const Bernstein &b,int n){bool positive=true,negative=true;for(int i=0;i<=n;i++){positive=positive && b[i].lo>0;negative=negative && b[i].hi<0;}return positive?1:(negative?-1:0);}
static int variation(const Bernstein &b,int n){int prior=0,count=0;for(int i=0;i<=n;i++){int s=b[i].lo>0?1:(b[i].hi<0?-1:0);if(!s)return -1;if(prior && s!=prior)count++;prior=s;}return count;}
struct Window {R lo,hi;};
struct Isolation {std::vector<Window> window;uint64_t panels=0,positive=0,negative=0,uncertainty_limited=0;bool budget_complete=true;R uncertain_length=0;};
static void split_at(const Bernstein &b,int n,I t,Bernstein &left,Bernstein &right){
    t.lo=std::max(0.L,t.lo);t.hi=std::min(1.L,t.hi);Bernstein work=b;left[0]=work[0];right[n]=work[n];for(int k=1;k<=n;k++){for(int i=0;i<=n-k;i++)work[i]=work[i]*(I(1)-t)+work[i+1]*t;left[k]=work[0];right[n-k]=work[n-k];}
}
static Bernstein restrict_to(const Bernstein &b,int n,R a,R z,R left,R right){
    Bernstein l={},r={},out={},unused={};if(left==a)r=b;else split_at(b,n,(I(left)-I(a))/(I(z)-I(a)),l,r);if(right==z)return r;split_at(r,n,(I(right)-I(left))/(I(z)-I(left)),out,unused);return out;
}
static void recurse(const Bernstein &b,const Bernstein &db,int n,R a,R z,R tolerance,int depth,int cap,uint64_t budget,Isolation &out){
    out.panels++;int s=sign(b,n);if(s){if(s>0)out.positive++;else out.negative++;return;}
    if(z-a<=tolerance || depth==cap || out.panels>=budget){out.window.push_back({a,z});if((depth==cap && z-a>tolerance) || out.panels>=budget)out.budget_complete=false;return;}
    if(n==0){out.window.push_back({a,z});out.uncertainty_limited++;return;}
    /* Interval Newton contracts the set of roots of every coefficient
     * polynomial in the enclosure. Once coefficient uncertainty controls
     * that set, retain the whole band instead of tiling it into millions of
     * tiny intervals. This preserves uncertainty; it does not call it zero. */
    if(n>0){
        /* Carry derivative coefficients independently. Differencing two
         * uncertain local Bernstein values and dividing by a shrinking
         * width would invent an uncertainty proportional to1/width. */
        I derivative=db[0];for(int k=1;k<n;k++){derivative.lo=std::min(derivative.lo,db[k].lo);derivative.hi=std::max(derivative.hi,db[k].hi);}
        Bernstein l={},r={};split_at(b,n,I(.5L),l,r);I midpoint=l[n];R center=(a+z)/2;
        if(derivative.lo>0 || derivative.hi<0){
            I next=I(center)-midpoint/derivative;R lo=std::max(a,next.lo),hi=std::min(z,next.hi);if(lo>hi)return;
            if(hi-lo<.8L*(z-a) && lo<hi){recurse(restrict_to(b,n,a,z,lo,hi),restrict_to(db,n-1,a,z,lo,hi),n,lo,hi,tolerance,depth+1,cap,budget,out);return;}
            R smallest=std::min(fabsl(derivative.lo),fabsl(derivative.hi));if(midpoint.lo<=0 && midpoint.hi>=0 && absmax(derivative)<4*smallest){out.window.push_back({a,z});out.uncertainty_limited++;return;}
        }
        if(midpoint.lo<=0 && midpoint.hi>=0 && absmax(derivative)*(z-a)<=4*(midpoint.hi-midpoint.lo)){out.window.push_back({a,z});out.uncertainty_limited++;return;}
    }
    Bernstein work=b,left={},right={};left[0]=work[0];right[n]=work[n];for(int k=1;k<=n;k++){for(int i=0;i<=n-k;i++)work[i]=(work[i]+work[i+1])/I(2);left[k]=work[0];right[n-k]=work[n-k];}R m=(a+z)/2;
    Bernstein dl={},dr={};split_at(db,n-1,I(.5L),dl,dr);recurse(left,dl,n,a,m,tolerance,depth+1,cap,budget,out);recurse(right,dr,n,m,z,tolerance,depth+1,cap,budget,out);
}
static Isolation isolate(const P &p,R a,R b,R tolerance=1e-11L,int cap=42,uint64_t budget=16384){
    need(a<b && tolerance>0,"finite angular root chart");Isolation out;auto coeff=bernstein(p,a,b);P derivative(std::max(0,p.degree-1));for(int i=1;i<=p.degree;i++)derivative.c[i-1]=p.c[i]*I(i);auto dcoeff=bernstein(derivative,a,b);recurse(coeff,dcoeff,p.degree,a,b,tolerance,0,cap,budget,out);std::vector<Window> merged;
    for(auto w:out.window){if(!merged.empty() && w.lo<=merged.back().hi)merged.back().hi=std::max(merged.back().hi,w.hi);else merged.push_back(w);}out.window.swap(merged);for(auto w:out.window)out.uncertain_length+=w.hi-w.lo;return out;
}
struct Ray {P denominator,position[2],velocity[2],bary[3];};
/* The physical reference domain is the triangle of the stored vertices.
 * Its predicates use an interval inverse of those edges. The stored double
 * gradients remain the map's basis coefficients; using them as exact domain
 * predicates would silently move an edge by their construction roundoff. */
static void geometric_barycentric(const curved::Record &r,const P &X,const P &Y,const P &D,P out[3]){
    I ax(r.triangle.geometry.x[0][0]),ay(r.triangle.geometry.x[0][1]),ux=I(r.triangle.geometry.x[1][0])-ax,uy=I(r.triangle.geometry.x[1][1])-ay,vx=I(r.triangle.geometry.x[2][0])-ax,vy=I(r.triangle.geometry.x[2][1])-ay,det=ux*vy-uy*vx;need(det.lo>0,"positive exact-vertex reference domain");P dx=sub(X,scale(D,ax)),dy=sub(Y,scale(D,ay));out[1]=scale(sub(scale(dx,vy),scale(dy,vx)),I(1)/det);out[2]=scale(sub(scale(dy,ux),scale(dx,uy)),I(1)/det);out[0]=sub(sub(D,out[1]),out[2]);
}
/* Two exact charts, centered at0 andpi, cover the circle with u in[-1,1].
 * Trigonometric constants are exactly0,+1,-1. Interval radius permits later
 * radial-panel enclosures, including coefficient arithmetic uncertainty. */
static Ray ray_direct(const curved::Record &r,I radius,int chart,R eps){
    need(chart==1 || chart==-1,"exact half-circle chart");Ray out;P D(2);D.c[0]=I(1);D.c[2]=I(1);out.denominator=D;P D2=mul(D,D),D3=mul(D2,D),X(2),Y(2);X.c[0]=radius*I(chart);X.c[2]=-radius*I(chart);Y.c[1]=radius*I(2*chart);
    I ax(r.triangle.geometry.x[0][0]),ay(r.triangle.geometry.x[0][1]);I za=curved_bounds::square(ax)+curved_bounds::square(ay)+curved_bounds::square(I(eps));P dx=sub(X,scale(D,ax)),dy=sub(Y,scale(D,ay));out.bary[1]=add(scale(dx,I(r.triangle.geometry.grad[1][0])),scale(dy,I(r.triangle.geometry.grad[1][1])));out.bary[2]=add(scale(dx,I(r.triangle.geometry.grad[2][0])),scale(dy,I(r.triangle.geometry.grad[2][1])));out.bary[0]=sub(sub(D,out.bary[1]),out.bary[2]);P shape[10];
    for(int a=0;a<10;a++){shape[a]=constant(I(1));int index[3];gfs_material_p3_index(a,index);for(int k=0;k<3;k++)for(int j=0;j<index[k];j++)shape[a]=mul(shape[a],scale(sub(scale(out.bary[k],I(3)),scale(D,I(j))),I(1)/I(index[k]-j)));}
    P L=scale(add(scale(dx,I(2)*ax),scale(dy,I(2)*ay)),I(1)/za),Q=scale(add(mul(dx,dx),mul(dy,dy)),I(1)/za),L2=mul(L,L),L3=mul(L2,L);
    P A3=add(add(add(D3,scale(mul(L,D2),I(-.75L))),scale(mul(Q,D),I(-.75L))),add(add(scale(mul(L2,D),I(.65625L)),scale(mul(L,Q),I(1.3125L))),scale(L3,I(-.6015625L))));
    P A2=add(add(D2,scale(mul(L,D),I(-.75L))),add(scale(Q,I(-.75L)),scale(L2,I(.65625L))));I w0=curved_bounds::root(curved_bounds::root(za)),omega0=I(1)/(w0*w0*w0);
    P T3[2]={scale(add(scale(A3,-ay),mul(scale(dy,I(-1)),A2)),omega0),scale(add(scale(A3,ax),mul(dx,A2)),omega0)};
    I z=curved_bounds::square(radius)+curved_bounds::square(I(eps)),w=curved_bounds::root(curved_bounds::root(z)),omega=I(1)/(w*w*w);P v0[2]={scale(mul(scale(Y,I(-1)),D2),omega),scale(mul(X,D2),omega)},residual[2]={sub(v0[0],T3[0]),sub(v0[1],T3[1])};
    for(int k=0;k<2;k++)for(int a=0;a<10;a++)residual[k]=sub(residual[k],scale(shape[a],I(r.triangle.remainder[a][k])));
    for(int k=0;k<2;k++){
        out.position[k]=scale(D3,I(r.origin[k]));out.velocity[k]=constant(I(0));for(int a=0;a<10;a++){out.position[k]=add(out.position[k],scale(shape[a],I(r.map.x[k][a])));out.velocity[k]=add(out.velocity[k],scale(shape[a],I(r.map.v[k][a])));}
        for(int j=0;j<2;j++){out.position[k]=add(out.position[k],scale(residual[j],I(r.map.x[k][10+j])/I(r.scale[j])));out.velocity[k]=add(out.velocity[k],scale(residual[j],I(r.map.v[k][10+j])/I(r.scale[j])));}
    }geometric_barycentric(r,X,Y,D,out.bary);return out;
}
static P pow_poly(const P &p,int n){P q=constant(I(1));for(int j=0;j<n;j++)q=mul(q,p);return q;}
static Ray ray(const curved::Record &r,I radius,int chart,R eps,const curved_compiled::Map *cached=nullptr){
    need(chart==1 || chart==-1,"compiled exact half-circle chart");curved_compiled::Map local;if(!cached){local=curved_compiled::compile(r,eps);cached=&local;}const auto &map=*cached;Ray out;P D(2);D.c[0]=I(1);D.c[2]=I(1);out.denominator=D;P D2=mul(D,D),D3=mul(D2,D),X(2),Y(2);X.c[0]=radius*I(chart);X.c[2]=-radius*I(chart);Y.c[1]=radius*I(2*chart);P dx=sub(X,scale(D,I(map.anchor[0]))),dy=sub(Y,scale(D,I(map.anchor[1]))),u=scale(dx,I(1)/I(map.h)),v=scale(dy,I(1)/I(map.h));
    out.bary[1]=add(scale(dx,I(r.triangle.geometry.grad[1][0])),scale(dy,I(r.triangle.geometry.grad[1][1])));out.bary[2]=add(scale(dx,I(r.triangle.geometry.grad[2][0])),scale(dy,I(r.triangle.geometry.grad[2][1])));out.bary[0]=sub(sub(D,out.bary[1]),out.bary[2]);
    auto numerator=[&](const curved_compiled::C &p){P sum=constant(I(0));for(int a=0;a<10;a++)if(!curved_compiled::zero(p.a[a])){auto e=curved_compiled::exponents(a);sum=add(sum,scale(mul(mul(pow_poly(u,e[0]),pow_poly(v,e[1])),pow_poly(D,3-e[0]-e[1])),p.a[a]));}return sum;};
    I z=curved_bounds::square(radius)+curved_bounds::square(I(eps)),w=curved_bounds::root(curved_bounds::root(z)),omega=I(1)/(w*w*w);P initial[2]={scale(mul(scale(Y,I(-1)),D2),omega),scale(mul(X,D2),omega)};
    for(int k=0;k<2;k++){out.position[k]=numerator(map.position[k]);out.velocity[k]=numerator(map.velocity[k]);for(int j=0;j<2;j++){out.position[k]=add(out.position[k],scale(initial[j],map.position_B[k][j]));out.velocity[k]=add(out.velocity[k],scale(initial[j],map.velocity_B[k][j]));}}geometric_barycentric(r,X,Y,D,out.bary);return out;
}
static P plane(const Ray &r,I ax,I ay,I b,R physical_position_uncertainty=0){
    P D3=mul(mul(r.denominator,r.denominator),r.denominator);P out=sub(add(scale(r.position[0],ax),scale(r.position[1],ay)),scale(D3,b));
    if(physical_position_uncertainty){R delta=physical_position_uncertainty*(absmax(ax)+absmax(ay));out=add(out,scale(D3,I(-delta,delta)));}return out;
}
static std::vector<Window> source_intervals(const Ray &r,R &uncertain_length,uint64_t &panels){
    std::vector<Window> uncertain;for(int k=0;k<3;k++){auto roots=isolate(r.bary[k],-1,1);panels+=roots.panels;need(roots.budget_complete,"source angular root isolation fixed budget");uncertain.insert(uncertain.end(),roots.window.begin(),roots.window.end());}
    std::sort(uncertain.begin(),uncertain.end(),[](Window a,Window b){return a.lo<b.lo;});std::vector<Window> merged;
    for(auto w:uncertain){if(!merged.empty() && w.lo<=merged.back().hi)merged.back().hi=std::max(merged.back().hi,w.hi);else merged.push_back(w);}
    std::vector<Window> result;R left=-1;auto consider=[&](R a,R b){if(b<=a)return;R mid=(a+b)/2;bool inside=true;for(int k=0;k<3;k++){I v=value(r.bary[k],I(mid));need(v.lo>0 || v.hi<0,"source interval classification separated from roots");inside=inside && v.lo>0;}if(inside)result.push_back({a,b});};
    for(auto w:merged){consider(left,w.lo);uncertain_length+=w.hi-w.lo;left=w.hi;}consider(left,1);return result;
}
}
#endif
