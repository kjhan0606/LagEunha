#ifndef MATERIAL_CURVED_EVENTS_H
#define MATERIAL_CURVED_EVENTS_H
#include "material_curved_partition.h"
/* Event isolation for the preserved current map.  All exclusion tests are
 * interval tests. Root boxes are retained, including uncertainty-limited
 * boxes; no small overlap is discarded. This is an engineering interval
 * certificate with the same explicit native-evaluation model as bounds.h. */
namespace curved_events {
using R=long double;using I=curved_bounds::I;using Pt=std::array<R,2>;
using curved_bounds::absmax;
static R mid(I a){return a.lo+(a.hi-a.lo)/2;}
static R width(I a){return a.hi-a.lo;}
static bool excludes(I a){return a.lo>0 || a.hi<0;}
static bool intersect(I a,I b,I &out){R l=std::max(a.lo,b.lo),h=std::min(a.hi,b.hi);if(l>h)return false;out=I(l,h);return true;}
struct Box {I x[2];int depth=0;};
struct Jet {I v,d[2],dd[2][2];};
static I polynomial_derivative(const curved_compiled::C &p,I x,I y,int dx,int dy){
    I z(0);for(int a=0;a<10;a++)if(!curved_compiled::zero(p.a[a])){auto n=curved_compiled::exponents(a);if(n[0]<dx || n[1]<dy)continue;int factor=1;for(int j=0;j<dx;j++)factor*=n[0]-j;for(int j=0;j<dy;j++)factor*=n[1]-j;z=z+p.a[a]*I(factor)*curved_compiled::power(x,n[0]-dx)*curved_compiled::power(y,n[1]-dy);}return z;
}
static std::array<Jet,2> evaluate(const curved_compiled::Map &map,const Box &b,R eps){
    I X=b.x[0],Y=b.x[1],x=(X-I(map.anchor[0]))/I(map.h),y=(Y-I(map.anchor[1]))/I(map.h),coord[2]={X,Y};
    I z=curved_bounds::square(X)+curved_bounds::square(Y)+curved_bounds::square(I(eps)),s=curved_bounds::root(curved_bounds::root(z)),o=I(1)/(s*s*s),od[2],odd[2][2];
    for(int i=0;i<2;i++){od[i]=I(-1.5L)*coord[i]*o/z;for(int j=0;j<2;j++)odd[i][j]=I(-1.5L*(i==j))*o/z+I(5.25L)*coord[i]*coord[j]*o/(z*z);}
    I rot[2]={-Y,X};const int E[2][2]={{0,-1},{1,0}};Jet initial[2];
    for(int k=0;k<2;k++){initial[k].v=o*rot[k];for(int i=0;i<2;i++){initial[k].d[i]=od[i]*rot[k]+o*I(E[k][i]);for(int j=0;j<2;j++)initial[k].dd[i][j]=odd[i][j]*rot[k]+od[i]*I(E[k][j])+od[j]*I(E[k][i]);}}
    std::array<Jet,2> out;for(int k=0;k<2;k++){
        auto &q=out[k];q.v=curved_compiled::value(map.position[k],x,y);for(int i=0;i<2;i++){q.d[i]=polynomial_derivative(map.position[k],x,y,i==0,i==1)/I(map.h);for(int j=0;j<2;j++)q.dd[i][j]=polynomial_derivative(map.position[k],x,y,(i==0)+(j==0),(i==1)+(j==1))/I(map.h*map.h);}
        for(int a=0;a<2;a++)if(!curved_compiled::zero(map.position_B[k][a])){q.v=q.v+map.position_B[k][a]*initial[a].v;for(int i=0;i<2;i++){q.d[i]=q.d[i]+map.position_B[k][a]*initial[a].d[i];for(int j=0;j<2;j++)q.dd[i][j]=q.dd[i][j]+map.position_B[k][a]*initial[a].dd[i][j];}}
    }return out;
}
struct Plane {I n[2],c;int i=-1,j=-1;};
static Plane plane(const std::vector<curved_partition::Site> &site,const std::vector<R> &w,int i,int j){
    Plane p;p.i=i;p.j=j;p.n[0]=I(site[j].x)-I(site[i].x);p.n[1]=I(site[j].y)-I(site[i].y);p.c=(curved_bounds::square(I(site[j].x))+curved_bounds::square(I(site[j].y))-curved_bounds::square(I(site[i].x))-curved_bounds::square(I(site[i].y))+I(w[i])-I(w[j]))/I(2);return p;
}
static Jet projected(const std::array<Jet,2> &m,const Plane &p,R position_error){
    Jet z;z.v=curved_bounds::widen(p.n[0]*m[0].v+p.n[1]*m[1].v-p.c,position_error*(absmax(p.n[0])+absmax(p.n[1])));for(int i=0;i<2;i++){z.d[i]=p.n[0]*m[0].d[i]+p.n[1]*m[1].d[i];for(int j=0;j<2;j++)z.dd[i][j]=p.n[0]*m[0].dd[i][j]+p.n[1]*m[1].dd[i][j];}return z;
}
static bool source_possible(const curved::Record &r,const Box &b){
    I ax(r.triangle.geometry.x[0][0]),ay(r.triangle.geometry.x[0][1]),dx=b.x[0]-ax,dy=b.x[1]-ay,ux=I(r.triangle.geometry.x[1][0])-ax,uy=I(r.triangle.geometry.x[1][1])-ay,vx=I(r.triangle.geometry.x[2][0])-ax,vy=I(r.triangle.geometry.x[2][1])-ay,det=ux*vy-uy*vx;need(det.lo>0,"positive event reference triangle");I a=(dx*vy-dy*vx)/det,c=(ux*dy-uy*dx)/det;return a.hi>=0 && c.hi>=0 && (I(1)-a-c).hi>=0;
}
struct Context {
    const curved::Record &record;const std::vector<curved_partition::Site> &site;const std::vector<R> &weight;
    curved_compiled::Map map;Box domain;R eps,position_error,gradient_error;std::vector<int> candidate;
    Context(const curved::Record &r,R e,const std::vector<curved_partition::Site> &s,const std::vector<R> &w,R pe,R ge):record(r),site(s),weight(w),map(curved_compiled::compile(r,e)),eps(e),position_error(pe),gradient_error(ge){
        for(int k=0;k<2;k++){R lo=r.triangle.geometry.x[0][k],hi=lo;for(int a=1;a<3;a++){lo=std::min(lo,(R)r.triangle.geometry.x[a][k]);hi=std::max(hi,(R)r.triangle.geometry.x[a][k]);}domain.x[k]=I(lo,hi);}
        int own=-1;for(int j=0;j<(int)s.size();j++)if(s[j].id==r.cell)own=j;need(own>=0,"event source site");candidate.push_back(own);auto values=evaluate(map,domain,eps);for(int j=0;j<(int)s.size();j++)if(j!=own && projected(values,plane(s,w,own,j),pe).v.hi>=0)candidate.push_back(j);
    }
    bool active(const Box &b,int i)const{auto m=evaluate(map,b,eps);for(int k:candidate)if(k!=i && projected(m,plane(site,weight,i,k),position_error).v.lo>0)return false;return true;}
};
enum Kind {SOURCE_VERTEX,SOURCE_RADIAL_EXTREMUM,EDGE_CROSSING,POWER_VERTEX,RADIAL_TANGENCY,CHART_CROSSING};
struct Event {Box box;I radius;Kind kind;int i=-1,j=-1,k=-1;};
struct Isolation {std::vector<Event> event;uint64_t evaluations=0,excluded=0,contractions=0,uncertainty_limited=0,budget_failures=0;bool complete=true;};
static I radial(const Box &b){I z=curved_bounds::square(b.x[0])+curved_bounds::square(b.x[1]);z.lo=std::max(0.L,z.lo);return curved_bounds::root(z);}
static void retain(Isolation &out,const Box &b,Kind kind,int i,int j,int k=-1){out.event.push_back({b,radial(b),kind,i,j,k});}
struct System {Plane a,b;bool tangent=false;};
static void system_value(const Context &c,const System &s,const Box &b,I f[2],I J[2][2]){
    auto m=evaluate(c.map,b,c.eps);auto a=projected(m,s.a,c.position_error);f[0]=a.v;for(int k=0;k<2;k++)J[0][k]=a.d[k];
    if(!s.tangent){auto q=projected(m,s.b,c.position_error);f[1]=q.v;for(int k=0;k<2;k++)J[1][k]=q.d[k];}
    else{f[1]=-b.x[1]*a.d[0]+b.x[0]*a.d[1];f[1]=curved_bounds::widen(f[1],c.gradient_error*(absmax(s.a.n[0])+absmax(s.a.n[1]))*(absmax(b.x[0])+absmax(b.x[1])));J[1][0]=-b.x[1]*a.dd[0][0]+a.d[1]+b.x[0]*a.dd[1][0];J[1][1]=-a.d[0]-b.x[1]*a.dd[0][1]+b.x[0]*a.dd[1][1];}
}
static void isolate_system(const Context &c,const System &s,Kind kind,Isolation &out,int third=-1){
    std::vector<Box> stack={c.domain};uint64_t evaluations=0;const uint64_t cap=4096;const int depth_cap=48;
    while(!stack.empty()){
        Box b=stack.back();stack.pop_back();if(!source_possible(c.record,b)){out.excluded++;continue;}
        if(evaluations++>=cap){out.complete=false;out.budget_failures++;retain(out,b,kind,s.a.i,s.a.j,third);for(auto q:stack)retain(out,q,kind,s.a.i,s.a.j,third);return;}
        out.evaluations++;I f[2],J[2][2];system_value(c,s,b,f,J);if(excludes(f[0]) || excludes(f[1]) || !c.active(b,s.a.i)){out.excluded++;continue;}
        Box center=b;for(int k=0;k<2;k++)center.x[k]=I(mid(b.x[k]));I fc[2],unused[2][2];system_value(c,s,center,fc,unused);R a=mid(J[0][0]),d=mid(J[1][1]),u=mid(J[0][1]),v=mid(J[1][0]),det=a*d-u*v;bool contracted=false;
        if(det!=0 && std::isfinite(det)){
            R C[2][2]={{d/det,-u/det},{-v/det,a/det}};Box next=b;bool empty=false;
            for(int i=0;i<2;i++){I N=center.x[i];for(int k=0;k<2;k++)N=N-I(C[i][k])*fc[k];for(int j=0;j<2;j++){I E(i==j);for(int k=0;k<2;k++)E=E-I(C[i][k])*J[k][j];N=N+E*(b.x[j]-center.x[j]);}if(!intersect(b.x[i],N,next.x[i])){empty=true;break;}}
            if(empty){out.excluded++;continue;}
            R prior=std::max(width(b.x[0]),width(b.x[1])),after=std::max(width(next.x[0]),width(next.x[1]));
            if(after<.8L*prior){out.contractions++;next.depth=b.depth+1;if(next.depth<depth_cap){stack.push_back(next);contracted=true;}else{out.complete=false;out.budget_failures++;retain(out,next,kind,s.a.i,s.a.j,third);continue;}}
            else{
                /* With an invertible interval derivative, the retained set
                 * may be controlled by coefficient/native uncertainty.
                 * Keep that set as an event box; its physical error is not
                 * discarded and must meet the downstream mass/K budget. */
                R norm=0;for(int i=0;i<2;i++){R row=0;for(int j=0;j<2;j++){I E(i==j);for(int k=0;k<2;k++)E=E-I(C[i][k])*J[k][j];row+=absmax(E);}norm=std::max(norm,row);}
                if(norm<.5L && fc[0].lo<=0 && fc[0].hi>=0 && fc[1].lo<=0 && fc[1].hi>=0){out.uncertainty_limited++;retain(out,next,kind,s.a.i,s.a.j,third);continue;}
            }
        }
        if(contracted)continue;
        R tolerance=128*LDBL_EPSILON*(1+c.map.h+absmax(b.x[0])+absmax(b.x[1]));int axis=width(b.x[0])>=width(b.x[1])?0:1;
        if(width(b.x[axis])<=tolerance){retain(out,b,kind,s.a.i,s.a.j,third);continue;}
        if(b.depth>=depth_cap){out.complete=false;out.budget_failures++;retain(out,b,kind,s.a.i,s.a.j,third);continue;}
        R half=mid(b.x[axis]);Box left=b,right=b;left.x[axis].hi=half;right.x[axis].lo=half;left.depth=right.depth=b.depth+1;stack.push_back(right);stack.push_back(left);
    }
}
static Box segment_box(Pt a,Pt b,I t){Box q;for(int k=0;k<2;k++)q.x[k]=I(a[k])+t*(I(b[k])-I(a[k]));return q;}
static void isolate_edge(const Context &c,const Plane &p,Pt a,Pt b,Kind kind,Isolation &out){
    struct Segment {I t;int depth;};std::vector<Segment> stack={{I(0,1),0}};uint64_t evaluations=0;const uint64_t cap=1024;
    while(!stack.empty()){
        auto segment=stack.back();stack.pop_back();auto box=segment_box(a,b,segment.t);if(evaluations++>=cap){out.complete=false;out.budget_failures++;retain(out,box,kind,p.i,p.j);for(auto q:stack)retain(out,segment_box(a,b,q.t),kind,p.i,p.j);return;}
        out.evaluations++;auto m=evaluate(c.map,box,c.eps);auto f=projected(m,p,c.position_error);if(excludes(f.v) || !c.active(box,p.i)){out.excluded++;continue;}I derivative=f.d[0]*I(b[0]-a[0])+f.d[1]*I(b[1]-a[1]);
        if(excludes(derivative)){
            R center=mid(segment.t);auto fm=projected(evaluate(c.map,segment_box(a,b,I(center)),c.eps),p,c.position_error).v;I next;if(!intersect(segment.t,I(center)-fm/derivative,next)){out.excluded++;continue;}
            if(width(next)<.8L*width(segment.t) && segment.depth<48){out.contractions++;stack.push_back({next,segment.depth+1});continue;}
            if(fm.lo<=0 && fm.hi>=0 && absmax(derivative)<4*std::min(fabsl(derivative.lo),fabsl(derivative.hi))){out.uncertainty_limited++;retain(out,segment_box(a,b,next),kind,p.i,p.j);continue;}
        }
        if(width(segment.t)<=128*LDBL_EPSILON){retain(out,box,kind,p.i,p.j);continue;}if(segment.depth>=48){out.complete=false;out.budget_failures++;retain(out,box,kind,p.i,p.j);continue;}R center=mid(segment.t);stack.push_back({I(center,segment.t.hi),segment.depth+1});stack.push_back({I(segment.t.lo,center),segment.depth+1});
    }
}
static Isolation isolate(const Context &c){
    Isolation out;const auto &g=c.record.triangle.geometry;
    for(int a=0;a<3;a++){int b=(a+1)%3;Pt A={g.x[a][0],g.x[a][1]},B={g.x[b][0],g.x[b][1]};retain(out,segment_box(A,A,I(0)),SOURCE_VERTEX,-1,-1);I dx=I(B[0])-I(A[0]),dy=I(B[1])-I(A[1]),t=-(I(A[0])*dx+I(A[1])*dy)/(dx*dx+dy*dy),q;if(intersect(t,I(0,1),q))retain(out,segment_box(A,B,q),SOURCE_RADIAL_EXTREMUM,-1,-1);}
    /* Active power vertices require three sites sharing the same score.
     * Enumerating all candidate triples includes a target cell entirely
     * inside the source; source-edge intersections alone cannot do so. */
    for(size_t a=0;a<c.candidate.size();a++)for(size_t b=a+1;b<c.candidate.size();b++)for(size_t d=b+1;d<c.candidate.size();d++){
        int i=c.candidate[a],j=c.candidate[b],k=c.candidate[d];System s={plane(c.site,c.weight,i,j),plane(c.site,c.weight,i,k),false};
        I determinant=s.a.n[0]*s.b.n[1]-s.a.n[1]*s.b.n[0];if(excludes(determinant)){
            I x=(s.a.c*s.b.n[1]-s.a.n[1]*s.b.c)/determinant,y=(s.a.n[0]*s.b.c-s.a.c*s.b.n[0])/determinant;auto m=evaluate(c.map,c.domain,c.eps);I unused;if(!intersect(curved_bounds::widen(m[0].v,c.position_error),x,unused) || !intersect(curved_bounds::widen(m[1].v,c.position_error),y,unused)){out.excluded++;continue;}
            bool active=true;for(int n:c.candidate)if(n!=i){auto f=plane(c.site,c.weight,i,n);if((f.n[0]*x+f.n[1]*y-f.c).lo>0)active=false;}if(!active){out.excluded++;continue;}
        }
        isolate_system(c,s,POWER_VERTEX,out,k);
    }
    std::vector<Pt> chart;for(int a=0;a<3;a++){int b=(a+1)%3;Pt A={g.x[a][0],g.x[a][1]},B={g.x[b][0],g.x[b][1]};if(A[0]==0)chart.push_back(A);if((A[0]<0 && B[0]>0)||(A[0]>0 && B[0]<0)){R t=-A[0]/(B[0]-A[0]);chart.push_back({0,A[1]+t*(B[1]-A[1])});}}
    for(size_t a=0;a<c.candidate.size();a++)for(size_t b=a+1;b<c.candidate.size();b++){
        Plane p=plane(c.site,c.weight,c.candidate[a],c.candidate[b]);for(int edge=0;edge<3;edge++){int next=(edge+1)%3;isolate_edge(c,p,{g.x[edge][0],g.x[edge][1]},{g.x[next][0],g.x[next][1]},EDGE_CROSSING,out);}
        if(chart.size()>=2){auto cmp=[](Pt a,Pt b){return a[1]<b[1];};isolate_edge(c,p,*std::min_element(chart.begin(),chart.end(),cmp),*std::max_element(chart.begin(),chart.end(),cmp),CHART_CROSSING,out);}
        isolate_system(c,{p,{},true},RADIAL_TANGENCY,out);
    }return out;
}
static std::vector<curved_roots::Window> bands(const Isolation &r){std::vector<curved_roots::Window> v;for(auto e:r.event)v.push_back({e.radius.lo,e.radius.hi});return curved_partition::merge(v);}
}
#endif
