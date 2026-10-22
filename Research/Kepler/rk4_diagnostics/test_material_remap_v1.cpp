#include <cfloat>
#include <array>
#include <algorithm>
#define GFS_ENRICH_TEST_LIBRARY
#include "test_enriched_euler_v1.cpp"
using Real=long double;
struct P {Real x,y;};
using Poly=std::vector<P>;
static P sub(P a,P b){return {a.x-b.x,a.y-b.y};}
static Real dot(P a,P b){return a.x*b.x+a.y*b.y;}
static Real cross(P a,P b){return a.x*b.y-a.y*b.x;}
static bool equal(P a,P b){return a.x==b.x && a.y==b.y;}
static Poly clean(Poly p){Poly out;for(auto x:p)if(out.empty() || !equal(x,out.back()))out.push_back(x);if(out.size()>1 && equal(out.front(),out.back()))out.pop_back();if(out.size()<3)out.clear();return out;}
static Poly clip(Poly p,P a,P origin,Real bound){
    Poly out;if(p.empty())return out;
    for(size_t i=0;i<p.size();i++){
        P x=p[i],y=p[(i+1)%p.size()];Real dx=dot(a,sub(x,origin))-bound,dy=dot(a,sub(y,origin))-bound;bool ix=dx<=0,iy=dy<=0;
        if(ix)out.push_back(x);
        if(ix!=iy){Real t=dx/(dx-dy);P z=t==0?x:t==1?y:P{x.x+t*(y.x-x.x),x.y+t*(y.y-x.y)};out.push_back(z);}
    }return clean(out);
}
struct Moment {Real a=0,x=0,y=0,xx=0,xy=0,yy=0;};
static Moment moments(const Poly &p){
    Moment out;if(p.empty())return out;
    Real signed_area=0;for(size_t i=1;i+1<p.size();i++)signed_area+=cross(sub(p[i],p[0]),sub(p[i+1],p[0]))/2;
    /* A face/vertex-only intersection has exactly zero area, not an invalid
     * material cell. No small positive polygon is removed by a threshold. */
    if(signed_area==0)return out;need(signed_area>0,"positive oriented overlap area");
    P c={0,0};for(auto q:p){c.x+=q.x/p.size();c.y+=q.y/p.size();}
    for(size_t e=0;e<p.size();e++){
        P b=p[e],d=p[(e+1)%p.size()];if(equal(b,d))continue;Real area=cross(sub(b,c),sub(d,c))/2;
        if(!(area>0)){fprintf(stderr,"nonpositive common-refinement fan area %.21Lg\n",area);need(false,"positive overlap polygon without sliver cutoff");}
        Real sx=c.x+b.x+d.x,sy=c.y+b.y+d.y;out.a+=area;out.x+=area*sx/3;out.y+=area*sy/3;
        out.xx+=area*(c.x*c.x+b.x*b.x+d.x*d.x+c.x*b.x+c.x*d.x+b.x*d.x)/6;
        out.yy+=area*(c.y*c.y+b.y*b.y+d.y*d.y+c.y*b.y+c.y*d.y+b.y*d.y)/6;
        out.xy+=area*(sx*sy+c.x*c.y+b.x*b.y+d.x*d.y)/12;
    }return out;
}
struct Fragment {int source,target;Poly p;};
struct Domain {std::vector<P> q;std::vector<Fragment> pieces;std::vector<Real> m0,mu;};
struct Euler {Real m=0,px=0,py=0,E=0,L=0,U=0,K=0,area=0;};
static Euler plus(Euler a,Euler b){a.m+=b.m;a.px+=b.px;a.py+=b.py;a.E+=b.E;a.L+=b.L;a.U+=b.U;a.K+=b.K;a.area+=b.area;return a;}
static Euler integrate(const Poly &p){
    auto z=moments(p);const Real rho=1.3L,shear=1.2L,P0=1.e-6L,gamma=(double)1.666667f;
    Euler a;a.area=z.a;a.m=rho*z.a;a.px=rho*shear*z.y;a.K=.5L*rho*shear*shear*z.yy;
    a.U=P0*z.a/(gamma-1);a.E=a.K+a.U;a.L=-rho*shear*z.yy;return a;
}
static std::vector<Poly> target_cells(const Domain &domain,const std::vector<Real> &W){
    Real xmin=LDBL_MAX,xmax=-LDBL_MAX,ymin=LDBL_MAX,ymax=-LDBL_MAX;
    for(auto &f:domain.pieces)for(auto p:f.p){xmin=std::min(xmin,p.x);xmax=std::max(xmax,p.x);ymin=std::min(ymin,p.y);ymax=std::max(ymax,p.y);}
    Real pad=std::max(xmax-xmin,ymax-ymin);std::vector<Poly> out;
    for(size_t i=0;i<domain.q.size();i++){
        Poly p={{xmin-pad,ymin-pad},{xmax+pad,ymin-pad},{xmax+pad,ymax+pad},{xmin-pad,ymax+pad}};
        for(size_t j=0;j<domain.q.size() && !p.empty();j++)if(j!=i){P a=sub(domain.q[j],domain.q[i]);p=clip(p,a,domain.q[i],.5L*(dot(a,a)+W[i]-W[j]));}
        out.push_back(p);
    }return out;
}
static std::vector<Fragment> common_refinement(const Domain &domain,const std::vector<Real> &W,std::vector<Euler> &cell){
    std::vector<Fragment> out;cell.assign(domain.q.size(),Euler{});
    for(size_t i=0;i<domain.q.size();i++)for(auto &f:domain.pieces){
        Poly p=f.p;
        /* Clip with original power supports, not normals rebuilt from rounded
         * target vertices. Reconstructing those normals creates fake slivers
         * when source and target share an exactly coincident interface. */
        for(size_t j=0;j<domain.q.size() && !p.empty();j++)if(i!=j){
            P a=sub(domain.q[j],domain.q[i]);p=clip(p,a,domain.q[i],.5L*(dot(a,a)+W[i]-W[j]));
        }
        if(p.empty())continue;
        auto u=integrate(p);if(!(u.area>0))continue;cell[i]=plus(cell[i],u);out.push_back({f.source,(int)i,p});
    }return out;
}
static bool interval_clip(Real alpha,Real beta,Real &lo,Real &hi){
    if(alpha>0)hi=std::min(hi,beta/alpha);else if(alpha<0)lo=std::max(lo,beta/alpha);else if(beta<0)return false;return hi>lo;
}
/* Exact active power-face length inside the common physical domain. Source
 * fragments may share boundaries; union intervals before integrating so a
 * coincident face is counted once. No negative overlap/area fallback. */
static Real face_length(const Domain &domain,const std::vector<Real> &W,int i,int j){
    P a=sub(domain.q[j],domain.q[i]);Real d=sqrtl(dot(a,a));need(d>0,"distinct Laguerre sites");P n={a.x/d,a.y/d},t={-n.y,n.x};
    P x0={.5L*(domain.q[i].x+domain.q[j].x)+(W[i]-W[j])*n.x/(2*d),.5L*(domain.q[i].y+domain.q[j].y)+(W[i]-W[j])*n.y/(2*d)};
    Real lo=LDBL_MAX,hi=-LDBL_MAX;
    for(auto &f:domain.pieces)for(auto x:f.p){Real s=dot(sub(x,x0),t);lo=std::min(lo,s);hi=std::max(hi,s);}
    for(size_t k=0;k<domain.q.size();k++)if((int)k!=i && (int)k!=j){
        P c=sub(domain.q[k],domain.q[i]);Real bound=.5L*(dot(c,c)+W[i]-W[k])-dot(c,sub(x0,domain.q[i]));
        if(!interval_clip(dot(c,t),bound,lo,hi))return 0;
    }
    std::vector<std::pair<Real,Real>> intervals;
    for(auto &f:domain.pieces){
        Real l=lo,h=hi;bool ok=true;
        for(size_t e=0;e<f.p.size();e++){
            P p=f.p[e],q=f.p[(e+1)%f.p.size()],c={q.y-p.y,p.x-q.x};
            if(!interval_clip(dot(c,t),dot(c,sub(p,x0)),l,h)){ok=false;break;}
        }if(ok && h>l)intervals.push_back({l,h});
    }
    std::sort(intervals.begin(),intervals.end());Real total=0,last_lo=0,last_hi=0;bool first=true;
    for(auto z:intervals){if(first){last_lo=z.first;last_hi=z.second;first=false;}else if(z.first<=last_hi)last_hi=std::max(last_hi,z.second);else{total+=last_hi-last_lo;last_lo=z.first;last_hi=z.second;}}
    if(!first)total+=last_hi-last_lo;return total;
}
static std::vector<std::vector<Real>> jacobian(const Domain &domain,const std::vector<Real> &W){
    int n=W.size();std::vector<std::vector<Real>> J(n,std::vector<Real>(n));for(int i=0;i<n;i++)J[i][i]=domain.mu[i];
    for(int i=0;i<n;i++)for(int j=i+1;j<n;j++){
        P d=sub(domain.q[j],domain.q[i]);Real a=1.3L*face_length(domain,W,i,j)/(2*sqrtl(dot(d,d)));
        J[i][i]+=a;J[j][j]+=a;J[i][j]-=a;J[j][i]-=a;
    }return J;
}
static std::vector<Real> solve(const std::vector<std::vector<Real>> &J,std::vector<Real> b){
    int n=b.size();std::vector<std::vector<Real>> L(n,std::vector<Real>(n));
    for(int i=0;i<n;i++)for(int j=0;j<=i;j++){
        Real s=J[i][j];for(int k=0;k<j;k++)s-=L[i][k]*L[j][k];need(i!=j || s>0,"finite screened Jacobian SPD");L[i][j]=i==j?sqrtl(s):s/L[j][j];
    }
    for(int i=0;i<n;i++){for(int j=0;j<i;j++)b[i]-=L[i][j]*b[j];b[i]/=L[i][i];}
    for(int i=n-1;i>=0;i--){for(int j=i+1;j<n;j++)b[i]-=L[j][i]*b[j];b[i]/=L[i][i];}return b;
}
static Real residual(const Domain &d,const std::vector<Real> &W,std::vector<Euler> &cell,std::vector<Fragment> &fragments,std::vector<Real> &F){
    fragments=common_refinement(d,W,cell);F.resize(W.size());Real error=0;
    for(size_t i=0;i<W.size();i++){F[i]=cell[i].m+d.mu[i]*W[i]-d.m0[i];if(!(cell[i].m>0))return LDBL_MAX;error=std::max(error,fabsl(F[i]/d.m0[i]));}return error;
}
struct Remap {std::vector<Real> W;std::vector<Euler> cells;std::vector<Fragment> fragments;int steps=0,backtracks=0;Real error=0,jac_error=0,jac_steps[3]={0},jac_errors[3]={0};};
static Remap remap(const Domain &d,std::vector<Real> W){
    Remap result;std::vector<Real> F;Real error=residual(d,W,result.cells,result.fragments,F);
    for(int it=0;error>1.e-12L && it<30;it++){
        auto J=jacobian(d,W);std::vector<Real> rhs=F;for(auto &v:rhs)v=-v;auto step=solve(J,rhs);bool accepted=false;
        for(int line=0;line<20;line++){
            Real alpha=ldexpl(1.L,-line);auto trial=W;for(size_t i=0;i<W.size();i++)trial[i]+=alpha*step[i];
            std::vector<Euler> cell;std::vector<Fragment> fragments;std::vector<Real> tf;Real e=residual(d,trial,cell,fragments,tf);
            if(e<error){W=trial;error=e;F=tf;result.cells=cell;result.fragments=fragments;result.backtracks+=line;accepted=true;break;}
        }need(accepted,"finite remap Newton descent, no mass-transfer deletion");result.steps++;
    }need(error<1.e-12L,"finite remap mass-weight convergence");result.error=error;result.W=W;
    auto J=jacobian(d,W);std::vector<Real> direction(W.size());Real scale=LDBL_MAX;
    for(size_t i=0;i<W.size();i++){direction[i]=sinl((Real)i+.37L);scale=std::min(scale,d.m0[i]/d.mu[i]);}
    for(int test=0;test<3;test++){
        Real h=ldexpl(1.e-5L*scale,-test);auto wp=W,wm=W;for(size_t i=0;i<W.size();i++){wp[i]+=h*direction[i];wm[i]-=h*direction[i];}
        std::vector<Euler> cp,cm;std::vector<Fragment> fp,fm;std::vector<Real> rp,rm;residual(d,wp,cp,fp,rp);residual(d,wm,cm,fm,rm);
        for(size_t i=0;i<W.size();i++){Real exact=0;for(size_t j=0;j<W.size();j++)exact+=J[i][j]*direction[j];result.jac_errors[test]=std::max(result.jac_errors[test],fabsl((rp[i]-rm[i])/(2*h)-exact)/d.mu[i]);}
        result.jac_steps[test]=h;
    }result.jac_error=result.jac_errors[0];
    return result;
}
static Domain initial(const std::map<uint64_t,Cell> &cells){
    Domain d;int i=0;for(auto &entry:cells){
        auto &c=entry.second;Poly p;for(size_t k=0;k<c.x.size();k++)p.push_back({(Real)c.gx-12+c.x[k],(Real)c.gy-12+c.y[k]});p=clean(p);
        d.q.push_back({(Real)c.gx-12,(Real)c.gy-12});d.pieces.push_back({i,i,p});auto u=integrate(p);d.m0.push_back(u.m);d.mu.push_back(u.m/u.area);i++;
    }return d;
}
static Domain deform(Domain d,Real strain){
    for(auto &q:d.q)q.x+=strain*q.y;for(auto &f:d.pieces)for(auto &x:f.p)x.x+=strain*x.y;return d;
}
static Euler total(const std::vector<Fragment> &p){Euler sum;for(auto &f:p)sum=plus(sum,integrate(f.p));return sum;}
static void report(const char *name,const Remap &r,const Euler &base){
    auto now=total(r.fragments);Real maxW=0,gauge=0,masschange=0,Kloss=0,meanL=0;
    for(size_t i=0;i<r.W.size();i++){
        maxW=std::max(maxW,fabsl(r.W[i]));gauge+=1.3L*r.W[i];
        auto &u=r.cells[i];Kloss+=u.K-.5L*(u.px*u.px+u.py*u.py)/u.m;
        /* Uniform rho, affine shear: exact intrinsic angular contribution
         * removed by replacing this field by only its mean momentum. */
        meanL+=u.L+2*(u.K-.5L*u.px*u.px/u.m)/1.2L;
    }
    printf("{\"name\":\"%s\",\"Newton_steps\":%d,\"backtracks\":%d,\"fragments\":%zu,\"max_mass_weight_relative_residual\":%.17Lg,\"max_abs_W\":%.17Lg,\"sum_mu_W\":%.17Lg,\"finite_Jacobian_FD_error_over_mu\":%.17Lg,\"mass_change\":%.17Lg,\"px_change\":%.17Lg,\"energy_change_over_thermal\":%.17Lg,\"angular_change\":%.17Lg,\"thermal_change_relative\":%.17Lg,\"kinetic_change_over_thermal\":%.17Lg,\"hypothetical_mean_velocity_projection_K_loss_over_thermal\":%.17Lg,\"hypothetical_mean_only_angular_error\":%.17Lg,\"actual_velocity_compression\":\"none: affine field retained exactly on common-refinement fragments\",\"Jacobian_FD\":[",
        name,r.steps,r.backtracks,r.fragments.size(),r.error,maxW,gauge,r.jac_error,now.m-base.m,now.px-base.px,(now.E-base.E)/base.U,now.L-base.L,(now.U-base.U)/base.U,(now.K-base.K)/base.U,Kloss/base.U,meanL-now.L);
    for(int i=0;i<3;i++){if(i)printf(",");printf("{\"step\":%.17Lg,\"error_over_mu\":%.17Lg}",r.jac_steps[i],r.jac_errors[i]);}
    printf("],\"weights\":[");for(size_t i=0;i<r.W.size();i++){if(i)printf(",");printf("%.21Lg",r.W[i]);}
    printf("],\"remapped_cell_mass\":[");for(size_t i=0;i<r.cells.size();i++){if(i)printf(",");printf("%.21Lg",r.cells[i].m);}printf("]}");
}
int main(int argc,char **argv){
    need(argc==2,"actual native polygon packets");std::set<uint64_t> ids={17894,20502,39878,44945,67927,67934,67935,67942,67943,67992,67993,68000,68001,68002,68525,68526,68527,68534,68535,68584,68585,68592,68593,68600};
    auto cells=load(argv[1],ids);auto start=initial(cells);auto base=total(start.pieces);std::vector<Real> W(cells.size());
    auto identity=remap(start,W);auto material=deform(start,.125L);auto forward=remap(material,W);
    Domain returned=material;returned.pieces=forward.fragments;returned=deform(returned,-.125L);auto backward=remap(returned,forward.W);
    printf("{\"scope\":\"24 actual polygons, manufactured uniform-density constant-P affine shear; compatible exact material control plus finite mass-weight common-refinement remap and inverse deformation; no curved-Kepler remap or runtime compression\",\"rho\":1.3,\"P\":1e-6,\"velocity\":\"(1.2*y,0)\",\"finite_shear_strain\":0.125,\"reports\":[");
    report("identity_W0",identity,base);printf(",");report("finite_shear_to_Laguerre",forward,base);printf(",");report("inverse_deformation_return",backward,base);printf("]}\n");return 0;
}
