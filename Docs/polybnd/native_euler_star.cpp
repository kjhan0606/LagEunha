// Exact Euler pressure/speed batch using the installed SciPy Brent C API.
// Compile without fast math or contraction. No alternate solver/fallback.
#include <cmath>
#include <cstddef>
#include <limits>
#include <algorithm>

struct scipy_output { int funcalls; int iterations; int error_num; double root; };
using callback = double (*)(double, void*);
using brent_api = double (*)(callback,double,double,void*,double,double,int,scipy_output*);
struct problem {
    double rl,rr,pl,pr,ul,ur,gamma,cl,cr,scale,lo,hi,flo,fhi;
    int calls;
};
static double wave(double p,double rho,double pk,double c,double gamma) {
    if(p>pk) return (p-pk)*std::sqrt(2/((gamma+1)*rho)/(p+(gamma-1)/(gamma+1)*pk));
    return 2*c/(gamma-1)*(std::pow(p/pk,(gamma-1)/(2*gamma))-1);
}
static double residual(double q,void* opaque) {
    auto& s=*static_cast<problem*>(opaque);
    double p=q*s.scale;
    double f=wave(p,s.rl,s.pl,s.cl,s.gamma)+wave(p,s.rr,s.pr,s.cr,s.gamma)+s.ur-s.ul;
    ++s.calls;
    if(f<=0 && q>=s.lo) {s.lo=q;s.flo=f;}
    if(f>=0 && q<=s.hi) {s.hi=q;s.fhi=f;}
    return f;
}
static double derivative(double p,double rho,double pk,double c,double g) {
    if(p>pk) {
        double B=(g-1)/(g+1)*pk;
        return std::sqrt(2/((g+1)*rho)/(p+B))*(1-(p-pk)/(2*(p+B)));
    }
    return c/(g*pk)*std::pow(p/pk,-(g+1)/(2*g));
}
extern "C" void euler_star_batch(std::size_t count,const double* input,double* output,void* brent_pointer) {
    auto brent=reinterpret_cast<brent_api>(brent_pointer);
    const double eps=std::numeric_limits<double>::epsilon();
    for(std::size_t k=0;k<count;++k) {
        const double* a=input+7*k;double* out=output+13*k;
        std::fill(out,out+13,std::numeric_limits<double>::quiet_NaN());out[11]=0;
        problem s{};s.rl=a[0];s.ul=a[1];s.pl=a[2];s.rr=a[3];s.ur=a[4];s.pr=a[5];s.gamma=a[6];
        if(!(s.rl>0 && s.rr>0 && s.pl>0 && s.pr>0 && s.gamma>1) ||
           !std::isfinite(s.rl) || !std::isfinite(s.rr) || !std::isfinite(s.pl) || !std::isfinite(s.pr) || !std::isfinite(s.gamma) || !std::isfinite(s.ul+s.ur)) {out[11]=1;continue;}
        s.cl=std::sqrt(s.gamma*s.pl/s.rl);s.cr=std::sqrt(s.gamma*s.pr/s.rr);
        if(s.ur-s.ul>=2*(s.cl+s.cr)/(s.gamma-1)) {out[11]=2;continue;}
        if(!std::isfinite(s.cl+s.cr)) {out[11]=6;continue;}
        s.scale=std::max(s.pl,s.pr);s.lo=0;s.hi=std::numeric_limits<double>::infinity();s.flo=-INFINITY;s.fhi=INFINITY;
        double upper=1;int expansions=0;double f=residual(upper,&s);
        while(f<0) {
            upper*=2;++expansions;
            if(!std::isfinite(upper*s.scale)) break;
            f=residual(upper,&s);
        }
        if(!std::isfinite(upper*s.scale)) {out[11]=3;continue;}
        if(!std::isfinite(f)) {out[11]=6;continue;}
        scipy_output stats{};
        double q=brent(residual,0.,upper,&s,std::numeric_limits<double>::min(),8*eps,200,&stats);
        if(stats.error_num!=0) {out[11]=4;out[9]=stats.iterations;out[10]=s.calls;continue;}
        double p=q*s.scale;
        if(!(p>0)) {out[11]=5;continue;}
        double fl=wave(p,s.rl,s.pl,s.cl,s.gamma),fr=wave(p,s.rr,s.pr,s.cr,s.gamma);
        double speed=.5*(s.ul+s.ur+fr-fl);
        if(!std::isfinite(p+speed)) {out[11]=6;continue;}
        double pres=s.lo*s.scale,phi=s.hi*s.scale;
        double slope=derivative(phi,s.rl,s.pl,s.cl,s.gamma)+derivative(phi,s.rr,s.pr,s.cr,s.gamma);
        double round=64*eps*(std::abs(fl)+std::abs(fr)+std::abs(s.ul)+std::abs(s.ur)+s.cl+s.cr);
        double perror=(s.hi-s.lo)*s.scale+round/slope+16*eps*std::abs(p);
        double maxslope=derivative(pres,s.rl,s.pl,s.cl,s.gamma)+derivative(pres,s.rr,s.pr,s.cr,s.gamma);
        double serror=.5*(wave(phi,s.rl,s.pl,s.cl,s.gamma)-wave(pres,s.rl,s.pl,s.cl,s.gamma)+wave(phi,s.rr,s.pr,s.cr,s.gamma)-wave(pres,s.rr,s.pr,s.cr,s.gamma))+.5*maxslope*round/slope+round;
        out[0]=p;out[1]=speed;out[2]=s.lo;out[3]=s.hi;out[4]=s.flo;out[5]=s.fhi;
        out[6]=fl+fr+s.ur-s.ul;out[7]=perror;out[8]=serror;out[9]=stats.iterations;out[10]=s.calls;out[12]=expansions;
    }
}
