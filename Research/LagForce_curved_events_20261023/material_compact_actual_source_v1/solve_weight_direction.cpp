// Diagnostic Newton direction for the actual step5 power-area residual.
// The RHS uses paired-order Gauss areas without a certified integration
// error. No weight in this program is accepted for a physical remap.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using R=long double;
struct Edge {size_t a,b;R c;};
static void need(bool good,const char*why){if(!good)throw std::runtime_error(why);}
static std::vector<std::string> fields(const std::string&line){
    std::vector<std::string> out;std::stringstream stream(line);std::string s;
    while(std::getline(stream,s,','))out.push_back(s);return out;
}
static R number(const std::string&s){
    char*end=nullptr;R value=strtold(s.c_str(),&end);
    need(end&&*end==0&&std::isfinite(value),"finite complete CSV real");return value;
}
static R trial_weight(uint64_t id){
    uint64_t z=id+0x9e3779b97f4a7c15ULL;
    z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL;
    z=(z^(z>>27))*0x94d049bb133111ebULL;
    z^=z>>31;
    return (R(z)/R(UINT64_MAX)-R(.5L))*R(1e-18L);
}
static R dot(const std::vector<R>&a,const std::vector<R>&b){
    need(a.size()==b.size(),"same CG vector length");R sum=0;
    for(size_t i=0;i<a.size();i++)sum+=a[i]*b[i];return sum;
}
static void apply(const std::vector<Edge>&edge,const std::vector<R>&x,
        std::vector<R>&out){
    std::fill(out.begin(),out.end(),R(0));
    for(const auto&e:edge){R flux=e.c*(x[e.a]-x[e.b]);
        out[e.a]+=flux;out[e.b]-=flux;}
}
int main(int argc,char**argv){
    try{
        need(argc==4||argc==5,
            "usage: solve_weight_direction area.csv edges.csv output.csv [previous_weights.csv]");
        std::ifstream area(argv[1]);need(bool(area),"area ledger readable");
        std::string line;need(bool(std::getline(area,line)),"area header readable");
        std::vector<R>source,target,b;
        while(std::getline(area,line)){
            auto f=fields(line);need(f.size()==6,"six complete area ledger fields");
            size_t id=std::stoull(f[0]);need(id==source.size(),"sorted contiguous actual cell identities");
            R s=number(f[1]),a=number(f[3]);need(s>0&&a>0,"positive actual source and target cell areas");
            source.push_back(s);target.push_back(a);b.push_back(s-a);
        }
        need(source.size()==71080,"full actual step5 area ledger");
        std::vector<R>current(source.size());
        for(size_t i=0;i<current.size();i++)current[i]=trial_weight(i);
        if(argc==5){
            std::ifstream previous(argv[4]);need(bool(previous),"previous weight ledger readable");
            need(bool(std::getline(previous,line)),"previous weight header readable");
            for(size_t i=0;i<current.size();i++){
                need(bool(std::getline(previous,line)),"one previous weight per actual cell");
                auto f=fields(line);need(f.size()==6&&std::stoull(f[0])==i,
                    "sorted complete previous weight ledger");
                need(number(f[1])==trial_weight(i),
                    "previous weight ledger starts at same canonical trial");
                current[i]=number(f[3]);
            }
            need(!std::getline(previous,line),"no extra previous weight entry");
        }
        std::ifstream adjacency(argv[2]);need(bool(adjacency),"power adjacency ledger readable");
        need(bool(std::getline(adjacency,line)),"power adjacency header readable");
        std::vector<Edge>edge;std::vector<R>diag(source.size());
        while(std::getline(adjacency,line)){
            auto f=fields(line);need(f.size()==7,"seven complete adjacency fields");
            size_t a=std::stoull(f[0]),z=std::stoull(f[1]);R c=number(f[4]);
            need(a<z&&z<source.size()&&c>0,"positive sorted power-edge conductance");
            edge.push_back({a,z,c});diag[a]+=c;diag[z]+=c;
        }
        need(edge.size()==212217,"full actual step5 power adjacency ledger");
        for(R d:diag)need(d>0&&std::isfinite(d),"positive finite power Jacobian diagonal");
        R sum=0;for(R v:b)sum+=v;R mean=sum/R(b.size());
        for(R&v:b)v-=mean; // constant-weight gauge/nullspace compatibility
        R norm0=sqrtl(dot(b,b));need(norm0>0,"nonzero actual target-area residual");
        std::vector<R>x(b.size()),r=b,z(b.size()),p(b.size()),Ap(b.size());
        for(size_t i=0;i<r.size();i++)z[i]=r[i]/diag[i];p=z;
        R rho=dot(r,z);need(rho>0,"positive initial preconditioned residual");
        size_t iteration=0;R residual=norm0;
        for(;iteration<10000;iteration++){
            apply(edge,p,Ap);R denom=dot(p,Ap);
            need(denom>0&&std::isfinite(denom),"positive graph Laplacian CG denominator");
            R alpha=rho/denom;
            for(size_t i=0;i<x.size();i++){x[i]+=alpha*p[i];r[i]-=alpha*Ap[i];}
            residual=sqrtl(dot(r,r));
            if(residual<=norm0*R(1e-12L)){iteration++;break;}
            for(size_t i=0;i<r.size();i++)z[i]=r[i]/diag[i];
            R next=dot(r,z);need(next>0&&std::isfinite(next),
                "finite positive next preconditioned residual");
            R beta=next/rho;
            for(size_t i=0;i<p.size();i++)p[i]=z[i]+beta*p[i];
            rho=next;
        }
        need(iteration<10000,"finite power-area Newton-direction CG budget");
        R gauge=0;for(R v:x)gauge+=v;gauge/=R(x.size());
        for(R&v:x)v-=gauge;
        apply(edge,x,Ap);R true_norm2=0,max_equation_error=0;
        for(size_t i=0;i<x.size();i++){
            R err=Ap[i]-b[i];true_norm2+=err*err;
            max_equation_error=std::max(max_equation_error,fabsl(err));
        }
        R true_relative=sqrtl(true_norm2)/norm0;
        need(true_relative<1e-10L,"direct recomputed Newton equation residual");
        std::ofstream out(argv[3]);need(bool(out),"Newton direction output writable");
        out<<std::setprecision(25)
            <<"cell,initial_trial_weight,cumulative_newton_delta,unit_step_candidate_weight,source_area,target_gauss16\n";
        R max_delta=0,max_new=0;
        for(size_t i=0;i<x.size();i++){
            R initial=trial_weight(i),next=current[i]+x[i];
            R cumulative=next-initial;
            out<<i<<','<<initial<<','<<cumulative<<','<<next<<','<<source[i]<<','<<target[i]<<'\n';
            max_delta=std::max(max_delta,fabsl(x[i]));
            max_new=std::max(max_new,fabsl(next));
        }
        out.close();need(bool(out),"complete Newton direction output flush");
        std::cout<<std::setprecision(21)
            <<"ACTUAL_STEP5_NEWTON_DIRECTION cells="<<x.size()
            <<" edges="<<edge.size()
            <<" iterations="<<iteration
            <<" initial_rhs_l2="<<norm0
            <<" recomputed_relative_residual="<<true_relative
            <<" max_equation_error="<<max_equation_error
            <<" max_delta_weight="<<max_delta
            <<" max_unit_step_weight="<<max_new
            <<" mean_rhs_removed="<<mean
            <<" weights_accepted=0"<<std::endl;
        return 0;
    }catch(const std::exception&e){std::cerr<<"NEWTON_DIRECTION_FAIL "<<e.what()<<'\n';return 2;}
}
