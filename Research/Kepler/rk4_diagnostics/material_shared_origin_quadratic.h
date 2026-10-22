#ifndef MATERIAL_SHARED_ORIGIN_QUADRATIC_H
#define MATERIAL_SHARED_ORIGIN_QUADRATIC_H
#include "material_shared_collected_radius.h"
namespace shared_target {
struct OriginBranchProof {uint64_t method=0;S coefficient[3];Q root_distance_lower=0,root_gap=0,sqrt_c_lower=0,derivative_upper=0;};
// P(t)=|U+tV|^2=a+b*t+c*t^2 (the holomorphic polynomial, no conjugate t).
// If the t disk misses both conjugate roots, its positive-real-center
// square-root germ exists on the entire simply connected disk. Bound its
// derivative there; no principal-cut or positive-half-plane assumption
// on the *enclosure disk of P* is necessary in that certified case.
static S origin_quadratic_radius(const sd::Compiled&m,const sd::Chart&chart,S s,S t,const sd::Point&point,OriginBranchProof&proof,S*ray=nullptr,S*argument=nullptr){
 bool origin=true;for(int k=0;k<2;k++)origin=origin&&m.base.vertex[chart.apex][k].n==0&&m.base.vertex[chart.apex][k].e==0;
 if(!origin)return reference_radius(m,chart,s,t,point);
 auto actual=frame(m,chart,t);S square=actual.D[0]*actual.D[0]+actual.D[1]*actual.D[1];if(ray)for(int k=0;k<2;k++)ray[k]=actual.D[k];if(argument)*argument=square;
 if(sf::lower(square)>0)return s*complex_sqrt(square);
 proof.method=1;auto center_frame=frame(m,chart,S(0));S a,b,c;for(int k=0;k<2;k++){a=a+center_frame.D[k]*center_frame.D[k];b=b+S(2)*center_frame.D[k]*center_frame.E[k];c=c+center_frame.E[k]*center_frame.E[k];}proof.coefficient[0]=a;proof.coefficient[1]=b;proof.coefficient[2]=c;
 need(sf::lower(c)>0,"origin ray quadratic has positive leading coefficient");S discriminant=S(4)*a*c-b*b;need(sf::lower(discriminant)>0,"origin ray quadratic conjugate roots separated from real line");S alpha=-b/(S(2)*c),beta2=discriminant/(S(4)*c*c),offset=S(t.n)-alpha;need(sf::lower(beta2)>0,"origin ray conjugate root imaginary distance positive");S distance=mf::root(offset*offset+beta2);proof.root_distance_lower=sf::lower(distance);proof.root_gap=mf::FP<Q>::down(proof.root_distance_lower-t.e);need(proof.root_gap>0,"angular complex disk excludes both exact origin-ray branch roots");proof.sqrt_c_lower=sf::lower(mf::root(c));proof.derivative_upper=sf::add(sf::mag(b+S(2)*c*S(t.n)),sf::mul(sf::mul(Q(2),sf::mag(c)),t.e));Q denominator=mf::FP<Q>::down(mf::FP<Q>::down(Q(2)*proof.sqrt_c_lower)*proof.root_gap);need(denominator>0,"origin ray square-root derivative denominator certified");auto fixed=frame(m,chart,S(t.n));S center=mf::root(fixed.D[0]*fixed.D[0]+fixed.D[1]*fixed.D[1]);Q variation=sf::quotient(sf::mul(t.e,proof.derivative_upper),denominator);center.e=sf::add(center.e,variation);return s*center;
}
}
#endif
