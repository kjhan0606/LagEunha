#ifndef MATERIAL_SHARED_LOGISTIC_DENSITY_H
#define MATERIAL_SHARED_LOGISTIC_DENSITY_H
#include "material_shared_atlas.h"
namespace shared_target {
struct LogisticDisk {S argument,exponential,denominator,value;uint64_t reflected=0,increasing=0;};
// For z=c+w, |w|<=R, use E=exp(-z) when c>0, exp(z) otherwise.
// E is enclosed by the original holomorphic_exp and native libm premise.
// lower(1+E)>0 excludes denominator zeros throughout this same disk.
// exp(z)/(1+exp(z)) is E/(1+E) or 1/(1+E), respectively.
// 1/(1+exp(z)) uses the complementary identity. No density-profile change.
static LogisticDisk stable_logistic_disk(S z,bool increasing){
 LogisticDisk out;out.argument=z;out.reflected=z.n>0;out.increasing=increasing;
 out.exponential=holomorphic_exp(out.reflected?-z:z);out.denominator=S(1)+out.exponential;
 need(sf::lower(out.denominator)>0,"same logistic exponential disk denominator excludes zero");
 bool numerator_exponential=bool(increasing)!=bool(out.reflected);
 out.value=numerator_exponential?out.exponential/out.denominator:S(1)/out.denominator;
 return out;
}
struct DensityDisk {S radius,radial_denominator;LogisticDisk inner,outer;S value;};
static DensityDisk stable_reference_density(const curved::Packet&p,S radius){
 DensityDisk out;out.radius=radius;out.radial_denominator=radius+sf::exact(p.parameters[7]);
 need(sf::lower(out.radial_denominator)>0,"original density r+epsilon complex pole not excluded");
 out.inner=stable_logistic_disk(S(2)*(radius-sf::exact(p.parameters[4]))/sf::exact(p.parameters[8]),true);
 out.outer=stable_logistic_disk(S(2)*(radius-sf::exact(p.parameters[5]))/sf::exact(p.parameters[8]),false);
 out.value=sf::exact(p.parameters[6])+sf::exact(p.parameters[9])/out.radial_denominator*out.inner.value*out.outer.value;
 return out;
}
static Q stable_complex_density_bound(const curved::Packet&p,S radius){return sf::mag(stable_reference_density(p,radius).value);}
}
#endif
