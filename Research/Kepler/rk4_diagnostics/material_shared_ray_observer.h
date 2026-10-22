#ifndef MATERIAL_SHARED_RAY_OBSERVER_H
#define MATERIAL_SHARED_RAY_OBSERVER_H
namespace ray_observer {
struct Trace {uint64_t slot,target,chart;shared_function::S free,eta,t,width;shared_duffy::AffineAngular lo,hi,mixed;shared_function::S lo_value,hi_value;shared_function::Q lo_contraction,hi_contraction,lo_invariance,hi_invariance,lo_radius,hi_radius,lo_prior,hi_prior;};
static thread_local std::vector<Trace>* records=nullptr;static thread_local uint64_t slot=0,target=0;
template<class Lower,class Upper>static void observe(const shared_duffy::Compiled&m,const shared_duffy::Chart&chart,shared_function::S free,shared_function::S eta,const Lower&lo,const Upper&hi,const shared_duffy::AffineAngular&mixed,shared_function::S width,shared_function::S t){if(!records)return;Trace z={};z.slot=slot;z.target=target;z.chart=2*chart.apex+chart.part;z.free=free;z.eta=eta;z.t=t;z.width=width;z.lo=lo.coordinate;z.hi=hi.coordinate;z.mixed=mixed;z.lo_value=lo.certificate.value;z.hi_value=hi.certificate.value;z.lo_contraction=lo.certificate.contraction;z.hi_contraction=hi.certificate.contraction;z.lo_invariance=lo.certificate.invariance;z.hi_invariance=hi.certificate.invariance;z.lo_radius=lo.certificate.correction_radius;z.hi_radius=hi.certificate.correction_radius;z.lo_prior=lo.certificate.prior_radius;z.hi_prior=hi.certificate.prior_radius;records->push_back(z);}
}
#endif
