#ifndef MATERIAL_SHARED_REAL_GRAPH_VALUE_H
#define MATERIAL_SHARED_REAL_GRAPH_VALUE_H
// Same already selected implicit branch, evaluated at one actual real Gauss
// free-coordinate ball. The saved complex-domain certificate supplies the
// initial root enclosure. No new canonical/event root or arrangement search.
namespace real_graph_value {
static thread_local uint64_t observed_face_evaluations=0;
struct Result {S value,derivative;uint64_t iterations=0,evaluations=0;bool family_floor=false;};
template<class Evaluate>static Result refine(S initial,Evaluate evaluate){
 Result out;out.value=initial;
 for(uint64_t it=0;it<256;it++){
  out.iterations=it+1;
  observed_face_evaluations++;auto box=evaluate(out.value);out.evaluations++;
  out.derivative=box.second;
  need(sf::lower(out.derivative)>0||sf::upper(out.derivative)<0,"saved real graph-value derivative separated");
  observed_face_evaluations++;auto center=evaluate(S(out.value.n));out.evaluations++;
  S inclusion=S(out.value.n)-center.first/out.derivative;
  Q lo=std::max(sf::lower(out.value),sf::lower(inclusion)),hi=std::min(sf::upper(out.value),sf::upper(inclusion));
  need(lo<=hi,"same selected real graph-value interval Newton inclusion");
  S next=sf::range(lo,hi);
  if(hi-lo<1e-27Q){out.value=next;return out;}
  if(!(next.e<out.value.e)){out.family_floor=true;return out;}
  out.value=next;
 }
 throw std::runtime_error("original fixed256 real graph-value work budget");
}
static Result angular(const target::Face&face,const sd::Chart&chart,S free,sd::AffineAngular&coordinate){
 need(!coordinate.interpolated&&!coordinate.oblique_collected,"saved regular boundary coordinate, not a different interpolation branch");
 const auto base=coordinate;
 auto out=refine(base.correction,[&](S correction){auto c=base;c.correction=correction;auto jet=target::face_q(face,chart,free,c);return std::make_pair(jet.v,jet.g[1]);});
 coordinate.correction=out.value;return out;
}
static Result radial(const target::Face&face,const sd::Chart&chart,S free,S initial){
 return refine(initial,[&](S radial){auto jet=target::in_axis(face,chart,0,free,radial);return std::make_pair(jet.v,jet.g[0]);});
}
}
#endif
