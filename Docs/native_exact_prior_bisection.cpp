/* Fixed 64-step Python-equivalent binary64 bisection on a rational polygon.
 * Shape dispatch is restricted to polygons wholly inside the circular wall.
 * Rational clip areas are compared after the same correctly rounded binary64
 * conversion used by float(Fraction). No geometric predicate is approximated.
 */
#include <gmpxx.h>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <vector>
#include <stdexcept>
#include <limits>

using Q=mpq_class;
struct Point { Q x,y; };
static bool same(const Point&a,const Point&b){return a.x==b.x&&a.y==b.y;}
static Q area(const std::vector<Point>&p,const Q&nx,const Q&ny,const Q&cut){
 std::vector<Point> out;
 auto push=[&](const Point&p){if(out.empty()||!same(out.back(),p))out.push_back(p);};
 for(size_t k=0;k<p.size();++k){
  const auto&a=p[k];const auto&b=p[(k+1)%p.size()];
  Q sa=nx*a.x+ny*a.y-cut,sb=nx*b.x+ny*b.y-cut;
  bool ai=sa<=0,bi=sb<=0;
  if(ai)push(a);
  if(ai!=bi){Q t=sa/(sa-sb);push({a.x+t*(b.x-a.x),a.y+t*(b.y-a.y)});}
 }
 if(out.size()>1&&same(out.front(),out.back()))out.pop_back();
 if(out.size()<3)return Q(0);
 Q twice=0;
 for(size_t k=0;k<out.size();++k){auto&a=out[k];auto&b=out[(k+1)%out.size()];twice+=a.x*b.y-a.y*b.x;}
 if(twice<=0)throw std::runtime_error("nonpositive exact clipped area");
 return twice/2;
}
extern "C" int exact_prior_bisect64(int n,const char*const*xy,double nx,double ny,
 double lo,double hi,double target,double radius,double*out,uint64_t*witness){
 try{
  if(n<3||!std::isfinite(nx)||!std::isfinite(ny)||!std::isfinite(lo)||!std::isfinite(hi)
   ||!(lo<hi)||!std::isfinite(target)||target<0||!std::isfinite(radius)||radius<=0)return 1;
  std::vector<Point> p;Q R(radius);R*=R;
  for(int k=0;k<n;++k){Q x(xy[2*k],10),y(xy[2*k+1],10);x.canonicalize();y.canonicalize();
   if(x*x+y*y>R)return 2;p.push_back({x,y});}
  // Python clips with (-normal,-mid), retaining normal.x >= mid.
  Q qx(-nx),qy(-ny);
  double next=std::nextafter(target,std::numeric_limits<double>::infinity());
  if(!std::isfinite(next))return 3;
  Q threshold=(Q(target)+Q(next))/2;
  uint64_t bits;std::memcpy(&bits,&target,sizeof bits);
  bool target_is_even=(bits&1)==0;
  uint64_t branches=0;
  for(int k=0;k<64;++k){
   double mid=(lo+hi)/2;
   Q a=area(p,qx,qy,Q(-mid));
   // rounded(area)>target iff area exceeds the adjacent rounding midpoint;
   // at an exact tie the even significand wins.
   bool greater=a>threshold||(a==threshold&&!target_is_even);
   if(greater){lo=mid;branches|=(uint64_t(1)<<k);}else hi=mid;
  }
  *out=(lo+hi)/2;*witness=branches;return 0;
 }catch(...){return 4;}
}
