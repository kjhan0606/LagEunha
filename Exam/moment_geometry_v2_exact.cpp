/* Explicit native rational resolution. This host stage does not silently
 * replace an unresolved CUDA call; the caller requests it by this API. */
#include "moment_geometry_v2.h"
#include "moment_geometry_interval.h"
#include <gmpxx.h>
#include <vector>
#include <map>
#include <set>
#include <algorithm>
#include <fstream>
#include <cmath>
#include <cstring>
#include <cfenv>
#include <stdexcept>
#include <omp.h>
using Q=mpq_class;
struct PointQ{Q x,y;};
struct PlaneQ{Q x,y,b;};
struct WitnessQ{int64_t added;PointQ global;Q violation;};
struct PredicateQ{int64_t owner,a,b,plane;int exact_sign,filter_sign;};
struct CellQ{
 int64_t row=0;std::vector<PointQ> p;std::vector<int64_t> tag,candidates;
 std::vector<WitnessQ> witnesses;std::vector<PredicateQ> queries;uint64_t filtered=0,exact=0;
};
struct RegistryQ{
 int n;const int64_t *id;const double *xy,*w,*v,*z,*bounds;const int *kind;
 std::vector<Q> x,y,W;
 RegistryQ(int nr,const int64_t *ids,const double *pos,const double *weights,
 const double *vel,const double *rates,const int *k,const double *b):n(nr),id(ids),xy(pos),w(weights),v(vel),z(rates),bounds(b){
  kind=k;x.reserve(nr);y.reserve(nr);W.reserve(nr);
  for(int i=0;i<nr;i++){x.emplace_back(xy[2*i]);y.emplace_back(xy[2*i+1]);W.emplace_back(w[i]);}
 }
 PlaneQ plane(int64_t i,int64_t j)const{
  if(j==-1)return {Q(0),Q(-1),y[i]-Q(bounds[2])};
  if(j==-2)return {Q(1),Q(0),Q(bounds[1])-x[i]};
  if(j==-3)return {Q(0),Q(1),Q(bounds[3])-y[i]};
  if(j==-4)return {Q(-1),Q(0),x[i]-Q(bounds[0])};
  Q dx=x[j]-x[i],dy=y[j]-y[i];
  if(dx==0&&dy==0)throw 6;
  return {dx,dy,(dx*dx+dy*dy+W[i]-W[j])/2};
 }
 int64_t label(int64_t j)const{return j<0?j:id[j];}
};
static PointQ intersection(const PlaneQ &a,const PlaneQ &b){
 Q D=a.x*b.y-a.y*b.x;if(D==0)throw 6;
 return {(a.b*b.y-a.y*b.b)/D,(a.x*b.b-a.b*b.x)/D};
}
static bool same(const PointQ&a,const PointQ&b){return a.x==b.x&&a.y==b.y;}
static uint64_t bits(double x){uint64_t a;std::memcpy(&a,&x,8);return a;}
/* GMP get_d truncates. Choose nearest adjacent binary64 using exact midpoints. */
static double nearest(const Q &q){
 if(q==0)return 0.;bool neg=q<0;Q a=neg?-q:q;double lo=a.get_d();
 if(!std::isfinite(lo))throw 10;
 if(Q(lo)>a)throw 10;
 if(Q(lo)==a)return neg?-lo:lo;
 double hi=std::nextafter(lo,INFINITY);if(!std::isfinite(hi))throw 10;
 if(Q(hi)<a)throw 10;
 Q mid=(Q(lo)+Q(hi))/2;double r=a<mid?lo:(a>mid?hi:((bits(lo)&1)?hi:lo));
 return neg?-r:r;
}
static MG2Interval enclosure(const Q&q){
 double d=nearest(q);Q exact(d);
 if(exact==q)return {d,d};
 return exact<q?MG2Interval{d,std::nextafter(d,INFINITY)}:MG2Interval{std::nextafter(d,-INFINITY),d};
}
struct SqrtQ{double value,lo,hi;};
static SqrtQ sqrt_round(const Q&q){
 if(q<=0)throw 10;
 mpf_class f(q,256);mpf_sqrt(f.get_mpf_t(),f.get_mpf_t());double lo=f.get_d();
 if(!std::isfinite(lo)||lo<0)throw 10;
 int steps=0;while(Q(lo)*Q(lo)>q){lo=std::nextafter(lo,0.);if(++steps>8)throw 10;}
 double hi=std::nextafter(lo,INFINITY);if(!std::isfinite(hi))throw 10;
 while(Q(hi)*Q(hi)<q){lo=hi;hi=std::nextafter(hi,INFINITY);if(!std::isfinite(hi))throw 10;if(++steps>8)throw 10;}
 if(!std::isfinite(hi))throw 10;
 if(Q(lo)*Q(lo)==q)return {lo,lo,lo};
 if(Q(hi)*Q(hi)==q)return {hi,hi,hi};
 Q mid=(Q(lo)+Q(hi))/2,square=mid*mid;
 double value=q<square?lo:(q>square?hi:((bits(lo)&1)?hi:lo));
 if(!(value>0.))throw 10;
 return {value,lo,hi};
}
static int classify(CellQ &c,const RegistryQ&r,int e,int64_t j){
 int n=c.p.size();int64_t a=c.tag[(e+n-1)%n],b=c.tag[e];
 int s=mg2_support_sign(mg2_plane(c.row,a,r.xy,r.w,r.bounds),mg2_plane(c.row,b,r.xy,r.w,r.bounds),mg2_plane(c.row,j,r.xy,r.w,r.bounds));
 PlaneQ plane=r.plane(c.row,j);Q value=plane.x*c.p[e].x+plane.y*c.p[e].y-plane.b;
 int exact=sgn(value);c.queries.push_back({c.row,a,b,j,exact,s});
 if(s){if(s!=exact)throw 13;c.filtered++;return s;}
 c.exact++;return exact;
}
static void clip(CellQ&c,const RegistryQ&r,int64_t j){
 int n=c.p.size();std::vector<int> signs(n);for(int e=0;e<n;e++)signs[e]=classify(c,r,e,j);
 std::vector<PointQ> p;std::vector<int64_t>tag;
 for(int e=0;e<n;e++){
  int en=(e+1)%n;bool ai=signs[e]<=0,bi=signs[en]<=0;
  if(ai){p.push_back(c.p[e]);tag.push_back(c.tag[e]);}
  if(ai!=bi){p.push_back(intersection(r.plane(c.row,c.tag[e]),r.plane(c.row,j)));tag.push_back(ai?j:c.tag[e]);}
 }
 std::vector<PointQ> clean;std::vector<int64_t>ct;
 for(size_t e=0;e<p.size();e++){
  if(!clean.empty()&&same(p[e],clean.back()))ct.back()=tag[e];
  else{clean.push_back(p[e]);ct.push_back(tag[e]);}
 }
 if(clean.size()>1&&same(clean.back(),clean.front())){clean.pop_back();ct.pop_back();}
 if(clean.size()<3)throw 4;if(clean.size()>MG2_CAP)throw 2;
 c.p.swap(clean);c.tag.swap(ct);
}
static void construct(CellQ&c,const RegistryQ&r){
 Q xi=r.x[c.row],yi=r.y[c.row];
 c.p={{Q(r.bounds[0])-xi,Q(r.bounds[2])-yi},{Q(r.bounds[1])-xi,Q(r.bounds[2])-yi},
      {Q(r.bounds[1])-xi,Q(r.bounds[3])-yi},{Q(r.bounds[0])-xi,Q(r.bounds[3])-yi}};
 c.tag={-1,-2,-3,-4};for(auto j:c.candidates)clip(c,r,j);
}
static void complete(CellQ&c,const RegistryQ&r){
 for(int iteration=0;iteration<32;iteration++){
  construct(c,r);std::set<int64_t> add;
  for(const auto &p:c.p){
   MG2Interval x=enclosure(p.x),y=enclosure(p.y);
   for(int64_t j=0;j<r.n;j++){
    if(j==c.row)continue;
    auto plane=mg2_plane(c.row,j,r.xy,r.w,r.bounds);
    auto residual=mg2_sub(mg2_add(mg2_mul(plane.x,x),mg2_mul(plane.y,y)),plane.b);
    int sign=mg2_sign(residual);Q exact;
    if(sign)c.filtered++;
    else{c.exact++;PlaneQ q=r.plane(c.row,j);exact=q.x*p.x+q.y*p.y-q.b;sign=sgn(exact);}
    if(sign>0){
     if(std::find(c.candidates.begin(),c.candidates.end(),j)!=c.candidates.end())throw 6;
     if(add.insert(j).second){PlaneQ q=r.plane(c.row,j);exact=q.x*p.x+q.y*p.y-q.b;c.witnesses.push_back({r.id[j],{p.x+r.x[c.row],p.y+r.y[c.row]},exact});}
    }
   }
  }
  if(add.empty())return;
  c.candidates.insert(c.candidates.end(),add.begin(),add.end());
  std::sort(c.candidates.begin(),c.candidates.end(),[&](int64_t a,int64_t b){return r.id[a]<r.id[b];});
 }
 throw 6;
}
static void export_cell(const CellQ&c,const RegistryQ&r,uint64_t epoch,MomentGeometryV2&o){
 o={};o.epoch=epoch;o.count=c.p.size();o.filtered_checks=c.filtered;o.exact_checks=c.exact;o.added_candidates=c.witnesses.size();
 Q twice=0,qx=0,qy=0,xx=0,yy=0,xy=0;int n=c.p.size();
 for(int e=0;e<n;e++){
  const PointQ&a=c.p[e],&b=c.p[(e+1)%n];Q cross=a.x*b.y-b.x*a.y;
  twice+=cross;qx+=(a.x+b.x)*cross;qy+=(a.y+b.y)*cross;
  xx+=(a.x*a.x+a.x*b.x+b.x*b.x)*cross;yy+=(a.y*a.y+a.y*b.y+b.y*b.y)*cross;
  xy+=(2*a.x*a.y+a.x*b.y+b.x*a.y+2*b.x*b.y)*cross;
 }
 if(twice<=0)throw 6;
 o.area=nearest(twice/2);o.centroid[0]=nearest(qx/(3*twice));o.centroid[1]=nearest(qy/(3*twice));
 o.second[0]=nearest(xx/12);o.second[1]=nearest(xy/24);o.second[2]=o.second[1];o.second[3]=nearest(yy/12);
 if(!(o.area>0.))throw 10;
 for(int e=0;e<n;e++){
  int en=(e+1)%n;int64_t j=c.tag[e];if(j<0)throw 5;
  const PointQ&a=c.p[e],&b=c.p[en];PointQ mid{(a.x+b.x)/2,(a.y+b.y)/2};Q ex=b.x-a.x,ey=b.y-a.y;
  Q dx=r.x[j]-r.x[c.row],dy=r.y[j]-r.y[c.row];
  if(ey*dx-ex*dy<=0||ex*dx+ey*dy!=0)throw 6;
  auto len=sqrt_round(ex*ex+ey*ey),dist=sqrt_round(dx*dx+dy*dy);
  o.x[e]=nearest(a.x);o.y[e]=nearest(a.y);o.midpoint_x[e]=nearest(mid.x);o.midpoint_y[e]=nearest(mid.y);
  o.edge_dx[e]=nearest(ex);o.edge_dy[e]=nearest(ey);o.length[e]=len.value;o.length_lo[e]=len.lo;o.length_hi[e]=len.hi;
  o.support_a[e]=r.label(c.tag[(e+n-1)%n]);o.support_b[e]=r.label(j);o.neighbor[e]=r.id[j];o.neighbor_kind[e]=r.kind[j];
  o.nx[e]=nearest(dx)/dist.value;o.ny[e]=nearest(dy)/dist.value;
  Q vxj(r.v[2*j]),vyj(r.v[2*j+1]),vxi(r.v[2*c.row]),vyi(r.v[2*c.row+1]);
  Q numerator=dx*vxj+dy*vyj-(vxj-vxi)*mid.x-(vyj-vyi)*mid.y+(Q(r.z[c.row])-Q(r.z[j]))/2;
  o.geometry_wn[e]=nearest(numerator)/dist.value;
  if(!std::isfinite(o.geometry_wn[e])||(!o.edge_dx[e]&&!o.edge_dy[e]))throw 10;
  double sx=nearest(b.x)-o.x[e],sy=nearest(b.y)-o.y[e];
  if(sx!=o.edge_dx[e]||sy!=o.edge_dy[e])o.endpoint_rounding_difference++;
  if(sx==0.&&sy==0.)o.endpoint_collapse++;
  if(Q(sy)*dx-Q(sx)*dy<0)o.endpoint_reversal++;
  double relative=std::hypot(sx-o.edge_dx[e],sy-o.edge_dy[e])/len.value;
  o.max_relative_vector_error=std::max(o.max_relative_vector_error,relative);
 }
}
static void shared(const std::vector<CellQ>&cells,const RegistryQ&r){
 std::map<int64_t,int> selected;for(size_t k=0;k<cells.size();k++)if(!selected.emplace(r.id[cells[k].row],k).second)throw 1;
 for(const auto &c:cells)for(size_t e=0;e<c.p.size();e++){
  int64_t j=c.tag[e];if(j<0)throw 5;auto other=selected.find(r.id[j]);if(other==selected.end())continue;
  const auto &d=cells[other->second];int found=-1;
  for(size_t k=0;k<d.p.size();k++)if(d.tag[k]==c.row){if(found>=0)throw 9;found=k;}
  if(found<0)throw 9;
  PointQ a{c.p[e].x+r.x[c.row],c.p[e].y+r.y[c.row]},b{c.p[(e+1)%c.p.size()].x+r.x[c.row],c.p[(e+1)%c.p.size()].y+r.y[c.row]};
  PointQ da{d.p[found].x+r.x[d.row],d.p[found].y+r.y[d.row]},db{d.p[(found+1)%d.p.size()].x+r.x[d.row],d.p[(found+1)%d.p.size()].y+r.y[d.row]};
  if(!same(a,db)||!same(b,da))throw 9;
 }
}
static void write_point(std::ostream&s,const PointQ&p){s<<"[\""<<p.x.get_str()<<"\",\""<<p.y.get_str()<<"\"]";}
static bool audit(const char*path,const std::vector<CellQ>&cells,const RegistryQ&r){
 if(!path||!*path)return true;std::ifstream exists(path);if(exists.good())return false;
 std::ofstream out(path);if(!out)return false;
 out<<"{\"stage\":\"explicit_native_host_exact_resolution\",\"all_registry_count\":"<<r.n<<",\"shared_global_endpoints_exact_reversed_equal\":true,\"cells\":[";
 for(size_t i=0;i<cells.size();i++){
  const auto &c=cells[i];if(i)out<<",";
  out<<"{\"owner\":"<<r.id[c.row]<<",\"filtered_checks\":"<<c.filtered<<",\"exact_checks\":"<<c.exact<<",\"witnesses\":[";
  for(size_t k=0;k<c.witnesses.size();k++){
   if(k)out<<",";const auto&w=c.witnesses[k];out<<"{\"added_id\":"<<w.added<<",\"vertex_global\":";write_point(out,w.global);out<<",\"positive_violation\":\""<<w.violation.get_str()<<"\"}";
  }
  out<<"],\"final_candidate_ids\":[";
  for(size_t k=0;k<c.candidates.size();k++){if(k)out<<",";out<<r.id[c.candidates[k]];}
  out<<"],\"clipping_predicate_queries\":[";
  for(size_t k=0;k<c.queries.size();k++){if(k)out<<",";const auto&q=c.queries[k];out<<"["<<q.owner<<","<<q.a<<","<<q.b<<","<<q.plane<<","<<q.exact_sign<<","<<q.filter_sign<<"]";}
  out<<"],\"edges\":[";
  for(size_t e=0;e<c.p.size();e++){
   if(e)out<<",";const auto&a=c.p[e];const auto&b=c.p[(e+1)%c.p.size()];Q dx=b.x-a.x,dy=b.y-a.y,length2=dx*dx+dy*dy;
   out<<"{\"neighbor\":"<<r.label(c.tag[e])<<",\"support\":["<<r.label(c.tag[(e+c.p.size()-1)%c.p.size()])<<","<<r.label(c.tag[e])<<"],\"A_global\":";write_point(out,{a.x+r.x[c.row],a.y+r.y[c.row]});out<<",\"B_global\":";write_point(out,{b.x+r.x[c.row],b.y+r.y[c.row]});out<<",\"length_squared\":\""<<length2.get_str()<<"\"}";
  }
  out<<"]}";
 }
 out<<"]}\n";out.flush();return out.good();
}
extern "C" size_t moment_geometry_v2_size(void){return sizeof(MomentGeometryV2);}
extern "C" int moment_geometry_v2_resolve(int threads,int n,int nr,const int64_t*ids,
 const double*xy,const double*W,const double*v,const double*z,const int*kind,
 const int64_t*selected,const int64_t*off,const int64_t*cand,const double*bounds,
 uint64_t epoch,MomentGeometryV2*output,int*status,const char*audit_path){
 if(threads<1||n<1||nr<1||!ids||!xy||!W||!v||!z||!kind||!selected||!off||!cand||!bounds||!output||!status||std::fegetround()!=FE_TONEAREST)return 1;
 if(off[0]!=0||off[n]<0)return 1;
 for(int k=0;k<4;k++)if(!std::isfinite(bounds[k]))return 1;
 if(bounds[0]>=bounds[1]||bounds[2]>=bounds[3])return 1;
 try{
  std::set<int64_t> unique;
  for(int i=0;i<nr;i++){
   if(ids[i]<0||!unique.insert(ids[i]).second||!std::isfinite(xy[2*i])||!std::isfinite(xy[2*i+1])||!std::isfinite(W[i])||!std::isfinite(v[2*i])||!std::isfinite(v[2*i+1])||!std::isfinite(z[i])||kind[i]<0||kind[i]>2)return 1;
  }
  for(int i=0;i<n;i++){
   if(selected[i]<0||selected[i]>=nr||off[i]>off[i+1])return 1;
   int64_t prev=-1;for(int64_t k=off[i];k<off[i+1];k++){int64_t j=cand[k];if(j<0||j>=nr||j==selected[i]||ids[j]<=prev)return 1;prev=ids[j];}
  }
  RegistryQ r(nr,ids,xy,W,v,z,kind,bounds);std::vector<CellQ> cells(n);std::vector<MomentGeometryV2> tmp(n);std::vector<int> st(n,0);
  #pragma omp parallel for num_threads(threads)
  for(int i=0;i<n;i++)try{cells[i].row=selected[i];cells[i].candidates.assign(cand+off[i],cand+off[i+1]);complete(cells[i],r);export_cell(cells[i],r,epoch,tmp[i]);}catch(int e){st[i]=e;}catch(...){st[i]=8;}
  std::copy(st.begin(),st.end(),status);for(int s:st)if(s)return s;
  shared(cells,r);if(!audit(audit_path,cells,r))return 11;
  std::copy(tmp.begin(),tmp.end(),output);return 0;
 }catch(int e){return e;}catch(...){return 8;}
}

/* Independent exact sign oracle for recorded support-plane queries.
 * Each row is [owner registry row, support-a row/box, support-b, tested row/box].
 * exact sign2 means the two supports have D=0 (no unique intersection). */
extern "C" int moment_geometry_v2_exact_queries(int nr,const double*xy,const double*W,
 const double*bounds,int n,const int64_t*queries,int*output){
 if(nr<1||n<0||!xy||!W||!bounds||!queries||!output)return 1;
 for(int i=0;i<nr;i++)if(!std::isfinite(xy[2*i])||!std::isfinite(xy[2*i+1])||!std::isfinite(W[i]))return 1;
 for(int k=0;k<4;k++)if(!std::isfinite(bounds[k]))return 1;
 if(bounds[0]>=bounds[1]||bounds[2]>=bounds[3])return 1;
 try{
  std::vector<int64_t>ids(nr);std::vector<double>v(2*nr,0.),z(nr,0.);std::vector<int>kind(nr,0),tmp(n);
  for(int i=0;i<nr;i++)ids[i]=i;
  RegistryQ r(nr,ids.data(),xy,W,v.data(),z.data(),kind.data(),bounds);
  for(int k=0;k<n;k++){
   const int64_t*q=queries+4*k;int64_t i=q[0];if(i<0||i>=nr)return 1;
   for(int t=1;t<4;t++)if(q[t]<-4||q[t]>=nr||q[t]==i)return 1;
   PlaneQ a=r.plane(i,q[1]),b=r.plane(i,q[2]),c=r.plane(i,q[3]);Q D=a.x*b.y-a.y*b.x;
   if(D==0){tmp[k]=2;continue;}
   Q X=a.b*b.y-a.y*b.b,Y=a.x*b.b-a.b*b.x,N=c.x*X+c.y*Y-c.b*D;
   tmp[k]=sgn(N)*sgn(D);
  }
  std::copy(tmp.begin(),tmp.end(),output);return 0;
 }catch(int e){return e;}catch(...){return 8;}
}
