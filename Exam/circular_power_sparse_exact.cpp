/* Sparse weighted-power completion of EXACT polygon supports before disk clip.
 * Reuses the audited MG2 rational intersections and outward interval filters.
 * A spatial bounding tree certifies all omitted sites: box distance^2-max(W)
 * is an outward lower bound on any site's power distance at each exact vertex.
 * No Qhull completeness assumption, dropped short edge, or approximate sign.
 */
#include "moment_geometry_v2_exact.cpp"
#include <numeric>
#include <iomanip>

struct PowerTreeNode {
 double xmin,xmax,ymin,ymax,wmax;
 int begin,end,left=-1,right=-1;
};
struct PowerSearchStats {uint64_t visits=0,pruned=0,leaf_checks=0,exact_comparisons=0;};
class PowerTree {
 const RegistryQ &r;
 std::vector<int> order;
 std::vector<PowerTreeNode> nodes;
 int build(int begin,int end){
  PowerTreeNode node{INFINITY,-INFINITY,INFINITY,-INFINITY,-INFINITY,begin,end};
  for(int k=begin;k<end;k++){
   int j=order[k];node.xmin=std::min(node.xmin,r.xy[2*j]);node.xmax=std::max(node.xmax,r.xy[2*j]);
   node.ymin=std::min(node.ymin,r.xy[2*j+1]);node.ymax=std::max(node.ymax,r.xy[2*j+1]);node.wmax=std::max(node.wmax,r.w[j]);
  }
  int index=nodes.size();nodes.push_back(node);
  if(end-begin>12){
   int axis=(node.xmax-node.xmin>=node.ymax-node.ymin)?0:1,mid=(begin+end)/2;
   std::nth_element(order.begin()+begin,order.begin()+mid,order.begin()+end,[&](int a,int b){return r.xy[2*a+axis]<r.xy[2*b+axis]||(r.xy[2*a+axis]==r.xy[2*b+axis]&&a<b);});
   int l=build(begin,mid),rr=build(mid,end);nodes[index].left=l;nodes[index].right=rr;
  }
  return index;
 }
 static MG2Interval distance_axis(MG2Interval x,double lo,double hi){
  if(x.hi<lo)return mg2_sub(mg2_point(lo),x);
  if(x.lo>hi)return mg2_sub(x,mg2_point(hi));
  return mg2_point(0.);
 }
 MG2Interval lower(const PowerTreeNode &node,MG2Interval px,MG2Interval py)const{
  auto dx=distance_axis(px,node.xmin,node.xmax),dy=distance_axis(py,node.ymin,node.ymax);
  auto xx=mg2_mul(dx,dx),yy=mg2_mul(dy,dy);xx.lo=std::max(0.,xx.lo);yy.lo=std::max(0.,yy.lo);
  return mg2_sub(mg2_add(xx,yy),mg2_point(node.wmax));
 }
 void visit(int index,const PointQ &p,MG2Interval px,MG2Interval py,int &best,Q &power,MG2Interval &bound,PowerSearchStats &stats)const{
  const auto &node=nodes[index];stats.visits++;
  if(lower(node,px,py).lo>=bound.hi){stats.pruned++;return;}
  if(node.left<0){
   for(int k=node.begin;k<node.end;k++){
    int j=order[k];if(j==best)continue;stats.leaf_checks++;
    auto dx=mg2_sub(px,mg2_point(r.xy[2*j])),dy=mg2_sub(py,mg2_point(r.xy[2*j+1]));
    auto value=mg2_sub(mg2_add(mg2_mul(dx,dx),mg2_mul(dy,dy)),mg2_point(r.w[j]));
    if(value.lo>bound.hi)continue;
    Q qdx=p.x-r.x[j],qdy=p.y-r.y[j];Q actual=qdx*qdx+qdy*qdy-r.W[j];stats.exact_comparisons++;
    if(actual<power){best=j;power=actual;bound=enclosure(power);}
   }
  }else{
   int a=node.left,b=node.right;
   if(lower(nodes[b],px,py).lo<lower(nodes[a],px,py).lo)std::swap(a,b);
   visit(a,p,px,py,best,power,bound,stats);visit(b,p,px,py,best,power,bound,stats);
  }
 }
 public:
 PowerTree(const RegistryQ &registry):r(registry),order(registry.n){std::iota(order.begin(),order.end(),0);nodes.reserve(2*r.n);build(0,r.n);}
 int nearest_power(const PointQ &p,int owner,PowerSearchStats &stats)const{
  Q dx=p.x-r.x[owner],dy=p.y-r.y[owner];Q power=dx*dx+dy*dy-r.W[owner];int best=owner;
  auto bound=enclosure(power);visit(0,p,enclosure(p.x),enclosure(p.y),best,power,bound,stats);return best;
 }
};

static void complete_sparse(CellQ &c,const RegistryQ &r,const PowerTree &tree,PowerSearchStats &stats){
 for(int iteration=0;iteration<64;iteration++){
  construct(c,r);std::set<int64_t> add;
  for(const auto &p:c.p){
   PointQ global{p.x+r.x[c.row],p.y+r.y[c.row]};int j=tree.nearest_power(global,c.row,stats);
   if(j==c.row)continue;
   PlaneQ plane=r.plane(c.row,j);Q violation=plane.x*p.x+plane.y*p.y-plane.b;
   if(violation<=0)throw 13;
   if(std::find(c.candidates.begin(),c.candidates.end(),j)!=c.candidates.end())throw 6;
   if(add.insert(j).second)c.witnesses.push_back({r.id[j],global,violation});
  }
  if(add.empty())return;
  c.candidates.insert(c.candidates.end(),add.begin(),add.end());
  std::sort(c.candidates.begin(),c.candidates.end(),[&](int64_t a,int64_t b){return r.id[a]<r.id[b];});
 }
 throw 6;
}

static void shared_with_box(const std::vector<CellQ>&cells,const RegistryQ&r){
 for(const auto &c:cells)for(size_t e=0;e<c.p.size();e++){
  int64_t j=c.tag[e];if(j<0)continue;
  const auto &d=cells.at(j);int found=-1;
  for(size_t k=0;k<d.p.size();k++)if(d.tag[k]==c.row){if(found>=0)throw 9;found=k;}
  if(found<0)throw 9;
  PointQ a{c.p[e].x+r.x[c.row],c.p[e].y+r.y[c.row]},b{c.p[(e+1)%c.p.size()].x+r.x[c.row],c.p[(e+1)%c.p.size()].y+r.y[c.row]};
  PointQ da{d.p[found].x+r.x[d.row],d.p[found].y+r.y[d.row]},db{d.p[(found+1)%d.p.size()].x+r.x[d.row],d.p[(found+1)%d.p.size()].y+r.y[d.row]};
  if(!same(a,db)||!same(b,da))throw 9;
 }
}
extern "C" int circular_power_sparse_resolve(int threads,int n,const int64_t*ids,
 const double*xy,const double*W,const int64_t*offsets,const int64_t*candidates,
 double radius,int*status,const char*path){
 if(threads<1||n<1||!ids||!xy||!W||!offsets||!candidates||!status||!path||!std::isfinite(radius)||radius<=0||std::fegetround()!=FE_TONEAREST)return 1;
 try{
  std::set<int64_t> unique;
  for(int i=0;i<n;i++)if(ids[i]<0||!unique.insert(ids[i]).second||!std::isfinite(xy[2*i])||!std::isfinite(xy[2*i+1])||!std::isfinite(W[i])||offsets[i]>offsets[i+1])return 1;
  if(offsets[0]!=0)return 1;
  for(int i=0;i<n;i++){
   int64_t previous=-1;
   for(int64_t k=offsets[i];k<offsets[i+1];k++){auto j=candidates[k];if(j<0||j>=n||j==i||ids[j]<=previous)return 1;previous=ids[j];}
  }
  std::vector<double>v(2*n,0.),z(n,0.);std::vector<int>kind(n,0);double bounds[4]={-radius,radius,-radius,radius};
  RegistryQ r(n,ids,xy,W,v.data(),z.data(),kind.data(),bounds);PowerTree tree(r);std::vector<CellQ>cells(n);std::vector<PowerSearchStats>stats(n);std::vector<int>st(n,0);
  #pragma omp parallel for num_threads(threads)
  for(int i=0;i<n;i++)try{cells[i].row=i;cells[i].candidates.assign(candidates+offsets[i],candidates+offsets[i+1]);complete_sparse(cells[i],r,tree,stats[i]);}catch(int e){st[i]=e;}catch(...){st[i]=8;}
  std::copy(st.begin(),st.end(),status);for(int s:st)if(s)return s;
  shared_with_box(cells,r);
  std::ifstream exists(path);if(exists.good())return 11;std::ofstream out(path);if(!out)return 11;
  out<<"{\"stage\":\"sparse_native_exact_polygon_before_disk_clip\",\"registry_count\":"<<n<<",\"shared_exact_reversed\":true,\"cells\":[";
  for(int i=0;i<n;i++){
   const auto &c=cells[i];const auto &s=stats[i];if(i)out<<",";
   out<<"{\"owner\":"<<ids[i]<<",\"tree_visits\":"<<s.visits<<",\"pruned_nodes\":"<<s.pruned<<",\"leaf_checks\":"<<s.leaf_checks<<",\"exact_power_comparisons\":"<<s.exact_comparisons<<",\"filtered_checks\":"<<c.filtered<<",\"exact_checks\":"<<c.exact<<",\"added_candidates\":"<<c.witnesses.size()<<",\"candidate_ids\":[";
   for(size_t k=0;k<c.candidates.size();k++){if(k)out<<",";out<<ids[c.candidates[k]];}out<<"],\"edges\":[";
   for(size_t e=0;e<c.p.size();e++){
    if(e)out<<",";const auto &a=c.p[e];const auto &b=c.p[(e+1)%c.p.size()];
    out<<"{\"neighbor\":"<<r.label(c.tag[e])<<",\"A\":";write_point(out,{a.x+r.x[i],a.y+r.y[i]});out<<",\"B\":";write_point(out,{b.x+r.x[i],b.y+r.y[i]});out<<"}";
   }out<<"]}";
  }
  out<<"]}\n";out.flush();return out.good()?0:11;
 }catch(int e){return e;}catch(...){return 8;}
}
