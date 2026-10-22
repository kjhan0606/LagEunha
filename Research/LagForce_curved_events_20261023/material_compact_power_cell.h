#ifndef MATERIAL_COMPACT_POWER_CELL_H
#define MATERIAL_COMPACT_POWER_CELL_H
// Certified half-plane construction in a source-enclosing rectangle. These
// are BOUNDED PLANE power cells, NOT T_i=Omega* intersect PowerCell_i. Actual
// curved boundary events/components/coverage must still be assembled.
#include "material_compact_power_curve.h"
#include <cstdio>
#include <queue>
#include <tuple>
namespace compact_power {
struct LineKey {
    uint64_t kind=0,a=0,b=0; // 0 enclosing box (axis,side);1 sorted cell pair
    bool operator<(const LineKey&x)const{return std::tie(kind,a,b)<std::tie(x.kind,x.a,x.b);}
    bool operator==(const LineKey&x)const{return kind==x.kind&&a==x.a&&b==x.b;}
};
struct Constraint {LineKey key;Plane plane;int side=1;}; // side*f<=0 inside
struct CornerKey {
    uint64_t kind=0;std::array<uint64_t,3>id={};
    bool operator<(const CornerKey&x)const{return std::tie(kind,id)<std::tie(x.kind,x.id);}
};
struct Corner {
    Constraint line[2];CornerKey key;std::array<S,2>x;
    // x is an enclosure of the canonical TWO-LINE INTERSECTION, not a
    // rounded vertex chosen as new exact geometry. Keep its defining keys.
};
inline Corner crossing(const compact_target::Registry&registry,const std::vector<R>&weights,Constraint a,Constraint b){
    require(!(a.key==b.key),"distinct canonical lines at a power corner");if(b.key<a.key)std::swap(a,b);
    Corner c;c.line[0]=a;c.line[1]=b;
    if(a.key.kind==1&&b.key.kind==1){
        std::set<uint64_t>ids={a.key.a,a.key.b,b.key.a,b.key.b};
        require(ids.size()==3,"one owning site and two neighbors at a simple power vertex");c.key.kind=2;std::copy(ids.begin(),ids.end(),c.key.id.begin());
        // All three incident cells use the SAME smallest-ID site and two
        // SAME canonical planes, even when their local edge pairs differ.
        auto at=[&](uint64_t id){auto it=std::lower_bound(registry.site.begin(),registry.site.end(),id,[](const compact_source::Site&s,uint64_t x){return s.cell<x;});require(it!=registry.site.end()&&it->cell==id,"canonical corner generator identity");return size_t(it-registry.site.begin());};
        size_t i=at(c.key.id[0]),j=at(c.key.id[1]),k=at(c.key.id[2]);
        a.plane=plane(registry.site[i],weights[i],registry.site[j],weights[j]);b.plane=plane(registry.site[i],weights[i],registry.site[k],weights[k]);
    }else if(a.key.kind==0&&b.key.kind==1){c.key.kind=1;c.key.id={2*a.key.a+a.key.b,b.key.a,b.key.b};}
    else{require(a.key.kind==0&&b.key.kind==0,"known canonical bounding or power lines");c.key.kind=0;c.key.id={2*a.key.a+a.key.b,2*b.key.a+b.key.b,0};}
    const auto&A=a.plane;const auto&B=b.plane;S det=A.normal[0]*B.normal[1]-A.normal[1]*B.normal[0];
    require(strict_sign(det)!=0,"certified nonparallel corner lines, no tolerance merge");
    S rhs=B.offset;for(int k=0;k<2;k++)rhs=rhs-B.normal[k]*(A.origin[k]-B.origin[k]);
    S ux=(A.offset*B.normal[1]-A.normal[1]*rhs)/det;
    S uy=(A.normal[0]*rhs-A.offset*B.normal[0])/det;
    c.x={A.origin[0]+ux,A.origin[1]+uy};return c;
}
inline int side(const Constraint&line,const Corner&point){
    if(point.line[0].key==line.key||point.line[1].key==line.key)return 0; // algebraic identity
    std::array<D,2>x={D(point.x[0]),D(point.x[1])};int sign=strict_sign(value(line.plane,x).v);
    if(!sign){
        auto v=value(line.plane,x).v;
        char message[320];std::snprintf(message,sizeof(message),
            "unresolved third-line corner incidence corner_kind=%llu ids=%llu,%llu,%llu line_kind=%llu pair=%llu,%llu value=[%.21Le,%.21Le]",
            static_cast<unsigned long long>(point.key.kind),
            static_cast<unsigned long long>(point.key.id[0]),
            static_cast<unsigned long long>(point.key.id[1]),
            static_cast<unsigned long long>(point.key.id[2]),
            static_cast<unsigned long long>(line.key.kind),
            static_cast<unsigned long long>(line.key.a),
            static_cast<unsigned long long>(line.key.b),
            static_cast<long double>(sf::lower(v)),
            static_cast<long double>(sf::upper(v)));
        throw std::runtime_error(message);
    }
    return sign*line.side;
}
inline Constraint common_line(const Corner&a,const Corner&b){
    int matches=0;Constraint out;
    for(const auto&A:a.line)for(const auto&B:b.line)if(A.key==B.key){out=A;matches++;}
    require(matches==1,"one canonical support line per nondegenerate polygon edge");return out;
}
inline void clip(const compact_target::Registry&registry,const std::vector<R>&weights,std::vector<Corner>&polygon,const Constraint&line,size_t maximum_corners){
    if(polygon.empty())return;require(polygon.size()>=3,"nondegenerate current bounded power cell");
    std::vector<int>sign(polygon.size());for(size_t j=0;j<polygon.size();j++)sign[j]=side(line,polygon[j]);
    std::vector<Corner>out;out.reserve(polygon.size()+1);
    for(size_t j=0;j<polygon.size();j++){
        size_t next=(j+1)%polygon.size();const auto&a=polygon[j];const auto&b=polygon[next];int sa=sign[j],sb=sign[next];
        if(sb<=0){if(sa>0&&sb<0)out.push_back(crossing(registry,weights,common_line(a,b),line));out.push_back(b);}
        else if(sa<0)out.push_back(crossing(registry,weights,common_line(a,b),line));
    }
    require(out.empty()||out.size()>=3,"lower-dimensional bounded power cell is not an accepted material cell");
    require(out.size()<=maximum_corners,"finite bounded power cell corner budget");polygon=std::move(out);
}
struct Node {
    Q lower[2]={},upper[2]={},maximum_weight=0;
    size_t begin=0,end=0,left=SIZE_MAX,right=SIZE_MAX;
};
struct Spatial {
    const compact_target::Registry&registry;
    const std::vector<R>weight; // OWNED immutable copy for this construction
    uint64_t weight_hash=1469598103934665603ULL;
    std::vector<size_t>order;std::vector<Node>node;
    double lower[2]={},upper[2]={};
    Spatial(const compact_target::Registry&r,const std::vector<R>&w):registry(r),weight(w){
        require(!r.site.empty()&&r.site.size()==w.size()&&!r.exterior.empty(),"actual small geometry and one fixed weight per site");
        for(double&x:lower)x=INFINITY;for(double&x:upper)x=-INFINITY;
        for(const auto&e:r.exterior)for(int k=0;k<2;k++){require(std::isfinite(e.lower[k])&&std::isfinite(e.upper[k])&&e.lower[k]<=e.upper[k],"same source exterior enclosure");lower[k]=std::min(lower[k],e.lower[k]);upper[k]=std::max(upper[k],e.upper[k]);}
        for(int k=0;k<2;k++)require(lower[k]<upper[k],"positive source-enclosing rectangle, never original-square clipping");
        for(size_t j=0;j<w.size();j++){
            require(std::isfinite(w[j])&&(j==0||r.site[j-1].cell<r.site[j].cell),"fixed finite weights and unique sorted source generators");
            double hi=double(w[j]),lo=double(w[j]-R(hi));require(std::isfinite(hi)&&std::isfinite(lo)&&R(hi)+R(lo)==w[j],"finite lossless two-double weight seal");
            for(double v:{hi,lo}){uint64_t bits=compact_target::bits(v);for(int k=0;k<8;k++){weight_hash^=(bits>>(8*k))&255;weight_hash*=1099511628211ULL;}}
            for(int k=0;k<2;k++)require(std::isfinite(r.site[j].position[k])&&std::isfinite(r.site[j].position_error[k])&&r.site[j].position_error[k]>=0,"finite actual canonical site bounds");
        }
        order.resize(r.site.size());std::iota(order.begin(),order.end(),size_t(0));node.reserve(2*order.size());build(0,order.size());
    }
    Spatial(const Spatial&)=delete;Spatial&operator=(const Spatial&)=delete;
    size_t build(size_t begin,size_t end){
        size_t id=node.size();node.emplace_back();auto&n=node[id];n.begin=begin;n.end=end;n.maximum_weight=Q(weight[order[begin]]);
        for(int k=0;k<2;k++){n.lower[k]=HUGE_VALQ;n.upper[k]=-HUGE_VALQ;}
        for(size_t p=begin;p<end;p++){
            size_t j=order[p];n.maximum_weight=std::max(n.maximum_weight,Q(weight[j]));
            for(int k=0;k<2;k++){S q=compact_source::canonical_position(registry.site[j],k);n.lower[k]=std::min(n.lower[k],sf::lower(q));n.upper[k]=std::max(n.upper[k],sf::upper(q));}
        }
        if(end-begin>8){
            int axis=n.upper[1]-n.lower[1]>n.upper[0]-n.lower[0]?1:0;size_t mid=begin+(end-begin)/2;
            std::nth_element(order.begin()+begin,order.begin()+mid,order.begin()+end,[&](size_t a,size_t b){double x=registry.site[a].position[axis],y=registry.site[b].position[axis];return x<y||(x==y&&registry.site[a].cell<registry.site[b].cell);});
            size_t left=build(begin,mid),right=build(mid,end);node[id].left=left;node[id].right=right;
        }return id;
    }
};
struct CellBudget {size_t maximum_nodes=0,maximum_planes=0,maximum_corners=0;};
struct Exclusion {size_t node=0;Q lower_margin=0;};
struct BoundedCell {
    uint64_t cell=0,source_registry_hash=0,weight_hash=0;std::vector<Corner>corner;std::vector<size_t>tested_sites;
    std::vector<Exclusion>excluded;size_t nodes=0;
    // Even a nonempty cell here may have empty/disconnected intersection with
    // Omega*. No star-shapedness, collar or actual mass is inferred here.
};
inline std::vector<Corner> rectangle(const Spatial&tree){
    Constraint line[4];for(int k=0;k<2;k++)for(int s=0;s<2;s++){
        auto&c=line[2*k+s];c.key={0,uint64_t(k),uint64_t(s)};c.plane.normal[k]=S(1);c.plane.offset=S(Q(s?tree.upper[k]:tree.lower[k]));c.side=s?1:-1;
    }
    return {crossing(tree.registry,tree.weight,line[0],line[2]),crossing(tree.registry,tree.weight,line[1],line[2]),crossing(tree.registry,tree.weight,line[1],line[3]),crossing(tree.registry,tree.weight,line[0],line[3])};
}
inline Q exclusion_margin(const Spatial&tree,size_t site,const std::vector<Corner>&polygon,const Node&node){
    Q lower[2]={HUGE_VALQ,HUGE_VALQ},upper[2]={-HUGE_VALQ,-HUGE_VALQ},radius2=0;
    const auto&q=tree.registry.site[site];
    for(const auto&v:polygon){
        S norm;for(int k=0;k<2;k++){
            lower[k]=std::min(lower[k],sf::lower(v.x[k]));upper[k]=std::max(upper[k],sf::upper(v.x[k]));
            S d=v.x[k]-compact_source::canonical_position(q,k);norm=norm+d*d;
        }radius2=std::max(radius2,sf::upper(norm));
    }
    // For EVERY x in this convex polygon and EVERY generator in the node:
    // power_j(x)-power_i(x) >= dist(box(P),box(q_node))² - R_i² + Wi-Wmax.
    // Future clipping only shrinks P, so this strict exclusion remains valid.
    S distance2;for(int k=0;k<2;k++){
        S gap;if(node.lower[k]>upper[k])gap=S(node.lower[k])-S(upper[k]);
        else if(lower[k]>node.upper[k])gap=S(lower[k])-S(node.upper[k]);
        distance2=distance2+gap*gap;
    }
    return sf::lower(distance2-S(radius2)+S(Q(tree.weight[site]))-S(node.maximum_weight));
}
inline BoundedCell bounded_cell(const Spatial&tree,size_t index,const CellBudget&budget){
    require(index<tree.registry.site.size()&&budget.maximum_nodes>0&&budget.maximum_nodes<=262144&&budget.maximum_planes>0&&budget.maximum_planes<=71080&&budget.maximum_corners>=4&&budget.maximum_corners<=4096,"explicit finite power-construction budgets");
    BoundedCell out;out.cell=tree.registry.site[index].cell;out.source_registry_hash=tree.registry.wire_hash;out.weight_hash=tree.weight_hash;out.corner=rectangle(tree);
    auto distance=[&](const Node&node){
        // Traversal priority ONLY. This approximate distance never prunes.
        Q d=0;for(int k=0;k<2;k++){Q x=Q(tree.registry.site[index].position[k]);Q delta=std::max({node.lower[k]-x,x-node.upper[k],Q(0)});d+=delta*delta;}return d;
    };
    struct Visit {Q priority;size_t node;};
    auto farther=[](const Visit&a,const Visit&b){
        return a.priority>b.priority||(a.priority==b.priority&&a.node>b.node);
    };
    std::priority_queue<Visit,std::vector<Visit>,decltype(farther)>pending(farther);
    pending.push({distance(tree.node[0]),0});
    while(!pending.empty()&&!out.corner.empty()){
        size_t id=pending.top().node;pending.pop();require(out.nodes++<budget.maximum_nodes,"bounded spatial power-cell traversal");const auto&node=tree.node[id];
        Q margin=exclusion_margin(tree,index,out.corner,node);if(margin>0){out.excluded.push_back({id,margin});continue;}
        if(node.left!=SIZE_MAX){
            pending.push({distance(tree.node[node.left]),node.left});
            pending.push({distance(tree.node[node.right]),node.right});
        }else for(size_t p=node.begin;p<node.end&&!out.corner.empty();p++){
            size_t j=tree.order[p];if(j==index)continue;require(out.tested_sites.size()<budget.maximum_planes,"bounded actual power-plane intersections");out.tested_sites.push_back(j);
            Constraint line;line.plane=plane(tree.registry.site[index],tree.weight[index],tree.registry.site[j],tree.weight[j]);line.key={1,line.plane.cell[0],line.plane.cell[1]};line.side=out.cell==line.plane.cell[0]?1:-1;
            clip(tree.registry,tree.weight,out.corner,line,budget.maximum_corners);
        }
    }
    if(!out.corner.empty()){
        S twice_area;const auto&origin=out.corner[0].x;
        for(size_t j=1;j+1<out.corner.size();j++){
            const auto&a=out.corner[j].x;const auto&b=out.corner[j+1].x;
            twice_area=twice_area+(a[0]-origin[0])*(b[1]-origin[1])-(a[1]-origin[1])*(b[0]-origin[0]);
        }
        require(sf::lower(twice_area)>0,"certified positive oriented bounded power cell area");
    }return out;
}
} // namespace compact_power
#endif
