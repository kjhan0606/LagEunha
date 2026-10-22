#ifndef MATERIAL_COMPACT_CHART_TOPOLOGY_H
#define MATERIAL_COMPACT_CHART_TOPOLOGY_H
// Astra's ONE fixed connectivity candidate. All coordinates are enclosures
// of canonical functions; this does NOT yet certify curved J or coverage.
#include "material_compact_domain.h"
#include <cstdio>
namespace compact_chart {
namespace cd=compact_domain;namespace cp=compact_power;namespace sf=shared_function;
using Q=sf::Q;using S=sf::S;
inline void require(bool b,const char*s){compact_p3::require(b,s);}
using Key=std::vector<uint64_t>;
struct Node {Key key;std::array<S,2>x;};
struct CurveEdge {
    bool boundary=false,curved=false;cp::LineKey line;size_t source=SIZE_MAX;
    S parameter[2];Key endpoint[2];
    // Parameters retain their defining endpoint/root/midpoint functions in
    // endpoint keys; centers of these balls are NOT new exact geometry.
};
struct Patch {
    Node vertex[3];CurveEdge edge[3];size_t component=0;
};
struct Candidate {
    uint64_t cell=0,source_registry_hash=0,weight_hash=0,source_epoch=0,source_steps=0;double source_time=0;
    std::vector<Patch>patch;
    size_t components=0,curved_midpoints=0,ears=0,collar_fans=0,bulk_fans=0,star_fans=0;
    // Positive chord orientations/ear connectivity only. Physical curved
    // Jacobian, C0/degree-one coverage and mode rank remain mandatory guards.
};
struct Budget {size_t maximum_vertices=0,maximum_ears=0,maximum_patches=0;};
inline Key graph_key(const cd::Key&key){Key out={0,key.kind};out.insert(out.end(),key.id.begin(),key.id.end());return out;}
inline Node graph_node(const cd::Vertex&v){return {graph_key(v.key),v.x};}
inline S orientation(const Node&a,const Node&b,const Node&c){return (b.x[0]-a.x[0])*(c.x[1]-a.x[1])-(b.x[1]-a.x[1])*(c.x[0]-a.x[0]);}
inline bool disjoint(const Node&a,const Node&b,const Node&c,const Node&d){
    // Strict separation of either coordinate projection proves the closed
    // chord segments disjoint, including exact or near-collinear cases.
    // Canonical endpoint balls are used directly; no coordinate snap or
    // tolerance is introduced.
    for(int axis=0;axis<2;axis++){
        Q ab_lo=std::min(sf::lower(a.x[axis]),sf::lower(b.x[axis]));
        Q ab_hi=std::max(sf::upper(a.x[axis]),sf::upper(b.x[axis]));
        Q cd_lo=std::min(sf::lower(c.x[axis]),sf::lower(d.x[axis]));
        Q cd_hi=std::max(sf::upper(c.x[axis]),sf::upper(d.x[axis]));
        if(ab_hi<cd_lo||cd_hi<ab_lo)return true;
    }
    int ab_c=cp::strict_sign(orientation(a,b,c)),ab_d=cp::strict_sign(orientation(a,b,d));
    int cd_a=cp::strict_sign(orientation(c,d,a)),cd_b=cp::strict_sign(orientation(c,d,b));
    if(ab_c&&ab_d&&ab_c==ab_d)return true;if(cd_a&&cd_b&&cd_a==cd_b)return true;
    // Collinear disjointness is deliberately not guessed by a tolerance.
    return false;
}
inline Node centroid(uint64_t cell,const std::vector<Node>&nodes,uint64_t kind){
    std::vector<const Node*>sorted;for(const auto&v:nodes)sorted.push_back(&v);std::sort(sorted.begin(),sorted.end(),[](const Node*a,const Node*b){return a->key<b->key;});
    Node out;out.key={kind,cell,uint64_t(sorted.size())};for(const auto*v:sorted){out.key.push_back(v->key.size());out.key.insert(out.key.end(),v->key.begin(),v->key.end());for(int k=0;k<2;k++)out.x[k]=out.x[k]+v->x[k];}
    for(int k=0;k<2;k++)out.x[k]=out.x[k]/S(Q(sorted.size()));return out;
}
inline void polygon_guard(const std::vector<Node>&v){
    require(v.size()>=3,"positive chord polygon after finite curved midpoint split");std::set<Key>ids;S area;
    for(size_t j=0;j<v.size();j++){
        require(ids.insert(v[j].key).second,"unique canonical node on a component cycle");const auto&a=v[j];const auto&b=v[(j+1)%v.size()];
        S length;for(int k=0;k<2;k++){S d=b.x[k]-a.x[k];length=length+d*d;}require(sf::lower(length)>0,"nonzero certified chord edge");
        if(j&&j+1<v.size())area=area+orientation(v[0],v[j],v[j+1]);
        for(size_t k=j+1;k<v.size();k++)if(k!=(j+1)%v.size()&&j!=(k+1)%v.size())require(disjoint(v[j],v[(j+1)%v.size()],v[k],v[(k+1)%v.size()]),"simple chord polygon, no crossing or unresolved collinear boundary segments");
    }require(sf::lower(area)>0,"positive oriented component chord area; holes are not silently filled");
}
inline std::vector<std::array<size_t,3>> triangulate(const std::vector<Node>&v,const Budget&budget){
    std::vector<size_t>ring(v.size());std::iota(ring.begin(),ring.end(),size_t(0));std::vector<std::array<size_t,3>>out;
    while(ring.size()>3){
        require(out.size()<budget.maximum_ears,"finite certified ear budget");std::vector<size_t>order(ring.size());std::iota(order.begin(),order.end(),size_t(0));std::sort(order.begin(),order.end(),[&](size_t a,size_t b){return v[ring[a]].key<v[ring[b]].key;});bool found=false;
        for(size_t at:order){
            size_t before=(at+ring.size()-1)%ring.size(),after=(at+1)%ring.size();size_t a=ring[before],b=ring[at],c=ring[after];
            if(cp::strict_sign(orientation(v[a],v[b],v[c]))!=1)continue;
            bool ear=true;
            for(size_t p:ring)if(p!=a&&p!=b&&p!=c){
                int s[3]={cp::strict_sign(orientation(v[a],v[b],v[p])),cp::strict_sign(orientation(v[b],v[c],v[p])),cp::strict_sign(orientation(v[c],v[a],v[p]))};
                // Require exclusion of the CLOSED triangle; an uncertain or
                // exact boundary equality is not broken by an epsilon.
                if(s[0]>=0&&s[1]>=0&&s[2]>=0){ear=false;break;}
            }
            if(!ear)continue;
            for(size_t j=0;j<ring.size();j++){
                size_t p=ring[j],q=ring[(j+1)%ring.size()];if(p==a||p==c||q==a||q==c)continue;
                if(!disjoint(v[a],v[c],v[p],v[q])){ear=false;break;}
            }
            // For a simple oriented polygon, convex local wedge + no ear
            // vertex/diagonal crossing puts the diagonal in the interior.
            if(ear){out.push_back({a,b,c});ring.erase(ring.begin()+at);found=true;break;}
        }
        require(found,"no certified ear in fixed canonical-ID order; no refinement/perturbation retry");
    }
    require(out.size()<budget.maximum_ears&&cp::strict_sign(orientation(v[ring[0]],v[ring[1]],v[ring[2]]))==1,"positive last ear, not a collinear or reversed triangle");out.push_back({ring[0],ring[1],ring[2]});return out;
}
inline Candidate connect(const cd::Domain&domain,const cd::Graph&graph,const Budget&budget){
    require(budget.maximum_vertices>=3&&budget.maximum_vertices<=65536&&budget.maximum_ears>0&&budget.maximum_ears<=65536&&budget.maximum_patches>0&&budget.maximum_patches<=262144,"explicit finite chart connectivity budgets");
    require(graph.source_registry_hash==domain.source.wire_hash&&graph.source_epoch==domain.source.source_epoch&&graph.source_steps==domain.source.source_steps&&graph.source_time==domain.source.source_time,"same sealed source for chart graph");
    Candidate out;out.cell=graph.cell;out.source_registry_hash=graph.source_registry_hash;out.weight_hash=graph.weight_hash;out.source_epoch=graph.source_epoch;out.source_steps=graph.source_steps;out.source_time=graph.source_time;out.components=graph.component.size();
    for(size_t component=0;component<graph.component.size();component++){
        const auto&cycle=graph.component[component];std::vector<Node>polygon;std::vector<CurveEdge>edge;
        for(size_t id:cycle){
            const auto&a=graph.arc.at(id);const auto&A=graph.vertex.at(a.from);const auto&B=graph.vertex.at(a.to);polygon.push_back(graph_node(A));
            if(!a.curved){CurveEdge straight;straight.boundary=true;straight.line=a.line;edge.push_back(straight);}
            else{
                require(a.source_curve<domain.curve.size(),"same source curve for chart midpoint");
                S pa=sf::range(a.parameter_lower[0],a.parameter_upper[0]),pb=sf::range(a.parameter_lower[1],a.parameter_upper[1]),mid=(pa+pb)/S(2);
                Node m;Key ka=graph_key(A.key),kb=graph_key(B.key);if(a.parameter_lower[1]<a.parameter_lower[0])std::swap(ka,kb);
                m.key={1,domain.source.exterior[a.source_curve].edge,uint64_t(ka.size())};m.key.insert(m.key.end(),ka.begin(),ka.end());m.key.push_back(kb.size());m.key.insert(m.key.end(),kb.begin(),kb.end());
                auto point=cp::position(domain.curve[a.source_curve],mid);m.x={point[0].v,point[1].v};
                CurveEdge first;first.boundary=true;first.curved=true;first.line=a.line;first.source=a.source_curve;first.parameter[0]=pa;first.parameter[1]=mid;first.endpoint[0]=graph_key(A.key);first.endpoint[1]=m.key;edge.push_back(first);
                CurveEdge second=first;second.parameter[0]=mid;second.parameter[1]=pb;second.endpoint[0]=m.key;second.endpoint[1]=graph_key(B.key);polygon.push_back(m);edge.push_back(second);out.curved_midpoints++;
            }
            require(polygon.size()<=budget.maximum_vertices,"finite one-midpoint-per-arc boundary node budget");
        }
        require(polygon.size()==edge.size(),"one reference edge per canonical component node");polygon_guard(polygon);
        std::map<std::pair<Key,Key>,CurveEdge>boundary;bool curved=false;
        for(size_t j=0;j<polygon.size();j++){
            curved|=edge[j].curved;
            require(edge[j].boundary&&boundary.emplace(std::make_pair(polygon[j].key,polygon[(j+1)%polygon.size()].key),edge[j]).second,"unique complete curved or canonical straight physical boundary edge");
        }
        auto emit=[&](Node a,Node b,Node c,const CurveEdge*physical_edge,
                const S*algebraic_area=nullptr){
            require(out.patch.size()<budget.maximum_patches,"finite patch budget");
            S evaluated=orientation(a,b,c);
            S area=algebraic_area?*algebraic_area:evaluated;
            if(algebraic_area)require(
                std::max(sf::lower(area),sf::lower(evaluated))<=
                    std::min(sf::upper(area),sf::upper(evaluated)),
                "centroid fan algebraic area and direct endpoint enclosure agree");
            if(cp::strict_sign(area)!=1){
                char message[256];std::snprintf(message,sizeof(message),
                    "unresolved positive chord patch cell=%llu component=%zu patch=%zu curved=%d area=[%.21Le,%.21Le]",
                    static_cast<unsigned long long>(graph.cell),component,out.patch.size(),
                    physical_edge&&physical_edge[1].curved?1:0,
                    static_cast<long double>(sf::lower(area)),
                    static_cast<long double>(sf::upper(area)));
                throw std::runtime_error(message);
            }
            Patch p;p.vertex[0]=std::move(a);p.vertex[1]=std::move(b);p.vertex[2]=std::move(c);if(physical_edge)for(int j=0;j<3;j++)p.edge[j]=physical_edge[j];p.component=component;out.patch.push_back(std::move(p));
        };
        if(!curved){
            for(size_t j=0;j<polygon.size();j++)require(cp::strict_sign(orientation(polygon[j],polygon[(j+1)%polygon.size()],polygon[(j+2)%polygon.size()]))==1,"convex bulk polygon before original centroid fan");
            Node center=centroid(graph.cell,polygon,2);for(size_t j=0;j<polygon.size();j++){CurveEdge edge[3]={};edge[1]=boundary.at({polygon[j].key,polygon[(j+1)%polygon.size()].key});emit(center,polygon[j],polygon[(j+1)%polygon.size()],edge);}out.bulk_fans++;continue;
        }
        // A strictly positive orientation from one canonical interior point
        // to EVERY boundary chord proves it lies in the kernel of this simple
        // polygon. Use that whole-component fan before ear decomposition; it
        // avoids almost flat ear centroids at nearly coincident source nodes.
        // Curved images still require the separate whole-chart J/degree gate.
        Node star=centroid(graph.cell,polygon,4);bool star_visible=true;
        for(size_t j=0;j<polygon.size();j++)
            star_visible&=cp::strict_sign(orientation(
                star,polygon[j],polygon[(j+1)%polygon.size()]))==1;
        if(star_visible){
            for(size_t j=0;j<polygon.size();j++){
                CurveEdge fan_edge[3]={};
                fan_edge[1]=boundary.at({polygon[j].key,polygon[(j+1)%polygon.size()].key});
                emit(star,polygon[j],polygon[(j+1)%polygon.size()],fan_edge);
            }
            out.star_fans++;continue;
        }
        auto ears=triangulate(polygon,budget);out.ears+=ears.size();
        for(const auto&t:ears){
            Node v[3]={polygon[t[0]],polygon[t[1]],polygon[t[2]]};CurveEdge physical_edge[3];bool collar=false;
            for(int j=0;j<3;j++){auto found=boundary.find({v[j].key,v[(j+1)%3].key});if(found!=boundary.end()){physical_edge[j]=found->second;collar|=physical_edge[j].curved;}}
            if(!collar)emit(v[0],v[1],v[2],physical_edge);
            else{
                // The canonical center is (v0+v1+v2)/3. Every oriented
                // centroid fan patch therefore has exactly one third of the
                // parent ear's signed chord area. Use this identity when
                // independent vertex balls lose a thin positive lower bound.
                S one_third=orientation(v[0],v[1],v[2])/S(3);
                require(cp::strict_sign(one_third)==1,
                    "positive parent ear gives positive canonical centroid fan area");
                Node center=centroid(graph.cell,{v[0],v[1],v[2]},3);
                for(int j=0;j<3;j++){
                    CurveEdge fan_edge[3]={};fan_edge[1]=physical_edge[j];
                    emit(center,v[j],v[(j+1)%3],fan_edge,&one_third);
                }out.collar_fans++;
            }
        }
    }return out;
}
} // namespace compact_chart
#endif
