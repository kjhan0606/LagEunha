#ifndef MATERIAL_COMPACT_CHART_CERTIFICATE_H
#define MATERIAL_COMPACT_CHART_CERTIFICATE_H
// Global chart connectivity and exact shared-edge provenance. This establishes
// a conforming connected chart complex; curved-J positivity and degree-one
// physical coverage remain separate required certificates.
#include "material_compact_chart_topology.h"

namespace compact_chart {
inline bool same_ball(const S&a,const S&b){return a.n==b.n&&a.e==b.e;}
inline bool contains_exact(const S&x,Q value){return sf::lower(x)<=value&&sf::upper(x)>=value;}
struct PatchRef {size_t candidate=0,patch=0;uint64_t cell=0;};
struct EdgeUse {
    size_t flat_patch=0,local_edge=0;uint64_t cell=0;Key from,to;
    CurveEdge geometry;
};
struct Certificate {
    uint64_t source_registry_hash=0,weight_hash=0,source_epoch=0,source_steps=0,law_hash=0;
    uint64_t chart_hash=0;
    double source_time=0;
    std::vector<PatchRef>patch;
    std::vector<uint64_t>global_component;
    size_t straight_shared_edges=0,internal_edges=0,curved_boundary_segments=0,global_components=0;
    size_t unique_nodes=0,unique_edges=0,boundary_edges=0,vertex_stars=0;
    int64_t euler_characteristic=0;
    size_t geometry_modes=0;
    bool exact_C0_connectivity=false,oriented_source_boundary_complete=false,manifold_vertex_stars=false,disk_topology=false;
    bool curved_J_positive=false,degree_one_coverage=false;
};
struct DisjointSet {
    std::vector<size_t>parent,rank;
    explicit DisjointSet(size_t n):parent(n),rank(n){std::iota(parent.begin(),parent.end(),size_t(0));}
    size_t find(size_t x){if(parent[x]!=x)parent[x]=find(parent[x]);return parent[x];}
    void join(size_t a,size_t b){a=find(a);b=find(b);if(a==b)return;if(rank[a]<rank[b])std::swap(a,b);parent[b]=a;if(rank[a]==rank[b])rank[a]++;}
};
struct EdgeKey {Key a,b;bool operator<(const EdgeKey&x)const{return std::tie(a,b)<std::tie(x.a,x.b);}};
struct BoundaryUse {size_t use=0;Key from,to;CurveEdge geometry;};

inline void chart_hash_value(uint64_t&h,const void*data,size_t bytes){h=compact_p3::bytes_hash(data,bytes,h);}
template<class T> inline void chart_hash_value(uint64_t&h,const T&value){chart_hash_value(h,&value,sizeof(value));}
inline void chart_hash_key(uint64_t&h,const Key&key){chart_hash_value(h,key.size());for(uint64_t word:key)chart_hash_value(h,word);}
inline void chart_hash_ball(uint64_t&h,const S&value){chart_hash_value(h,value.n);chart_hash_value(h,value.e);}
inline uint64_t chart_collection_hash(const std::vector<Candidate>&candidate){
    uint64_t h=14695981039346656037ULL;chart_hash_value(h,candidate.size());
    for(const auto&c:candidate){
        chart_hash_value(h,c.cell);chart_hash_value(h,c.source_registry_hash);chart_hash_value(h,c.weight_hash);
        chart_hash_value(h,c.source_epoch);chart_hash_value(h,c.source_steps);chart_hash_value(h,c.source_time);
        chart_hash_value(h,c.components);chart_hash_value(h,c.patch.size());chart_hash_value(h,c.curved_midpoints);
        for(const auto&p:c.patch){
            chart_hash_value(h,p.component);
            for(const auto&v:p.vertex){chart_hash_key(h,v.key);for(const auto&x:v.x)chart_hash_ball(h,x);}
            for(const auto&e:p.edge){
                chart_hash_value(h,e.boundary);chart_hash_value(h,e.curved);chart_hash_value(h,e.line.kind);
                chart_hash_value(h,e.line.a);chart_hash_value(h,e.line.b);chart_hash_value(h,e.source);
                for(const auto&x:e.parameter)chart_hash_ball(h,x);
                chart_hash_key(h,e.endpoint[0]);chart_hash_key(h,e.endpoint[1]);
            }
        }
    }
    return h;
}

inline Certificate certify_global_connectivity(const compact_target::Registry&registry,
        const std::vector<Candidate>&candidate,uint64_t law_hash,size_t maximum_patches){
    require(candidate.size()==registry.site.size()&&!candidate.empty()&&maximum_patches>0,"one chart candidate per complete global target-site registry");
    Certificate out;out.source_registry_hash=registry.wire_hash;out.weight_hash=candidate.front().weight_hash;
    out.source_epoch=registry.source_epoch;out.source_steps=registry.source_steps;out.source_time=registry.source_time;out.law_hash=law_hash;
    std::map<uint64_t,size_t>candidate_cell;
    size_t total_patch=0;
    for(size_t j=0;j<candidate.size();j++){
        const auto&c=candidate[j];require(c.source_registry_hash==out.source_registry_hash&&c.weight_hash==out.weight_hash&&
            c.source_epoch==out.source_epoch&&c.source_steps==out.source_steps&&c.source_time==out.source_time&&
            c.cell==registry.site[j].cell&&(j==0||candidate[j-1].cell<c.cell),"sorted complete chart candidates share the exact source and immutable W epoch");
        require(candidate_cell.emplace(c.cell,j).second&&!c.patch.empty(),"one nonempty chart candidate per transported generator");
        for(const auto&p:c.patch){
            require(p.component<c.components,"local clipped-graph component retained by its own chart patches");
            require(++total_patch<=maximum_patches,"finite global chart-patch budget");
            for(int k=0;k<3;k++)require(p.vertex[k].key!=p.vertex[(k+1)%3].key,"distinct implicit chart node identities");
        }
    }
    out.patch.reserve(total_patch);out.global_component.resize(total_patch);
    std::map<Key,std::array<S,2>>node_function;
    std::map<EdgeKey,std::vector<EdgeUse>>edge;
    std::map<size_t,std::vector<BoundaryUse>>curved;
    std::map<Key,std::map<Key,std::vector<Key>>>vertex_link;
    for(size_t c=0;c<candidate.size();c++)for(size_t p=0;p<candidate[c].patch.size();p++){
        const auto&patch=candidate[c].patch[p];size_t flat=out.patch.size();out.patch.push_back({c,p,candidate[c].cell});
        for(const auto&node:patch.vertex){
            auto inserted=node_function.emplace(node.key,node.x);
            if(!inserted.second)for(int k=0;k<2;k++)require(same_ball(inserted.first->second[k],node.x[k]),"shared global chart node retains one canonical function enclosure");
        }
        for(size_t k=0;k<3;k++){
            const Key&center=patch.vertex[k].key;const Key&left=patch.vertex[(k+1)%3].key;const Key&right=patch.vertex[(k+2)%3].key;
            vertex_link[center][left].push_back(right);vertex_link[center][right].push_back(left);
        }
        size_t curved_count=0;
        for(size_t k=0;k<3;k++){
            const auto&description=patch.edge[k];const Key&from=patch.vertex[k].key;const Key&to=patch.vertex[(k+1)%3].key;
            require(!description.curved||description.boundary,"curved trace must be a physical boundary edge");
            if(description.curved){
                curved_count++;require(description.source<registry.exterior.size()&&description.endpoint[0]==from&&description.endpoint[1]==to&&
                    cp::strict_sign(description.parameter[1]-description.parameter[0])!=0,
                    "curved chart edge keeps its exact source function, oriented endpoint IDs and nonzero parameter span");
                BoundaryUse b;b.use=flat;b.from=from;b.to=to;b.geometry=description;curved[description.source].push_back(b);
                out.curved_boundary_segments++;
            }else if(description.boundary){
                require(description.line.kind==1&&(description.line.a==candidate[c].cell||description.line.b==candidate[c].cell),"straight physical boundary uses this target cell's canonical sorted neighbor plane");
            }else require(description.source==SIZE_MAX,"internal chart edge carries no exterior curve source");
            Key a=from,b=to;if(b<a)std::swap(a,b);edge[{a,b}].push_back({flat,k,candidate[c].cell,from,to,description});
        }
        require(curved_count<=1,"selected compact chart has at most one physical curved boundary edge");
    }

    DisjointSet sets(total_patch);
    std::map<Key,size_t>boundary_degree;
    for(const auto&item:edge){
        const auto&uses=item.second;require(uses.size()==1||uses.size()==2,"chart edge has one exterior use or exactly one conforming pair");
        if(uses.size()==1){require(uses[0].geometry.boundary&&uses[0].geometry.curved,"unpaired straight edge would leave a coverage gap");boundary_degree[uses[0].from]++;boundary_degree[uses[0].to]++;out.boundary_edges++;continue;}
        const auto&a=uses[0];const auto&b=uses[1];
        require(a.from==b.to&&a.to==b.from&&a.geometry.curved==b.geometry.curved&&a.geometry.boundary==b.geometry.boundary,
            "shared chart traces have identical implicit endpoints and opposite orientation");
        require(!a.geometry.curved,"curved source boundary is owned exactly once");
        if(a.geometry.boundary){
            require(a.cell!=b.cell&&a.geometry.line.kind==1&&a.geometry.line==b.geometry.line,
                "two target charts share the same exact sorted-pair power-plane trace");out.straight_shared_edges++;
        }else out.internal_edges++;
        sets.join(a.flat_patch,b.flat_patch);
    }

    require(curved.size()==registry.exterior.size()&&registry.boundary_cycle.size()==registry.exterior.size(),"every original exterior curve participates in one oriented chart-boundary cycle");
    for(size_t j=0;j<registry.exterior.size();j++){
        auto found=curved.find(j);require(found!=curved.end()&&!found->second.empty(),"complete retained curved source-edge coverage");
        const auto&e=registry.exterior[j];
        std::map<Key,size_t>next;std::map<Key,size_t>previous;
        for(size_t k=0;k<found->second.size();k++){
            const auto&segment=found->second[k];const auto&descriptor=segment.geometry;
            const int direction=e.orientation==1?1:-1;int pdir=cp::strict_sign(descriptor.parameter[1]-descriptor.parameter[0]);
            require((e.orientation==1||e.orientation==-1)&&pdir==direction,"curve parameter direction matches the original oriented exterior edge");
            require(next.emplace(segment.from,k).second&&previous.emplace(segment.to,k).second,"one incoming/outgoing curved segment at each source-boundary event");
            bool first=direction>0?contains_exact(descriptor.parameter[0],0):contains_exact(descriptor.parameter[0],1);
            bool last=direction>0?contains_exact(descriptor.parameter[1],1):contains_exact(descriptor.parameter[1],0);
            if(segment.from==graph_key(cd::endpoint_key(e,direction>0?0:1)))require(first,"original curve start enclosure contains its exact endpoint parameter");
            if(segment.to==graph_key(cd::endpoint_key(e,direction>0?1:0)))require(last,"original curve end enclosure contains its exact endpoint parameter");
        }
        Key start=graph_key(cd::endpoint_key(e,e.orientation==1?0:1)),finish=graph_key(cd::endpoint_key(e,e.orientation==1?1:0));
        require(next.count(start)==1&&previous.count(finish)==1,"curve chart chain starts and ends at the same canonical source vertices");
        size_t visited=0;Key current=start;S parameter=e.orientation==1?S(0):S(1);
        while(current!=finish){auto at=next.find(current);require(at!=next.end()&&visited<found->second.size(),"no gap or branch in the exact curved boundary chain");
            const auto&segment=found->second[at->second];const auto&description=segment.geometry;
            if(visited==0)require(contains_exact(description.parameter[0],e.orientation==1?0:1),"source curve chain starts at the exact oriented endpoint");
            else require(same_ball(description.parameter[0],parameter),"successive source-curve intervals meet at the same exact implicit parameter function");
            parameter=description.parameter[1];current=segment.to;visited++;
        }
        require(visited==found->second.size()&&contains_exact(parameter,e.orientation==1?1:0),"whole original exterior curve is covered exactly once without a discarded interval");
    }
    for(size_t j=0;j<registry.boundary_cycle.size();j++){
        const auto&e=registry.exterior[registry.boundary_cycle[j]];const auto&n=registry.exterior[registry.boundary_cycle[(j+1)%registry.boundary_cycle.size()]];
        uint64_t end=e.vertex[e.orientation==1?1:0],start=n.vertex[n.orientation==1?0:1];require(end==start,"successive source curves share the exact accepted boundary vertex");
    }
    std::map<size_t,uint64_t>component_id;uint64_t next_component=0;
    for(size_t j=0;j<total_patch;j++){size_t root=sets.find(j);auto q=component_id.emplace(root,next_component);if(q.second)next_component++;out.global_component[j]=q.first->second;}
    out.global_components=next_component;require(out.global_components==1,"the complete accepted square image has one edge-connected global chart component");
    out.unique_nodes=node_function.size();out.unique_edges=edge.size();
    out.euler_characteristic=int64_t(out.unique_nodes)-int64_t(out.unique_edges)+int64_t(total_patch);
    for(const auto&node:node_function){
        auto star=vertex_link.find(node.first);require(star!=vertex_link.end()&&!star->second.empty(),"every chart node has a nonempty incident triangle star");
        bool boundary=boundary_degree[node.first]!=0;
        require(!boundary||boundary_degree[node.first]==2,"each boundary node has exactly two distinct exterior edges");
        size_t degree_one=0;for(const auto&neighbor:star->second){size_t degree=neighbor.second.size();
            if(degree==1)degree_one++;else require(degree==2,"vertex link has only path/cycle degrees");
        }
        require((boundary&&degree_one==2&&star->second.size()>=2)||(!boundary&&degree_one==0&&star->second.size()>=3),
            "interior vertex link is a cycle and boundary vertex link is a path");
        std::set<Key>visited;std::vector<Key>stack={star->second.begin()->first};
        while(!stack.empty()){Key current=stack.back();stack.pop_back();if(!visited.insert(current).second)continue;
            for(const auto&next:star->second.at(current))stack.push_back(next);
        }
        require(visited.size()==star->second.size(),"each vertex star has one connected link, not a pinch or multiple fan");out.vertex_stars++;
    }
    out.manifold_vertex_stars=out.vertex_stars==out.unique_nodes;
    require(out.manifold_vertex_stars&&out.boundary_edges>0&&out.euler_characteristic==1,
        "connected orientable chart complex is a topological disk with one complete exterior boundary");
    out.disk_topology=true;
    out.geometry_modes=2*out.global_components;out.chart_hash=chart_collection_hash(candidate);
    out.exact_C0_connectivity=true;out.oriented_source_boundary_complete=true;
    return out;
}
inline void require_same_chart_complex(const compact_target::Registry&registry,
        const std::vector<Candidate>&candidate,const Certificate&certificate){
    require(certificate.source_registry_hash==registry.wire_hash&&certificate.source_epoch==registry.source_epoch&&
        certificate.source_steps==registry.source_steps&&certificate.source_time==registry.source_time&&
        candidate.size()==registry.site.size()&&!candidate.empty()&&candidate.front().weight_hash==certificate.weight_hash&&
        certificate.chart_hash==chart_collection_hash(candidate),
        "same immutable source, W epoch and byte-sealed chart candidate set as the connectivity certificate");
}
}
#endif
