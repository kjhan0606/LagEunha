#ifndef MATERIAL_COMPACT_SOURCE_SPATIAL_H
#define MATERIAL_COMPACT_SOURCE_SPATIAL_H
// Conservative broad-phase only: boxes may return false positives, never
// certify that a target patch is covered by any source triangle. Exact
// source inverse/edge routing and overlap integration remain separate gates.
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace compact_source_spatial {
using R=long double;
inline void require(bool good,const char*message){if(!good)throw std::runtime_error(message);}
struct Box {
    R lo[2]={},hi[2]={};
};
inline bool valid(const Box&b){
    for(int d=0;d<2;d++)if(!std::isfinite(b.lo[d])||!std::isfinite(b.hi[d])||b.lo[d]>b.hi[d])return false;
    return true;
}
inline bool intersects(const Box&a,const Box&b){
    return a.lo[0]<=b.hi[0]&&b.lo[0]<=a.hi[0]&&a.lo[1]<=b.hi[1]&&b.lo[1]<=a.hi[1];
}
inline Box hull(const Box&a,const Box&b){
    Box out;for(int d=0;d<2;d++){out.lo[d]=std::min(a.lo[d],b.lo[d]);out.hi[d]=std::max(a.hi[d],b.hi[d]);}
    return out;
}
struct Entry {
    uint64_t original_triangle=UINT64_MAX,local_triangle=UINT64_MAX;
    int owner=-1;
    Box physical;
};
struct Budget {
    size_t maximum_entries=0,maximum_nodes=0,maximum_results_per_query=0;
    size_t leaf_size=8,maximum_depth=64;
    int ranks=0;
};
class Index {
    struct Node {Box box;size_t begin=0,end=0,left=SIZE_MAX,right=SIZE_MAX;};
    std::vector<Entry>entry_;
    std::vector<size_t>order_;
    std::vector<Node>node_;
    Budget budget_;
    size_t build(size_t begin,size_t end,size_t depth){
        require(begin<end&&depth<=budget_.maximum_depth&&node_.size()<budget_.maximum_nodes,
            "finite conservative source BVH node/depth budget");
        const size_t at=node_.size();node_.push_back(Node{});
        Box box=entry_[order_[begin]].physical;
        for(size_t j=begin+1;j<end;j++)box=hull(box,entry_[order_[j]].physical);
        node_[at].box=box;node_[at].begin=begin;node_[at].end=end;
        if(end-begin<=budget_.leaf_size)return at;
        int axis=(box.hi[1]-box.lo[1]>box.hi[0]-box.lo[0])?1:0;
        size_t middle=begin+(end-begin)/2;
        std::nth_element(order_.begin()+begin,order_.begin()+middle,order_.begin()+end,
            [&](size_t a,size_t b){
                R ca=entry_[a].physical.lo[axis]/2+entry_[a].physical.hi[axis]/2;
                R cb=entry_[b].physical.lo[axis]/2+entry_[b].physical.hi[axis]/2;
                return ca<cb||(ca==cb&&entry_[a].original_triangle<entry_[b].original_triangle);
            });
        size_t left=build(begin,middle,depth+1),right=build(middle,end,depth+1);
        node_[at].left=left;node_[at].right=right;return at;
    }
public:
    Index(std::vector<Entry>entry,const Budget&budget):entry_(std::move(entry)),budget_(budget){
        require(budget_.maximum_entries>0&&budget_.maximum_nodes>0&&budget_.maximum_results_per_query>0&&
            budget_.leaf_size>0&&budget_.maximum_depth>0&&budget_.maximum_depth<=64&&budget_.ranks>0&&
            !entry_.empty()&&entry_.size()<=budget_.maximum_entries,"finite source broad-phase budgets and nonempty registry");
        std::sort(entry_.begin(),entry_.end(),[](const Entry&a,const Entry&b){return a.original_triangle<b.original_triangle;});
        for(size_t j=0;j<entry_.size();j++){
            const auto&e=entry_[j];require(e.original_triangle!=UINT64_MAX&&e.local_triangle!=UINT64_MAX&&
                e.owner>=0&&e.owner<budget_.ranks&&valid(e.physical),"finite identified source physical enclosure");
            require(j==0||entry_[j-1].original_triangle<e.original_triangle,"unique original source triangle across MPI owners");
        }
        order_.resize(entry_.size());std::iota(order_.begin(),order_.end(),size_t(0));
        (void)build(0,entry_.size(),0);
    }
    size_t size()const{return entry_.size();}
    size_t nodes()const{return node_.size();}
    std::vector<Entry> query(const Box&target)const{
        require(valid(target),"finite target physical box for source broad phase");
        std::vector<Entry>out;std::vector<size_t>stack={0};
        while(!stack.empty()){
            size_t at=stack.back();stack.pop_back();const auto&node=node_[at];
            if(!intersects(node.box,target))continue;
            if(node.left==SIZE_MAX){
                for(size_t j=node.begin;j<node.end;j++){
                    const auto&e=entry_[order_[j]];
                    if(intersects(e.physical,target)){
                        require(out.size()<budget_.maximum_results_per_query,
                            "finite source candidates per target-box query; no silent truncation");
                        out.push_back(e);
                    }
                }
            }else{stack.push_back(node.right);stack.push_back(node.left);}
        }
        std::sort(out.begin(),out.end(),[](const Entry&a,const Entry&b){return a.original_triangle<b.original_triangle;});
        return out;
    }
};
}
#endif
