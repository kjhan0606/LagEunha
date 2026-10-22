#include "material_compact_source_spatial.h"
#include <iostream>
#include <string>
using namespace compact_source_spatial;
template<class F>static void rejects(F f,const char*name){
    try{f();throw std::runtime_error(std::string("unexpected pass: ")+name);}
    catch(const std::runtime_error&e){if(std::string(e.what()).find("unexpected pass:")!=std::string::npos)throw;}
}
static Box box(R x0,R y0,R x1,R y1){Box b;b.lo[0]=x0;b.lo[1]=y0;b.hi[0]=x1;b.hi[1]=y1;return b;}
int main(){
    std::vector<Entry>source={{107,2,1,box(1,1,3,2)},{101,0,0,box(-1,-1,1,1)},
        {111,3,1,box(3,0,5,3)},{103,1,0,box(0,0,2,2)}};
    Budget budget;budget.maximum_entries=8;budget.maximum_nodes=15;
    budget.maximum_results_per_query=4;budget.leaf_size=1;budget.maximum_depth=8;budget.ranks=2;
    Index index(source,budget);
    require(index.size()==4&&index.nodes()==7,"four source boxes and complete binary BVH");
    uint64_t checks=0;
    for(int a=0;a<64;a++)for(int b=0;b<64;b++){
        R x=R(a)/10-1.5L,y=R(b)/10-1.5L;
        Box q=box(x,y,x+.37L,y+.29L);
        std::vector<uint64_t>direct,found;
        for(const auto&e:source)if(intersects(e.physical,q))direct.push_back(e.original_triangle);
        std::sort(direct.begin(),direct.end());
        for(const auto&e:index.query(q))found.push_back(e.original_triangle);
        require(found==direct,"BVH broad phase has no false negative against all source AABBs");
        checks++;
    }
    auto bad=source;bad[1].original_triangle=bad[0].original_triangle;
    rejects([&]{Index forbidden(bad,budget);},"duplicate original source ID");
    bad=source;bad[0].physical=box(1,1,0,2);
    rejects([&]{Index forbidden(bad,budget);},"invalid source enclosure");
    Budget small=budget;small.maximum_results_per_query=1;
    Index limited(source,small);
    rejects([&]{(void)limited.query(box(0,0,2,2));},"candidate budget cannot silently truncate");
    small=budget;small.maximum_entries=3;
    rejects([&]{Index forbidden(source,small);},"finite global entry budget");
    small=budget;small.maximum_nodes=6;
    rejects([&]{Index forbidden(source,small);},"finite BVH node budget");
    std::cout<<"SOURCE_PHYSICAL_AABB_BVH_BROAD_PHASE_PASS_SYNTHETIC_ONLY\n"
        <<"entries="<<index.size()<<" nodes="<<index.nodes()<<" queries="<<checks
        <<" actual_source_imported=false overlap_constructed=false\n";
}
