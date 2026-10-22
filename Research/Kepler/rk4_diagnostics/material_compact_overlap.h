#ifndef MATERIAL_COMPACT_OVERLAP_H
#define MATERIAL_COMPACT_OVERLAP_H
// Provenance and positivity contract for a future conservative source/target
// overlap quadrature. This file does not construct geometric intersections,
// project fields, or certify continuous physical coverage.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace compact_overlap {
using R=long double;
inline void require(bool ok,const char*message){if(!ok)throw std::runtime_error(message);}

// One positive physical-area quadrature point shared by the outgoing source
// and incoming target operators. Both coordinate records must enclose the
// same physical location. Barycentrics are local to the source triangle and
// target chart patch, respectively.
struct Point {
    uint64_t id=UINT64_MAX,source_triangle=UINT64_MAX,target_cell=UINT64_MAX,target_patch=UINT64_MAX;
    int source_owner=-1,target_owner=-1;
    double source_l[3]={},target_l[3]={};
    double source_x[2]={},target_x[2]={};
    double source_x_error[2]={},target_x_error[2]={};
    R physical_area_weight=0,physical_area_error=0;
};

struct Seal {
    uint64_t source_coefficient_hash=0,source_registry_hash=0,source_epoch=0,source_steps=0;
    double source_time=0;
    uint64_t target_chart_hash=0,weight_hash=0,mass_law_hash=0;
    uint64_t source_ids_hash=0,target_ids_hash=0,patch_cells_hash=0;
    uint64_t source_triangles=0,target_cells=0,target_patches=0;
    int ranks=0;
};

// Stable byte-order-independent identity of the exact ordered ID list. These
// hashes are transaction consistency guards, not collision-resistant digests.
inline uint64_t id_hash(const std::vector<uint64_t>&ids){
    uint64_t h=14695981039346656037ULL;
    auto word=[&](uint64_t v){for(int k=0;k<8;k++){h^=(v>>(8*k))&255U;h*=1099511628211ULL;}};
    word(ids.size());for(uint64_t id:ids)word(id);return h;
}

struct Budget {
    uint64_t maximum_points=0;
    R maximum_total_area_error=0;
    R maximum_coordinate_gap=0;
};

struct Marginals {
    std::vector<R>source_triangle_area,target_patch_area;
    std::vector<R>source_triangle_error,target_patch_error;
    R total_area=0,total_error=0,maximum_coordinate_gap=0;
    uint64_t points=0;
    bool positive_packet_measure=false;
    bool source_target_coordinates_enclose_same_points=false;
    bool source_target_area_marginals_balance=false;
    // A finite packet measure cannot prove a continuous partition by itself.
    bool continuous_physical_coverage=false;
};

inline bool barycentric(const double l[3]){
    R sum=0;
    for(int k=0;k<3;k++){if(!std::isfinite(l[k])||l[k]<0||l[k]>1)return false;sum+=R(l[k]);}
    return fabsl(sum-1)<=64*std::numeric_limits<double>::epsilon();
}

inline void require_same_seal(const Seal&expected,const Seal&actual){
    require(expected.source_coefficient_hash==actual.source_coefficient_hash&&
        expected.source_registry_hash==actual.source_registry_hash&&expected.source_epoch==actual.source_epoch&&
        expected.source_steps==actual.source_steps&&expected.source_time==actual.source_time&&
        expected.target_chart_hash==actual.target_chart_hash&&expected.weight_hash==actual.weight_hash&&
        expected.mass_law_hash==actual.mass_law_hash&&
        expected.source_ids_hash==actual.source_ids_hash&&expected.target_ids_hash==actual.target_ids_hash&&
        expected.patch_cells_hash==actual.patch_cells_hash&&expected.source_triangles==actual.source_triangles&&
        expected.target_cells==actual.target_cells&&expected.target_patches==actual.target_patches&&
        expected.ranks==actual.ranks,"packet source/target/mass-law provenance matches the independent transaction expectation");
}

class Validator {
public:
    struct Indices {size_t source=0,target=0;};
private:
    Seal seal_;
    Budget budget_;
    std::vector<uint64_t>source_ids_,target_ids_,patch_cells_;
    std::vector<R>source_area_,source_error_,target_area_,target_error_;
    Marginals out_;
    uint64_t previous_=0;
    bool finished_=false;
public:
    Validator(const Seal&expected_seal,const Seal&seal,
            const std::vector<uint64_t>&source_triangle_ids,const std::vector<uint64_t>&target_cell_ids,
            const std::vector<uint64_t>&target_patch_cell_ids,
            const std::vector<R>&expected_source_area,const std::vector<R>&expected_source_area_error,
            const std::vector<R>&expected_target_area,const std::vector<R>&expected_target_area_error,
            const Budget&budget):seal_(seal),budget_(budget),source_ids_(source_triangle_ids),
            target_ids_(target_cell_ids),patch_cells_(target_patch_cell_ids),
            source_area_(expected_source_area),source_error_(expected_source_area_error),
            target_area_(expected_target_area),target_error_(expected_target_area_error){
        require_same_seal(expected_seal,seal_);
        require(seal_.source_coefficient_hash&&seal_.source_registry_hash&&seal_.source_epoch&&seal_.source_steps&&
            seal_.source_time>0&&std::isfinite(seal_.source_time),"complete accepted-source identity and physical epoch");
        require(seal_.target_chart_hash&&seal_.weight_hash&&seal_.mass_law_hash,"sealed target chart, weight and fixed-mass law");
        require(seal_.source_triangles>0&&seal_.target_cells>0&&seal_.target_patches>0&&seal_.ranks>0,
            "nonempty bounded source/target ownership");
        require(source_area_.size()==seal_.source_triangles&&source_error_.size()==seal_.source_triangles&&
            target_area_.size()==seal_.target_patches&&target_error_.size()==seal_.target_patches,
            "one independent expected area and error enclosure per source triangle and target chart patch");
        auto valid_ids=[](const std::vector<uint64_t>&ids){return !ids.empty()&&ids.back()!=UINT64_MAX&&
            std::adjacent_find(ids.begin(),ids.end(),[](uint64_t a,uint64_t b){return a>=b;})==ids.end();};
        require(source_ids_.size()==seal_.source_triangles&&valid_ids(source_ids_)&&
            target_ids_.size()==seal_.target_cells&&valid_ids(target_ids_)&&
            patch_cells_.size()==seal_.target_patches,
            "independent sorted unique original source/target IDs and one target cell per chart patch");
        for(uint64_t cell:patch_cells_)require(std::binary_search(target_ids_.begin(),target_ids_.end(),cell),
            "every target chart patch belongs to a registered original cell");
        require(seal_.source_ids_hash&&seal_.target_ids_hash&&seal_.patch_cells_hash&&
            seal_.source_ids_hash==id_hash(source_ids_)&&seal_.target_ids_hash==id_hash(target_ids_)&&
            seal_.patch_cells_hash==id_hash(patch_cells_),
            "exact source/target original ID lists and patch-to-cell map match the transaction seal");
        require(budget_.maximum_points>0&&budget_.maximum_total_area_error>=0&&
            std::isfinite(budget_.maximum_total_area_error)&&budget_.maximum_coordinate_gap>=0&&
            std::isfinite(budget_.maximum_coordinate_gap),"finite nonnegative overlap proof budgets");
        out_.source_triangle_area.assign(source_ids_.size(),0);out_.source_triangle_error.assign(source_ids_.size(),0);
        out_.target_patch_area.assign(patch_cells_.size(),0);out_.target_patch_error.assign(patch_cells_.size(),0);
    }
    Indices append(const Point&p){
        require(!finished_&&out_.points<budget_.maximum_points,"streaming overlap packet remains open and under its finite point budget");
        require(p.id!=UINT64_MAX&&(!out_.points||previous_<p.id),"strictly ordered unique packet point IDs");
        auto source_it=std::lower_bound(source_ids_.begin(),source_ids_.end(),p.source_triangle);
        auto target_it=std::lower_bound(target_ids_.begin(),target_ids_.end(),p.target_cell);
        require(source_it!=source_ids_.end()&&*source_it==p.source_triangle&&
            target_it!=target_ids_.end()&&*target_it==p.target_cell&&p.target_patch<seal_.target_patches,
            "packet original source/target identities and target chart patch are registered");
        require(patch_cells_[p.target_patch]==p.target_cell,"packet target chart patch belongs to its original target cell");
        require(p.source_owner>=0&&p.source_owner<seal_.ranks&&p.target_owner>=0&&p.target_owner<seal_.ranks,
            "packet source/target owners are valid ranks");
        require(barycentric(p.source_l)&&barycentric(p.target_l),"source and target local coordinates lie in their closed chart triangles");
        require(p.physical_area_weight>0&&std::isfinite(p.physical_area_weight)&&
            p.physical_area_error>=0&&std::isfinite(p.physical_area_error)&&p.physical_area_error<p.physical_area_weight,
            "strictly positive physical-area packet weight with a finite smaller error bound");
        R gap2=0;
        for(int d=0;d<2;d++){
            require(std::isfinite(p.source_x[d])&&std::isfinite(p.target_x[d])&&
                p.source_x_error[d]>=0&&std::isfinite(p.source_x_error[d])&&
                p.target_x_error[d]>=0&&std::isfinite(p.target_x_error[d]),"finite source/target physical-point enclosures");
            require(R(p.source_x_error[d])+p.target_x_error[d]<=budget_.maximum_coordinate_gap,
                "source and target coordinate uncertainty fits the finite shared-point budget");
            R gap=fabsl(R(p.source_x[d])-p.target_x[d])-R(p.source_x_error[d])-p.target_x_error[d];
            if(gap>0)gap2+=gap*gap;
        }
        R gap=sqrtl(gap2);require(gap<=budget_.maximum_coordinate_gap,
            "source and target inverse maps enclose the same physical quadrature point");
        const size_t s=size_t(source_it-source_ids_.begin()),t=size_t(target_it-target_ids_.begin());
        out_.maximum_coordinate_gap=std::max(out_.maximum_coordinate_gap,gap);
        out_.source_triangle_area[s]+=p.physical_area_weight;out_.target_patch_area[p.target_patch]+=p.physical_area_weight;
        out_.source_triangle_error[s]+=p.physical_area_error;out_.target_patch_error[p.target_patch]+=p.physical_area_error;
        out_.total_area+=p.physical_area_weight;out_.total_error+=p.physical_area_error;
        out_.points++;previous_=p.id;
        require(std::isfinite(out_.total_area)&&std::isfinite(out_.total_error),"finite streaming overlap area ledger");
        require(out_.total_error<=budget_.maximum_total_area_error,"streaming overlap quadrature error remains inside its finite budget");
        return {s,t};
    }
    Marginals finish(){
        require(!finished_&&out_.points>0,"nonempty overlap stream finalized exactly once");finished_=true;
        auto close=[](R got,R error,R expected){return std::isfinite(got)&&std::isfinite(error)&&std::isfinite(expected)&&expected>0&&
            fabsl(got-expected)<=error+64*std::numeric_limits<R>::epsilon()*std::max(fabsl(got),fabsl(expected));};
        R source_total=0,target_total=0,source_error=0,target_error=0;bool source_ok=true,target_ok=true;
        for(size_t i=0;i<source_area_.size();i++){
            require(source_area_[i]>0&&std::isfinite(source_area_[i])&&source_error_[i]>=0&&
                std::isfinite(source_error_[i]),"positive independently assembled source-triangle area and finite error enclosure");
            source_total+=source_area_[i];source_error+=out_.source_triangle_error[i]+source_error_[i];
            source_ok&=close(out_.source_triangle_area[i],out_.source_triangle_error[i]+source_error_[i],source_area_[i]);
        }
        for(size_t i=0;i<target_area_.size();i++){
            require(target_area_[i]>0&&std::isfinite(target_area_[i])&&target_error_[i]>=0&&
                std::isfinite(target_error_[i]),"positive independently assembled target-chart area and finite error enclosure");
            target_total+=target_area_[i];target_error+=out_.target_patch_error[i]+target_error_[i];
            target_ok&=close(out_.target_patch_area[i],out_.target_patch_error[i]+target_error_[i],target_area_[i]);
        }
        require(source_error+target_error-out_.total_error<=budget_.maximum_total_area_error,
            "packet and independent source/target area errors fit the finite transaction budget");
        R common_error=out_.total_error+source_error+target_error;
        out_.source_target_area_marginals_balance=source_ok&&target_ok&&close(source_total,common_error,target_total);
        require(out_.source_target_area_marginals_balance,"both packet marginals match independent source/target areas within recorded integration errors");
        out_.positive_packet_measure=true;out_.source_target_coordinates_enclose_same_points=true;
        // A finite packet measure does not prove continuous domain coverage.
        out_.continuous_physical_coverage=false;
        return out_;
    }
};

inline Marginals validate(const Seal&expected_seal,const Seal&seal,const std::vector<Point>&points,
        const std::vector<uint64_t>&source_triangle_ids,const std::vector<uint64_t>&target_cell_ids,
        const std::vector<uint64_t>&target_patch_cell_ids,
        const std::vector<R>&expected_source_area,const std::vector<R>&expected_source_area_error,
        const std::vector<R>&expected_target_area,const std::vector<R>&expected_target_area_error,
        const Budget&budget){
    Validator stream(expected_seal,seal,source_triangle_ids,target_cell_ids,target_patch_cell_ids,
        expected_source_area,expected_source_area_error,expected_target_area,expected_target_area_error,budget);
    for(const Point&p:points)(void)stream.append(p);
    return stream.finish();
}
}
#endif
