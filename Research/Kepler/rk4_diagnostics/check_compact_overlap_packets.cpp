#include "material_compact_moments.h"
#include <iostream>

using namespace compact_overlap;
template<class F>static void expect_reject(F f,const char*name){
    try{f();throw std::runtime_error(std::string("negative fixture unexpectedly passed: ")+name);}
    catch(const std::runtime_error&e){if(std::string(e.what()).find("negative fixture unexpectedly passed")!=std::string::npos)throw;}
}
int main(){
    const std::vector<uint64_t>source_ids={101,103},target_ids={10,20},patch_cells={10,20};
    Seal s;s.source_coefficient_hash=11;s.source_registry_hash=12;s.source_epoch=3;s.source_steps=5;s.source_time=1e-4;
    s.target_chart_hash=21;s.weight_hash=22;s.mass_law_hash=23;s.source_triangles=2;s.target_cells=2;s.target_patches=2;s.ranks=2;
    s.source_ids_hash=id_hash(source_ids);s.target_ids_hash=id_hash(target_ids);s.patch_cells_hash=id_hash(patch_cells);
    Budget budget;budget.maximum_points=8;budget.maximum_total_area_error=1e-10L;budget.maximum_coordinate_gap=1e-14L;
    const std::vector<R>source_area={.5L,.5L},target_area={.5L,.5L};
    const std::vector<R>source_area_error={1e-12L,1e-12L},target_area_error={1e-12L,1e-12L};
    std::vector<Point>p(4);
    const R w[4]={.3L,.2L,.2L,.3L};const size_t st[4]={0,0,1,1},tp[4]={0,1,0,1};
    for(size_t i=0;i<p.size();i++){
        p[i].id=i+1;p[i].source_triangle=source_ids[st[i]];p[i].target_cell=target_ids[tp[i]];p[i].target_patch=tp[i];
        p[i].source_owner=int(st[i]);p[i].target_owner=int(tp[i]);
        for(int k=0;k<3;k++)p[i].source_l[k]=p[i].target_l[k]=1.0/3;
        p[i].source_x[0]=p[i].target_x[0]=double(i);p[i].source_x[1]=p[i].target_x[1]=double(i*i);
        for(int d=0;d<2;d++)p[i].source_x_error[d]=p[i].target_x_error[d]=1e-15;
        p[i].physical_area_weight=w[i];p[i].physical_area_error=1e-12L;
    }
    FieldStream f;f.seal=s;
    for(size_t i=0;i<p.size();i++){
        PhysicalSample a;a.point_id=p[i].id;a.rho=i+1;a.internal_energy_density=2*(i+1);
        a.vx=i+1;a.vy=2*(i+1);a.phi=-1;f.sample.push_back(a);
    }
    auto check=[&](const Seal&expected,const Seal&actual,const std::vector<Point>&points,
            const FieldStream&field,const std::vector<uint64_t>&sources,
            const std::vector<uint64_t>&targets,const std::vector<uint64_t>&patches){
        return accumulate(expected,actual,points,field,sources,targets,patches,
            source_area,source_area_error,target_area,target_area_error,budget);
    };
    auto area=validate(s,s,p,source_ids,target_ids,patch_cells,source_area,
        source_area_error,target_area,target_area_error,budget);
    auto m=check(s,s,p,f,source_ids,target_ids,patch_cells);
    auto make_stream=[&](const Budget&b){return StreamAccumulator(s,s,f.seal,source_ids,target_ids,
        patch_cells,source_area,source_area_error,target_area,target_area_error,b);};
    auto stream=make_stream(budget);
    for(size_t begin=0;begin<p.size();begin+=2)
        for(size_t i=begin;i<std::min(begin+2,p.size());i++)stream.append(p[i],f.sample[i]);
    auto streamed=stream.finish();
    for(size_t j=0;j<NMOMENTS;j++)if(fabsl(m.packet_total[j]-streamed.packet_total[j])>1e-18L)
        throw std::runtime_error("buffered and two-chunk streaming moments disagree");
    if(!area.positive_packet_measure||!area.source_target_coordinates_enclose_same_points||
        !area.source_target_area_marginals_balance||area.continuous_physical_coverage||area.points!=4||
        fabsl(area.total_area-1)>1e-18L||!m.source_target_algebraic_balance||
        !streamed.source_target_algebraic_balance||streamed.area.points!=4||
        m.physical_quadrature_error_certified||fabsl(m.packet_total[M]-2.5L)>1e-18L||
        fabsl(m.packet_total[U]-5)>1e-18L||fabsl(m.packet_total[PHI]+2.5L)>1e-18L||
        fabsl(m.packet_total[K]-66.25L)>1e-16L)
        throw std::runtime_error("synthetic common packet moments fixture contract");
    auto bad=p;bad[0].physical_area_weight=-.3L;
    expect_reject([&]{check(s,s,bad,f,source_ids,target_ids,patch_cells);},"negative physical weight");
    bad=p;bad[1].id=bad[0].id;
    expect_reject([&]{check(s,s,bad,f,source_ids,target_ids,patch_cells);},"duplicate packet ID");
    bad=p;bad[2].target_x[0]+=1e-6;
    expect_reject([&]{check(s,s,bad,f,source_ids,target_ids,patch_cells);},"inconsistent inverse locations");
    bad=p;bad[2].target_x[0]+=1e-6;bad[2].source_x_error[0]=1e-5;
    expect_reject([&]{check(s,s,bad,f,source_ids,target_ids,patch_cells);},"inflated coordinate enclosure");
    bad=p;bad[3].target_l[2]=1.1;
    expect_reject([&]{check(s,s,bad,f,source_ids,target_ids,patch_cells);},"outside target chart");
    bad=p;bad[3].target_patch=4;
    expect_reject([&]{check(s,s,bad,f,source_ids,target_ids,patch_cells);},"out-of-range target patch");
    bad=p;bad[3].physical_area_weight=.1L;
    expect_reject([&]{check(s,s,bad,f,source_ids,target_ids,patch_cells);},"missing packet area");
    bad=p;bad[0].source_triangle=0;
    expect_reject([&]{check(s,s,bad,f,source_ids,target_ids,patch_cells);},"unregistered original source triangle ID");
    bad=p;bad[0].target_cell=0;
    expect_reject([&]{check(s,s,bad,f,source_ids,target_ids,patch_cells);},"unregistered original target cell ID");
    bad=p;bad[0].target_cell=20;
    expect_reject([&]{check(s,s,bad,f,source_ids,target_ids,patch_cells);},"target chart belongs to different cell");
    auto wrong_patches=patch_cells;wrong_patches[1]=30;
    expect_reject([&]{check(s,s,p,f,source_ids,target_ids,wrong_patches);},"unregistered patch owner cell");
    auto wrong_sources=source_ids;wrong_sources[1]=wrong_sources[0];
    expect_reject([&]{check(s,s,p,f,wrong_sources,target_ids,patch_cells);},"duplicate registered source ID");
    wrong_sources=source_ids;wrong_sources[1]=104;bad=p;bad[2].source_triangle=104;bad[3].source_triangle=104;
    expect_reject([&]{check(s,s,bad,f,wrong_sources,target_ids,patch_cells);},"changed source ID list under same seal");
    auto loose_area=source_area_error;loose_area[0]=1;
    expect_reject([&]{accumulate(s,s,p,f,source_ids,target_ids,patch_cells,
        source_area,loose_area,target_area,target_area_error,budget);},"inflated independent area error");
    Seal stale=s;stale.source_epoch++;
    expect_reject([&]{check(s,stale,p,f,source_ids,target_ids,patch_cells);},"stale packet source epoch");
    auto wrong_field=f;wrong_field.seal.source_coefficient_hash++;
    expect_reject([&]{check(s,s,p,wrong_field,source_ids,target_ids,patch_cells);},"stale source field seal");
    wrong_field=f;wrong_field.sample[0].point_id++;
    expect_reject([&]{check(s,s,p,wrong_field,source_ids,target_ids,patch_cells);},"field packet ID mismatch");
    wrong_field=f;wrong_field.sample[0].rho=-1;
    expect_reject([&]{check(s,s,p,wrong_field,source_ids,target_ids,patch_cells);},"negative physical density");
    wrong_field=f;wrong_field.sample[0].internal_energy_density=-1;
    expect_reject([&]{check(s,s,p,wrong_field,source_ids,target_ids,patch_cells);},"negative physical heat");
    expect_reject([&]{stream.finish();},"stream finalized twice");
    auto short_stream=make_stream(budget);short_stream.append(p[0],f.sample[0]);
    expect_reject([&]{short_stream.finish();},"premature stream without complete marginals");
    auto duplicate_stream=make_stream(budget);duplicate_stream.append(p[0],f.sample[0]);
    expect_reject([&]{duplicate_stream.append(p[0],f.sample[0]);},"duplicate ID across stream boundary");
    Budget small_budget=budget;small_budget.maximum_points=3;
    auto bounded_stream=make_stream(small_budget);
    for(size_t i=0;i<3;i++)bounded_stream.append(p[i],f.sample[i]);
    expect_reject([&]{bounded_stream.append(p[3],f.sample[3]);},"global stream point budget");
    std::cout<<"OVERLAP_PACKET_ORIGINAL_IDS_AND_COMMON_MOMENTS_PASS_SYNTHETIC_ONLY\n"
        <<"points="<<area.points<<" area="<<double(area.total_area)
        <<" mass="<<double(m.packet_total[M])<<" heat="<<double(m.packet_total[U])
        <<" kinetic="<<double(m.packet_total[K])<<" max_roundoff_gap="<<double(m.maximum_roundoff_gap)
        <<" continuous_coverage=false quadrature_error_certified=false\n";
}
