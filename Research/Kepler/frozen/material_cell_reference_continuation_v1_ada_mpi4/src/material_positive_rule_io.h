#ifndef MATERIAL_POSITIVE_RULE_IO_H
#define MATERIAL_POSITIVE_RULE_IO_H
#include "material_full_reduced_rule.h"
// Included after native hash/json helpers, inside full_evolution.
static void load_positive_rule(Material&m,const md::Ownership&layout,const std::string&directory,bool require_anchor=true){
    auto&r=m.positive;need(!r.loaded,"one immutable positive-rule generation per runtime");
    std::string report=textfile(directory+"/contract.json");
    need(report.find("\"candidate_ready_for_actual_RK_oracle\":true")!=std::string::npos&&report.find("\"all_original_triangles_processed\":true")!=std::string::npos,"complete full-domain candidate, not a live/partial prefix");
    need(uint64_t(json_number(report,"completed_triangles"))==m.owned.global_triangles&&uint64_t(json_number(report,"original_positive_points"))==2040256512ULL,"same complete original reference problem");
    uint64_t expected_points=uint64_t(json_number(report,"retained_points_with_original_rule_fallbacks"));need(expected_points>0&&expected_points*4<=2040256512ULL,"full candidate point-reduction gate");
    std::ifstream exit_file(directory+"/native_exit_code.txt");int exit_code=-1;exit_file>>exit_code;need(exit_code==0,"actual completed native preprocessing");
    FILE*f=fopen((directory+"/run/positive_rule_r"+std::to_string(layout.rank)+".bin").c_str(),"rb");need(f,"complete native positive-rule shard");
    uint64_t h[12];material_full::read(f,h,12);need(h[0]==0x4d46504f53525531ULL&&h[1]==1&&h[2]==uint64_t(layout.rank)&&h[3]==uint64_t(layout.ranks)&&h[4]==m.owned.f.triangle.size()&&h[5]==m.owned.global_triangles&&h[6]==sizeof(R)&&h[8]==hash(m.immutable.data(),m.immutable.size())&&h[9]==MF_EVOL_METRICS&&h[10]==16,"native rule reference/owner/moment ABI");
    material_full::read(f,&r.anchor_time);r.anchor_state=h[7];if(require_anchor)need(r.anchor_time==m.ledger.time&&r.anchor_state==hash(m.owned.f.coefficient.data(),m.owned.f.coefficient.size()*sizeof(STCoefficient)),"candidate belongs to this exact retained accepted state");
    size_t nt=h[4];r.offset.resize(nt+1);r.row_offset.resize(nt+1);r.count.resize(nt);r.compressed.resize(nt);r.diagnostic.resize(nt);
    uint64_t local_original=0,local_evaluated=0,stream=hash(h,sizeof(h));stream=hash(&r.anchor_time,sizeof(r.anchor_time),stream);
    for(size_t s=0;s<nt;s++){
        uint64_t record[6];material_full::read(f,record,6);size_t rows=record[4],nodes=record[2];int dofs=m.owned.f.triangle[s].n;
        need(record[0]==m.owned.source_triangle[s]&&record[1]>0&&record[1]<=m.cache_capacity&&nodes<=record[1]&&rows==size_t(dofs*dofs+4*dofs+4)&&record[5]==0&&record[3]<=1,"ordered complete rule records with every original constraint");
        bool compressed=record[3];need((compressed&&nodes>0&&4*nodes<=record[1])||(!compressed&&nodes==0),"strict positive subset or explicit original-rule fallback");
        r.offset[s+1]=r.offset[s]+nodes;r.row_offset[s+1]=r.row_offset[s]+rows;r.count[s]=compressed?nodes:record[1];r.compressed[s]=compressed;
        material_full::read(f,r.diagnostic[s].data(),16);size_t at=r.row_offset[s];r.residual.resize(at+rows);r.rounding.resize(at+rows);material_full::read(f,r.residual.data()+at,rows);material_full::read(f,r.rounding.data()+at,rows);
        size_t point_at=r.index.size();r.index.resize(point_at+nodes);r.weight.resize(point_at+nodes);material_full::read(f,r.index.data()+point_at,nodes);material_full::read(f,r.weight.data()+point_at,nodes);uint64_t got;material_full::read(f,&got);
        uint64_t checksum=hash(record,sizeof(record));checksum=hash(r.diagnostic[s].data(),16*sizeof(R),checksum);checksum=hash(r.residual.data()+at,rows*sizeof(R),checksum);checksum=hash(r.rounding.data()+at,rows*sizeof(R),checksum);checksum=hash(r.index.data()+point_at,nodes*sizeof(uint32_t),checksum);checksum=hash(r.weight.data()+point_at,nodes*sizeof(double),checksum);need(checksum==got,"each native rule/residual/roundoff payload checksum");
        for(size_t j=0;j<nodes;j++)need(r.index[point_at+j]<record[1]&&(!j||r.index[point_at+j-1]<r.index[point_at+j])&&r.weight[point_at+j]>0&&std::isfinite(r.weight[point_at+j]),"unique original positive nodes and strictly positive finite stored weights");
        for(size_t j=0;j<rows;j++)need(std::isfinite(r.residual[at+j])&&std::isfinite(r.rounding[at+j])&&r.rounding[at+j]>=0,"finite signed all-row discrepancies and rounding bounds");
        if(compressed){for(R value:r.diagnostic[s])need(std::isfinite(value),"finite native adopted candidate contractions");const auto&d=r.diagnostic[s];need(d[0]>=0&&d[8]>=0&&d[9]>=0&&d[11]>=0&&d[14]>0&&fabsl(d[1])<1e-8L&&fabsl(d[2])<1e-8L&&fabsl(d[3])<1e-11L&&fabsl(d[4])<1e-8L&&fabsl(d[5])<1e-8L&&fabsl(d[15])<1e-8L&&d[8]<1e-10L&&fabsl(d[7])<1e-8L&&d[9]<1e-11L&&fabsl(d[10])<1e-11L&&d[11]<1e-11L,"every compressed record retains the original native cold-heat/metric/volume/force gates");}
        if(!compressed)for(size_t j=0;j<rows;j++)need(r.residual[at+j]==0&&r.rounding[at+j]==0,"original-rule fallback has no invented discrepancy");
        stream=hash(&checksum,sizeof(checksum),stream);local_original+=record[1];local_evaluated+=r.count[s];
    }
    need(fgetc(f)==EOF&&fclose(f)==0,"no incomplete/trailing candidate records");uint64_t counts[2]={local_original,local_evaluated};MPI_Allreduce(MPI_IN_PLACE,counts,2,MPI_UINT64_T,MPI_SUM,layout.comm);need(counts[0]==2040256512ULL&&counts[1]==expected_points,"all native source/fallback point counts agree with producer");
    std::vector<uint64_t>generation(layout.ranks);MPI_Allgather(&stream,1,MPI_UINT64_T,generation.data(),1,MPI_UINT64_T,layout.comm);r.generation=hash(generation.data(),generation.size()*sizeof(uint64_t));r.global_points=counts[1];r.directory=directory;r.binding="POSITIVE_ORIGINAL_POINT_SUBSET_V1\n"+directory+"\n"+std::to_string(r.generation);r.loaded=true;m.numerical.generation=r.generation;
}
static void upload_positive_rule(Material&m,const md::Ownership&layout,Device&gpu){auto&r=m.positive;need(r.loaded,"loaded immutable rule before device upload");int bad=mf_evolved_rule_upload(gpu.ptr,r.offset.data(),m.owned.f.triangle.size(),r.index.data(),r.weight.data(),r.index.size())!=0,any;MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,layout.comm);need(!any,"native positive subset/weights upload");}
static void use_positive_rule(Material&m,const md::Ownership&layout,Device&gpu,bool active){need(!active||m.positive.loaded,"cannot activate an unverified/incomplete rule");int bad=mf_evolved_rule_use(gpu.ptr,active?1:0)!=0,any;MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,layout.comm);need(!any,"native full/reduced same-state rule switch");m.positive.active=active;}
#endif
