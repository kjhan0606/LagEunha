#ifndef MATERIAL_COMPACT_MASS_LAW_H
#define MATERIAL_COMPACT_MASS_LAW_H
// Fixed-M reference masses and immutable parameters for the finite screened
// Laguerre weight law. This is deliberately separate from selected point mass.
#include <climits>
#include <cfloat>
#include <cstdio>
#include <array>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <memory>
#include <map>
#include <numeric>
#include <string>
#include <vector>

namespace compact_mass_law {
using R=long double;
inline void require(bool b,const char*s){compact_p3::require(b,s);}

struct FixedMassPayload {
    std::vector<R>cell_mass;
    uint64_t payload_hash=14695981039346656037ULL;
    uint64_t triangles=0;
    R total_mass=0;
};

inline uint64_t hash_bytes(uint64_t h,const void*data,size_t n){
    const auto*p=static_cast<const unsigned char*>(data);
    for(size_t j=0;j<n;j++){h^=p[j];h*=1099511628211ULL;}return h;
}
inline void compensated_add(R x,R&sum,R&correction){
    R next=sum+x;
    correction+=fabsl(sum)>=fabsl(x)?(sum-next)+x:(x-next)+sum;
    sum=next;
}

// `triangle_cell` is the original per-triangle local cell index. The payload
// is streamed one fixed-M block at a time; it is never duplicated in memory.
// For reference density 1, each triangle's mass is 1^T M_triangle 1 over the
// first ten P3 basis functions. Enriched residual rows are intentionally zero.
inline FixedMassPayload read_fixed_mass(const std::string&path,
        const std::vector<uint64_t>&triangle_cell,size_t cell_count){
    static_assert(sizeof(R)==16&&LDBL_MANT_DIG==64,"original binary80 fixed-M ABI");
    require(cell_count>0&&triangle_cell.size()>0,"nonempty original fixed-M ownership shard");
    require(triangle_cell.size()<=UINT64_MAX/(uint64_t(ST_MAX)*ST_MAX),"fixed-M block count does not overflow");
    std::unique_ptr<FILE,int(*)(FILE*)>file(fopen(path.c_str(),"rb"),fclose);
    require(bool(file),"original immutable fixed-M payload readable");
    uint64_t header[4];require(fread(header,sizeof(uint64_t),4,file.get())==4,"complete original fixed-M header");
    uint64_t blocks=uint64_t(triangle_cell.size())*ST_MAX*ST_MAX;
    require(header[0]==0x4d464d4153533031ULL&&header[1]==1&&header[2]==blocks&&header[3]==sizeof(R),"original fixed-M magic/version/count/binary80 ABI");
    require(blocks<=(UINT64_MAX-sizeof(header))/sizeof(R),"fixed-M payload byte count does not overflow");
    off_t expected=off_t(sizeof(header)+blocks*sizeof(R));
    require(fseeko(file.get(),0,SEEK_END)==0&&ftello(file.get())==expected&&fseeko(file.get(),sizeof(header),SEEK_SET)==0,"exact complete original fixed-M file extent");

    FixedMassPayload out;out.triangles=triangle_cell.size();out.cell_mass.assign(cell_count,0);
    std::vector<R>cell_correction(cell_count,0);std::array<R,size_t(ST_MAX)*ST_MAX>block{};
    out.payload_hash=hash_bytes(out.payload_hash,header,sizeof(header));
    for(size_t s=0;s<triangle_cell.size();s++){
        uint64_t c=triangle_cell[s];require(c<cell_count,"fixed-M triangle maps to its original source cell");
        require(fread(block.data(),sizeof(R),block.size(),file.get())==block.size(),"complete fixed-M triangle block");
        out.payload_hash=hash_bytes(out.payload_hash,block.data(),sizeof(block));
        for(R value:block)require(std::isfinite(value),"finite original fixed-M block entries");
        R triangle_sum=0,triangle_correction=0;
        for(int a=0;a<10;a++)for(int b=0;b<10;b++)compensated_add(block[a*ST_MAX+b],triangle_sum,triangle_correction);
        triangle_sum+=triangle_correction;
        require(std::isfinite(triangle_sum)&&triangle_sum>0,"positive fixed-M mass of each original reference triangle");
        compensated_add(triangle_sum,out.cell_mass[c],cell_correction[c]);
    }
    require(fgetc(file.get())==EOF&&fclose(file.release())==0,"complete unchanged fixed-M payload and clean close");
    for(size_t c=0;c<cell_count;c++){
        out.cell_mass[c]+=cell_correction[c];
        require(std::isfinite(out.cell_mass[c])&&out.cell_mass[c]>0,"positive fixed-M mass for every original material cell");
        out.total_mass+=out.cell_mass[c];
    }
    require(std::isfinite(out.total_mass)&&out.total_mass>0,"finite positive total fixed-M reference mass");
    return out;
}

struct Entry {
    uint64_t cell=0;
    int owner=-1;
    R mref=0,Vsource=0,mu=0,point_mass=0;
    R point_to_fixed_relative=0;
};
inline Entry make_entry(uint64_t cell,int owner,R mref,R Vsource,R point_mass){
    require(owner>=0&&std::isfinite(mref)&&mref>0&&std::isfinite(Vsource)&&Vsource>0&&
        std::isfinite(point_mass)&&point_mass>0,"positive fixed-M mass, canonical source area and separate point mass");
    Entry e;e.cell=cell;e.owner=owner;e.mref=mref;e.Vsource=Vsource;e.mu=mref/Vsource;e.point_mass=point_mass;
    e.point_to_fixed_relative=fabsl(point_mass-mref)/mref;
    require(std::isfinite(e.mu)&&e.mu>0&&std::isfinite(e.point_to_fixed_relative),"finite immutable screened-law coefficient and retained point/M discrepancy");
    return e;
}

class BoundLaw {
public:
    const std::vector<Entry>&entries()const{return entries_;}
    const std::array<uint64_t,4>&fixed_mass_payload_hashes()const{return payload_hashes_;}
    uint64_t source_registry_hash()const{return source_registry_hash_;}
    uint64_t manifest_content_hash()const{return manifest_content_hash_;}
    uint64_t law_hash()const{return law_hash_;}
    uint64_t source_epoch()const{return source_epoch_;}
    uint64_t source_steps()const{return source_steps_;}
    double source_time()const{return source_time_;}
    R coefficient()const{return coefficient_;}
    R total_fixed_mass()const{return total_fixed_mass_;}
    R total_point_mass()const{return total_point_mass_;}
    R total_mass_relative_difference()const{return total_mass_relative_difference_;}
    R fixed_to_checkpoint_mass_relative_difference()const{return fixed_to_checkpoint_mass_relative_difference_;}
    R max_cell_mass_relative_difference()const{return max_cell_mass_relative_difference_;}
    R inherited_mass_bridge()const{return inherited_mass_bridge_;}
private:
    friend BoundLaw bind_actual(const compact_source::Inputs&,const compact_source::Model&,
        const compact_source::GeometryBound&,const compact_source::CurrentContext&,
        const compact_target::Registry&,MPI_Comm);
    std::vector<Entry>entries_;
    std::array<uint64_t,4>payload_hashes_{};
    uint64_t source_registry_hash_=0,manifest_content_hash_=0,law_hash_=0,source_epoch_=0,source_steps_=0;
    double source_time_=0;
    R coefficient_=1,total_fixed_mass_=0,total_point_mass_=0,checkpoint_mass_=0;
    R total_mass_relative_difference_=0,fixed_to_checkpoint_mass_relative_difference_=0;
    R max_cell_mass_relative_difference_=0,inherited_mass_bridge_=0;
};

struct Pair {double hi=0,lo=0;};
inline Pair lossless_pair(R value){
    Pair out;out.hi=double(value);out.lo=double(value-R(out.hi));
    require(std::isfinite(out.hi)&&std::isfinite(out.lo)&&R(out.hi)+R(out.lo)==value,"lossless binary80 to two-double wire pair");return out;
}
inline R from_pair(const double*p){return R(p[0])+R(p[1]);}
inline void add_u64(uint64_t&h,uint64_t value){for(int k=0;k<8;k++){h^=(value>>(8*k))&255;h*=1099511628211ULL;}}
inline void add_real(uint64_t&h,R value){auto p=lossless_pair(value);for(double x:{p.hi,p.lo}){uint64_t bits=compact_target::bits(x);add_u64(h,bits);}}

// Called only after actual source geometry/context and the global target
// registry are bound. The FNV law seal is reproducible metadata; the full
// transaction's prelaunch SHA256 manifest remains the file-integrity gate.
inline BoundLaw bind_actual(const compact_source::Inputs&input,const compact_source::Model&m,
        const compact_source::GeometryBound&geometry,const compact_source::CurrentContext&context,
        const compact_target::Registry&registry,MPI_Comm comm){
    int rank=0,ranks=0;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&ranks);
    uint64_t source_hash=compact_p3::bytes_hash(m.owned.f.coefficient.data(),m.owned.f.coefficient.size()*sizeof(STCoefficient));
    compact_mass::collective_require(ranks==4&&m.layout&&m.layout->rank==rank&&m.layout->ranks==ranks&&
        m.original_reference_receipt_bound&&m.owned.global_triangles==1138432&&m.owned.global_dofs==10816441&&
        registry.site.size()==71080&&registry.source_epoch==geometry.source_epoch&&
        registry.source_steps==geometry.source_steps&&registry.source_time==geometry.source_time&&
        context.coefficient_hash==geometry.source_coefficient_hash&&source_hash==geometry.source_coefficient_hash&&
        registry.source_shard_hash[rank]==source_hash&&
        geometry.source_epoch==m.accepted.physical.epoch&&geometry.source_steps==m.accepted.physical.steps&&
        geometry.source_time==m.accepted.physical.time&&context.baseline.size()==m.owned.f.cell.size(),
        comm,"fixed-M law requires one full actual accepted source/geometry/context/registry epoch");

    FixedMassPayload local;bool read_ok=true;
    try{
        std::vector<uint64_t>triangle_cell(m.owned.f.triangle.size());
        for(size_t s=0;s<triangle_cell.size();s++){
            const auto&t=m.owned.f.triangle[s];
            require(t.cell>=0&&size_t(t.cell)<m.owned.f.cell.size()&&m.owned.source_cell[s]==m.owned.f.cell[t.cell],"fixed-M triangle preserves original source-cell identity");
            triangle_cell[s]=uint64_t(t.cell);
        }
        local=read_fixed_mass(m.fixed_mass_path,triangle_cell,m.owned.f.cell.size());
    }catch(const std::exception&){read_ok=false;}
    compact_mass::collective_require(read_ok,comm,"all owners read complete finite original fixed-M shards");

    std::map<uint64_t,size_t>registry_index;
    bool local_ok=true;
    for(size_t j=0;j<registry.site.size();j++)local_ok&=registry_index.emplace(registry.site[j].cell,j).second;
    std::vector<uint64_t>ids(m.owned.f.cell.size());std::vector<double>values(ids.size()*6);
    try{
        for(size_t c=0;c<ids.size();c++){
            uint64_t cell=m.owned.f.cell[c];ids[c]=cell;
            const auto&base=context.baseline[c];auto found=registry_index.find(cell);
            require(base.cell==cell&&base.volume>0&&base.mass>0&&found!=registry_index.end(),"saved canonical source area and endpoint point mass for the same fixed-M cell");
            const auto&site=registry.site[found->second];
            require(site.owner==rank&&site.source_mass==base.mass,"point-functional mass remains separately identified from fixed-M reference mass");
            Pair mass=lossless_pair(local.cell_mass[c]),volume=lossless_pair(base.volume),point=lossless_pair(base.mass);
            double*dst=values.data()+6*c;dst[0]=mass.hi;dst[1]=mass.lo;dst[2]=volume.hi;dst[3]=volume.lo;dst[4]=point.hi;dst[5]=point.lo;
        }
    }catch(const std::exception&){local_ok=false;}
    compact_mass::collective_require(local_ok,comm,"fixed-M cell masses bind to their exact owner, saved Vsource and separate point mass");

    uint64_t local_count=ids.size();std::vector<uint64_t>counts(ranks);MPI_Allgather(&local_count,1,MPI_UINT64_T,counts.data(),1,MPI_UINT64_T,comm);
    std::vector<int>id_count(ranks),id_disp(ranks),value_count(ranks),value_disp(ranks);uint64_t total=0;
    bool bounded=true;
    for(int r=0;r<ranks;r++){
        bounded&=counts[r]<=uint64_t(INT_MAX/6)&&total<=uint64_t(INT_MAX/6);
        id_count[r]=int(counts[r]);id_disp[r]=int(total);value_count[r]=int(6*counts[r]);value_disp[r]=int(6*total);total+=counts[r];
    }
    compact_mass::collective_require(bounded&&total==registry.site.size(),comm,"bounded complete globally owned fixed-M cell records");
    std::vector<uint64_t>all_ids(total);std::vector<double>all_values(6*total);
    MPI_Allgatherv(ids.data(),int(ids.size()),MPI_UINT64_T,all_ids.data(),id_count.data(),id_disp.data(),MPI_UINT64_T,comm);
    MPI_Allgatherv(values.data(),int(values.size()),MPI_DOUBLE,all_values.data(),value_count.data(),value_disp.data(),MPI_DOUBLE,comm);
    uint64_t payload_hash=local.payload_hash;BoundLaw law;MPI_Allgather(&payload_hash,1,MPI_UINT64_T,law.payload_hashes_.data(),1,MPI_UINT64_T,comm);

    std::map<uint64_t,Entry>global;
    for(int owner=0;owner<ranks;owner++)for(int j=0;j<id_count[owner];j++){
        size_t q=size_t(id_disp[owner]+j);const double*p=all_values.data()+6*q;Entry e;e.cell=all_ids[q];e.owner=owner;
        try{e=make_entry(e.cell,owner,from_pair(p),from_pair(p+2),from_pair(p+4));}
        catch(const std::exception&){local_ok=false;}
        local_ok&=global.emplace(e.cell,e).second;
    }
    local_ok&=global.size()==registry.site.size();
    for(size_t j=0;j<registry.site.size();j++){
        const auto&site=registry.site[j];auto it=global.find(site.cell);
        local_ok&=it!=global.end()&&it->second.owner==site.owner&&(j==0||registry.site[j-1].cell<site.cell);
    }
    compact_mass::collective_require(local_ok,comm,"one positive fixed-M mass and saved physical area per canonical source site");

    law.entries_.reserve(registry.site.size());
    for(const auto&site:registry.site){const auto&e=global.at(site.cell);law.entries_.push_back(e);law.total_fixed_mass_+=e.mref;law.total_point_mass_+=e.point_mass;law.max_cell_mass_relative_difference_=std::max(law.max_cell_mass_relative_difference_,e.point_to_fixed_relative);}
    require(law.total_fixed_mass_>0&&law.total_point_mass_>0&&std::isfinite(law.total_fixed_mass_)&&std::isfinite(law.total_point_mass_),"positive global fixed-M and separately retained point-functional mass totals");
    law.total_mass_relative_difference_=fabsl(law.total_point_mass_-law.total_fixed_mass_)/law.total_fixed_mass_;
    law.checkpoint_mass_=m.accepted.current.mass;
    require(std::isfinite(law.checkpoint_mass_)&&law.checkpoint_mass_>0,"positive inherited accepted checkpoint mass");
    law.fixed_to_checkpoint_mass_relative_difference_=fabsl(law.total_fixed_mass_-law.checkpoint_mass_)/law.checkpoint_mass_;
    law.inherited_mass_bridge_=m.accepted.numerical.mass_bridge;
    law.source_registry_hash_=registry.wire_hash;law.source_epoch_=registry.source_epoch;law.source_steps_=registry.source_steps;law.source_time_=registry.source_time;
    law.coefficient_=R(1); // mu = 1 * mref / Vsource; never reset during W Newton.
    std::string manifest;bool manifest_ok=true;
    try{manifest=compact_source::small_text(input.immutable_manifest);}catch(const std::exception&){manifest_ok=false;}
    compact_mass::collective_require(manifest_ok,comm,"same immutable transaction manifest is readable on every owner");
    law.manifest_content_hash_=compact_p3::bytes_hash(manifest.data(),manifest.size());
    uint64_t manifest_hash_min=law.manifest_content_hash_,manifest_hash_max=law.manifest_content_hash_;
    MPI_Allreduce(MPI_IN_PLACE,&manifest_hash_min,1,MPI_UINT64_T,MPI_MIN,comm);
    MPI_Allreduce(MPI_IN_PLACE,&manifest_hash_max,1,MPI_UINT64_T,MPI_MAX,comm);
    compact_mass::collective_require(manifest_hash_min==manifest_hash_max,comm,"one exact immutable transaction manifest on all owners");
    uint64_t h=14695981039346656037ULL;add_u64(h,law.source_registry_hash_);add_u64(h,law.source_epoch_);add_u64(h,law.source_steps_);
    add_u64(h,compact_target::bits(law.source_time_));add_u64(h,law.manifest_content_hash_);add_real(h,law.coefficient_);
    for(int r=0;r<ranks;r++){add_u64(h,registry.source_shard_hash[r]);add_u64(h,law.payload_hashes_[r]);}
    for(const auto&e:law.entries_){add_u64(h,e.cell);add_u64(h,uint64_t(e.owner));add_real(h,e.mref);add_real(h,e.Vsource);add_real(h,e.mu);add_real(h,e.point_mass);}
    law.law_hash_=h;
    return law;
}
BoundLaw bind_accepted_reference_law(const compact_source::Inputs&,const compact_source::Model&,
    const compact_source::GeometryBound&,const compact_source::CurrentContext&,
    const compact_target::Registry&,MPI_Comm);
} // namespace compact_mass_law
#endif
