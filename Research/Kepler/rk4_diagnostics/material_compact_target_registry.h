#ifndef MATERIAL_COMPACT_TARGET_REGISTRY_H
#define MATERIAL_COMPACT_TARGET_REGISTRY_H
// Include after compact_source_target. Only small geometry is replicated:
// transported generators and frozen exterior functions, never field DOFs/M.
#include <climits>
#include <numeric>
namespace compact_target {
using compact_source::Site;
using compact_source::Exterior;
using compact_source::TargetSource;
using R=long double;
inline void require(bool b,const char*s){compact_p3::require(b,s);}
inline uint64_t bits(double x){uint64_t out;std::memcpy(&out,&x,sizeof out);return out;}
struct Registry {
    std::vector<Site>site;             // sorted by ORIGINAL cell identity
    std::vector<Exterior>exterior;     // sorted by ORIGINAL exterior edge ID
    std::vector<size_t>boundary_cycle; // positive source orientation
    std::array<uint64_t,4>source_shard_hash={};
    uint64_t source_steps=0,source_epoch=0,wire_hash=0;
    double source_time=0;
    R maximum_position_error=0,maximum_velocity_error=0;
    // Connectivity of the SOURCE boundary is checked below. This is NOT a
    // power partition, intersection-root proof, target coverage or transfer.
};
constexpr int site_tags=5,site_values=20,edge_tags=12,edge_values=53;
struct Wire {std::vector<uint64_t>tag;std::vector<double>value;};
inline Wire pack_sites(const std::vector<Site>&sites,int rank){
    Wire out;out.tag.reserve(site_tags*sites.size());out.value.reserve(site_values*sites.size());
    for(const auto&s:sites){
        require(s.owner==rank&&s.chart>=0&&s.chart<=2,"local transported site owner/chart");
        require(s.cell!=UINT64_MAX&&s.source_triangle!=UINT64_MAX,"source site identities");
        out.tag.insert(out.tag.end(),{s.cell,s.source_triangle,uint64_t(s.owner),uint64_t(s.chart),uint64_t(s.canonical_bound)});
        for(const double*a:{s.reference,s.position,s.velocity,s.position_error,s.velocity_error})for(int k=0;k<2;k++){
            require(std::isfinite(a[k]),"finite transported site fields");out.value.push_back(a[k]);
        }
        for(int k=0;k<2;k++)require(s.position_error[k]>=0&&s.velocity_error[k]>=0,"nonnegative source transport enclosures");
        require(s.source_mass>0&&std::isfinite(s.source_mass),"actual positive source mass, not a law parameter");
        double hi=double(s.source_mass),lo=double(s.source_mass-R(hi));
        require(std::isfinite(hi)&&std::isfinite(lo)&&R(hi)+R(lo)==s.source_mass,"lossless binary80 source mass transport");
        out.value.push_back(hi);out.value.push_back(lo);
        for(const double*a:{s.position_low,s.velocity_low,s.canonical_position_error,s.canonical_velocity_error})for(int k=0;k<2;k++){
            require(std::isfinite(a[k]),"finite refined canonical site wire");out.value.push_back(a[k]);
        }
        for(int k=0;k<2;k++){
            require(s.canonical_position_error[k]>=0&&s.canonical_velocity_error[k]>=0,"positive canonical tail enclosures");
            if(s.canonical_bound)require(R(s.position_error[k])>=fabsl(R(s.position_low[k]))+R(s.canonical_position_error[k])&&R(s.velocity_error[k])>=fabsl(R(s.velocity_low[k]))+R(s.canonical_velocity_error[k]),"refined canonical source ball contained in retained native cache ball");
            else require(s.position_low[k]==0&&s.velocity_low[k]==0&&s.canonical_position_error[k]==0&&s.canonical_velocity_error[k]==0,"unbound site must retain its original uncertainty, no implicit precision claim");
        }
    }return out;
}
inline Wire pack_edges(const std::vector<Exterior>&edges,int rank,uint64_t hash,uint64_t epoch){
    Wire out;out.tag.reserve(edge_tags*edges.size());out.value.reserve(edge_values*edges.size());
    for(const auto&e:edges){
        const auto&f=e.frozen;const auto&q=f.function;
        require(e.owner==rank&&(e.orientation==1||e.orientation==-1)&&e.local_edge>=0&&e.local_edge<3,"owned oriented actual exterior edge");
        require(e.edge==f.edge&&f.source_hash==hash&&f.source_epoch==epoch,"same immutable source exterior function");
        require(e.vertex[0]!=e.vertex[1]&&e.vertex[0]!=UINT64_MAX&&e.vertex[1]!=UINT64_MAX&&e.edge!=UINT64_MAX&&e.cell!=UINT64_MAX&&e.source_triangle!=UINT64_MAX,"exterior identities");
        require(q.rank>=0&&q.rank<=2,"complete frozen source edge rank");
        for(int j=0;j<4;j++)require(q.poly[j]==j,"canonical frozen P3 indices");
        for(int j=0;j<q.rank;j++)require(q.id[j]==4+j,"canonical frozen scalar indices");
        out.tag.insert(out.tag.end(),{e.edge,e.vertex[0],e.vertex[1],e.cell,e.source_triangle,uint64_t(e.owner),uint64_t(e.orientation==1),uint64_t(e.local_edge),f.edge,f.source_hash,f.source_epoch,uint64_t(q.rank)});
        auto append=[&](double x){require(std::isfinite(x),"finite frozen exterior wire values");out.value.push_back(x);};
        for(double x:q.A)append(x);for(double x:q.D)append(x);
        for(const auto&row:q.qnode)for(double x:row)append(x);append(q.eps);
        for(const auto&row:q.encode)for(double x:row)append(x);
        for(const auto&row:q.decode)for(double x:row)append(x);
        for(const auto&c:f.coefficient){for(double x:c.x)append(x);for(double v:c.v)append(v);}
        for(const auto&x:e.reference)for(double v:x)append(v);
        for(double x:e.lower)append(x);for(double x:e.upper)append(x);
        require(q.eps>0,"original frozen edge softening");
        for(int k=0;k<2;k++){
            require(bits(q.A[k])==bits(e.reference[0][k])&&bits(q.D[k])==bits(e.reference[1][k]-e.reference[0][k]),"frozen function and original oriented reference endpoints");
            require(e.lower[k]<=e.upper[k],"source curve bounding box ordering");
            for(int j:{0,3})require(f.coefficient[j].x[k]>=e.lower[k]&&f.coefficient[j].x[k]<=e.upper[k],"actual exterior vertices inside source box");
        }
    }return out;
}
inline Site unpack_site(const uint64_t*t,const double*v){
    Site s;s.cell=t[0];s.source_triangle=t[1];s.owner=int(t[2]);s.chart=int(t[3]);s.canonical_bound=t[4]!=0;size_t at=0;
    for(double*a:{s.reference,s.position,s.velocity,s.position_error,s.velocity_error})for(int k=0;k<2;k++)a[k]=v[at++];
    s.source_mass=R(v[10])+R(v[11]);at=12;
    for(double*a:{s.position_low,s.velocity_low,s.canonical_position_error,s.canonical_velocity_error})for(int k=0;k<2;k++)a[k]=v[at++];return s;
}
inline Exterior unpack_edge(const uint64_t*t,const double*v){
    Exterior e;e.edge=t[0];e.vertex[0]=t[1];e.vertex[1]=t[2];e.cell=t[3];e.source_triangle=t[4];e.owner=int(t[5]);e.orientation=t[6]?1:-1;e.local_edge=int(t[7]);
    auto&f=e.frozen;f.edge=t[8];f.source_hash=t[9];f.source_epoch=t[10];auto&q=f.function;q.rank=int(t[11]);size_t at=0;
    for(int j=0;j<4;j++)q.poly[j]=j;for(int j=0;j<q.rank;j++)q.id[j]=4+j;
    for(double&x:q.A)x=v[at++];for(double&x:q.D)x=v[at++];
    for(auto&row:q.qnode)for(double&x:row)x=v[at++];q.eps=v[at++];
    for(auto&row:q.encode)for(double&x:row)x=v[at++];
    for(auto&row:q.decode)for(double&x:row)x=v[at++];
    for(auto&c:f.coefficient){for(double&x:c.x)x=v[at++];for(double&x:c.v)x=v[at++];}
    for(auto&row:e.reference)for(double&x:row)x=v[at++];
    for(double&x:e.lower)x=v[at++];for(double&x:e.upper)x=v[at++];
    require(at==edge_values,"complete fieldwise frozen-curve wire, no struct padding");return e;
}
inline Wire gather_wire(const Wire&local,const std::vector<uint64_t>&count,int nt,int nv,MPI_Comm comm){
    std::vector<int>tc(count.size()),td(count.size()),vc(count.size()),vd(count.size());uint64_t n=0;
    for(size_t j=0;j<count.size();j++){
        require(count[j]<=uint64_t(INT_MAX/std::max(nt,nv))&&n<=uint64_t(INT_MAX/std::max(nt,nv))-count[j],"bounded MPI geometry registry counts");
        tc[j]=int(count[j]*nt);vc[j]=int(count[j]*nv);td[j]=int(n*nt);vd[j]=int(n*nv);n+=count[j];
    }
    Wire out;out.tag.resize(n*nt);out.value.resize(n*nv);
    MPI_Allgatherv(local.tag.data(),int(local.tag.size()),MPI_UINT64_T,out.tag.data(),tc.data(),td.data(),MPI_UINT64_T,comm);
    MPI_Allgatherv(local.value.data(),int(local.value.size()),MPI_DOUBLE,out.value.data(),vc.data(),vd.data(),MPI_DOUBLE,comm);
    return out;
}
inline void verify_boundary(Registry&out){
    require(out.exterior.size()>=3,"nonempty closed actual source boundary");
    struct Vertex {double reference[2],position[2],velocity[2];size_t incoming=SIZE_MAX,outgoing=SIZE_MAX;};
    std::map<uint64_t,Vertex>vertex;
    for(size_t j=0;j<out.exterior.size();j++){
        const auto&e=out.exterior[j];
        auto site=std::lower_bound(out.site.begin(),out.site.end(),e.cell,[](const Site&s,uint64_t c){return s.cell<c;});
        require(site!=out.site.end()&&site->cell==e.cell&&site->owner==e.owner,"exterior belongs to a registered original material cell/owner");
        for(int a=0;a<2;a++){
            const auto&coefficient=e.frozen.coefficient[a?3:0];
            auto p=vertex.emplace(e.vertex[a],Vertex{});auto&v=p.first->second;
            for(int k=0;k<2;k++){
                if(p.second){v.reference[k]=e.reference[a][k];v.position[k]=coefficient.x[k];v.velocity[k]=coefficient.v[k];}
                else require(bits(v.reference[k])==bits(e.reference[a][k])&&bits(v.position[k])==bits(coefficient.x[k])&&bits(v.velocity[k])==bits(coefficient.v[k]),"exact shared reference and accepted x/v vertex bits across source owners");
            }
            bool start=(a==0)==(e.orientation==1);size_t&slot=start?v.outgoing:v.incoming;
            require(slot==SIZE_MAX,"one incoming/outgoing exterior curve per original boundary vertex");slot=j;
        }
    }
    require(vertex.size()==out.exterior.size(),"closed source boundary edge/vertex count");
    for(const auto&entry:vertex)require(entry.second.incoming!=SIZE_MAX&&entry.second.outgoing!=SIZE_MAX,"no open source boundary vertex");
    size_t current=0;std::vector<bool>seen(out.exterior.size());
    for(size_t j=0;j<out.exterior.size();j++){
        require(!seen[current],"one source boundary component, no disconnected closed extra cycle");seen[current]=true;out.boundary_cycle.push_back(current);
        const auto&e=out.exterior[current];current=vertex.at(e.vertex[e.orientation==1?1:0]).outgoing;
    }
    require(current==0,"exactly one oriented exterior cycle for original square source domain");
    // This graph check deliberately does not approximate Gamma, resolve new
    // power-plane roots, prove target chart Jacobians or certify Ω coverage.
}
inline Registry gather(const TargetSource&local,MPI_Comm comm){
    int rank=0,ranks=0;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&ranks);
    // Candidate-specific bounded small geometry. Small synthetic registries
    // may test this adapter; the guarded importer separately requires71080.
    bool valid=ranks>=1&&ranks<=4&&local.global_sites>0&&local.global_sites<=71080&&local.global_exterior>=3&&local.global_exterior<=2048;
    valid&=std::isfinite(local.source_time)&&local.source_time>=0&&local.site.size()<=local.global_sites&&local.exterior.size()<=local.global_exterior;
    compact_mass::collective_require(valid,comm,"bounded source registry metadata on every rank");
    uint64_t header[9]={local.source_steps,local.source_epoch,bits(local.source_time),local.global_sites,local.global_exterior,local.source_shard_hash[0],local.source_shard_hash[1],local.source_shard_hash[2],local.source_shard_hash[3]};
    std::vector<uint64_t>headers(9*ranks);MPI_Allgather(header,9,MPI_UINT64_T,headers.data(),9,MPI_UINT64_T,comm);
    valid=true;for(int j=0;j<ranks;j++)for(int k=0;k<9;k++)valid&=header[k]==headers[9*j+k];
    compact_mass::collective_require(valid,comm,"same exact source epoch/time/hash set and declared geometry counts");
    Wire sw,ew;valid=true;try{sw=pack_sites(local.site,rank);ew=pack_edges(local.exterior,rank,local.source_shard_hash[rank],local.source_epoch);}catch(const std::exception&){valid=false;}
    compact_mass::collective_require(valid,comm,"valid fieldwise source geometry packets on every owner");
    uint64_t counts[2]={local.site.size(),local.exterior.size()};std::vector<uint64_t>all_counts(2*ranks),sc(ranks),ec(ranks);
    MPI_Allgather(counts,2,MPI_UINT64_T,all_counts.data(),2,MPI_UINT64_T,comm);
    for(int j=0;j<ranks;j++){sc[j]=all_counts[2*j];ec[j]=all_counts[2*j+1];}
    require(std::accumulate(sc.begin(),sc.end(),uint64_t(0))==local.global_sites&&std::accumulate(ec.begin(),ec.end(),uint64_t(0))==local.global_exterior,"complete source site/exterior global counts");
    auto sites=gather_wire(sw,sc,site_tags,site_values,comm);auto edges=gather_wire(ew,ec,edge_tags,edge_values,comm);
    Registry out;out.source_shard_hash=local.source_shard_hash;out.source_steps=local.source_steps;out.source_epoch=local.source_epoch;out.source_time=local.source_time;
    size_t sat=0,eat=0;for(int owner=0;owner<ranks;owner++){
        for(uint64_t j=0;j<sc[owner];j++,sat++){
            auto s=unpack_site(sites.tag.data()+sat*site_tags,sites.value.data()+sat*site_values);require(s.owner==owner,"site sender is actual owner");out.site.push_back(s);
            for(int k=0;k<2;k++){out.maximum_position_error=std::max(out.maximum_position_error,R(s.position_error[k]));out.maximum_velocity_error=std::max(out.maximum_velocity_error,R(s.velocity_error[k]));}
        }
        for(uint64_t j=0;j<ec[owner];j++,eat++){
            auto e=unpack_edge(edges.tag.data()+eat*edge_tags,edges.value.data()+eat*edge_values);require(e.owner==owner,"exterior sender is actual owner");out.exterior.push_back(e);
        }
    }
    std::sort(out.site.begin(),out.site.end(),[](const Site&a,const Site&b){return a.cell<b.cell;});
    std::sort(out.exterior.begin(),out.exterior.end(),[](const Exterior&a,const Exterior&b){return a.edge<b.edge;});
    for(size_t j=1;j<out.site.size();j++)require(out.site[j-1].cell!=out.site[j].cell,"unique global transported generator per original cell");
    for(size_t j=1;j<out.exterior.size();j++)require(out.exterior[j-1].edge!=out.exterior[j].edge,"unique actual source exterior edge");
    verify_boundary(out);
    // Hash fieldwise wire data; never hash native long-double/struct padding.
    uint64_t hash=1469598103934665603ULL;auto add=[&](uint64_t v){for(int k=0;k<8;k++){hash^=(v>>(8*k))&255;hash*=1099511628211ULL;}};
    for(uint64_t x:header)add(x);for(const Wire*w:{&sites,&edges}){for(uint64_t t:w->tag)add(t);for(double v:w->value)add(bits(v));}out.wire_hash=hash;
    return out;
}
} // namespace compact_target
#endif
