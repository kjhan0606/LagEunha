#ifndef MATERIAL_COMPACT_SOURCE_TARGET_H
#define MATERIAL_COMPACT_SOURCE_TARGET_H
// Included after the actual source, whole-map geometry and current-context
// adapters. Inputs for the ONE actual compact transaction, not a partial run.
#include "material_compact_geometry.h"
#include "material_compact_predicate.h"
#include <filesystem>
namespace compact_source {
struct Site {
    uint64_t cell=0,source_triangle=0;int owner=0,chart=0;
    double reference[2]={},position[2]={},velocity[2]={};
    double position_error[2]={},velocity_error[2]={};
    // Native high values above are retained. Geometry can use this tighter
    // SAME canonical field enclosure instead of a binary64 cache-error floor.
    double position_low[2]={},velocity_low[2]={};
    double canonical_position_error[2]={},canonical_velocity_error[2]={};
    bool canonical_bound=false;
    // Actual accepted positive-functional source mass, NOT legacy generator
    // mass and NOT an unilaterally chosen screened-law reference mass or mu.
    R source_mass=0;
};
inline sf::S canonical_position(const Site&s,int k){
    if(!s.canonical_bound)return sf::S(sf::Q(s.position[k]),sf::Q(s.position_error[k]));
    return sf::exact(s.position[k])+sf::exact(s.position_low[k])+sf::S(0,sf::Q(s.canonical_position_error[k]));
}
inline sf::S canonical_velocity(const Site&s,int k){
    if(!s.canonical_bound)return sf::S(sf::Q(s.velocity[k]),sf::Q(s.velocity_error[k]));
    return sf::exact(s.velocity[k])+sf::exact(s.velocity_low[k])+sf::S(0,sf::Q(s.canonical_velocity_error[k]));
}
struct Exterior {
    uint64_t edge=0,vertex[2]={},cell=0,source_triangle=0;
    int owner=0,orientation=0,local_edge=0;
    double reference[2][2]={},lower[2]={},upper[2]={};
    compact_geometry::FrozenTrace frozen;
};
struct TargetSource {
    std::vector<Site>site;
    std::vector<Exterior>exterior;
    std::array<uint64_t,4>source_shard_hash={};
    uint64_t source_steps=0,source_epoch=0,global_sites=0,global_exterior=0;
    double source_time=0;
    R maximum_position_error=0,maximum_velocity_error=0;
    // No target partition, weight solve, coverage proof, projection or commit.
};
inline std::vector<std::array<double,2>> raw_generators(const Inputs&input,const Model&m){
    auto path=std::filesystem::path(input.reference_receipt).parent_path()/("material_native_geometry_r"+std::to_string(m.layout->rank)+".bin");
    // This exact raw geometry shard must be in the full transaction manifest.
    std::unique_ptr<FILE,int(*)(FILE*)>file(fopen(path.c_str(),"rb"),fclose);require(bool(file),"original native generator coordinates");
    uint64_t h[8];material_full::read(file.get(),h,8);
    require(h[0]==0x4d46524157474531ULL&&h[1]==1&&h[2]==m.owned.f.cell.size()&&h[3]==sizeof(MFCell)&&h[4]==sizeof(MFFace)&&h[5]==uint64_t(m.layout->rank)&&h[6]==0&&h[7]==0,"same original retained generator owner/ABI");
    std::vector<std::array<double,2>>out(m.owned.f.cell.size());std::vector<bool>seen(out.size());
    for(uint64_t j=0;j<h[2];j++){
        MFCell cell;material_full::read(file.get(),&cell);auto it=std::lower_bound(m.owned.f.cell.begin(),m.owned.f.cell.end(),cell.id);
        require(it!=m.owned.f.cell.end()&&*it==cell.id&&cell.n>=3&&cell.n<=128,"same original physical cell generator identity");
        size_t at=it-m.owned.f.cell.begin();require(!seen[at],"one original generator per owned material cell");seen[at]=true;
        double x[2]={cell.x,cell.y};
        for(int k=0;k<2;k++){
            require(std::isfinite(x[k])&&x[k]>0&&x[k]<24,"same original unwrapped box24 native generator");
            out[at][k]=x[k]-12;require(sf::Q(out[at][k])==sf::Q(x[k])-12,"exact original coordinate centering, no new orbit");
        }
        require(fseeko(file.get(),off_t(cell.n*sizeof(MFFace)),SEEK_CUR)==0,"skip already manifested raw face data");
    }
    require(fgetc(file.get())==EOF&&fclose(file.release())==0,"complete original native generator shard");return out;
}
struct Location {size_t triangle=SIZE_MAX;int chart=0,edge=-1,vertex=-1;};
inline Location locate_generator(const Model&m,const std::vector<size_t>&triangles,const double*X){
    Location out;uint64_t selected=UINT64_MAX;
    for(size_t s:triangles){
        const auto&t=m.owned.f.triangle[s];int sign[3],zeros=0;bool inside=true;
        for(int k=0;k<3;k++){sign[k]=compact_predicate::orientation(t.vertex[k],t.vertex[(k+1)%3],X);inside&=sign[k]>=0;zeros+=sign[k]==0;}
        if(!inside||m.owned.source_triangle[s]>=selected)continue;
        selected=m.owned.source_triangle[s];out={s,zeros?1:0,-1,-1};
        for(int k=0;k<3;k++)if(sign[k]==0)out.edge=k;
        if(zeros>=2){out.chart=2;for(int k=0;k<3;k++)if(X[0]==t.vertex[k][0]&&X[1]==t.vertex[k][1])out.vertex=k;require(out.vertex>=0,"exact original vertex chart, no tolerance snapping");}
    }
    require(out.triangle!=SIZE_MAX,"original generator lies in its original canonical material cell");return out;
}
inline Site transported_site(const Model&m,const std::vector<size_t>&triangles,const double*X,R mass){
    auto at=locate_generator(m,triangles,X);const auto&t=m.owned.f.triangle[at.triangle];
    Site out;out.cell=m.owned.f.cell[t.cell];out.source_triangle=m.owned.source_triangle[at.triangle];out.owner=m.layout->rank;out.chart=at.chart;out.source_mass=mass;
    for(int k=0;k<2;k++)out.reference[k]=X[k];
    sf::S position[2],velocity[2];
    if(at.chart==2){
        int local=-1;for(int a=0;a<10;a++){int index[3];gfs_material_p3_index(a,index);if(index[at.vertex]==3)local=t.id[a];}
        require(local>=0,"actual shared source vertex coefficient");
        for(int k=0;k<2;k++){out.position[k]=m.owned.f.coefficient[local].x[k];out.velocity[k]=m.owned.f.coefficient[local].v[k];position[k]=sf::exact(out.position[k]);velocity[k]=sf::exact(out.velocity[k]);}
    }else if(at.chart==1){
        const auto&e=t.edge[at.edge];auto canonical=sf::edge(e,t,at.edge);
        int dimension=sf::mag(canonical.Dref[0])>=sf::mag(canonical.Dref[1])?0:1;
        require(sf::lower(canonical.Dref[dimension])>0||sf::upper(canonical.Dref[dimension])<0,"nonzero original canonical edge direction");
        auto parameter=(sf::exact(X[dimension])-canonical.A[dimension])/canonical.Dref[dimension];
        require(sf::lower(parameter)>0&&sf::upper(parameter)<1,"exact nonvertex original edge chart");
        double u=double(parameter.n);require(st_trace(e,m.owned.f.coefficient.data(),u,out.position,out.velocity),"actual source trace transports generator");
        auto raw=sf::edge_raw(canonical,sf::D(parameter));
        for(int mode=0;mode<2;mode++)for(int k=0;k<2;k++){
            auto coefficient=[&](int a){const auto&q=m.owned.f.coefficient[a];return sf::exact(mode?q.v[k]:q.x[k]);};
            auto origin=coefficient(e.poly[0]);sf::S value=origin;
            for(int j=0;j<4;j++)value=value+(coefficient(e.poly[j])-origin)*sf::lagrange(j,sf::D(parameter)).v;
            for(int j=0;j<e.rank;j++){sf::S N;for(int c=0;c<2;c++)N=N+sf::exact(e.encode[j][c])*raw[c].v;value=value+coefficient(e.id[j])*N;}
            if(mode)velocity[k]=value;else position[k]=value;
        }
    }else{
        STPoint q={};q.X=X[0];q.Y=X[1];q.boundary=-1;
        sf::Q ux=sf::Q(t.vertex[1][0])-t.vertex[0][0],uy=sf::Q(t.vertex[1][1])-t.vertex[0][1];
        sf::Q vx=sf::Q(t.vertex[2][0])-t.vertex[0][0],vy=sf::Q(t.vertex[2][1])-t.vertex[0][1];
        sf::Q dx=sf::Q(X[0])-t.vertex[0][0],dy=sf::Q(X[1])-t.vertex[0][1],det=ux*vy-uy*vx;
        require(det>0,"same positive original reference triangle");q.l[1]=double((dx*vy-dy*vx)/det);q.l[2]=double((ux*dy-uy*dx)/det);q.l[0]=1-q.l[1]-q.l[2];
        STJet N[ST_MAX];GfsMaterialPoint12 point;
        require(st_shape(t,q,N)&&st_point(t,m.owned.f.coefficient.data(),N,&point),"actual full source transports generator, no IC orbit reset");
        auto map=sf::compile(m.owned.f,at.triangle);auto canonical=sf::evaluate(map,sf::exact(X[0]),sf::exact(X[1]));
        for(int k=0;k<2;k++){out.position[k]=point.position[k];out.velocity[k]=point.velocity[k];position[k]=canonical.position[k];velocity[k]=canonical.velocity[k];}
    }
    for(int k=0;k<2;k++){
        require(std::isfinite(out.position[k])&&std::isfinite(out.velocity[k]),"finite actual transported source generator");
        out.position_error[k]=sf::cache(out.position[k],position[k]).e;out.velocity_error[k]=sf::cache(out.velocity[k],velocity[k]).e;
        auto split=[&](double high,sf::S canonical,double&low,double&error,double&cache_error){
            low=double(canonical.n-sf::Q(high));
            sf::Q tail=sf::mf::plus_up(fabsq((sf::Q(high)+sf::Q(low))-canonical.n),canonical.e);
            error=nextafter(double(tail),INFINITY);
            cache_error=std::max(cache_error,nextafter(fabs(low)+error,INFINITY));
            require(std::isfinite(low)&&std::isfinite(error),"finite hi/lo canonical source transport, no position reset");
        };
        split(out.position[k],position[k],out.position_low[k],out.canonical_position_error[k],out.position_error[k]);
        split(out.velocity[k],velocity[k],out.velocity_low[k],out.canonical_velocity_error[k],out.velocity_error[k]);
        require(std::isfinite(out.position_error[k])&&std::isfinite(out.velocity_error[k]),"finite source-generator canonical transport enclosure");
    }out.canonical_bound=true;return out;
}
inline TargetSource target_source_inputs(const Model&m,const Inputs&input,const GeometryBound&g,const CurrentContext&context){
    require(m.layout&&m.layout->ranks==4&&g.canonical_global_injective&&g.exact_shared_trace&&m.current_context_bound,"actual accepted source geometry/context before target construction");
    uint64_t hash=compact_p3::bytes_hash(m.owned.f.coefficient.data(),m.owned.f.coefficient.size()*sizeof(STCoefficient));
    require(hash==g.source_coefficient_hash&&hash==context.coefficient_hash&&g.source_epoch==m.accepted.physical.epoch&&g.source_steps==m.accepted.physical.steps&&g.source_time==m.accepted.physical.time,"one immutable actual source epoch for target inputs");
    require(context.baseline.size()==m.owned.f.cell.size()&&g.triangle.size()==m.owned.f.triangle.size(),"same full source cells and geometry enclosures");
    TargetSource out;out.source_steps=g.source_steps;out.source_epoch=g.source_epoch;out.source_time=g.source_time;
    MPI_Allgather(&hash,1,MPI_UINT64_T,out.source_shard_hash.data(),1,MPI_UINT64_T,m.layout->comm);
    auto initial=raw_generators(input,m);std::vector<std::vector<size_t>>triangles(m.owned.f.cell.size());
    for(size_t s=0;s<m.owned.f.triangle.size();s++)triangles.at(m.owned.f.triangle[s].cell).push_back(s);
    out.site.resize(initial.size());bool valid=true;
    #pragma omp parallel for schedule(dynamic,32) reduction(&&:valid)
    for(size_t c=0;c<initial.size();c++)try{
        require(context.baseline[c].cell==m.owned.f.cell[c]&&context.baseline[c].mass>0,"actual accepted source mass identity, no legacy mass substitution");
        out.site[c]=transported_site(m,triangles[c],initial[c].data(),context.baseline[c].mass);
    }catch(const std::exception&){valid=false;}
    compact_mass::collective_require(valid,m.layout->comm,"all original source generators transported with canonical enclosures");
    std::set<uint64_t>edges;
    for(size_t s=0;s<m.owned.f.triangle.size();s++)for(int k=0;k<3;k++){
        const auto&t=m.owned.f.triangle[s];uint64_t global=m.owned.triangle_edge.at(3*s+k);if(!m.owned.external_edge.at(global))continue;
        require(edges.insert(global).second,"one original exterior curve per actual source owner");
        Exterior curve;curve.edge=global;curve.cell=m.owned.f.cell[t.cell];curve.source_triangle=m.owned.source_triangle[s];curve.owner=m.layout->rank;curve.local_edge=k;
        curve.orientation=t.endpoint[k][0]==k?1:-1;
        const auto&e=t.edge[k];curve.vertex[0]=m.owned.global[e.poly[0]];curve.vertex[1]=m.owned.global[e.poly[3]];
        curve.frozen=compact_geometry::freeze_trace(e,m.owned.f.coefficient,global,hash,g.source_epoch);
        for(int d=0;d<2;d++){
            for(int a=0;a<2;a++)curve.reference[a][d]=t.vertex[t.endpoint[k][a]][d];
            sf::Q displacement=g.triangle[s].displacement[d];
            sf::Q lower=sf::Q(std::min(curve.reference[0][d],curve.reference[1][d]))-displacement;
            sf::Q upper=sf::Q(std::max(curve.reference[0][d],curve.reference[1][d]))+displacement;
            curve.lower[d]=nextafter(double(lower),-INFINITY);curve.upper[d]=nextafter(double(upper),INFINITY);
            require(std::isfinite(curve.lower[d])&&std::isfinite(curve.upper[d]),"same-epoch whole-source enclosure of the actual exterior curve");
        }out.exterior.push_back(curve);
    }
    uint64_t counts[2]={out.site.size(),out.exterior.size()};MPI_Allreduce(MPI_IN_PLACE,counts,2,MPI_UINT64_T,MPI_SUM,m.layout->comm);
    require(counts[0]==71080&&counts[1]==g.external_edges,"complete current generator/exterior ownership, no source strips discarded");out.global_sites=counts[0];out.global_exterior=counts[1];
    for(const auto&s:out.site)for(int k=0;k<2;k++){out.maximum_position_error=std::max(out.maximum_position_error,R(s.position_error[k]));out.maximum_velocity_error=std::max(out.maximum_velocity_error,R(s.velocity_error[k]));}
    MPI_Allreduce(MPI_IN_PLACE,&out.maximum_position_error,1,MPI_LONG_DOUBLE,MPI_MAX,m.layout->comm);MPI_Allreduce(MPI_IN_PLACE,&out.maximum_velocity_error,1,MPI_LONG_DOUBLE,MPI_MAX,m.layout->comm);
    return out;
}
} // namespace compact_source
#endif
