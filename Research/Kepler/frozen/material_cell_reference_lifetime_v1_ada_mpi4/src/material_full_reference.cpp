// Native full-domain geometry/canonical reference front end. No Python physics.
#include <mpi.h>
#include <cuda_runtime.h>
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <fstream>
#include <map>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#include <unistd.h>
#include "exam_gpu.h"
#include "material_full_reference.h"
using R=long double;using Q=__float128;using V=std::array<double,2>;
static void check(bool p,const char*s){if(!p)throw std::runtime_error(s);}
template<class T>static void put(FILE*f,const T*p,size_t n=1){check(fwrite(p,sizeof(T),n,f)==n,"native reference write");}
template<class T>static void get(FILE*f,T*p,size_t n=1){check(fread(p,sizeof(T),n,f)==n,"native reference read");}
static void close_durable(FILE*f){check(fflush(f)==0&&fsync(fileno(f))==0&&fclose(f)==0,"durable reference close");}
struct Cell {MFCell s;std::vector<MFFace>face;uint64_t offset=0;};
struct DSU {
    std::vector<uint64_t>p;
    explicit DSU(size_t n):p(n){std::iota(p.begin(),p.end(),0);}
    uint64_t root(uint64_t i){while(p[i]!=i){p[i]=p[p[i]];i=p[i];}return i;}
    void join(uint64_t a,uint64_t b){a=root(a);b=root(b);if(a!=b)p[std::max(a,b)]=std::min(a,b);}
};
struct Corner {Q X,Y,D;}; // exact rational binary128 numerators, not coordinates
struct Triangle {uint64_t cell,v[3];};
struct Mesh {std::vector<V>vertex;std::vector<Triangle>triangle;};
struct Counts {
    uint64_t cells=0,raw_faces=0,exact_zero_edges=0,paired=0,external=0;
    uint64_t vertices=0,triangles=0,edges=0,rank1=0,rank2=0,dofs=0;
    uint64_t quadrature_slots=0,positive_quadrature_slots=0;
    R area=0,legacy_mass=0,legacy_U=0,legacy_K=0,min_area=LDBL_MAX;
    R raw_vertex_shift=0,max_area_change=0,max_input_volume_relative_difference=0;
};
static int wrap(double x,double canonical,double box){
    double n=nearbyint((x-canonical)/box);check(std::isfinite(n)&&fabs(n)<=1&&x==canonical+n*box,"exact native periodic image identity");return int(n);
}
static void exact_input(double x){
    // Existing IC coordinates lie on 3/512 lattice. At |x|<64 all required
    // degree<=5 determinant numerators fit within113 significand bits.
    // This is an algebraic predicate for this retained IC, not near merging.
    check(std::isfinite(x)&&fabs(x)<64&&x*512==nearbyint(x*512),"retained dyadic IC required by exact topology predicate");
}
static void export_native(GPUContext*ctx,int np,int rank){
    check(np>=0&&np<=ctx->h_parts.n_total,"native owned geometry extent");
    std::vector<PolygonOut>poly(np);check(cudaMemcpy(poly.data(),ctx->d_polygons,np*sizeof(PolygonOut),cudaMemcpyDeviceToHost)==cudaSuccess,"download actual raw CUDA polygons");
    std::string path="material_native_geometry_r"+std::to_string(rank)+".bin";FILE*f=fopen(path.c_str(),"wbx");check(f,"fresh owned reference shard");
    uint64_t header[8]={0x4d46524157474531ULL,1,uint64_t(np),sizeof(MFCell),sizeof(MFFace),uint64_t(rank),0,0};put(f,header,8);
    for(int i=0;i<np;i++){
        const auto&p=poly[i];auto&h=ctx->h_parts;check(p.n_corners>=3&&p.n_corners<=GPU_MAX_CORNERS&&!p.halo_clipped,"complete unclipped native raw cell");
        MFCell c={uint64_t(h.indx[i]),uint64_t(p.n_corners),h.x[i],h.y[i],h.mass[i],h.volume[i],h.vx[i],h.vy[i],h.ie[i],h.pressure[i]};
        check(h.w2[i]==0,"initial W0 native geometry");put(f,&c);
        for(int k=0;k<p.n_corners;k++){
            int j=p.edge_neigh[k],next=(k+1)%p.n_corners;check(j>=0&&j<h.n_total&&uint64_t(h.indx[j])<72057594037927935ULL&&h.w2[j]==0,"complete real support plane including periodic image");
            MFFace e={uint64_t(h.indx[j]),h.x[j],h.y[j],p.cx[k],p.cy[k],p.cx[next],p.cy[next]};put(f,&e);
        }
    }close_durable(f);
}
static std::vector<Cell> read_all(int ranks){
    std::map<uint64_t,Cell>all;
    for(int rank=0;rank<ranks;rank++){
        std::string path="material_native_geometry_r"+std::to_string(rank)+".bin";FILE*f=fopen(path.c_str(),"rb");check(f,"actual native geometry shard");uint64_t h[8];get(f,h,8);
        check(h[0]==0x4d46524157474531ULL&&h[1]==1&&h[3]==sizeof(MFCell)&&h[4]==sizeof(MFFace)&&h[5]==uint64_t(rank),"owned reference ABI");
        for(uint64_t n=0;n<h[2];n++){Cell c;get(f,&c.s);check(c.s.n>=3&&c.s.n<=GPU_MAX_CORNERS&&c.s.mass>0&&c.s.U>0,"physical native IC cell");c.face.resize(c.s.n);get(f,c.face.data(),c.face.size());check(all.emplace(c.s.id,std::move(c)).second,"unique real global cell owner");}
        check(fgetc(f)==EOF&&fclose(f)==0,"complete native reference shard");
    }
    check(all.size()==71080,"retained full fine_integrated IC count");std::vector<Cell>out;out.reserve(all.size());
    for(auto&v:all){check(v.first==out.size(),"complete initial global cell IDs");out.push_back(std::move(v.second));}return out;
}
static Mesh canonical_mesh(std::vector<Cell>&cell,double box,Counts&count){
    count.cells=cell.size();uint64_t total=0;
    for(auto&c:cell){c.offset=total;total+=c.face.size();exact_input(c.s.x);exact_input(c.s.y);count.legacy_mass+=c.s.mass;count.legacy_U+=c.s.U;count.legacy_K+=.5L*c.s.mass*(R(c.s.vx)*c.s.vx+R(c.s.vy)*c.s.vy);}
    count.raw_faces=total;DSU dsu(total);std::vector<Corner>corner(total);std::vector<bool>zero(total);
    for(auto&c:cell){
        size_t n=c.face.size();for(size_t k=0;k<n;k++){
            const auto&a=c.face[(k+n-1)%n];const auto&b=c.face[k];check(a.neighbor<cell.size()&&b.neighbor<cell.size(),"all support cell IDs present");exact_input(a.nx);exact_input(a.ny);exact_input(b.nx);exact_input(b.ny);
            Q ax=Q(a.nx)-c.s.x,ay=Q(a.ny)-c.s.y,bx=Q(b.nx)-c.s.x,by=Q(b.ny)-c.s.y;
            Q ar=(ax*ax+ay*ay)/2,br=(bx*bx+by*by)/2,D=ax*by-ay*bx;check(D!=0,"distinct nonparallel actual vertex support planes");
            Q X=ar*by-ay*br,Y=ax*br-ar*bx;corner[c.offset+k]={X+Q(c.s.x)*D,Y+Q(c.s.y)*D,D};
        }
        for(size_t k=0;k<n;k++){auto&a=corner[c.offset+k];auto&b=corner[c.offset+(k+1)%n];
            if(a.X*b.D==b.X*a.D&&a.Y*b.D==b.Y*a.D){zero[c.offset+k]=true;dsu.join(c.offset+k,c.offset+(k+1)%n);count.exact_zero_edges++;}}
    }
    struct Half {uint64_t source,k;};std::map<std::pair<uint64_t,uint64_t>,std::vector<Half>>paired;
    for(auto&c:cell)for(size_t k=0;k<c.face.size();k++)if(!zero[c.offset+k]){
        const auto&e=c.face[k];const auto&other=cell[e.neighbor].s;
        int wx=wrap(e.nx,other.x,box),wy=wrap(e.ny,other.y,box);
        if(wx||wy){count.external++;continue;}
        check(c.s.id!=e.neighbor,"nonself interior face");paired[std::minmax(c.s.id,e.neighbor)].push_back({c.s.id,k});
    }
    for(auto&v:paired){check(v.second.size()==2,"each nonzero interior face has exactly two opposite owners");auto a=v.second[0],b=v.second[1];check(a.source!=b.source&&cell[a.source].face[a.k].neighbor==b.source&&cell[b.source].face[b.k].neighbor==a.source,"reciprocal native interior neighbor IDs");
        auto&ca=cell[a.source];auto&cb=cell[b.source];dsu.join(ca.offset+a.k,cb.offset+(b.k+1)%cb.face.size());dsu.join(ca.offset+(a.k+1)%ca.face.size(),cb.offset+b.k);count.paired++;
    }
    Mesh mesh;std::map<uint64_t,uint64_t>index;
    for(auto&c:cell)for(size_t k=0;k<c.face.size();k++){
        auto id=c.offset+k,root=dsu.root(id);const auto&a=corner[id];const auto&b=corner[root];
        check(a.X*b.D==b.X*a.D&&a.Y*b.D==b.Y*a.D,"shared vertex exact support-plane coincidence after topology union");
        if(!index.count(root)){uint64_t next=mesh.vertex.size();index[root]=next;mesh.vertex.push_back({double(b.X/b.D-Q(box)/2),double(b.Y/b.D-Q(box)/2)});}
        V x=mesh.vertex[index.at(root)];count.raw_vertex_shift=std::max(count.raw_vertex_shift,hypotl(R(x[0])+box/2-c.s.x-c.face[k].ax,R(x[1])+box/2-c.s.y-c.face[k].ay));
    }
    FILE*physical=fopen("material_reference_cells.csv","wx");check(physical,"fresh full reference cell ledger");fprintf(physical,"id,native_mass,native_U,native_K,native_float_volume,canonical_area,canonical_minus_raw_area\n");
    for(auto&c:cell){
        std::vector<uint64_t>nodes;for(size_t k=0;k<c.face.size();k++)if(!zero[c.offset+k])nodes.push_back(index.at(dsu.root(c.offset+k)));
        check(nodes.size()>=3&&std::set<uint64_t>(nodes.begin(),nodes.end()).size()==nodes.size(),"simple positive-area material polygon nodes");
        R area2=0,cx=0,cy=0,raw_area=0;R ox=R(c.s.x)-box/2,oy=R(c.s.y)-box/2;
        for(size_t k=0;k<nodes.size();k++){auto a=mesh.vertex[nodes[k]],b=mesh.vertex[nodes[(k+1)%nodes.size()]];R ax=R(a[0])-ox,ay=R(a[1])-oy,bx=R(b[0])-ox,by=R(b[1])-oy,z=ax*by-ay*bx;area2+=z;cx+=(ax+bx)*z;cy+=(ay+by)*z;}
        for(const auto&e:c.face)raw_area+=R(e.ax)*e.by-R(e.ay)*e.bx;
        check(area2>0,"positive canonical material polygon orientation");R area=area2/2;count.area+=area;count.min_area=std::min(count.min_area,area);count.max_area_change=std::max(count.max_area_change,fabsl((area2-raw_area)/area2));count.max_input_volume_relative_difference=std::max(count.max_input_volume_relative_difference,fabsl(R(c.s.volume)/area-1));
        V center={double(ox+cx/(3*area2)),double(oy+cy/(3*area2))};uint64_t apex=mesh.vertex.size();mesh.vertex.push_back(center);
        for(size_t k=0;k<nodes.size();k++){auto a=mesh.vertex[nodes[k]],b=mesh.vertex[nodes[(k+1)%nodes.size()]];check((R(a[0])-center[0])*(R(b[1])-center[1])-(R(a[1])-center[1])*(R(b[0])-center[0])>0,"positive centroid fan triangle");mesh.triangle.push_back({c.s.id,{apex,nodes[k],nodes[(k+1)%nodes.size()]}});}
        fprintf(physical,"%llu,%.21Lg,%.21Lg,%.21Lg,%.21Lg,%.21Lg,%.21Lg\n",(unsigned long long)c.s.id,R(c.s.mass),R(c.s.U),.5L*c.s.mass*(R(c.s.vx)*c.s.vx+R(c.s.vy)*c.s.vy),R(c.s.volume),area,(area2-raw_area)/2);
    }close_durable(physical);check(fabsl(count.area-R(box)*box)<1e-10L,"full canonical material reference covers original box area");
    // Same single h/2 conforming reference subdivision as the retained patch.
    std::map<std::pair<uint64_t,uint64_t>,uint64_t>mid;std::vector<Triangle>refined;refined.reserve(mesh.triangle.size()*4);
    auto midpoint=[&](uint64_t a,uint64_t b){auto key=std::minmax(a,b);auto it=mid.find(key);if(it!=mid.end())return it->second;auto x=mesh.vertex[a],y=mesh.vertex[b];uint64_t i=mesh.vertex.size();mesh.vertex.push_back({double((R(x[0])+y[0])/2),double((R(x[1])+y[1])/2)});mid[key]=i;return i;};
    for(auto&t:mesh.triangle){auto a=t.v[0],b=t.v[1],c=t.v[2],ab=midpoint(a,b),bc=midpoint(b,c),ca=midpoint(c,a);for(auto v:std::vector<std::array<uint64_t,3>>{{a,ab,ca},{ab,b,bc},{ca,bc,c},{ab,bc,ca}})refined.push_back({t.cell,{v[0],v[1],v[2]}});}
    mesh.triangle.swap(refined);count.vertices=mesh.vertex.size();count.triangles=mesh.triangle.size();return mesh;
}
static void write_mesh_plan(const Mesh&mesh,double box,Counts&c){
    std::map<std::pair<uint64_t,uint64_t>,unsigned>edge;
    for(auto&t:mesh.triangle)for(int k=0;k<3;k++)edge[std::minmax(t.v[k],t.v[(k+1)%3])]++;
    for(auto&v:edge){check(v.second==1||v.second==2,"conforming refined full reference manifold");auto a=mesh.vertex[v.first.first],b=mesh.vertex[v.first.second];Q cross=Q(a[0])*b[1]-Q(a[1])*b[0];if(cross==0)c.rank1++;else c.rank2++;}
    c.edges=edge.size();c.dofs=c.vertices+2*c.edges+c.triangles+c.rank1+2*c.rank2+2*c.triangles;
    check(c.vertices-c.edges+c.triangles==1,"full open material disk Euler characteristic");
    FILE*out=fopen("material_canonical_reference.bin","wbx");check(out,"fresh canonical reference mesh");uint64_t h[8]={0x4d4643414e4f4e31ULL,1,c.cells,c.vertices,c.triangles,c.edges,sizeof(V),sizeof(Triangle)};put(out,h,8);put(out,mesh.vertex.data(),mesh.vertex.size());put(out,mesh.triangle.data(),mesh.triangle.size());close_durable(out);
    // Persist compact immutable per-triangle quadrature plans, not billions
    // of STPoint or GPOutput records. No field quadrature is evaluated here.
    GfsEnrichConfig cfg={};cfg.na=64;cfg.nr=28;cfg.field.ic.kind=2;cfg.field.ic.cx=cfg.field.ic.cy=box/2;GfsEnrichBasis basis={};basis.generator[0]=basis.generator[1]=box/2;
    FILE*plans=fopen("material_reference_quadrature_plan.bin","wbx");check(plans,"compact fixed quadrature plan");uint64_t qh[6]={0x4d46515541443131ULL,1,c.triangles,64,28,sizeof(GfsEnrichQuadrature)};put(plans,qh,6);
    for(auto&t:mesh.triangle){double x[3],y[3];for(int k=0;k<3;k++){x[k]=mesh.vertex[t.v[k]][0];y[k]=mesh.vertex[t.v[k]][1];}GfsMomentPolygon p={x,y,3,0};GfsEnrichQuadrature q={};check(gfs_enrich_quadrature_init(p,&cfg,&basis,&q),"full material fixed quadrature plan");c.quadrature_slots+=gfs_enrich_nquad(p,&cfg,&q);
        for(int k=0;k<3;k++){int j=(k+1)%3;double bx=x[k]-q.anchor[0],by=y[k]-q.anchor[1],cx=x[j]-q.anchor[0],cy=y[j]-q.anchor[1];R jac=R(bx)*cy-R(by)*cx;if(k==q.anchor_edge||(bx==0&&by==0)||(cx==0&&cy==0))jac=0;check(jac>=0,"nonnegative full material quadrature support");if(jac>0)c.positive_quadrature_slots+=uint64_t(q.ns)*64*28;}
        put(plans,&q);
    }close_durable(plans);
}
extern "C" int gfs_material_prepare_reference(void*gpu,int np,int fcomm,double gamma,double pressure,double eps,double box){
    MPI_Comm comm=MPI_Comm_f2c(fcomm);int rank,ranks;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&ranks);double start=MPI_Wtime();int bad=0;
    try{check(eps>0&&gamma>1&&pressure>0&&box==24,"retained continuous physical parameters");export_native(static_cast<GPUContext*>(gpu),np,rank);
        int device=-1;cudaDeviceProp device_prop;check(cudaGetDevice(&device)==cudaSuccess&&cudaGetDeviceProperties(&device_prop,device)==cudaSuccess,"actual native GPU identity");char name[128];snprintf(name,sizeof(name),"material_runtime_r%d.json",rank);FILE*f=fopen(name,"wx");check(f,"runtime device record");fprintf(f,"{\"rank\":%d,\"MPI_ranks\":%d,\"visible_device_ordinal\":%d,\"device_UUID_hex\":\"",rank,ranks,device);for(unsigned char ch:device_prop.uuid.bytes)fprintf(f,"%02x",unsigned(ch));fprintf(f,"\",\"owned_cells\":%d,\"geometry_from_original_CUDA\":true}\n",np);close_durable(f);
    }catch(const std::exception&e){fprintf(stderr,"rank%d MATERIAL_REFERENCE_EXPORT_FAIL %s\n",rank,e.what());bad=1;}
    int any=0;MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,comm);if(any)return 2;MPI_Barrier(comm);
    if(rank==0)try{auto cells=read_all(ranks);Counts count;auto mesh=canonical_mesh(cells,box,count);write_mesh_plan(mesh,box,count);
        FILE*f=fopen("material_full_reference.json","wx");check(f,"full material preparation result");fprintf(f,"{\"scope\":\"actual original native full fine_integrated IC and raw CUDA reference geometry; canonical reference and streaming quadrature plans; no physical evolution\",\"MPI_ranks\":%d,\"cells\":%llu,\"raw_face_ends\":%llu,\"exact_zero_edges_collapsed\":%llu,\"paired_physical_faces\":%llu,\"external_physical_faces\":%llu,\"reference_vertices\":%llu,\"reference_triangles\":%llu,\"reference_edges\":%llu,\"canonical_edge_rank1\":%llu,\"canonical_edge_rank2\":%llu,\"prospective_scalar_dofs\":%llu,\"fixed_quadrature_slots\":%llu,\"positive_quadrature_slots_estimate\":%llu,\"total_area\":%.21Lg,\"min_cell_area\":%.21Lg,\"legacy_mass\":%.21Lg,\"legacy_U\":%.21Lg,\"legacy_point_generator_K\":%.21Lg,\"raw_to_canonical_vertex_max_shift\":%.21Lg,\"canonical_vs_raw_cell_area_max_relative_difference\":%.21Lg,\"stored_float_volume_vs_canonical_max_relative_difference\":%.21Lg,\"boundary\":\"moving_material_constant_Pext_P0\",\"gamma\":%.17g,\"P0\":%.17g,\"eps\":%.17g,\"exact_support_topology_no_proximity_merge\":true,\"reference_preparation_complete\":true,\"full_field_initialized\":false,\"distributed_mass_operator_implemented\":false,\"accepted_steps\":0,\"long_Kepler_validated\":false,\"seconds\":%.9g}\n",ranks,(unsigned long long)count.cells,(unsigned long long)count.raw_faces,(unsigned long long)count.exact_zero_edges,(unsigned long long)count.paired,(unsigned long long)count.external,(unsigned long long)count.vertices,(unsigned long long)count.triangles,(unsigned long long)count.edges,(unsigned long long)count.rank1,(unsigned long long)count.rank2,(unsigned long long)count.dofs,(unsigned long long)count.quadrature_slots,(unsigned long long)count.positive_quadrature_slots,count.area,count.min_area,count.legacy_mass,count.legacy_U,count.legacy_K,count.raw_vertex_shift,count.max_area_change,count.max_input_volume_relative_difference,gamma,pressure,eps,MPI_Wtime()-start);close_durable(f);
    }catch(const std::exception&e){fprintf(stderr,"MATERIAL_REFERENCE_BUILD_FAIL %s\n",e.what());FILE*f=fopen("material_reference_failure.txt","wx");if(f){fprintf(f,"%s\n",e.what());fclose(f);}bad=1;}
    MPI_Bcast(&bad,1,MPI_INT,0,comm);return bad?2:0;
}
