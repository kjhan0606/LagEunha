#ifndef MATERIAL_FULL_OWNED_H
#define MATERIAL_FULL_OWNED_H
#include <atomic>
#include <climits>
#include <unistd.h>
#include "material_shared_trace_factory.h"
#include "gfs_enriched_euler.h"
#include "material_full_reference.h"
namespace material_full {
using R=long double;using Pt=shared_trace::Pt;using Key=std::pair<uint64_t,uint64_t>;
#include "material_full_edge_basis.h"
struct RefTriangle {uint64_t cell,v[3];};
struct RefEdge {uint64_t a,b,prefix=0;int owner=INT_MAX,rank=0,count=0;};
struct Owned {
    shared_trace::Factory f;std::vector<uint64_t>global,source_triangle,source_cell,triangle_edge;
    std::vector<int>owner;std::vector<GfsEnrichQuadrature>rule;
    uint64_t global_vertices=0,global_edges=0,global_triangles=0,global_dofs=0;
    std::vector<uint8_t>external_edge;
};
template<class T>static void read(FILE*f,T*p,size_t n=1){need(fread(p,sizeof(T),n,f)==n,"full owned input read");}
template<class T>static void write(FILE*f,const T*p,size_t n=1){need(fwrite(p,sizeof(T),n,f)==n,"full owned state write");}
static void durable(FILE*f){need(fflush(f)==0&&fsync(fileno(f))==0&&fclose(f)==0,"durable full owned state close");}
static Owned construct(const std::string&input,MPI_Comm comm){
    int rank,ranks;MPI_Comm_rank(comm,&rank);MPI_Comm_size(comm,&ranks);Owned out;
    FILE*file=fopen((input+"/material_canonical_reference.bin").c_str(),"rb");need(file,"saved full canonical reference");uint64_t h[8];read(file,h,8);need(h[0]==0x4d4643414e4f4e31ULL&&h[1]==1&&h[6]==sizeof(Pt)&&h[7]==sizeof(RefTriangle),"full canonical reference ABI");
    uint64_t nc=h[2],nv=h[3],nt=h[4],ne=h[5];need(nc>0&&nv<INT_MAX&&nt<INT_MAX&&ne<INT_MAX,"runtime full reference capacities");
    std::vector<Pt>vertex(nv);std::vector<RefTriangle>tri(nt);read(file,vertex.data(),nv);read(file,tri.data(),nt);need(fgetc(file)==EOF&&fclose(file)==0,"complete saved full canonical reference");
    std::vector<int>cell_owner(nc,-1),vertex_owner(nv,INT_MAX);
    for(int p=0;p<ranks;p++){
        file=fopen((input+"/material_native_geometry_r"+std::to_string(p)+".bin").c_str(),"rb");need(file,"retained native cell ownership");uint64_t b[8];read(file,b,8);need(b[0]==0x4d46524157474531ULL&&b[1]==1&&b[3]==sizeof(MFCell)&&b[4]==sizeof(MFFace)&&b[5]==uint64_t(p),"native owner shard identity");
        for(uint64_t j=0;j<b[2];j++){MFCell c;read(file,&c);need(c.id<nc&&cell_owner[c.id]==-1,"unique source material owner");cell_owner[c.id]=p;need(fseek(file,c.n*sizeof(MFFace),SEEK_CUR)==0,"skip saved native face metadata");}
        need(fgetc(file)==EOF&&fclose(file)==0,"complete owner shard");
    }
    need(std::find(cell_owner.begin(),cell_owner.end(),-1)==cell_owner.end(),"all source cells owned");
    std::vector<int>local_cell(nc,-1);for(size_t c=0;c<nc;c++)if(cell_owner[c]==rank){local_cell[c]=out.f.cell.size();out.f.cell.push_back(c);}
    std::map<Key,size_t>index;
    for(auto&t:tri){need(t.cell<nc,"source context range");for(int a=0;a<3;a++){need(t.v[a]<nv,"canonical vertex ID");vertex_owner[t.v[a]]=std::min(vertex_owner[t.v[a]],cell_owner[t.cell]);index.emplace(std::minmax(t.v[a],t.v[(a+1)%3]),0);}}
    need(index.size()==ne,"saved canonical edge count");std::vector<RefEdge>edge(ne);size_t ei=0;uint64_t residuals=0;
    for(auto&entry:index){entry.second=ei;auto&a=vertex[entry.first.first];auto&b=vertex[entry.first.second];int r=__float128(a[0])*b[1]-__float128(a[1])*b[0]==0?1:2;edge[ei]={entry.first.first,entry.first.second,residuals,INT_MAX,r,0};residuals+=r;ei++;}
    std::vector<uint8_t>needed(ne);
    for(auto&t:tri)for(int a=0;a<3;a++){size_t e=index.at(std::minmax(t.v[a],t.v[(a+1)%3]));edge[e].owner=std::min(edge[e].owner,cell_owner[t.cell]);edge[e].count++;if(cell_owner[t.cell]==rank)needed[e]=1;}
    uint64_t npoly=nv+2*ne+nt,bubble=npoly+residuals,ndof=bubble+2*nt;need(ndof<INT_MAX,"full DOF integer kernel capacity");out.global_vertices=nv;out.global_edges=ne;out.global_triangles=nt;out.global_dofs=ndof;
    std::vector<int>global_owner(ndof,INT_MAX);for(size_t v=0;v<nv;v++)global_owner[v]=vertex_owner[v];
    for(size_t e=0;e<ne;e++){const auto&v=edge[e];global_owner[nv+2*e]=global_owner[nv+2*e+1]=v.owner;for(int j=0;j<v.rank;j++)global_owner[npoly+v.prefix+j]=v.owner;}
    for(size_t s=0;s<nt;s++){int p=cell_owner[tri[s].cell];global_owner[nv+2*ne+s]=p;global_owner[bubble+2*s]=global_owner[bubble+2*s+1]=p;}
    auto add=[&](uint64_t id){out.global.push_back(id);};
    for(size_t s=0;s<nt;s++)if(cell_owner[tri[s].cell]==rank){const auto&t=tri[s];out.source_triangle.push_back(s);out.source_cell.push_back(t.cell);for(int a=0;a<3;a++)add(t.v[a]);add(nv+2*ne+s);add(bubble+2*s);add(bubble+2*s+1);for(int a=0;a<3;a++){auto e=index.at(std::minmax(t.v[a],t.v[(a+1)%3]));out.triangle_edge.push_back(e);add(nv+2*e);add(nv+2*e+1);for(int j=0;j<edge[e].rank;j++)add(npoly+edge[e].prefix+j);}}
    std::sort(out.global.begin(),out.global.end());out.global.erase(std::unique(out.global.begin(),out.global.end()),out.global.end());need(out.global.size()<INT_MAX,"owned/ghost local DOF capacity");
    std::vector<int>local(ndof,-1);out.owner.resize(out.global.size());out.f.coefficient.resize(out.global.size());out.f.dof.resize(out.global.size());
    for(size_t i=0;i<out.global.size();i++){uint64_t g=out.global[i];local[g]=i;out.owner[i]=global_owner[g];int kind=g<nv?0:g<nv+2*ne?1:g<npoly?2:g<bubble?3:4;out.f.dof[i]={g,g,kind,out.owner[i]};}
    std::vector<size_t>owned_edge;for(size_t e=0;e<ne;e++)if(needed[e])owned_edge.push_back(e);
    std::vector<STEdge>cache(owned_edge.size());std::vector<int>cache_index(ne,-1);for(size_t i=0;i<owned_edge.size();i++)cache_index[owned_edge[i]]=i;
    std::vector<double>gx,gw;shared_trace::gauss(64,gx,gw);std::atomic<int>first_bad{-1};int failures=0;
    #pragma omp parallel for schedule(static) reduction(+:failures)
    for(size_t k=0;k<owned_edge.size();k++)try{size_t e=owned_edge[k];auto&z=cache[k];cached_edge_basis(z,vertex[edge[e].a],vertex[edge[e].b],.047,gx,gw);need(z.rank==edge[e].rank,"canonical edge rank identity");z.poly[0]=local[edge[e].a];z.poly[3]=local[edge[e].b];z.poly[1]=local[nv+2*e];z.poly[2]=local[nv+2*e+1];for(int j=0;j<z.rank;j++)z.id[j]=local[npoly+edge[e].prefix+j];}catch(...){failures++;int expected=-1;first_bad.compare_exchange_strong(expected,int(k));}
    int all_failures=0;MPI_Allreduce(&failures,&all_failures,1,MPI_INT,MPI_SUM,comm);
    if(failures){size_t e=owned_edge[first_bad];fprintf(stderr,"rank%d FULL_EDGE_BASIS_FAIL global_edge=%zu A=(%.17g,%.17g) B=(%.17g,%.17g) local_failures=%d\n",rank,e,vertex[edge[e].a][0],vertex[edge[e].a][1],vertex[edge[e].b][0],vertex[edge[e].b][1],failures);}need(!all_failures,"full native canonical edge basis construction");
    for(size_t i=0;i<out.global.size();i++)if(out.global[i]<npoly){uint64_t g=out.global[i];Pt x;
        if(g<nv)x=vertex[g];else if(g<nv+2*ne){size_t e=(g-nv)/2;int part=(g-nv)%2+1;for(int d=0;d<2;d++)x[d]=double(((R)(3-part)*vertex[edge[e].a][d]+R(part)*vertex[edge[e].b][d])/3);}else{const auto&t=tri[g-(nv+2*ne)];for(int d=0;d<2;d++)x[d]=double((R(vertex[t.v[0]][d])+vertex[t.v[1]][d]+vertex[t.v[2]][d])/3);}
        auto v=shared_trace::initial_velocity(x,.047);for(int d=0;d<2;d++){out.f.coefficient[i].x[d]=x[d];out.f.coefficient[i].v[d]=v[d];}
    }
    for(auto&z:cache)for(int j=0;j<z.rank;j++)for(int d=0;d<2;d++)out.f.coefficient[z.id[j]].v[d]=z.decode[d][j];
    out.f.triangle.resize(out.source_triangle.size());failures=0;first_bad=-1;
    #pragma omp parallel for schedule(static) reduction(+:failures)
    for(size_t s=0;s<out.source_triangle.size();s++)try{
        uint64_t source=out.source_triangle[s];const auto&r=tri[source];auto&t=out.f.triangle[s];t.eps=.047;t.cell=local_cell[r.cell];GfsMaterialTriangle g={};
        for(int a=0;a<3;a++)for(int d=0;d<2;d++)g.x[a][d]=t.vertex[a][d]=vertex[r.v[a]][d];need(gfs_material_triangle_init(&g),"full owned positive triangle");for(int a=0;a<3;a++)for(int d=0;d<2;d++)t.grad[a][d]=g.grad[a][d];
        for(int a=0;a<10;a++){int m[3];gfs_material_p3_index(a,m);int v=-1,z=-1;for(int j=0;j<3;j++){if(m[j]==3)v=j;if(m[j]==0)z=j;}if(v>=0)t.id[a]=local[r.v[v]];else if(z>=0){int i=(z+1)%3,j=(z+2)%3;size_t e=index.at(std::minmax(r.v[i],r.v[j]));const auto&E=cache[cache_index[e]];int part=r.v[i]<r.v[j]?m[j]:m[i];t.id[a]=E.poly[part];}else t.id[a]=local[nv+2*ne+source];
            double Dq[4];const auto&c=out.f.coefficient[t.id[a]];need(gfs_material_p3_radial_remainder(t.vertex[0],c.x[0],c.x[1],t.eps,t.qnode[a],Dq),"full initial fixed remainder cache");}
        t.n=10;for(int a=0;a<3;a++){int b=(a+1)%3;size_t e=out.triangle_edge[s*3+a];t.edge[a]=cache[cache_index[e]];t.endpoint[a][r.v[a]<r.v[b]?0:1]=a;t.endpoint[a][r.v[a]<r.v[b]?1:0]=b;t.offset[a]=t.n;for(int j=0;j<t.edge[a].rank;j++)t.id[t.n++]=t.edge[a].id[j];}
        t.bubble_offset=t.n;for(int j=0;j<2;j++){t.id[t.n++]=local[bubble+2*source+j];t.bubble_encode[j][j]=t.bubble_decode[j][j]=1;out.f.coefficient[t.id[t.bubble_offset+j]].v[j]=1;}
    }catch(...){failures++;int expected=-1;first_bad.compare_exchange_strong(expected,int(s));}
    MPI_Allreduce(&failures,&all_failures,1,MPI_INT,MPI_SUM,comm);if(failures)fprintf(stderr,"rank%d FULL_TRIANGLE_BASIS_FAIL global_triangle=%llu failures=%d\n",rank,(unsigned long long)out.source_triangle[first_bad],failures);need(!all_failures,"full native owned triangle basis construction");
    out.external_edge.resize(ne);for(size_t e=0;e<ne;e++){need(edge[e].count==1||edge[e].count==2,"retained global edge manifold");out.external_edge[e]=edge[e].count==1;}
    file=fopen((input+"/material_reference_quadrature_plan.bin").c_str(),"rb");need(file,"saved compact full quadrature plan");uint64_t qh[6];read(file,qh,6);need(qh[0]==0x4d46515541443131ULL&&qh[2]==nt&&qh[3]==64&&qh[4]==28&&qh[5]==sizeof(GfsEnrichQuadrature),"full compact rule identity");
    std::vector<GfsEnrichQuadrature>chunk(8192);size_t selected=0;out.rule.reserve(out.source_triangle.size());
    for(size_t at=0;at<nt;at+=chunk.size()){size_t n=std::min(chunk.size(),size_t(nt-at));read(file,chunk.data(),n);for(size_t j=0;j<n;j++)if(cell_owner[tri[at+j].cell]==rank){need(out.source_triangle[selected++]==at+j,"fixed material rule source order");out.rule.push_back(chunk[j]);}}
    need(fgetc(file)==EOF&&fclose(file)==0&&out.rule.size()==out.f.triangle.size(),"all owned compact rules restored");return out;
}
static void save_basis(const Owned&m,const std::string&name){
    FILE*f=fopen(name.c_str(),"wbx");need(f,"new owned basis state");uint64_t h[12]={0x4d464f574e454431ULL,1,m.global_vertices,m.global_edges,m.global_triangles,m.global_dofs,m.f.cell.size(),m.f.triangle.size(),m.global.size(),sizeof(STTriangle),sizeof(STCoefficient),sizeof(GfsEnrichQuadrature)};write(f,h,12);write(f,m.global.data(),m.global.size());write(f,m.owner.data(),m.owner.size());write(f,m.f.cell.data(),m.f.cell.size());write(f,m.source_triangle.data(),m.source_triangle.size());write(f,m.source_cell.data(),m.source_cell.size());write(f,m.triangle_edge.data(),m.triangle_edge.size());write(f,m.f.triangle.data(),m.f.triangle.size());write(f,m.f.coefficient.data(),m.f.coefficient.size());write(f,m.rule.data(),m.rule.size());write(f,m.external_edge.data(),m.external_edge.size());durable(f);
}
static Owned load_basis(const std::string&name){
    FILE*f=fopen(name.c_str(),"rb");need(f,"saved full owned field");uint64_t h[12];read(f,h,12);need(h[0]==0x4d464f574e454431ULL&&h[1]==1&&h[9]==sizeof(STTriangle)&&h[10]==sizeof(STCoefficient)&&h[11]==sizeof(GfsEnrichQuadrature),"saved owned field version/ABI");Owned m;
    m.global_vertices=h[2];m.global_edges=h[3];m.global_triangles=h[4];m.global_dofs=h[5];m.f.cell.resize(h[6]);m.source_triangle.resize(h[7]);m.source_cell.resize(h[7]);m.triangle_edge.resize(h[7]*3);m.f.triangle.resize(h[7]);m.f.coefficient.resize(h[8]);m.global.resize(h[8]);m.owner.resize(h[8]);m.rule.resize(h[7]);m.external_edge.resize(h[3]);
    read(f,m.global.data(),m.global.size());read(f,m.owner.data(),m.owner.size());read(f,m.f.cell.data(),m.f.cell.size());read(f,m.source_triangle.data(),m.source_triangle.size());read(f,m.source_cell.data(),m.source_cell.size());read(f,m.triangle_edge.data(),m.triangle_edge.size());read(f,m.f.triangle.data(),m.f.triangle.size());read(f,m.f.coefficient.data(),m.f.coefficient.size());read(f,m.rule.data(),m.rule.size());read(f,m.external_edge.data(),m.external_edge.size());need(fgetc(f)==EOF&&fclose(f)==0,"complete owned full field checkpoint");
    m.f.dof.resize(m.global.size());uint64_t npoly=m.global_vertices+2*m.global_edges+m.global_triangles,bubble=m.global_dofs-2*m.global_triangles;
    for(size_t i=0;i<m.global.size();i++){uint64_t g=m.global[i];need(g<m.global_dofs&&(i==0||m.global[i-1]<g),"saved unique sorted global DOF IDs");m.f.dof[i]={g,g,g<m.global_vertices?0:g<m.global_vertices+2*m.global_edges?1:g<npoly?2:g<bubble?3:4,m.owner[i]};}
    return m;
}
}
#endif
