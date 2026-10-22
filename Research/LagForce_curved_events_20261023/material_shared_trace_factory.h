#ifndef MATERIAL_SHARED_TRACE_FACTORY_H
#define MATERIAL_SHARED_TRACE_FACTORY_H
#include "material_shared_trace_batch.h"
#include <map>
#include <array>
#include <vector>
#include <algorithm>
#include <limits>
namespace shared_trace {
using R=long double;using Pt=std::array<double,2>;using Key=std::pair<int,int>;
struct Dof {uint64_t id,parent;int kind,owner;};
struct Factory {
    std::vector<Pt> vertex,polynomial;
    std::vector<STEdge> edge;
    std::map<Key,int> edge_index;
    std::vector<STTriangle> triangle;
    std::vector<Dof> dof;
    std::vector<STCoefficient> coefficient;
    std::vector<uint64_t> cell;
    size_t rank1=0,rank2=0;double min_edge_scale=DBL_MAX,min_bubble_scale=DBL_MAX;
};
static void gauss(int n,std::vector<double>&x,std::vector<double>&w){x.resize(n);w.resize(n);for(int i=0;i<(n+1)/2;i++){R z=cosl(acosl(-1.L)*(i+.75L)/(n+.5L)),d=0;for(int it=0;it<64;it++){R a=1,b=z;for(int k=2;k<=n;k++){R c=((2*k-1)*z*b-(k-1)*a)/k;a=b;b=c;}d=n*(z*b-a)/(z*z-1);R v=z-b/d;if(fabsl(v-z)<4*LDBL_EPSILON){z=v;break;}z=v;}R a=1,b=z;for(int k=2;k<=n;k++){R c=((2*k-1)*z*b-(k-1)*a)/k;a=b;b=c;}d=n*(z*b-a)/(z*z-1);x[i]=(1-z)/2;x[n-1-i]=(1+z)/2;w[i]=w[n-1-i]=1/((1-z*z)*d*d);}}
static Pt initial_velocity(Pt x,double eps){double omega=pow(x[0]*x[0]+x[1]*x[1]+eps*eps,-.75);return {-x[1]*omega,x[0]*omega};}
static void edge_basis(STEdge&e,Pt A,Pt B,double eps){
    e.eps=eps;for(int k=0;k<2;k++){e.A[k]=A[k];e.D[k]=B[k]-A[k];}
    for(int j=0;j<4;j++){double X=e.A[0]+(j/3.)*e.D[0],Y=e.A[1]+(j/3.)*e.D[1],Dq[4];need(gfs_material_p3_radial_remainder(e.A,X,Y,eps,e.qnode[j],Dq),"canonical edge original IC remainder");}
    // Exact algebraic rank, no sampled rank truncation. For a nonradial
    // straight edge the functions (a+b t)(r(t)^2+eps^2)^(-3/4) modulo P3
    // are independent. A radial edge has one scalar residual direction.
    __float128 det=__float128(A[0])*__float128(B[1])-__float128(A[1])*__float128(B[0]);e.rank=det==0?1:2;
    std::vector<double>x,w;gauss(64,x,w);R G[2][2]={};for(int j=0;j<64;j++){STJet r[2];need(st_edge_raw(e,STJet(x[j]),r),"finite canonical residual");for(int a=0;a<2;a++)for(int b=0;b<2;b++)G[a][b]+=R(w[j])*r[a].v*r[b].v;}
    int p=G[0][0]>=G[1][1]?0:1,o=1-p;need(G[p][p]>0,"nonzero canonical residual norm");double s0=sqrtl(G[p][p]);e.encode[0][p]=1/s0;e.decode[p][0]=s0;
    if(e.rank==1){Pt V=fabs(A[0])+fabs(A[1])>fabs(B[0])+fabs(B[1])?A:B;double component[2]={-V[1],V[0]};need(component[p]!=0,"radial canonical pivot");e.decode[o][0]=s0*(component[o]/component[p]);}
    else {double beta=G[p][o]/G[p][p];R variance=0;for(int j=0;j<64;j++){STJet r[2];need(st_edge_raw(e,STJet(x[j]),r),"canonical second residual");R z=R(r[o].v)-R(beta)*r[p].v;variance+=R(w[j])*z*z;}need(variance>0,"independent canonical second residual norm");double s1=sqrtl(variance);e.encode[1][p]=-beta/s1;e.encode[1][o]=1/s1;e.decode[o][0]=beta*s0;e.decode[o][1]=s1;}
}
// Both retained CONTROL and fresh full-domain geometry provide this view.
// Only cell IDs, reference vertex coordinates and eps are read here.
template<class Packet> static Factory topology(const Packet&p,int ranks){
    Factory f;std::map<Pt,int> vertex;
    for(const auto&r:p.record){f.cell.push_back(r.cell);for(int a=0;a<3;a++)vertex[{r.triangle.geometry.x[a][0],r.triangle.geometry.x[a][1]}]=0;}
    std::sort(f.cell.begin(),f.cell.end());f.cell.erase(std::unique(f.cell.begin(),f.cell.end()),f.cell.end());
    for(auto&v:vertex){v.second=f.vertex.size();f.vertex.push_back(v.first);}
    for(const auto&r:p.record)for(int a=0;a<3;a++){int x=vertex.at({r.triangle.geometry.x[a][0],r.triangle.geometry.x[a][1]}),y=vertex.at({r.triangle.geometry.x[(a+1)%3][0],r.triangle.geometry.x[(a+1)%3][1]});f.edge_index[std::minmax(x,y)]=0;}
    f.polynomial=f.vertex;
    for(auto&entry:f.edge_index){entry.second=f.edge.size();STEdge edge={};Pt A=f.vertex[entry.first.first],B=f.vertex[entry.first.second];edge_basis(edge,A,B,p.parameters[3]);edge.poly[0]=entry.first.first;edge.poly[3]=entry.first.second;
        for(int j=1;j<=2;j++){edge.poly[j]=f.polynomial.size();Pt node;for(int k=0;k<2;k++)node[k]=double(((R)(3-j)*A[k]+R(j)*B[k])/3);f.polynomial.push_back(node);}f.edge.push_back(edge);
    }
    for(const auto&r:p.record){Pt node;for(int k=0;k<2;k++)node[k]=double((R(r.triangle.geometry.x[0][k])+r.triangle.geometry.x[1][k]+r.triangle.geometry.x[2][k])/3);f.polynomial.push_back(node);}
    auto add=[&](int kind,uint64_t parent){uint64_t id=f.dof.size();f.dof.push_back({id,parent,kind,int(id%ranks)});return int(id);};
    for(size_t i=0;i<f.polynomial.size();i++)add(i<f.vertex.size()?0:i<f.vertex.size()+2*f.edge.size()?1:2,i);
    for(size_t e=0;e<f.edge.size();e++){auto&edge=f.edge[e];for(int j=0;j<edge.rank;j++)edge.id[j]=add(3,e);f.rank1+=edge.rank==1;f.rank2+=edge.rank==2;for(int j=0;j<edge.rank;j++){double scale=0;for(int k=0;k<2;k++)scale=hypot(scale,edge.decode[k][j]);f.min_edge_scale=std::min(f.min_edge_scale,scale);}}
    f.triangle.resize(p.record.size());
    for(size_t s=0;s<p.record.size();s++){
        const auto&r=p.record[s];auto&t=f.triangle[s];t.eps=p.parameters[3];t.cell=std::lower_bound(f.cell.begin(),f.cell.end(),r.cell)-f.cell.begin();GfsMaterialTriangle g={};int vi[3];for(int a=0;a<3;a++){for(int k=0;k<2;k++)g.x[a][k]=t.vertex[a][k]=r.triangle.geometry.x[a][k];vi[a]=vertex.at({g.x[a][0],g.x[a][1]});}need(gfs_material_triangle_init(&g),"same positive original reference geometry");for(int a=0;a<3;a++)for(int k=0;k<2;k++)t.grad[a][k]=g.grad[a][k];
        for(int a=0;a<10;a++){int m[3];gfs_material_p3_index(a,m);int v=-1,z=-1;for(int k=0;k<3;k++){if(m[k]==3)v=k;if(m[k]==0)z=k;}if(v>=0)t.id[a]=vi[v];else if(z>=0){int i=(z+1)%3,j=(z+2)%3;auto&e=f.edge[f.edge_index.at(std::minmax(vi[i],vi[j]))];int part=vi[i]<vi[j]?m[j]:m[i];t.id[a]=e.poly[part];}else t.id[a]=f.vertex.size()+2*f.edge.size()+s;double Dq[4];auto x=f.polynomial[t.id[a]];need(gfs_material_p3_radial_remainder(t.vertex[0],x[0],x[1],t.eps,t.qnode[a],Dq),"fresh local IC remainder cache");}
        t.n=10;for(int a=0;a<3;a++){int b=(a+1)%3;auto key=std::minmax(vi[a],vi[b]);t.edge[a]=f.edge[f.edge_index.at(key)];t.endpoint[a][vi[a]<vi[b]?0:1]=a;t.endpoint[a][vi[a]<vi[b]?1:0]=b;t.offset[a]=t.n;for(int j=0;j<t.edge[a].rank;j++)t.id[t.n++]=t.edge[a].id[j];}
        t.bubble_offset=t.n;for(int j=0;j<2;j++){t.id[t.n++]=add(4,s);t.bubble_encode[j][j]=t.bubble_decode[j][j]=1;}
    }
    f.coefficient.resize(f.dof.size());return f;
}
static void bubble_basis(STTriangle&t,const R*matrix,double mass){
    int b=t.bubble_offset;R G0=matrix[9*ST_MAX+9];need(G0>0&&mass>0,"positive local bubble mass metric");
    for(int j=0;j<2;j++)t.bubble_shift[j]=matrix[9*ST_MAX+b+j]/G0;
    R H[2][2];for(int i=0;i<2;i++)for(int j=0;j<2;j++)H[i][j]=matrix[(b+i)*ST_MAX+b+j]-R(t.bubble_shift[i])*matrix[9*ST_MAX+b+j]-R(t.bubble_shift[j])*matrix[(b+i)*ST_MAX+9]+R(t.bubble_shift[i])*t.bubble_shift[j]*G0;
    int p=H[0][0]>=H[1][1]?0:1,o=1-p;need(H[p][p]>0,"positive interior residual leading norm");double s0=sqrtl(H[p][p]/mass),beta=H[p][o]/H[p][p];R variance=H[o][o]-2*R(beta)*H[p][o]+R(beta)*beta*H[p][p];need(variance>0,"independent interior residual second norm");double s1=sqrtl(variance/mass);
    for(int i=0;i<2;i++)for(int j=0;j<2;j++)t.bubble_encode[i][j]=t.bubble_decode[i][j]=0;
    t.bubble_encode[0][p]=1/s0;t.bubble_encode[1][p]=-beta/s1;t.bubble_encode[1][o]=1/s1;t.bubble_decode[p][0]=s0;t.bubble_decode[o][0]=beta*s0;t.bubble_decode[o][1]=s1;
}
static void initialise(Factory&f,double eps,int rank,int ranks){
    std::vector<STCoefficient> generated(f.coefficient.size());
    for(size_t i=0;i<f.polynomial.size();i++){auto v=initial_velocity(f.polynomial[i],eps);for(int k=0;k<2;k++){generated[i].x[k]=f.polynomial[i][k];generated[i].v[k]=v[k];}}
    for(const auto&e:f.edge)for(int j=0;j<e.rank;j++)for(int k=0;k<2;k++)generated[e.id[j]].v[k]=e.decode[k][j];
    for(const auto&t:f.triangle){for(int k=0;k<2;k++)generated[t.id[9]].v[k]+=t.bubble_shift[k];for(int j=0;j<2;j++)for(int k=0;k<2;k++)generated[t.id[t.bubble_offset+j]].v[k]=t.bubble_decode[k][j];}
    for(size_t i=0;i<generated.size();i++)if(f.dof[i].owner==rank)f.coefficient[i]=generated[i];
    static_assert(sizeof(STCoefficient)==4*sizeof(double),"native shared coefficient wire layout");
    MPI_Allreduce(MPI_IN_PLACE,f.coefficient.data(),f.coefficient.size()*4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
}
}
#endif
