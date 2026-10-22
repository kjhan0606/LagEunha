#ifndef MATERIAL_COMPACT_P3_H
#define MATERIAL_COMPACT_P3_H
// Compact rezoning space: no initial-condition function or residual basis.
// Shared native coefficient ABI and original barycentric P3 convention.
#include "material_shared_trace_batch.h"
#include <cmath>
#include <cstdint>
#include <vector>
#include <array>
#include <map>
#include <set>
#include <stdexcept>
#include <limits>
#include <algorithm>
#ifdef __CUDACC__
#define COMPACT_P3_HD __host__ __device__
#else
#define COMPACT_P3_HD
#endif

namespace compact_p3 {
using R=long double;
using Point=std::array<double,2>;
using Edge=std::pair<uint64_t,uint64_t>;
inline void require(bool b,const char*s){if(!b)throw std::runtime_error(s);}
struct Triangle {
    double vertex[3][2]={},grad[3][2]={},area=0;
    uint64_t id[10]={},cell=0;
};
struct ThermalDensity {
    // Reference-volume densities. They are independent physical fields.
    // Nonnegative coefficients are candidates, not a completed positivity proof.
    double mass[10]={},internal_energy[10]={};
};
COMPACT_P3_HD static inline void shape(const double*l,const double grad[3][2],STJet*N){
    STJet L[3];for(int k=0;k<3;k++)L[k]={l[k],grad[k][0],grad[k][1]};
    for(int a=0;a<10;a++){
        int index[3];gfs_material_p3_index(a,index);N[a]=STJet(1);
        for(int k=0;k<3;k++)for(int j=0;j<index[k];j++)
            N[a]=N[a]*(STJet(3)*L[k]-STJet(j))/STJet(index[k]-j);
    }
}
COMPACT_P3_HD static inline void bernstein(const double*l,double*B){
    const int factorial[4]={1,1,2,6};
    for(int a=0;a<10;a++){
        int index[3];gfs_material_p3_index(a,index);
        B[a]=6.0/(factorial[index[0]]*factorial[index[1]]*factorial[index[2]]);
        for(int k=0;k<3;k++)for(int j=0;j<index[k];j++)B[a]*=l[k];
    }
}
COMPACT_P3_HD static inline double field(const double*control,const double*B){
    double out=0;for(int a=0;a<10;a++)out+=control[a]*B[a];return out;
}
COMPACT_P3_HD static inline int point(const Triangle&t,const STCoefficient*c,
                              const double*l,GfsMaterialPoint12*out){
    STJet N[10];shape(l,t.grad,N);*out={};double V[4]={};
    const auto&origin=c[t.id[0]];
    // Translation cancellation is identical in meaning to the original kernel.
    for(int a=0;a<10;a++)for(int k=0;k<2;k++){
        const auto&q=c[t.id[a]];double x=q.x[k]-origin.x[k],v=q.v[k]-origin.v[k];
        out->position[k]+=N[a].v*x;out->velocity[k]+=N[a].v*v;
        out->F[2*k]+=N[a].x*x;out->F[2*k+1]+=N[a].y*x;
        V[2*k]+=N[a].x*v;V[2*k+1]+=N[a].y*v;
    }
    for(int k=0;k<2;k++){out->position[k]+=origin.x[k];out->velocity[k]+=origin.v[k];}
    out->J=out->F[0]*out->F[3]-out->F[1]*out->F[2];
    if(!(out->J>0)||!isfinite(out->J))return 0;
    double A[4]={out->F[3]/out->J,-out->F[1]/out->J,-out->F[2]/out->J,out->F[0]/out->J};
    for(int k=0;k<2;k++){
        out->grad_v[2*k]=V[2*k]*A[0]+V[2*k+1]*A[2];
        out->grad_v[2*k+1]=V[2*k]*A[1]+V[2*k+1]*A[3];
        if(!isfinite(out->position[k])||!isfinite(out->velocity[k]))return 0;
    }
    for(int k=0;k<4;k++)if(!isfinite(out->F[k])||!isfinite(out->grad_v[k]))return 0;
    return 1;
}
inline void geometry(Triangle&t){
    R ux=R(t.vertex[1][0])-t.vertex[0][0],uy=R(t.vertex[1][1])-t.vertex[0][1];
    R vx=R(t.vertex[2][0])-t.vertex[0][0],vy=R(t.vertex[2][1])-t.vertex[0][1];
    R det=ux*vy-uy*vx;require(std::isfinite(det)&&det>0,"positive compact straight reference fan");
    t.area=double(det/2);require(t.area>0&&std::isfinite(t.area),"representable compact reference area");
    t.grad[1][0]=double(vy/det);t.grad[1][1]=double(-vx/det);
    t.grad[2][0]=double(-uy/det);t.grad[2][1]=double(ux/det);
    for(int k=0;k<2;k++)t.grad[0][k]=-t.grad[1][k]-t.grad[2][k];
}
struct QuadraturePoint {double l[3],weight;};
// Positive 6x5 Duffy rule. Total-degree9 becomes degree<=10 in u and<=9 in v
// after the Jacobian u: Gauss6/Gauss5 integrate it exactly in real arithmetic.
// This establishes polynomial mass exactness only, not EOS/gravity accuracy.
inline std::array<QuadraturePoint,30> mass_rule(const Triangle&t){
    static constexpr double u[6]={.033765242898423986,.16939530676686774,.38069040695840155,.61930959304159845,.83060469323313226,.96623475710157601};
    static constexpr double w[6]={.085662246189585173,.18038078652406930,.23395696728634552,.23395696728634552,.18038078652406930,.085662246189585173};
    static constexpr double v[5]={.046910077030668004,.23076534494715845,.5,.76923465505284155,.95308992296933200};
    static constexpr double z[5]={.11846344252809454,.23931433524968323,.28444444444444444,.23931433524968323,.11846344252809454};
    std::array<QuadraturePoint,30>q;int at=0;
    for(int i=0;i<6;i++)for(int j=0;j<5;j++)q[at++]={{1-u[i],u[i]*(1-v[j]),u[i]*v[j]},2*t.area*u[i]*w[i]*z[j]};
    return q;
}
inline bool nonnegative_controls(const double*c){
    R sum=0;for(int a=0;a<10;a++){if(!std::isfinite(c[a])||c[a]<0)return false;sum+=c[a];}return sum>0;
}
inline std::array<R,100> mass(const Triangle&t,const ThermalDensity&f){
    require(nonnegative_controls(f.mass),"candidate Bernstein density without negative clamping");
    std::array<R,100>M={};
    for(const auto&q:mass_rule(t)){
        double B[10];STJet N[10];bernstein(q.l,B);shape(q.l,t.grad,N);
        R rho=0;for(int a=0;a<10;a++)rho+=R(f.mass[a])*B[a];
        require(rho>0&&q.weight>0,"positive polynomial mass quadrature");
        for(int a=0;a<10;a++)for(int b=0;b<=a;b++)M[a*10+b]+=R(q.weight)*rho*R(N[a].v)*N[b].v;
    }
    for(int a=0;a<10;a++)for(int b=0;b<a;b++)M[b*10+a]=M[a*10+b];return M;
}
struct Polygon {
    uint64_t cell=0;Point centroid={};std::vector<uint64_t>vertex;
};
struct Mesh {
    std::vector<Point>node;std::vector<Triangle>triangle;
    std::vector<Edge>edge;std::vector<uint64_t>edge_incidence;
    uint64_t reference_vertices=0,external_edges=0;
};
inline Mesh topology(const std::vector<Point>&vertices,const std::vector<Polygon>&polygons){
    Mesh out;out.node=vertices;std::set<uint64_t>cells;
    // Input registry must come from canonical CURRENT power-cell geometry.
    // No tolerance merge, old-cell subdivision, or analytic-orbit substitution.
    for(const auto&p:polygons){
        require(p.vertex.size()>=3&&cells.insert(p.cell).second,"unique compact physical target polygon");
        require(std::isfinite(p.centroid[0])&&std::isfinite(p.centroid[1]),"actual target area centroid");
        uint64_t center=out.node.size();out.node.push_back(p.centroid);
        std::set<uint64_t>seen;
        for(auto id:p.vertex)require(id<vertices.size()&&seen.insert(id).second,"canonical unique polygon vertices");
        for(size_t j=0;j<p.vertex.size();j++){
            uint64_t vertex[3]={center,p.vertex[j],p.vertex[(j+1)%p.vertex.size()]};Triangle t;t.cell=p.cell;
            for(int k=0;k<3;k++)for(int d=0;d<2;d++)t.vertex[k][d]=out.node[vertex[k]][d];
            geometry(t);
            for(int a=0;a<10;a++){int index[3];gfs_material_p3_index(a,index);for(int k=0;k<3;k++)if(index[k]==3)t.id[a]=vertex[k];}
            out.triangle.push_back(t);
        }
    }
    out.reference_vertices=out.node.size();std::map<Edge,std::array<uint64_t,3>>edge;
    for(auto&t:out.triangle){
        uint64_t vertex[3]={};for(int a=0;a<10;a++){int index[3];gfs_material_p3_index(a,index);for(int k=0;k<3;k++)if(index[k]==3)vertex[k]=t.id[a];}
        for(int a=0;a<10;a++){
            int index[3];gfs_material_p3_index(a,index);int zero=-1;bool isvertex=false;
            for(int k=0;k<3;k++){if(index[k]==0)zero=k;isvertex|=index[k]==3;}if(isvertex)continue;
            if(zero>=0){
                int i=(zero+1)%3,j=(zero+2)%3;uint64_t lo=vertex[i],hi=vertex[j];int part=index[j];if(lo>hi){std::swap(lo,hi);part=3-part;}
                auto key=Edge{lo,hi};auto it=edge.find(key);
                if(it==edge.end()){
                    uint64_t first=out.node.size();for(int k=1;k<=2;k++){Point x;for(int d=0;d<2;d++)x[d]=double((R(3-k)*out.node[lo][d]+R(k)*out.node[hi][d])/3);out.node.push_back(x);}
                    it=edge.emplace(key,std::array<uint64_t,3>{first,first+1,0}).first;
                }
                require(part==1||part==2,"canonical compact edge parameter orientation");t.id[a]=it->second[part-1];
            }else{Point x;for(int d=0;d<2;d++)x[d]=double((R(t.vertex[0][d])+t.vertex[1][d]+t.vertex[2][d])/3);t.id[a]=out.node.size();out.node.push_back(x);}
        }
        for(int k=0;k<3;k++){auto p=std::minmax(vertex[k],vertex[(k+1)%3]);edge.at(Edge{p.first,p.second})[2]++;}
    }
    for(const auto&e:edge){require(e.second[2]==1||e.second[2]==2,"compact manifold edge incidence");out.edge.push_back(e.first);out.edge_incidence.push_back(e.second[2]);out.external_edges+=e.second[2]==1;}
    require(!out.triangle.empty()&&out.reference_vertices+out.triangle.size()==out.edge.size()+1,"compact connected disk Euler incidence prerequisite");
    require(out.node.size()==out.reference_vertices+2*out.edge.size()+out.triangle.size(),"compact shared P3 scalar extent");
    return out;
}
} // namespace compact_p3
#undef COMPACT_P3_HD
#endif
