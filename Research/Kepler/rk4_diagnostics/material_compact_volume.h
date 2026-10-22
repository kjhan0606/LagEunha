#ifndef MATERIAL_COMPACT_VOLUME_H
#define MATERIAL_COMPACT_VOLUME_H
// Canonical finite boundary-area functional and its analytic adjoint. Uses
// the completed geometry space, including its dynamic component modes.
// Topology, canonical event identity and curved quadrature accuracy are
// prerequisites of the complete transaction, not proved by this adapter.
#include "material_compact_geometry.h"
#include <map>
#include <vector>
#include <array>
namespace compact_volume {
using R=long double;using V2=std::array<R,2>;
inline V2 plus(V2 a,V2 b){return {a[0]+b[0],a[1]+b[1]};}
inline V2 minus(V2 a,V2 b){return {a[0]-b[0],a[1]-b[1]};}
inline V2 times(R a,V2 b){return {a*b[0],a*b[1]};}
inline V2 cw(V2 a){return {a[1],-a[0]};}
inline R cross(V2 a,V2 b){return a[0]*b[1]-a[1]*b[0];}
inline R dot(V2 a,V2 b){return a[0]*b[0]+a[1]*b[1];}
using compact_p3::require;
struct Rule {std::vector<double>point,weight;};
struct Edge {
    uint64_t event=UINT64_MAX,id[6]={};
    int n=4;
    compact_geometry::FrozenTrace trace;
    // node stores the SAME small nodal collar extension as the patch, not
    // absolute Gamma coordinates. chord endpoints keep their canonical order.
    double parameter[2]={},node[4][2]={},chord[2][2]={},scale[2]={1,1};
};
struct Side {uint64_t cell=0;int sign=1;};
struct CanonicalEdge {
    Edge basis;std::vector<Side>side;
    // A one-sided LOCAL edge may meet a remote MPI cell. Exterior status
    // comes from the actual domain registry, never local incidence alone.
    bool physical_exterior=false;
};
struct Value {
    R volume=0,volume_dot=0;
    V2 G[6]={},H[6]={},cofactor[6]={};
    R derivative_error=0,derivative_scale=0;
};
inline Edge from_patch(const compact_geometry::Patch&p,int local_edge,bool forward,uint64_t event){
    require(local_edge>=0&&local_edge<3&&event!=UINT64_MAX,"canonical actual target edge event");
    int a=local_edge,b=(local_edge+1)%3;if(!forward)std::swap(a,b);
    Edge out;out.event=event;int polynomial_index[4]={};
    for(int j=0;j<4;j++){
        int selected=-1;for(int k=0;k<10;k++){
            int index[3];gfs_material_p3_index(k,index);
            if(index[a]==3-j&&index[b]==j){bool match=true;for(int z=0;z<3;z++)if(z!=a&&z!=b)match&=index[z]==0;if(match)selected=k;}
        }
        require(selected>=0,"same canonical P3 trace node ordering");out.id[j]=p.polynomial.id[selected];polynomial_index[j]=selected;
    }
    if(p.curved&&local_edge==1){
        out.n=6;out.trace=p.trace;
        for(int k=0;k<2;k++){out.id[4+k]=p.modes.id[k];out.scale[k]=p.modes.scale[k];}
        out.parameter[0]=p.parameter[forward?0:1];out.parameter[1]=p.parameter[forward?1:0];
        for(int j=0;j<2;j++)for(int k=0;k<2;k++)out.chord[j][k]=p.polynomial.vertex[forward?1+j:2-j][k];
        for(int j=0;j<4;j++)for(int k=0;k<2;k++)out.node[j][k]=p.interpolation[polynomial_index[j]][k];
    }
    return out;
}
inline void shapes(const Edge&e,double q,R*beta,R*dbeta,R*N){
    require(q>=0&&q<=1&&std::isfinite(q)&&(e.n==4||e.n==6),"valid canonical edge sample");
    STJet t(q,1,0),L[4];for(int j=0;j<4;j++){L[j]=st_lagrange(j,t);N[j]=L[j].v;}
    beta[1]=L[1].v;beta[2]=L[2].v;beta[0]=-(2*beta[1]+beta[2])/3;beta[3]=-(beta[1]+2*beta[2])/3;
    dbeta[1]=L[1].x;dbeta[2]=L[2].x;dbeta[0]=-(2*dbeta[1]+dbeta[2])/3;dbeta[3]=-(dbeta[1]+2*dbeta[2])/3;
    if(e.n==6){
        STJet relative[2];require(compact_geometry::relative_trace(e.trace,t,e.parameter[0],e.parameter[1],e.chord[0],e.chord[1],relative),"stable accepted exterior relative trace and derivative");
        for(int k=0;k<2;k++){
            require(e.scale[k]>0&&std::isfinite(e.scale[k]),"finite actual component-mode normalization");
            R h=relative[k].v,dh=relative[k].x;
            for(int j=0;j<4;j++){h-=R(L[j].v)*e.node[j][k];dh-=R(L[j].x)*e.node[j][k];}
            beta[4+k]=N[4+k]=h/e.scale[k];dbeta[4+k]=dh/e.scale[k];
        }
    }
}
inline Value evaluate(const Edge&e,const std::vector<STCoefficient>&c,const Rule&rule){
    require((e.n==4||e.n==6)&&rule.point.size()==rule.weight.size()&&!rule.point.empty(),"explicit positive canonical boundary functional");
    Value out;V2 x[6],v[6];for(int a=0;a<e.n;a++){
        require(e.id[a]<c.size(),"same local coefficient table for area/force");
        x[a]={c[e.id[a]].x[0],c[e.id[a]].x[1]};v[a]={c[e.id[a]].v[0],c[e.id[a]].v[1]};
        for(int k=0;k<2;k++)require(std::isfinite(x[a][k])&&std::isfinite(v[a][k]),"finite dynamic boundary coefficients");
    }
    V2 A=x[0],B=x[3],D=minus(B,A);
    out.volume=cross(A,B)/2;out.G[0]=times(.5L,cw(B));out.G[3]=times(-.5L,cw(A));
    V2 curve1=minus(minus(x[1],A),times(1.L/3,D)),curve2=minus(minus(x[2],A),times(2.L/3,D));
    for(size_t q=0;q<rule.point.size();q++){
        R w=rule.weight[q];require(w>0&&std::isfinite(w),"positive actual boundary quadrature weight");
        R beta[6]={},dbeta[6]={},N[6]={};shapes(e,rule.point[q],beta,dbeta,N);
        V2 delta=plus(times(beta[1],curve1),times(beta[2],curve2));
        V2 prime=plus(times(dbeta[1],curve1),times(dbeta[2],curve2));
        for(int a=4;a<e.n;a++){delta=plus(delta,times(beta[a],x[a]));prime=plus(prime,times(dbeta[a],x[a]));}
        out.volume+=w*(cross(delta,D)+cross(delta,prime)/2);
        for(int a=0;a<e.n;a++){
            R endpoint=(a==3?1:0)-(a==0?1:0);
            V2 g=minus(times(beta[a],cw(D)),times(endpoint,cw(delta)));
            g=plus(g,minus(times(beta[a]/2,cw(prime)),times(dbeta[a]/2,cw(delta))));
            out.G[a]=plus(out.G[a],times(w,g));out.cofactor[a]=plus(out.cofactor[a],times(w*N[a],cw(plus(D,prime))));
        }
    }
    for(int a=0;a<e.n;a++)out.H[a]=out.G[a];
    // Same open-edge endpoint gauge as the original native operator. It
    // cancels on an exactly closed canonical cell; H is the physical face
    // force adjoint, while G differentiates the open-edge area functional.
    out.H[0]=minus(out.H[0],times(.5L,cw(A)));out.H[3]=plus(out.H[3],times(.5L,cw(B)));
    for(int a=0;a<e.n;a++){
        out.volume_dot+=dot(out.G[a],v[a]);
        for(int k=0;k<2;k++){
            out.derivative_error=std::max(out.derivative_error,fabsl(out.H[a][k]-out.cofactor[a][k]));
            out.derivative_scale=std::max(out.derivative_scale,fabsl(out.cofactor[a][k]));
        }
    }
    require(std::isfinite(out.volume)&&std::isfinite(out.volume_dot),"finite canonical curved area functional");return out;
}
struct Cell {
    R volume=0,volume_dot=0,translation_error=0,rotation_error=0;
    std::map<uint64_t,V2>G;
};
struct Geometry {
    std::vector<Value>edge;
    std::vector<Cell>cell;
    R derivative_error=0,derivative_scale=0;
    uint64_t paired=0,external=0,partition_boundary=0;
};
inline Geometry assemble(const std::vector<CanonicalEdge>&edges,const std::vector<STCoefficient>&c,
    size_t cells,const Rule&straight,const Rule&curved,const std::vector<uint8_t>&constant){
    require(constant.size()==c.size(),"explicit constant coefficients: P3 ones and geometry-mode zeros");
    Geometry out;out.edge.resize(edges.size());out.cell.resize(cells);std::set<uint64_t>events;
    for(size_t k=0;k<edges.size();k++){
        const auto&e=edges[k];require(events.insert(e.basis.event).second,"one canonical face evaluation per actual event");
        require(e.side.size()==1||e.side.size()==2,"actual target boundary manifold incidence");
        if(e.side.size()==2){require(!e.physical_exterior&&e.side[0].cell!=e.side[1].cell&&e.side[0].sign==-e.side[1].sign,"paired actual cell faces have opposite orientation");out.paired++;}
        else if(e.physical_exterior)out.external++;else out.partition_boundary++;
        auto&v=out.edge[k];v=evaluate(e.basis,c,e.basis.n==6?curved:straight);
        out.derivative_error=std::max(out.derivative_error,v.derivative_error);out.derivative_scale=std::max(out.derivative_scale,v.derivative_scale);
        for(const auto&side:e.side){
            require(side.cell<cells&&(side.sign==1||side.sign==-1),"actual compact cell and canonical orientation");auto&cell=out.cell[side.cell];
            cell.volume+=side.sign*v.volume;
            for(int a=0;a<e.basis.n;a++)cell.G[e.basis.id[a]]=plus(cell.G[e.basis.id[a]],times(side.sign,v.H[a]));
        }
    }
    for(auto&cell:out.cell){
        V2 translation={};R rotation=0;
        for(const auto&entry:cell.G){
            auto id=entry.first;require(id<c.size()&&(constant[id]==0||constant[id]==1),"same area and dynamic coefficient indices");
            if(constant[id])translation=plus(translation,entry.second);
            rotation+=cross({c[id].x[0],c[id].x[1]},entry.second);cell.volume_dot+=dot({c[id].v[0],c[id].v[1]},entry.second);
        }
        cell.translation_error=hypotl(translation[0],translation[1]);cell.rotation_error=fabsl(rotation);
        require(cell.volume>0&&std::isfinite(cell.volume)&&std::isfinite(cell.volume_dot),"positive actual cell boundary area");
    }
    return out;
}
// Shared mean-pressure forces, including the EXACT same finite boundary
// functional for external work. Difference from integrated physical traction
// is retained in Geometry.derivative_error, never silently exempted.
inline std::vector<V2> mean_pressure_force(const std::vector<CanonicalEdge>&edges,const Geometry&g,
    const std::vector<R>&mean_pressure,size_t coefficients,R external_pressure){
    require(edges.size()==g.edge.size()&&mean_pressure.size()==g.cell.size()&&external_pressure>0,"same stage volume/pressure/exterior law");
    std::vector<V2>force(coefficients);
    for(size_t k=0;k<edges.size();k++){
        const auto&e=edges[k];R p=0;for(const auto&s:e.side){require(s.cell<mean_pressure.size()&&mean_pressure[s.cell]>0&&std::isfinite(mean_pressure[s.cell]),"positive stage mean pressure");p+=s.sign*mean_pressure[s.cell];}
        if(e.physical_exterior){require(e.side.size()==1,"single physical exterior incidence");p-=e.side[0].sign*external_pressure;}
        for(int a=0;a<e.basis.n;a++){auto id=e.basis.id[a];require(id<force.size(),"same mean-pressure force DOF");force[id]=plus(force[id],times(p,g.edge[k].H[a]));}
    }
    return force;
}
} // namespace compact_volume
#endif
