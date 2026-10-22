#ifndef MATERIAL_GEOMETRIC_PRESSURE_GEOMETRY_H
#define MATERIAL_GEOMETRIC_PRESSURE_GEOMETRY_H
namespace geometric_pressure {
using Vec=shared_trace_rhs::Vec;using V2=std::array<R,2>;
static V2 operator+(V2 a,V2 b){return {a[0]+b[0],a[1]+b[1]};}
static V2 operator-(V2 a,V2 b){return {a[0]-b[0],a[1]-b[1]};}
static V2 operator*(R a,V2 b){return {a*b[0],a*b[1]};}
static R cross(V2 a,V2 b){return a[0]*b[1]-a[1]*b[0];}
static R inner(V2 a,V2 b){return a[0]*b[0]+a[1]*b[1];}
static V2 cw(V2 a){return {a[1],-a[0]};}
struct Side {int triangle,local,sign,cell;};
struct Edge {STEdge reference;std::vector<Side>side;int id[6]={},n=0;R volume=0,vdot=0;V2 G[6]={},H[6]={},cofactor[6]={};};
struct Geometry {std::vector<Edge>edge;std::vector<std::map<int,V2>> G;Vec volume,vdot;R derivative_error=0,derivative_scale=0,translation=0,rotation=0;size_t internal_material=0,paired=0,external=0;};
// Canonical endpoint-subtracted residual coefficients; evaluate the SAME
// chosen finite quadrature functional and differentiate it analytically.
static void edge_shapes(const STEdge&e,double q,R*beta,R*dbeta,R*N){
    STJet t(q,1,0),L[4],raw[2];for(int j=0;j<4;j++){L[j]=st_lagrange(j,t);N[j]=L[j].v;}
    beta[1]=L[1].v;beta[2]=L[2].v;beta[0]=-(2*beta[1]+beta[2])/3;beta[3]=-(beta[1]+2*beta[2])/3;
    dbeta[1]=L[1].x;dbeta[2]=L[2].x;dbeta[0]=-(2*dbeta[1]+dbeta[2])/3;dbeta[3]=-(dbeta[1]+2*dbeta[2])/3;
    need(st_edge_raw(e,t,raw),"canonical original edge functions and derivatives");
    for(int j=0;j<e.rank;j++){R value=0,derivative=0;for(int k=0;k<2;k++){value+=R(e.encode[j][k])*raw[k].v;derivative+=R(e.encode[j][k])*raw[k].x;}beta[j+4]=N[j+4]=value;dbeta[j+4]=derivative;}
}
static Geometry geometry(const st::Factory&f){
    Geometry out;out.G.resize(f.cell.size());out.volume.resize(f.cell.size());out.vdot.resize(f.cell.size());std::map<int,Edge> edges;
    for(size_t s=0;s<f.triangle.size();s++)for(int a=0;a<3;a++){const auto&t=f.triangle[s];auto&e=edges[t.edge[a].poly[1]];e.reference=t.edge[a];e.side.push_back({int(s),a,t.endpoint[a][0]==a?1:-1,t.cell});}
    std::vector<double>x,w;st::gauss(64,x,w);
    for(auto&entry:edges){auto e=entry.second;need(e.side.size()==1||e.side.size()==2,"canonical manifold reference edge");if(e.side.size()==2&&e.side[0].cell==e.side[1].cell){out.internal_material++;continue;}if(e.side.size()==1)out.external++;else{out.paired++;need(e.side[0].sign==-e.side[1].sign,"opposite canonical face scatter");}
        e.n=4+e.reference.rank;for(int j=0;j<e.n;j++)e.id[j]=j<4?e.reference.poly[j]:e.reference.id[j-4];V2 c[6],v[6];for(int j=0;j<e.n;j++){const auto&q=f.coefficient[e.id[j]];c[j]={q.x[0],q.x[1]};v[j]={q.v[0],q.v[1]};}V2 a=c[0],b=c[3],d=b-a;e.volume=cross(a,b)/2;e.G[0]=.5L*cw(b);e.G[3]=-.5L*cw(a);
        for(int k=0;k<64;k++){R beta[6]={},dbeta[6]={},N[6]={};edge_shapes(e.reference,x[k],beta,dbeta,N);V2 delta={},delta_prime={};
            // Subtract the endpoint linear interpolant before evaluation.
            V2 curvature1=(c[1]-a)-(1.L/3)*d,curvature2=(c[2]-a)-(2.L/3)*d;
            delta=beta[1]*curvature1+beta[2]*curvature2;delta_prime=dbeta[1]*curvature1+dbeta[2]*curvature2;
            for(int j=4;j<e.n;j++){delta=delta+beta[j]*c[j];delta_prime=delta_prime+dbeta[j]*c[j];}
            R weight=w[k];e.volume+=weight*(cross(delta,d)+.5L*cross(delta,delta_prime));
            for(int j=0;j<e.n;j++){R endpoint=(j==3?1:0)-(j==0?1:0);V2 derivative=beta[j]*cw(d)-endpoint*cw(delta)+.5L*beta[j]*cw(delta_prime)-.5L*dbeta[j]*cw(delta);e.G[j]=e.G[j]+weight*derivative;e.cofactor[j]=e.cofactor[j]+(weight*N[j])*cw(d+delta_prime);}
        }
        // The open-edge area derivative has endpoint gauge terms. Remove
        // them for a physical canonical face force; the gauges cancel
        // exactly around a closed cell, so sum signed H == dV_cell/dx.
        for(int j=0;j<e.n;j++)e.H[j]=e.G[j];e.H[0]=e.H[0]-.5L*cw(a);e.H[3]=e.H[3]+.5L*cw(b);
        for(int j=0;j<e.n;j++){e.vdot+=inner(e.G[j],v[j]);for(int k=0;k<2;k++){out.derivative_error=std::max(out.derivative_error,fabsl(e.H[j][k]-e.cofactor[j][k]));out.derivative_scale=std::max(out.derivative_scale,fabsl(e.cofactor[j][k]));}}
        for(const auto&s:e.side){out.volume[s.cell]+=s.sign*e.volume;for(int j=0;j<e.n;j++)out.G[s.cell][e.id[j]]=out.G[s.cell][e.id[j]]+R(s.sign)*e.H[j];}out.edge.push_back(e);
    }
    for(size_t c=0;c<f.cell.size();c++){V2 translation={};R rotation=0;for(auto&entry:out.G[c]){int id=entry.first;const auto&q=f.coefficient[id];if(f.dof[id].kind<3)translation=translation+entry.second;rotation+=cross({q.x[0],q.x[1]},entry.second);out.vdot[c]+=inner({q.v[0],q.v[1]},entry.second);}out.translation=std::max(out.translation,hypotl(translation[0],translation[1]));out.rotation=std::max(out.rotation,fabsl(rotation));need(out.volume[c]>0,"positive canonical cell boundary volume");}return out;
}
}
#endif
