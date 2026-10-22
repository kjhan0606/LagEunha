// Read-only topology/count preflight. No field integration or evolution.
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <vector>
#include <cfloat>

using U = uint64_t;
using R = long double;
using Point = std::array<double, 2>;
using Key = std::pair<U, U>;
struct Triangle { U cell, v[3]; };
struct Edge { U midpoint = U(-1); std::vector<U> parents; };
static void need(bool b, const char* msg) { if (!b) throw std::runtime_error(msg); }
// Exact dyadic polynomial predicates via error-free binary64 expansions.
// Reject overflow/subnormal products rather than silently losing a term.
struct Q {
    std::vector<double> e;
    Q(double x=0) {if(x)e.push_back(x);}
    void add(double b) {
        std::vector<double> out;
        for(double a:e) {
            double sum=a+b,z=sum-a;
            double error=(a-(sum-z))+(b-z);
            need(std::isfinite(sum),"finite expansion sum");
            if(error)out.push_back(error);b=sum;
        }
        if(b)out.push_back(b);e.swap(out);
    }
    int sign()const {return e.empty()?0:(e.back()>0?1:-1);}
    friend Q operator+(Q a,const Q&b){for(double x:b.e)a.add(x);return a;}
    friend Q operator-(Q a,const Q&b){for(double x:b.e)a.add(-x);return a;}
    friend Q operator-(const Q&a){return Q()-a;}
    friend Q operator*(const Q&a,const Q&b){
        Q out;for(double x:a.e)for(double y:b.e){
            double p=x*y;need(std::isfinite(p)&&fabs(p)>=DBL_MIN,"normal expansion product");
            double tail=std::fma(x,y,-p);
            need(tail==0||fabs(tail)>=DBL_MIN,"normal expansion tail");
            out.add(tail);out.add(p);
        }return out;
    }
    friend bool operator<(const Q&a,const Q&b){return (a-b).sign()<0;}
    friend bool operator<=(const Q&a,const Q&b){return (a-b).sign()<=0;}
    friend bool operator>=(const Q&a,const Q&b){return (a-b).sign()>=0;}
};
template<class T> static void read(std::ifstream& f, T* p, size_t n) {
    f.read(reinterpret_cast<char*>(p), n*sizeof(T)); need(bool(f), "complete input read");
}
static Q exact(double x) {
    need(std::isfinite(x), "finite canonical coordinate");
    return Q(x);
}
static bool core_exact(const Triangle& t, const std::vector<Point>& v) {
    std::array<std::array<Q,2>,3> a;
    for (int k=0;k<3;k++) for (int d=0;d<2;d++) a[k][d]=exact(v[t.v[k]][d]);
    bool pos=true,neg=true;
    for(int k=0;k<3;k++) {
        const auto& x=a[k]; const auto& y=a[(k+1)%3];
        Q cross=x[0]*y[1]-x[1]*y[0];pos&=cross>=0;neg&=cross<=0;
    }
    if(pos||neg) return true;
    const Q radius2=Q(.25);
    for(int k=0;k<3;k++) {
        const auto& x=a[k]; const auto& y=a[(k+1)%3];
        Q dx=y[0]-x[0],dy=y[1]-x[1],den=dx*dx+dy*dy;
        Q numer=-(x[0]*dx+x[1]*dy);
        if(numer<=0) { if(x[0]*x[0]+x[1]*x[1]<radius2) return true; }
        else if(numer>=den) { if(y[0]*y[0]+y[1]*y[1]<radius2) return true; }
        else { Q cross=x[0]*y[1]-x[1]*y[0];if(cross*cross<radius2*den) return true; }
    }
    return false;
}
int main(int argc,char**argv) try {
    need(argc==2,"usage: inspector saved_canonical_reference.bin");
    const auto begin=std::chrono::steady_clock::now();
    std::ifstream in(argv[1],std::ios::binary);need(bool(in),"canonical input open");
    U h[8];read(in,h,8);
    need(h[0]==0x4d4643414e4f4e31ULL&&h[1]==1&&h[6]==sizeof(Point)&&h[7]==sizeof(Triangle),"canonical version/ABI");
    const U nc=h[2],nv=h[3],nt=h[4];need(nt%4==0,"four-child source layout");
    std::vector<Point> v(nv);std::vector<Triangle> fine(nt);read(in,v.data(),v.size());read(in,fine.data(),fine.size());
    need(in.peek()==EOF,"no extra canonical bytes");
    std::vector<Triangle> parent(nt/4);std::map<Key,Edge> pe;
    std::vector<uint8_t> core(nc);U exact_calls=0;R closest_margin=1e99L;
    for(U p=0;p<nt/4;p++) {
        const Triangle* f=fine.data()+4*p;
        U a=f[0].v[0],ab=f[0].v[1],ca=f[0].v[2],b=f[1].v[1],bc=f[1].v[2],c=f[2].v[2];
        parent[p]={f[0].cell,{a,b,c}};
        need(f[0].cell<nc,"source cell range");
        for(int k=1;k<4;k++)need(f[k].cell==f[0].cell,"four children share source");
        need(f[1].v[0]==ab&&f[2].v[0]==ca&&f[2].v[1]==bc&&f[3].v[0]==ab&&f[3].v[1]==bc&&f[3].v[2]==ca,"stored child connectivity");
        U ends[3][3]={{a,b,ab},{b,c,bc},{c,a,ca}};
        for(auto& e:ends) {
            need(e[0]<nv&&e[1]<nv&&e[2]<nv,"vertex range");
            for(int d=0;d<2;d++)need(v[e[2]][d]==double((R(v[e[0]][d])+R(v[e[1]][d]))/2),"bitwise stored midpoint construction");
            auto& edge=pe[std::minmax(e[0],e[1])];
            need(edge.midpoint==U(-1)||edge.midpoint==e[2],"shared midpoint identity");
            edge.midpoint=e[2];edge.parents.push_back(p);
        }
        // Cheap radial bounds avoid rational predicates far from the core.
        // Ambiguous triangles use exact rational geometry of binary64 input.
        R minr=1e99L,maxr=0;
        for(U id:parent[p].v) { R r=hypotl(R(v[id][0]),R(v[id][1]));minr=std::min(minr,r);maxr=std::max(maxr,r); }
        R diameter=0;for(int i=0;i<3;i++)for(int j=0;j<i;j++)diameter=std::max(diameter,hypotl(R(v[parent[p].v[i]][0])-v[parent[p].v[j]][0],R(v[parent[p].v[i]][1])-v[parent[p].v[j]][1]));
        bool intersects;
        if(minr<.5L-1e-10L)intersects=true;
        else if(minr-diameter>.5L+1e-10L)intersects=false;
        else { intersects=core_exact(parent[p],v);exact_calls++; }
        closest_margin=std::min(closest_margin,fabsl(minr-.5L));
        core[f[0].cell]|=intersects;
    }
    // The canonical disk is connected; verify via parent adjacency.
    std::vector<U> root(parent.size());for(U i=0;i<root.size();i++)root[i]=i;
    auto find=[&](U a){while(root[a]!=a){root[a]=root[root[a]];a=root[a];}return a;};
    U hanging=0,rounded_midpoint_edges=0,noncollinear_midpoint_edges=0;
    R max_midpoint_rounding=0;
    for(auto& item:pe) {
        const auto& e=item.second;need(e.parents.size()==1||e.parents.size()==2,"parent edge manifold");
        const auto& a0=v[item.first.first];const auto& b0=v[item.first.second];const auto& m0=v[e.midpoint];
        bool rounded=false;Q a[2],b[2],m[2];
        for(int d=0;d<2;d++) {
            a[d]=exact(a0[d]);b[d]=exact(b0[d]);m[d]=exact(m0[d]);
            rounded|=(Q(2)*m[d]-a[d]-b[d]).sign()!=0;
            max_midpoint_rounding=std::max(max_midpoint_rounding,fabsl(R(m0[d])-(R(a0[d])+R(b0[d]))/2));
        }
        rounded_midpoint_edges+=rounded;
        noncollinear_midpoint_edges+=((b[0]-a[0])*(m[1]-a[1])-(b[1]-a[1])*(m[0]-a[0])).sign()!=0;
        if(e.parents.size()==2) {
            U a=e.parents[0],b=e.parents[1];root[find(a)]=find(b);
            hanging+=core[parent[a].cell]!=core[parent[b].cell];
        }
    }
    std::set<U> components;for(U i=0;i<root.size();i++)components.insert(find(i));
    std::set<U> active_vertices;std::set<Key> active_edges;
    U active_triangles=0,protected_parents=0;
    auto include=[&](const Triangle&t){
        active_triangles++;
        for(int k=0;k<3;k++){active_vertices.insert(t.v[k]);active_edges.insert(std::minmax(t.v[k],t.v[(k+1)%3]));}
    };
    for(U p=0;p<parent.size();p++) {
        if(core[parent[p].cell]){protected_parents++;for(int k=0;k<4;k++)include(fine[4*p+k]);}
        else include(parent[p]);
    }
    // Each nonconforming edge has 2 coarse interior masters and 5 fine
    // trace slaves (midpoint + 4 interiors). All slaves are constrained.
    U raw=active_vertices.size()+2*active_edges.size()+active_triangles;
    need(raw>=5*hanging,"hanging constraint count");
    U polynomial=raw-5*hanging,ndof=polynomial+2*components.size();
    U corecells=std::count(core.begin(),core.end(),uint8_t(1));
    double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count();
    std::cout<<std::setprecision(17)
      <<"{\n  \"scope\":\"saved-reference topology only; no state transfer, integration or evolution\",\n"
      <<"  \"all_parent_child_connectivity_and_stored_midpoints_verified\":true,\n"
      <<"  \"source_cells\":"<<nc<<",\n  \"original_fine_triangles\":"<<nt<<",\n"
      <<"  \"parent_triangles\":"<<parent.size()<<",\n  \"protected_core_cells\":"<<corecells<<",\n"
      <<"  \"protected_parent_triangles\":"<<protected_parents<<",\n  \"active_triangles\":"<<active_triangles<<",\n"
      <<"  \"active_vertices_before_hanging_constraints\":"<<active_vertices.size()<<",\n"
      <<"  \"active_edges_before_hanging_constraints\":"<<active_edges.size()<<",\n"
      <<"  \"nonconforming_interfaces\":"<<hanging<<",\n  \"hanging_slave_dofs\":"<<5*hanging<<",\n"
      <<"  \"polynomial_scalar_dofs\":"<<polynomial<<",\n  \"connected_components\":"<<components.size()<<",\n"
      <<"  \"compact_scalar_dofs\":"<<ndof<<",\n  \"scalar_dof_upper_limit\":2500000,\n"
      <<"  \"dof_scale_gate_pass\":"<<(ndof<=2500000?"true":"false")<<",\n"
      <<"  \"bordered_mass_not_assembled\":true,\n  \"owned_ghost_partition_not_assembled\":true,\n"
      <<"  \"core_ambiguous_radial_bounds_resolved_exactly\":"<<exact_calls<<",\n"
      <<"  \"parent_edges_with_rounded_midpoint\":"<<rounded_midpoint_edges<<",\n"
      <<"  \"parent_edges_with_noncollinear_midpoint\":"<<noncollinear_midpoint_edges<<",\n"
      <<"  \"maximum_midpoint_coordinate_rounding\":"<<double(max_midpoint_rounding)<<",\n"
      <<"  \"exact_geometric_nesting_proven\":"<<(noncollinear_midpoint_edges==0?"true":"false")<<",\n"
      <<"  \"physical_function_embedding_not_yet_verified\":true,\n"
      <<"  \"minimum_vertex_radius_margin_to_core_boundary\":"<<double(closest_margin)<<",\n"
      <<"  \"seconds\":"<<seconds<<",\n  \"new_accepted_step\":false,\n  \"long_Kepler_validated\":false\n}\n";
    return 0;
} catch(const std::exception&e) {std::cerr<<"COMPACT_TOPOLOGY_REJECTED: "<<e.what()<<"\n";return 2;}
