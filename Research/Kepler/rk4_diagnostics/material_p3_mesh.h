#ifndef GFS_MATERIAL_P3_REFERENCE_MESH_H
#define GFS_MATERIAL_P3_REFERENCE_MESH_H
#include <tuple>
#include "gfs_material_p3.h"
#include "gfs_material_p3_radial_input.h"
/* Reference-only C++ topology. Original physical cells and boundaries persist. */
struct P3Mesh {
    SharedMesh linear;
    std::vector<Pair> node;
    std::vector<GfsMaterialP3Triangle> triangle;
    std::vector<std::array<int,10>> connectivity;
    std::vector<std::vector<int>> cell_nodes;
};
static SharedMesh p3_refine_once(const SharedMesh &input){
    SharedMesh mesh;mesh.node=input.node;mesh.cell=input.cell;mesh.components=input.components;std::map<std::pair<int,int>,int> mid;
    auto midpoint=[&](int a,int b){if(a>b)std::swap(a,b);auto key=std::make_pair(a,b);auto it=mid.find(key);if(it!=mid.end())return it->second;int id=mesh.node.size();mesh.node.push_back({.5*(mesh.node[a][0]+mesh.node[b][0]),.5*(mesh.node[a][1]+mesh.node[b][1])});mid[key]=id;return id;};
    for(const auto &t:input.tri){
        int a=t.vertex[0],b=t.vertex[1],c=t.vertex[2],ab=midpoint(a,b),bc=midpoint(b,c),ca=midpoint(c,a);
        int vertices[4][3]={{a,ab,ca},{ab,b,bc},{ca,bc,c},{ab,bc,ca}};
        for(int k=0;k<4;k++){Tri u={};u.cell=t.cell;u.component=t.component;for(int j=0;j<3;j++){u.vertex[j]=vertices[k][j];for(int d=0;d<2;d++)u.t.x[j][d]=mesh.node[u.vertex[j]][d];}need(gfs_material_triangle_init(&u.t),"positive h/2 P3 reference triangle");mesh.tri.push_back(u);}
    }
    for(size_t i=0;i<mesh.tri.size();i++)for(int a=0;a<3;a++){int x=mesh.tri[i].vertex[a],y=mesh.tri[i].vertex[(a+1)%3];if(x>y)std::swap(x,y);mesh.edge[{x,y}].push_back(i);}
    return mesh;
}
static P3Mesh p3_build(const SharedMesh &linear,const GfsEnrichConfig &cfg){
    P3Mesh mesh;mesh.linear=linear;mesh.node=linear.node;mesh.cell_nodes.resize(linear.cell.size());std::map<std::tuple<int,int,int>,int> edge;
    for(const auto &t:linear.tri){
        GfsMaterialP3Triangle poly={};poly.geometry=t.t;std::array<int,10> ids;
        for(int a=0;a<10;a++){
            int m[3];gfs_material_p3_index(a,m);int vertex=-1,zero=-1;for(int j=0;j<3;j++){if(m[j]==3)vertex=j;if(m[j]==0)zero=j;}
            if(vertex>=0)ids[a]=t.vertex[vertex];
            else if(zero>=0){
                int i=(zero+1)%3,j=(zero+2)%3,lo=t.vertex[i],hi=t.vertex[j],part=m[j];if(lo>hi){std::swap(lo,hi);part=3-part;}auto key=std::make_tuple(lo,hi,part);auto it=edge.find(key);
                if(it!=edge.end())ids[a]=it->second;
                else{ids[a]=mesh.node.size();Pair x;for(int d=0;d<2;d++)x[d]=(double)(((Real)(3-part)*mesh.node[lo][d]+(Real)part*mesh.node[hi][d])/3);mesh.node.push_back(x);edge[key]=ids[a];}
            }else{ids[a]=mesh.node.size();Pair x;for(int d=0;d<2;d++)x[d]=(double)(((Real)t.t.x[0][d]+t.t.x[1][d]+t.t.x[2][d])/3);mesh.node.push_back(x);}
            double B[4];velocity(cfg,mesh.node[ids[a]][0]+12,mesh.node[ids[a]][1]+12,poly.v0[a],B);mesh.cell_nodes[t.cell].push_back(ids[a]);
            double Dq[4];need(gfs_material_p3_radial_remainder(poly.geometry.x[0],mesh.node[ids[a]][0],mesh.node[ids[a]][1],cfg.field.ic.eps,poly.remainder[a],Dq),"stable P3 input remainder validity");
        }
        mesh.triangle.push_back(poly);mesh.connectivity.push_back(ids);
    }
    for(auto &ids:mesh.cell_nodes){std::sort(ids.begin(),ids.end());ids.erase(std::unique(ids.begin(),ids.end()),ids.end());}
    return mesh;
}
static void p3_shape_stable(const GfsMaterialP3Triangle &triangle,const GfsEnrichConfig &cfg,
    double X,double Y,const double scale[2],double N[12],double Nx[12],double Ny[12]){
    double q[2],Dq[4];need(gfs_material_p3_radial_remainder(triangle.geometry.x[0],X,Y,cfg.field.ic.eps,q,Dq),"stable P3 input remainder at point");gfs_material_p3_shape_remainder(&triangle,X,Y,q,Dq,scale,N,Nx,Ny);
}
#endif
