#ifndef GFS_FULL_REFERENCE_GEOMETRY_H
#define GFS_FULL_REFERENCE_GEOMETRY_H
// Exact-coordinate incidence of the retained native full-domain polygons.
// This constructs reference geometry only; it never fits an evolved field.
namespace full_reference {
using Pt=std::array<double,2>;
struct Cell {uint64_t id;Pt generator;};
struct Triangle {uint64_t cell,vertex[3];};
struct Mesh {std::vector<Pt>node;std::vector<Cell>cell;std::vector<Triangle>triangle;uint64_t external_edges=0,paired_edges=0;};
struct Polygon {uint64_t id;std::vector<Pt>local;};
static Mesh read_native(const std::string&project,const std::string&polygon,double cx,double cy,double box){
    std::map<uint64_t,Pt>generator;FILE*f=fopen(project.c_str(),"rb");need(f,"retained full native project metadata");uint64_t h[8];
    need(fread(h,sizeof(h),1,f)==1&&h[0]==0x4746534447503031ULL&&h[1]==1&&h[6]>=24&&h[6]<10000,"full project version/stride");std::vector<char>row(h[6]);
    for(uint64_t i=0;i<h[2];i++){need(fread(row.data(),row.size(),1,f)==1,"full native project row");uint64_t id;Pt p;memcpy(&id,row.data(),8);memcpy(p.data(),row.data()+8,16);need(generator.emplace(id,p).second,"unique full source cell identity");}
    need(fgetc(f)==EOF,"full native project EOF");fclose(f);
    f=fopen(polygon.c_str(),"rb");need(f,"retained raw native polygons");uint64_t hp[4];need(fread(hp,sizeof(hp),1,f)==1&&hp[0]==0x4746534447504f31ULL&&hp[1]==1&&hp[2]==h[2],"full native polygon epoch/count");
    std::map<uint64_t,Polygon>polygons;for(uint64_t i=0;i<hp[2];i++){uint64_t id;int tag[2];need(fread(&id,8,1,f)==1&&fread(tag,8,1,f)==1&&tag[0]>=3&&tag[0]<=48&&tag[1]==1,"valid retained native polygon");Polygon p;p.id=id;p.local.resize(tag[0]);need(fread(p.local.data(),16,p.local.size(),f)==p.local.size(),"all native polygon vertices");need(generator.count(id)&&polygons.emplace(id,std::move(p)).second,"polygon/project source identity");}
    need(fgetc(f)==EOF,"full native polygon EOF");fclose(f);Mesh m;std::map<Pt,uint64_t>index;
    auto add=[&](Pt p){need(std::isfinite(p[0])&&std::isfinite(p[1]),"finite full reference vertex");auto it=index.find(p);if(it!=index.end())return it->second;uint64_t id=m.node.size();index.emplace(p,id);m.node.push_back(p);return id;};
    for(const auto&entry:polygons){const auto&p=entry.second;Pt g=generator.at(p.id);m.cell.push_back({p.id,g});
        std::vector<double>x,y;for(auto a:p.local){x.push_back(a[0]);y.push_back(a[1]);}GfsMomentPolygon q={x.data(),y.data(),int(x.size()),0};GfsMomentDensity geom={};need(gfs_moment_geometry(q,&geom),"positive exact native cell geometry");
        uint64_t middle=add({(g[0]+geom.gx)-cx,(g[1]+geom.gy)-cy});
        for(size_t j=0;j<p.local.size();j++){auto a=p.local[j],b=p.local[(j+1)%p.local.size()];if(a==b)continue;
            Triangle t={p.id,{middle,add({(g[0]+a[0])-cx,(g[1]+a[1])-cy}),add({(g[0]+b[0])-cx,(g[1]+b[1])-cy})}};
            need(t.vertex[0]!=t.vertex[1]&&t.vertex[0]!=t.vertex[2]&&t.vertex[1]!=t.vertex[2],"nondegenerate exact native fan");m.triangle.push_back(t);
        }
    }
    std::map<std::pair<uint64_t,uint64_t>,uint64_t>edges;
    for(auto&t:m.triangle)for(int j=0;j<3;j++)edges[std::minmax(t.vertex[j],t.vertex[(j+1)%3])]++;
    for(auto&e:edges){need(e.second==1||e.second==2,"manifold native full reference fan");if(e.second==2){m.paired_edges++;continue;}
        auto a=m.node[e.first.first],b=m.node[e.first.second];bool exterior=false;
        for(int d=0;d<2;d++)for(double side:{-box/2,box/2})if(a[d]==side&&b[d]==side)exterior=true;
        need(exterior,"unpaired interior edge: no tolerance-based vertex repair permitted");m.external_edges++;
    }
    need(m.node.size()+m.triangle.size()==edges.size()+1,"full connected square Euler characteristic");return m;
}
static Mesh refine_once(const Mesh&base){
    Mesh m;m.node=base.node;m.cell=base.cell;m.external_edges=2*base.external_edges;std::map<std::pair<uint64_t,uint64_t>,uint64_t>mid;
    auto midpoint=[&](uint64_t a,uint64_t b){auto key=std::minmax(a,b);auto it=mid.find(key);if(it!=mid.end())return it->second;uint64_t id=m.node.size();m.node.push_back({double((R(m.node[a][0])+m.node[b][0])/2),double((R(m.node[a][1])+m.node[b][1])/2)});mid.emplace(key,id);return id;};
    m.triangle.reserve(4*base.triangle.size());for(auto&t:base.triangle){uint64_t a=t.vertex[0],b=t.vertex[1],c=t.vertex[2],ab=midpoint(a,b),bc=midpoint(b,c),ca=midpoint(c,a);
        for(auto v:{std::array<uint64_t,3>{a,ab,ca},std::array<uint64_t,3>{ab,b,bc},std::array<uint64_t,3>{ca,bc,c},std::array<uint64_t,3>{ab,bc,ca}})m.triangle.push_back({t.cell,{v[0],v[1],v[2]}});
    }
    m.paired_edges=2*base.paired_edges+3*base.triangle.size();return m;
}
static void save(const Mesh&m,const std::string&path){
    FILE*f=fopen(path.c_str(),"wbx");need(f,"new full native reference artifact");uint64_t h[8]={0x4746534652454631ULL,1,m.node.size(),m.cell.size(),m.triangle.size(),sizeof(Pt),sizeof(Cell),sizeof(Triangle)};
    need(fwrite(h,sizeof(h),1,f)==1&&fwrite(m.node.data(),sizeof(Pt),m.node.size(),f)==m.node.size()&&fwrite(m.cell.data(),sizeof(Cell),m.cell.size(),f)==m.cell.size()&&fwrite(m.triangle.data(),sizeof(Triangle),m.triangle.size(),f)==m.triangle.size()&&fflush(f)==0&&fsync(fileno(f))==0,"durable native full reference");need(fclose(f)==0,"reference close");
}
}
#endif
