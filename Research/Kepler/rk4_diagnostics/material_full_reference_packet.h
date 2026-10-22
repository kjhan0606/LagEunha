#ifndef GFS_FULL_REFERENCE_PACKET_H
#define GFS_FULL_REFERENCE_PACKET_H
#include <climits>
#include <set>
#include "material_full_reference_geometry.h"
namespace full_reference {
// The topology constructor needs no stored per-triangle velocity map.
// Fresh velocity/thermal data are initialized once from the continuum IC;
// evolved states subsequently come from the actual restart coefficients.
struct Record {
    uint64_t cell;
    struct {struct {double x[3][2];}geometry;}triangle;
};
struct Packet {double parameters[12]={};std::vector<Record>record;};
static Mesh load(const std::string&path){
    FILE*f=fopen(path.c_str(),"rb");need(f,"retained full reference geometry input");uint64_t h[8];
    need(fread(h,sizeof(h),1,f)==1&&h[0]==0x4746534652454631ULL&&h[1]==1&&h[5]==sizeof(Pt)&&h[6]==sizeof(Cell)&&h[7]==sizeof(Triangle),"full reference version/ABI");
    need(h[2]>0&&h[3]>0&&h[4]>0&&h[2]<uint64_t(INT_MAX)&&h[4]<uint64_t(INT_MAX),"full reference native index capacity");
    Mesh m;m.node.resize(h[2]);m.cell.resize(h[3]);m.triangle.resize(h[4]);
    need(fread(m.node.data(),sizeof(Pt),m.node.size(),f)==m.node.size()&&fread(m.cell.data(),sizeof(Cell),m.cell.size(),f)==m.cell.size()&&fread(m.triangle.data(),sizeof(Triangle),m.triangle.size(),f)==m.triangle.size()&&fgetc(f)==EOF,"complete full reference artifact");fclose(f);
    std::set<uint64_t>cell;for(auto&c:m.cell)need(cell.insert(c.id).second,"unique full material cell ID");
    for(auto&t:m.triangle){need(cell.count(t.cell),"triangle source material context exists");for(auto id:t.vertex)need(id<m.node.size(),"full reference vertex ID range");}
    // This format is the already verified connected open-square triangulation.
    // Euler V-E+T=1 and 3T=2E-E_boundary fix both incidence counts.
    uint64_t edges=m.node.size()+m.triangle.size()-1;
    need(2*edges>=3*m.triangle.size(),"full reference open-boundary incidence count");
    m.external_edges=2*edges-3*m.triangle.size();m.paired_edges=edges-m.external_edges;return m;
}
static Packet packet(const Mesh&m,double gamma,double pressure,double eps,double rin,double rout){
    need(gamma>1&&pressure>0&&eps>0&&rin>0&&rout>rin,"fresh full material physical parameters");Packet p;
    p.parameters[1]=gamma;p.parameters[2]=pressure;p.parameters[3]=eps;p.parameters[4]=rin;p.parameters[5]=rout;p.record.reserve(m.triangle.size());
    for(auto&t:m.triangle){Record r={};r.cell=t.cell;for(int a=0;a<3;a++)for(int d=0;d<2;d++)r.triangle.geometry.x[a][d]=m.node[t.vertex[a]][d];p.record.push_back(r);}return p;
}
}
#endif
