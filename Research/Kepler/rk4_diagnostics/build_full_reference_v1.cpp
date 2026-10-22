#include <algorithm>
#include <array>
#include <vector>
#include <map>
#include <string>
#include <cstdio>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <cfloat>
#include <stdexcept>
#include <chrono>
#include <unistd.h>
#include <omp.h>
using R=long double;
static void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
#include "gfs_enriched_euler.h"
#include "material_full_reference_geometry.h"
int main(int argc,char**argv){try{
    need(argc==4,"full native project polygon NEW_reference");auto start=std::chrono::steady_clock::now();auto base=full_reference::read_native(argv[1],argv[2],12,12,24);auto m=full_reference::refine_once(base);
    uint64_t logical=0,positive=0,invalid=0;R area=0;double minimum=DBL_MAX;
    #pragma omp parallel for schedule(static) reduction(+:logical,positive,invalid,area) reduction(min:minimum)
    for(size_t i=0;i<m.triangle.size();i++){
        const auto&t=m.triangle[i];double x[3],y[3];for(int a=0;a<3;a++){x[a]=m.node[t.vertex[a]][0];y[a]=m.node[t.vertex[a]][1];}
        R twice=(R(x[1])-x[0])*(R(y[2])-y[0])-(R(y[1])-y[0])*(R(x[2])-x[0]);if(!(twice>0)){invalid++;continue;}area+=twice/2;minimum=std::min(minimum,double(twice/2));
        GfsMomentPolygon p={x,y,3,0};GfsEnrichConfig cfg={};cfg.na=64;cfg.nr=28;cfg.field.ic.cx=cfg.field.ic.cy=12;GfsEnrichBasis basis={};basis.generator[0]=basis.generator[1]=12;GfsEnrichQuadrature q;
        if(!gfs_enrich_quadrature_init(p,&cfg,&basis,&q)){invalid++;continue;}logical+=gfs_enrich_nquad(p,&cfg,&q);
        for(int j=0;j<3;j++){int k=(j+1)%3;double ax=x[j]-q.anchor[0],ay=y[j]-q.anchor[1],bx=x[k]-q.anchor[0],by=y[k]-q.anchor[1],jac=ax*by-ay*bx;if(j==q.anchor_edge||(ax==0&&ay==0)||(bx==0&&by==0))jac=0;if(jac<0)invalid++;else if(jac>0)positive+=uint64_t(q.ns)*64*28;}
    }
    need(!invalid&&fabsl(area-576)<1e-10L,"full refined reference and actual fixed-rule geometry");full_reference::save(m,argv[3]);
    uint64_t edges=m.paired_edges+m.external_edges,scalar_p3=m.node.size()+2*edges+m.triangle.size();
    printf("{\"scope\":\"full retained native reference geometry and exact fixed-rule counts; no field/mass integration or state evolution\",\"cells\":%zu,\"base_triangles\":%zu,\"level1_triangles\":%zu,\"level1_vertices\":%zu,\"level1_edges\":%llu,\"external_edges\":%llu,\"P3_scalar_dofs\":%llu,\"enriched_scalar_dofs_lower_bound\":%llu,\"enriched_scalar_dofs_upper_bound\":%llu,\"actual_fixed_rule_logical_queries\":%llu,\"actual_fixed_rule_positive_points\":%llu,\"geometric_invalid_count\":%llu,\"reference_area\":%.21Lg,\"minimum_triangle_area\":%.17g,\"reference_artifact_bytes\":%llu,\"seconds\":%.9g,\"native_geometry_construction_pass\":true,\"full_disk_evolution_performed\":false}\n",m.cell.size(),base.triangle.size(),m.triangle.size(),m.node.size(),(unsigned long long)edges,(unsigned long long)m.external_edges,(unsigned long long)scalar_p3,(unsigned long long)(scalar_p3+edges+2*m.triangle.size()),(unsigned long long)(scalar_p3+2*edges+2*m.triangle.size()),(unsigned long long)logical,(unsigned long long)positive,(unsigned long long)invalid,area,minimum,(unsigned long long)(64+m.node.size()*sizeof(full_reference::Pt)+m.cell.size()*sizeof(full_reference::Cell)+m.triangle.size()*sizeof(full_reference::Triangle)),std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());
    return 0;
}catch(const std::exception&e){fprintf(stderr,"FULL_REFERENCE_FAIL %s\n",e.what());return 2;}}
