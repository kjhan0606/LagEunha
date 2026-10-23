#ifndef MATERIAL_CELL_FACTORS_IO_H
#define MATERIAL_CELL_FACTORS_IO_H
#include <cmath>
#include "material_cell_partition_io.h"
namespace cell_condensation {
using R=long double;
inline bool exact_DD_encoding(R x){double hi=double(x),lo=double(x-R(hi));return std::isfinite(x)&&std::isfinite(hi)&&std::isfinite(lo)&&R(hi)+R(lo)==x;}
struct Factors {
    std::vector<U>L_offset,IG_offset;std::vector<R>L,IG;U cells=0;
    Factors(const Partition&p,const std::string&path,int rank,int ranks){
        Reader f(path);U h[9];f.get(h,9);
        need(h[0]==0x4d43464143543031ULL&&h[1]==1&&h[2]==U(rank)&&h[3]==U(ranks)&&h[4]==p.cell.size()&&h[5]==p.nf&&h[6]==p.ng&&h[7]==p.total&&h[8]==sizeof(R),"original source factors ABI/partition identity");
        L_offset.resize(p.cell.size()+1);IG_offset.resize(p.cell.size()+1);cells=p.cell.size();
        for(size_t c=0;c<cells;c++){const auto&cell=p.cell[c];L_offset[c+1]=L_offset[c]+cell.I.size()*(cell.I.size()+1)/2;IG_offset[c+1]=IG_offset[c]+cell.I.size()*cell.G.size();}
        L.resize(L_offset.back());IG.resize(IG_offset.back());
        for(size_t c=0;c<cells;c++){
            const auto&cell=p.cell[c];U d[5];f.get(d,5);need(d[0]==cell.source&&d[1]==cell.begin&&d[2]==cell.count&&d[3]==cell.I.size()&&d[4]==cell.G.size()&&d[3]<=463&&d[4]<=140,"same complete original cell factor record");
            f.get(L.data()+L_offset[c],L_offset[c+1]-L_offset[c]);f.get(IG.data()+IG_offset[c],IG_offset[c+1]-IG_offset[c]);
            for(size_t i=0;i<cell.I.size();i++)need(L[L_offset[c]+i*(i+1)/2+i]>0,"every retained source Cholesky diagonal positive");
        }
        f.finish();int bad=0;
        #pragma omp parallel for schedule(static) reduction(+:bad)
        for(size_t i=0;i<L.size();i++)if(!exact_DD_encoding(L[i]))bad++;
        #pragma omp parallel for schedule(static) reduction(+:bad)
        for(size_t i=0;i<IG.size();i++)if(!exact_DD_encoding(IG[i]))bad++;
        need(!bad,"all original binary80 L/IG entries finite and exactly representable as DD");
    }
};
struct SkeletonMass {
    std::vector<U>row;std::vector<uint32_t>column;std::vector<R>value,diagonal;
    SkeletonMass(const Partition&p,const std::string&path,int rank,int ranks){
        Reader f(path);U h[8];f.get(h,8);
        need(h[0]==0x4d43534353523031ULL&&h[1]==1&&h[2]==U(rank)&&h[3]==U(ranks)&&h[4]==p.ng&&h[5]==p.column.size()&&h[6]==p.total&&h[7]==sizeof(R),"original Schur CSR ABI/partition identity");
        row.resize(p.row.size());column.resize(p.column.size());value.resize(p.column.size());diagonal.resize(p.ng);
        f.get(row.data(),row.size());f.get(column.data(),column.size());f.get(value.data(),value.size());f.get(diagonal.data(),diagonal.size());f.finish();need(row==p.row&&column==p.column,"same full preallocated exact source Schur CSR pattern");
        int bad=0;
        #pragma omp parallel for schedule(static) reduction(+:bad)
        for(size_t i=0;i<value.size();i++)if(!exact_DD_encoding(value[i]))bad++;
        #pragma omp parallel for schedule(static) reduction(+:bad)
        for(size_t i=0;i<diagonal.size();i++)if(!(diagonal[i]>0&&exact_DD_encoding(diagonal[i])))bad++;
        need(!bad,"all original Schur/assembled diagonal entries finite, exact DD encoding and positive diagonal");
        #pragma omp parallel for schedule(static) reduction(+:bad)
        for(size_t i=0;i<p.ng;i++)for(U j=row[i];j<row[i+1];j++){
            size_t k=column[j];auto at=std::lower_bound(column.begin()+row[k],column.begin()+row[k+1],uint32_t(i));
            if(at==column.begin()+row[k+1]||*at!=i||value[at-column.begin()]!=value[j])bad++;
        }
        need(!bad,"all retained original partial Schur entries exactly symmetric");
    }
};
}
#endif
