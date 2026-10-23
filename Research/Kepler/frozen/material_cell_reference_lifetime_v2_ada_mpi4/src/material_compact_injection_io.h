#ifndef MATERIAL_COMPACT_INJECTION_IO_H
#define MATERIAL_COMPACT_INJECTION_IO_H
#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>
#include "material_shared_trace_batch.h"

namespace compact_io {
using R=long double;using U=uint64_t;
inline void require(bool b,const char*s){if(!b)throw std::runtime_error(s);}
inline U hash(const void*p,size_t n,U h){const auto*b=static_cast<const uint8_t*>(p);for(size_t k=0;k<n;k++){h^=b[k];h*=1099511628211ULL;}return h;}
struct Reader {
    FILE*f;U checksum=14695981039346656037ULL;
    explicit Reader(const std::string&path):f(fopen(path.c_str(),"rb")){require(f,"saved injection input");}
    template<class T>void get(T*p,size_t n){require(fread(p,sizeof(T),n,f)==n,"complete injection read");checksum=hash(p,n*sizeof(T),checksum);}
    void finish(){U stored;require(fread(&stored,sizeof(stored),1,f)==1&&stored==checksum,"injection payload checksum");require(fgetc(f)==EOF&&fclose(f)==0,"injection complete file");f=nullptr;}
    ~Reader(){if(f)fclose(f);}
};
struct Injection {
    U polynomial=0,total=0,old_total=0;int rank=0,ranks=0;
    std::vector<U>fine_global,row,global;
    std::vector<uint32_t>column;
    std::vector<double>weight;
    std::array<std::vector<R>,2>residual;
    std::vector<int>owner;
    std::vector<STCoefficient>initial;
    size_t border[2]{};
    explicit Injection(const std::string&path,int expected_rank,int expected_ranks){
        Reader in(path);U h[12];in.get(h,12);
        require(h[0]==0x4d43494e4a454331ULL&&h[1]==1&&h[2]==U(expected_rank)&&h[3]==U(expected_ranks),"injection owner/version");
        require(h[9]==sizeof(R)&&h[10]==sizeof(STCoefficient)&&h[7]==h[6]+2,"injection physical ABI");
        rank=h[2];ranks=h[3];polynomial=h[6];total=h[7];old_total=h[11];
        fine_global.resize(h[4]);row.resize(h[4]+1);column.resize(h[8]);weight.resize(h[8]);
        global.resize(h[5]);owner.resize(h[5]);initial.resize(h[5]);for(auto&r:residual)r.resize(h[4]);
        in.get(fine_global.data(),fine_global.size());in.get(row.data(),row.size());in.get(column.data(),column.size());in.get(weight.data(),weight.size());
        for(auto&r:residual)in.get(r.data(),r.size());in.get(global.data(),global.size());in.get(owner.data(),owner.size());in.get(initial.data(),initial.size());in.finish();
        require(row.front()==0&&row.back()==column.size()&&std::is_sorted(row.begin(),row.end()),"injection CSR extent");
        require(std::is_sorted(global.begin(),global.end())&&std::adjacent_find(global.begin(),global.end())==global.end(),"unique compact local identities");
        for(uint32_t c:column)require(c<global.size()&&global[c]<polynomial,"only polynomial CSR columns");
        for(int d=0;d<2;d++){auto at=std::lower_bound(global.begin(),global.end(),polynomial+d);require(at!=global.end()&&*at==polynomial+U(d),"both rank2 border columns present");border[d]=at-global.begin();}
    }
    // Caller must synchronize compact ghosts before forward evaluation.
    std::vector<R>forward(const std::vector<R>&compact)const{
        require(compact.size()==global.size(),"compact forward extent");std::vector<R>fine(fine_global.size());
        for(size_t i=0;i<fine.size();i++){R z=0;for(U k=row[i];k<row[i+1];k++)z+=R(weight[k])*compact[column[k]];for(int d=0;d<2;d++)z+=residual[d][i]*compact[border[d]];fine[i]=z;}return fine;
    }
    // Unassembled, local-element loads: include every local row, then sum the
    // compact result across MPI exactly once. Never pass a ghost-replicated
    // assembled fine load here: that would count shared fine rows twice.
    std::vector<R>restrict_element_load(const std::vector<R>&fine)const{
        require(fine.size()==fine_global.size(),"local fine-load extent");std::vector<R>out(global.size());
        for(size_t i=0;i<fine.size();i++){for(U k=row[i];k<row[i+1];k++)out[column[k]]+=R(weight[k])*fine[i];for(int d=0;d<2;d++)out[border[d]]+=residual[d][i]*fine[i];}return out;
    }
    // Assembled fine vector: ONLY its unique owner rows contribute to T^T.
    std::vector<R>restrict_assembled_owned_load(const std::vector<R>&fine,const std::vector<size_t>&owned)const{
        require(fine.size()==fine_global.size(),"assembled fine-load extent");std::vector<R>out(global.size());
        for(size_t i:owned){require(i<fine.size(),"fine owner row");for(U k=row[i];k<row[i+1];k++)out[column[k]]+=R(weight[k])*fine[i];for(int d=0;d<2;d++)out[border[d]]+=residual[d][i]*fine[i];}return out;
    }
};
}
#endif
