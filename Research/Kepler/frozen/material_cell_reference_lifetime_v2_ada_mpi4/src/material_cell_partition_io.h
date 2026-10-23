#ifndef MATERIAL_CELL_PARTITION_IO_H
#define MATERIAL_CELL_PARTITION_IO_H
#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <string>
#include <vector>
namespace cell_condensation {
using U=uint64_t;
struct Cell {U source=0,begin=0,count=0;std::vector<uint32_t>I,G;};
struct Reader {
    FILE*f;U hash=14695981039346656037ULL;
    explicit Reader(const std::string&p):f(fopen(p.c_str(),"rb")){need(f,"saved exact cell partition");}
    template<class T>void get(T*p,size_t n){need(fread(p,sizeof(T),n,f)==n,"complete partition payload");auto*b=reinterpret_cast<const unsigned char*>(p);for(size_t j=0;j<n*sizeof(T);j++){hash^=b[j];hash*=1099511628211ULL;}}
    void finish(){U got;need(fread(&got,sizeof(U),1,f)==1&&got==hash&&fgetc(f)==EOF,"partition checksum and EOF");need(fclose(f)==0,"partition close");f=nullptr;}
    ~Reader(){if(f)fclose(f);}
};
struct Partition {
    U nf=0,ng=0,ni=0,nt=0,old_total=0,total=0,points=0;
    std::vector<uint32_t>old_G,column;std::vector<U>global_G,row,point_offset;std::vector<int>owner_G;std::vector<Cell>cell;
    Partition(const std::string&path,int rank,int ranks){
        Reader f(path);U h[14];f.get(h,14);
        need(h[0]==0x4d43434f4e443031ULL&&h[1]==1&&h[2]==U(rank)&&h[3]==U(ranks)&&h[9]==10816441&&h[10]<=2500000&&h[12]==sizeof(long double),"exact full field partition ABI and skeleton bound");
        nf=h[4];ng=h[5];ni=h[7];nt=h[11];old_total=h[9];total=h[10];points=h[13];
        need(nf<UINT32_MAX&&ng<=nf&&ng+ni==nf&&h[6]<=71080&&h[8]<=1000000000&&nt<=1138432,"bounded partition dimensions");
        old_G.resize(ng);global_G.resize(ng);owner_G.resize(ng);cell.resize(h[6]);row.resize(ng+1);column.resize(h[8]);point_offset.resize(nt+1);
        f.get(old_G.data(),ng);f.get(global_G.data(),ng);f.get(owner_G.data(),ng);
        need(std::is_sorted(old_G.begin(),old_G.end())&&std::adjacent_find(old_G.begin(),old_G.end())==old_G.end()&&old_G.back()<nf,"unique original skeleton IDs");
        U end=0,interiors=0;
        for(auto&c:cell){U a[5];f.get(a,5);c.source=a[0];c.begin=a[1];c.count=a[2];need(c.begin==end&&c.count>0&&c.begin+c.count<=nt&&a[3]>0&&a[3]<=nf&&a[4]>0&&a[4]<=ng,"whole owned source partition");end+=c.count;c.I.resize(a[3]);c.G.resize(a[4]);f.get(c.I.data(),c.I.size());f.get(c.G.data(),c.G.size());
            need(std::is_sorted(c.I.begin(),c.I.end())&&std::is_sorted(c.G.begin(),c.G.end())&&std::adjacent_find(c.I.begin(),c.I.end())==c.I.end()&&std::adjacent_find(c.G.begin(),c.G.end())==c.G.end()&&c.I.back()<nf&&c.G.back()<ng,"unique valid source supports");interiors+=c.I.size();}
        need(end==nt&&interiors==ni,"complete disjoint source counts");
        f.get(row.data(),row.size());f.get(column.data(),column.size());f.get(point_offset.data(),point_offset.size());f.finish();
        need(row[0]==0&&row.back()==column.size()&&std::is_sorted(row.begin(),row.end())&&point_offset.front()==0&&point_offset.back()==points&&std::is_sorted(point_offset.begin(),point_offset.end()),"complete CSR and point offsets");
        for(size_t i=0;i<ng;i++){need(std::is_sorted(column.begin()+row[i],column.begin()+row[i+1]),"sorted skeleton CSR columns");for(U j=row[i];j<row[i+1];j++)need(column[j]<ng,"bounded skeleton CSR column");}
    }
};
}
#endif
