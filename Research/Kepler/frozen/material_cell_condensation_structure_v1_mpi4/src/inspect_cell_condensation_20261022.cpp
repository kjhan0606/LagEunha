// Original physical-cell support partition and allocated Schur CSR pattern.
// Reads saved topology and rule metadata only. No field/force/trajectory change.
#include <mpi.h>
#include <omp.h>
#include <algorithm>
#include <array>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>
#include <unistd.h>
using U=uint64_t;using R=long double;
static void need(bool b,const char*s){if(!b)throw std::runtime_error(s);}
#include "material_full_owned.h"
#include "material_distributed_mass.h"
namespace md=material_distributed;
struct Writer {
    FILE*f;U hash=14695981039346656037ULL;
    explicit Writer(const std::string&p):f(fopen(p.c_str(),"wbx")){need(f,"new partition output");}
    template<class T>void put(const T*p,size_t n){need(fwrite(p,sizeof(T),n,f)==n,"partition write");auto*b=reinterpret_cast<const unsigned char*>(p);for(size_t j=0;j<n*sizeof(T);j++){hash^=b[j];hash*=1099511628211ULL;}}
    void finish(){need(fwrite(&hash,sizeof(hash),1,f)==1,"partition checksum");material_full::durable(f);f=nullptr;}
    ~Writer(){if(f)fclose(f);}
};
struct Cell {U begin=0,count=0;std::vector<uint32_t>I,G;};
struct EdgeUse {U edge;size_t triangle;int slot;bool operator<(const EdgeUse&b)const{return edge<b.edge;}};
int main(int argc,char**argv){
    MPI_Init(&argc,&argv);int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    try{
        need(argc==4&&ranks==4,"MPI4 usage: structure original_initial_dir positive_rule_run new_output");
        double start=MPI_Wtime();std::string suffix=std::to_string(rank),dest=argv[3];
        auto m=material_full::load_basis(std::string(argv[1])+"/full_conditioned_basis_r"+suffix+".bin");
        need(m.global_dofs==10816441&&m.global_triangles==1138432&&m.global.size()<UINT32_MAX,"original full field identity");
        md::Ownership fine(m.global,m.owner,MPI_COMM_WORLD);
        size_t nf=m.global.size(),nt=m.f.triangle.size();std::vector<Cell>cell(m.f.cell.size());
        std::vector<int>lo(nf,INT_MAX),hi(nf,-1);md::Vec boundary(nf);
        for(size_t s=0;s<nt;s++){
            const auto&t=m.f.triangle[s];need(t.cell>=0&&size_t(t.cell)<cell.size()&&m.source_cell[s]==m.f.cell[t.cell],"exact source-cell provenance");
            auto&c=cell[t.cell];if(!c.count)c.begin=s;need(c.begin+c.count==s,"one contiguous whole source cell per rank");c.count++;
            need(t.n>=16&&t.n<=ST_MAX&&t.bubble_offset+2==t.n,"frozen enriched basis layout");
            for(int a=0;a<t.n;a++){int i=t.id[a];need(i>=0&&size_t(i)<nf,"original triangle support");lo[i]=std::min(lo[i],t.cell);hi[i]=std::max(hi[i],t.cell);}
        }
        U boundary_edges=0,internal_edges=0;
        for(size_t ci=0;ci<cell.size();ci++){
            const auto&c=cell[ci];need(c.count>0,"each owned source has triangles");std::vector<EdgeUse>edge;edge.reserve(c.count*3);
            for(size_t s=c.begin;s<c.begin+c.count;s++)for(int e=0;e<3;e++)edge.push_back({m.triangle_edge[3*s+e],s,e});
            std::sort(edge.begin(),edge.end());
            for(size_t j=0;j<edge.size();){size_t k=j+1;while(k<edge.size()&&edge[k].edge==edge[j].edge)k++;
                need(k-j==1||k-j==2,"source-cell edge incidence one or two");
                const auto&use=edge[j];const auto&t=m.f.triangle[use.triangle];const auto&e=t.edge[use.slot];need(e.rank==1||e.rank==2,"original canonical edge rank");
                if(k-j==2){const auto&other=edge[j+1];const auto&b=m.f.triangle[other.triangle].edge[other.slot];need(e.rank==b.rank,"same internal edge rank");for(int a=0;a<4;a++)need(e.poly[a]==b.poly[a],"same internal polynomial trace IDs");for(int a=0;a<e.rank;a++)need(e.id[a]==b.id[a],"same internal enriched trace IDs");internal_edges++;}
                else {
                    boundary_edges++;
                    // Analytic trace support: four P3 edge functions + edge residuals.
                    // P3 interior and conditioned bubbles vanish on every edge;
                    // other edge extensions vanish here by their endpoint identity.
                    for(int a=0;a<4;a++){need(e.poly[a]>=0&&size_t(e.poly[a])<nf,"edge P3 ID range");boundary[e.poly[a]]=1;}
                    for(int a=0;a<e.rank;a++){need(e.id[a]==t.id[t.offset[use.slot]+a],"canonical nonzero trace support");boundary[e.id[a]]=1;}
                    for(int a=0;a<10;a++){int index[3];gfs_material_p3_index(a,index);int opposite=3-t.endpoint[use.slot][0]-t.endpoint[use.slot][1];if(index[opposite]==0)need(std::find(e.poly,e.poly+4,t.id[a])!=e.poly+4,"all nonzero P3 traces marked");}
                }
                j=k;
            }
        }
        fine.sum(boundary);std::vector<uint32_t>old_G;std::vector<int>local_G(nf,-1),owner_G;std::vector<U>global_G;
        U interior=0;
        for(size_t i=0;i<nf;i++){
            need(hi[i]>=0,"all original DOFs supported");
            if(boundary[i]>0){local_G[i]=old_G.size();old_G.push_back(i);global_G.push_back(m.global[i]);owner_G.push_back(m.owner[i]);}
            else{need(lo[i]==hi[i]&&m.owner[i]==rank,"every interior DOF supported in exactly one owned physical cell");interior++;}
        }
        md::Ownership skeleton(global_G,owner_G,MPI_COMM_WORLD);
        need(skeleton.global_owned+fine.global_owned-skeleton.global_owned==m.global_dofs,"original unique ownership count");
        U local_cells=cell.size(),global_cells=0,global_I=0;MPI_Allreduce(&local_cells,&global_cells,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);MPI_Allreduce(&interior,&global_I,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);
        need(global_I+skeleton.global_owned==m.global_dofs&&global_cells==71080,"complete disjoint global interior/skeleton partition");
        U L_entries=0,IG_entries=0,GG_entries=0,maxI=0,maxG=0,pair_count=0,checkedI=0;
        FILE*csv=fopen((dest+"/cell_dimensions_r"+suffix+".csv").c_str(),"wx");need(csv,"new source dimensions");fprintf(csv,"source_cell,triangle_begin,triangles,interior,trace\n");
        for(size_t ci=0;ci<cell.size();ci++){
            auto&c=cell[ci];std::vector<uint32_t>ids;ids.reserve(c.count*ST_MAX);
            for(size_t s=c.begin;s<c.begin+c.count;s++){const auto&t=m.f.triangle[s];for(int a=0;a<t.n;a++)ids.push_back(t.id[a]);}
            std::sort(ids.begin(),ids.end());ids.erase(std::unique(ids.begin(),ids.end()),ids.end());
            for(auto i:ids)if(local_G[i]>=0)c.G.push_back(local_G[i]);else{need(lo[i]==int(ci)&&hi[i]==int(ci),"interior has no other-cell mass coupling");c.I.push_back(i);}
            U ni=c.I.size(),ng=c.G.size();need(ni>0&&ng>0&&std::is_sorted(c.G.begin(),c.G.end()),"nonempty exact source partition");
            checkedI+=ni;L_entries+=ni*(ni+1)/2;IG_entries+=ni*ng;GG_entries+=ng*ng;maxI=std::max(maxI,ni);maxG=std::max(maxG,ng);pair_count+=ng*(ng+1)/2;
            fprintf(csv,"%llu,%llu,%llu,%llu,%llu\n",(unsigned long long)m.f.cell[ci],(unsigned long long)c.begin,(unsigned long long)c.count,(unsigned long long)ni,(unsigned long long)ng);
        }
        material_full::durable(csv);need(checkedI==interior,"interior partition appears exactly once");
        std::vector<U>pairs;pairs.reserve(pair_count);
        for(const auto&c:cell)for(size_t a=0;a<c.G.size();a++)for(size_t b=a;b<c.G.size();b++)pairs.push_back((U(c.G[a])<<32)|c.G[b]);
        std::sort(pairs.begin(),pairs.end());pairs.erase(std::unique(pairs.begin(),pairs.end()),pairs.end());
        U ng=old_G.size(),nnz=2*pairs.size()-ng;std::vector<U>row(ng+1);
        for(U key:pairs){uint32_t a=key>>32,b=key;row[a+1]++;if(a!=b)row[b+1]++;}
        for(size_t i=1;i<row.size();i++)row[i]+=row[i-1];need(row.back()==nnz,"actual complete allocated partial Schur CSR pattern");
        std::vector<uint32_t>column(nnz);auto offset=row;
        for(U key:pairs){uint32_t a=key>>32,b=key;column[offset[a]++]=b;if(a!=b)column[offset[b]++]=a;}
        for(size_t i=0;i<ng;i++){std::sort(column.begin()+row[i],column.begin()+row[i+1]);need(std::binary_search(column.begin()+row[i],column.begin()+row[i+1],uint32_t(i)),"every skeleton diagonal present");}
        pairs.clear();pairs.shrink_to_fit();
        // Scan record metadata. Payload integrity remains checked by full producer
        // audits/manifests and must be rechecked by the eventual cache reader.
        FILE*f=fopen((std::string(argv[2])+"/positive_rule_r"+suffix+".bin").c_str(),"rb");need(f,"selected original-point rule metadata");U h[12];material_full::read(f,h,12);double anchor;material_full::read(f,&anchor);
        need(h[0]==0x4d46504f53525531ULL&&h[1]==1&&h[2]==U(rank)&&h[3]==U(ranks)&&h[4]==nt&&h[5]==m.global_triangles&&h[6]==sizeof(R)&&h[10]==16,"selected rule ABI/source identity");
        U selected=0,original=0,fixed_cache=0,rule_residual_bytes=0;std::vector<U>point_offset(nt+1);
        for(size_t s=0;s<nt;s++){
            U rec[6];material_full::read(f,rec,6);U rows=rec[4],nodes=rec[2];int n=m.f.triangle[s].n;
            need(rec[0]==m.source_triangle[s]&&rec[1]>0&&rec[3]==1&&nodes>0&&4*nodes<=rec[1]&&rows==U(n*n+4*n+4)&&rec[5]==0,"existing fully compressed immutable original-point rule");
            U skip=16*sizeof(R)+2*rows*sizeof(R)+nodes*(sizeof(uint32_t)+sizeof(double))+sizeof(U);need(skip<=U(INT64_MAX)&&fseeko(f,off_t(skip),SEEK_CUR)==0,"skip rule payload without re-evaluation");
            original+=rec[1];selected+=nodes;point_offset[s+1]=selected;fixed_cache+=nodes*(48+24*(n-10));rule_residual_bytes+=rows*2*sizeof(R);
        }
        need(fgetc(f)==EOF&&fclose(f)==0,"complete ordered selected-rule metadata");
        Writer out(dest+"/cell_partition_r"+suffix+".bin");
        U header[14]={0x4d43434f4e443031ULL,1,U(rank),U(ranks),U(nf),ng,U(cell.size()),interior,nnz,m.global_dofs,skeleton.global_owned,U(nt),sizeof(R),selected};out.put(header,14);
        out.put(old_G.data(),old_G.size());out.put(global_G.data(),global_G.size());out.put(owner_G.data(),owner_G.size());
        for(size_t ci=0;ci<cell.size();ci++){const auto&c=cell[ci];U desc[5]={m.f.cell[ci],c.begin,c.count,U(c.I.size()),U(c.G.size())};out.put(desc,5);out.put(c.I.data(),c.I.size());out.put(c.G.data(),c.G.size());}
        out.put(row.data(),row.size());out.put(column.data(),column.size());out.put(point_offset.data(),point_offset.size());out.finish();
        U factor_bytes=16*L_entries,IG_bytes=16*IG_entries,CSR_bytes=20*nnz+8*(ng+1);
        // Explicit proposed GPU layout; bytes are NOT an actual allocation proof.
        // Full original M, original triangles/rules, current+generating fields,
        // six full DD scratch vectors, six skeleton DD vectors, point IDs and
        // weights, original rule discrepancy ledger, plus fixed million-point tile.
        U full_mass=nt*ST_MAX*ST_MAX*16,context=nt*(sizeof(STTriangle)+sizeof(GfsEnrichQuadrature))+nf*64;
        U scratch=nf*96+ng*96,rule_bytes=selected*12+rule_residual_bytes+nt*(16*sizeof(R)+32);
        U tile=1048576ULL*(sizeof(GfsMaterialPoint12)+48+24*(ST_MAX-10));
        U proposed_resident=factor_bytes+IG_bytes+CSR_bytes+full_mass+context+scratch+rule_bytes+tile+fixed_cache;
        U proposed_tiled=proposed_resident-fixed_cache;
        // Counts per rank and SUM; I/G maxima reduced separately.
        enum{NC=22};U counts[NC]={U(nf),ng,U(skeleton.owned.size()),interior,U(cell.size()),U(nt),boundary_edges,internal_edges,nnz,L_entries,IG_entries,GG_entries,selected,original,fixed_cache,factor_bytes,IG_bytes,CSR_bytes,full_mass,proposed_resident,proposed_tiled,rule_bytes};
        U totals[NC]={};std::vector<U>byrank(rank==0?NC*ranks:0);MPI_Reduce(counts,totals,NC,MPI_UINT64_T,MPI_SUM,0,MPI_COMM_WORLD);MPI_Gather(counts,NC,MPI_UINT64_T,byrank.data(),NC,MPI_UINT64_T,0,MPI_COMM_WORLD);
        U maxima[2]={maxI,maxG};MPI_Allreduce(MPI_IN_PLACE,maxima,2,MPI_UINT64_T,MPI_MAX,MPI_COMM_WORLD);
        if(rank==0){need(totals[12]==249691129&&totals[13]==2040256512ULL,"same complete selected and original point counts");
            f=fopen((dest+"/structure_result.json").c_str(),"wx");need(f,"new exact partition result");
            const char*names[NC]={"local_fine_dofs_with_ghosts","local_skeleton_dofs_with_ghosts","owned_skeleton_dofs","interior_dofs","cells","triangles","source_boundary_edge_incidence","internal_edges","partial_CSR_nnz","packed_II_entries","IG_entries","cell_GG_entries","selected_points","original_points","immutable_cache_bytes","packed_II_DD_bytes","IG_DD_bytes","partial_CSR_DD_bytes","original_binary80_mass_bytes","proposed_resident_device_bytes","proposed_tiled_device_bytes","selected_rule_payload_bytes"};
            fprintf(f,"{\"scope\":\"original exact cell support partition, actual allocated CSR pattern and fixed-cache byte counts; no factorization, force or evolution\",\"original_scalar_dofs\":%llu,\"global_skeleton_dofs\":%llu,\"global_interior_dofs\":%llu,\"max_cell_interior\":%llu,\"max_cell_trace\":%llu,\"skeleton_dof_gate\":%s,\"exact_single_cell_interior_support_verified\":true,\"zero_boundary_trace_by_original_basis_identity\":true,\"rule_metadata_only\":true,\"device_capacity_gate_pending\":true,\"factorization_pending\":true,\"physical_trajectory_advanced\":false,\"seconds\":%.9g,\"totals\":{",(unsigned long long)m.global_dofs,(unsigned long long)skeleton.global_owned,(unsigned long long)global_I,(unsigned long long)maxima[0],(unsigned long long)maxima[1],skeleton.global_owned<=2500000?"true":"false",MPI_Wtime()-start);
            for(int j=0;j<NC;j++)fprintf(f,"%s\"%s\":%llu",j?",":"",names[j],(unsigned long long)totals[j]);fprintf(f,"},\"ranks\":[");
            for(int p=0;p<ranks;p++){fprintf(f,"%s{\"rank\":%d",p?",":"",p);for(int j=0;j<NC;j++)fprintf(f,",\"%s\":%llu",names[j],(unsigned long long)byrank[p*NC+j]);fprintf(f,"}");}fprintf(f,"]}\n");material_full::durable(f);
        }
    }catch(const std::exception&e){fprintf(stderr,"rank%d CELL_CONDENSATION_STRUCTURE_REJECTED %s\n",rank,e.what());MPI_Abort(MPI_COMM_WORLD,2);return 2;}
    MPI_Finalize();return 0;
}
