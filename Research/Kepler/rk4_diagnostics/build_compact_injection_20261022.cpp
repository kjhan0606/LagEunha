// Native preparation of the ONE approved compact trial space. No evolution.
// Reuse the frozen preflight's dyadic predicates without changing that source.
#define main preserved_topology_count_entry_not_invoked
#include "inspect_compact_reference_topology_20261022.cpp"
#undef main
#include <mpi.h>
#include <omp.h>
#include <climits>
#include <cstring>
#include <string>
#include <unistd.h>
#include "material_full_owned.h"
#include "material_distributed_mass.h"

namespace compact_injection {
using Index=uint32_t;
constexpr Index absent=UINT32_MAX;
struct Entry {Index column;double weight;};
struct Row {
    std::array<Entry,10> entry{};int n=0;
    void add(Index c,double w){
        if(w==0)return;need(c!=absent,"valid compact polynomial master");
        for(int k=0;k<n;k++)if(entry[k].column==c){entry[k].weight+=w;return;}
        need(n<int(entry.size()),"bounded polynomial injection row");entry[n++]={c,w};
    }
    void check()const{R sum=0;for(int k=0;k<n;k++)sum+=entry[k].weight;need(sum==1,"exact constant in polynomial injection");}
};
struct FineEdge {U first_triangle=U(-1),index=U(-1);int owner=INT_MAX;};
struct ParentEdge {U midpoint=U(-1);std::vector<U> parents;};
struct Topology {
    U cells=0,nv=0,nt=0,ne=0,fine_poly=0;
    std::vector<Point> vertex,master_position,master_velocity;
    std::vector<Triangle> fine,parent;
    std::vector<int> cell_owner,fine_vertex_owner,master_owner;
    std::vector<uint8_t> core;
    std::map<Key,FineEdge> fine_edge;
    std::vector<Key> fine_edge_key;
    std::map<Key,ParentEdge> parent_edge;
    std::map<U,Key> midpoint_parent,slave_vertex;
    std::map<Key,Key> slave_edge;
    std::vector<Index> vertex_master,fine_interior,parent_interior;
    std::map<Key,std::array<Index,2>> edge_master;
    U interfaces=0,protected_cells=0;

    Index add(Point x){
        need(master_position.size()<absent,"compact metadata index capacity");
        Index id=master_position.size();master_position.push_back(x);
        master_velocity.push_back(shared_trace::initial_velocity(x,.047));
        return id;
    }
    Point edge_point(Key e,int part)const{
        Point x;for(int d=0;d<2;d++)x[d]=double((R(3-part)*vertex[e.first][d]+R(part)*vertex[e.second][d])/3);return x;
    }
    Point center(const Triangle&t)const{
        Point x;for(int d=0;d<2;d++)x[d]=double((R(vertex[t.v[0]][d])+vertex[t.v[1]][d]+vertex[t.v[2]][d])/3);return x;
    }
    Topology(const std::string&reference,int ranks){
        std::ifstream input(reference+"/material_canonical_reference.bin",std::ios::binary);need(bool(input),"saved canonical input");
        U h[8];read(input,h,8);need(h[0]==0x4d4643414e4f4e31ULL&&h[1]==1&&h[6]==sizeof(Point)&&h[7]==sizeof(Triangle),"canonical ABI");
        cells=h[2];nv=h[3];nt=h[4];ne=h[5];need(nt%4==0,"single stored h/2 subdivision");fine_poly=nv+2*ne+nt;
        vertex.resize(nv);fine.resize(nt);read(input,vertex.data(),vertex.size());read(input,fine.data(),fine.size());need(input.peek()==EOF,"complete canonical input");
        cell_owner.assign(cells,-1);
        for(int rank=0;rank<ranks;rank++){
            FILE*f=fopen((reference+"/material_native_geometry_r"+std::to_string(rank)+".bin").c_str(),"rb");need(f,"saved owner shard");
            U b[8];material_full::read(f,b,8);need(b[0]==0x4d46524157474531ULL&&b[1]==1&&b[3]==sizeof(MFCell)&&b[4]==sizeof(MFFace)&&b[5]==U(rank),"native owner shard ABI");
            for(U k=0;k<b[2];k++){MFCell c;material_full::read(f,&c);need(c.id<cells&&cell_owner[c.id]==-1,"unique source owner");cell_owner[c.id]=rank;need(fseek(f,c.n*sizeof(MFFace),SEEK_CUR)==0,"saved native face skip");}
            need(fgetc(f)==EOF&&fclose(f)==0,"complete owner shard");
        }
        need(std::find(cell_owner.begin(),cell_owner.end(),-1)==cell_owner.end(),"all source owners");
        parent.resize(nt/4);core.resize(cells);fine_vertex_owner.assign(nv,INT_MAX);
        for(U p=0;p<parent.size();p++){
            const Triangle*f=fine.data()+4*p;U a=f[0].v[0],ab=f[0].v[1],ca=f[0].v[2],b=f[1].v[1],bc=f[1].v[2],c=f[2].v[2];
            parent[p]={f[0].cell,{a,b,c}};need(f[0].cell<cells,"source range");
            need(f[1].cell==f[0].cell&&f[2].cell==f[0].cell&&f[3].cell==f[0].cell,"child source identity");
            need(f[1].v[0]==ab&&f[2].v[0]==ca&&f[2].v[1]==bc&&f[3].v[0]==ab&&f[3].v[1]==bc&&f[3].v[2]==ca,"stored four-child connectivity");
            U ends[3][3]={{a,b,ab},{b,c,bc},{c,a,ca}};
            for(auto&e:ends){
                need(e[0]<nv&&e[1]<nv&&e[2]<nv,"vertex range");Key key=std::minmax(e[0],e[1]);auto&z=parent_edge[key];
                need(z.midpoint==U(-1)||z.midpoint==e[2],"canonical midpoint sharing");z.midpoint=e[2];z.parents.push_back(p);
                auto old=midpoint_parent.emplace(e[2],key);need(old.second||old.first->second==key,"unique midpoint parent edge");
                for(int d=0;d<2;d++)need(vertex[e[2]][d]==double((R(vertex[e[0]][d])+vertex[e[1]][d])/2),"stored rounded midpoint construction");
            }
            R minr=1e99L,diameter=0;
            for(U id:parent[p].v)minr=std::min(minr,hypotl(vertex[id][0],vertex[id][1]));
            for(int i=0;i<3;i++)for(int j=0;j<i;j++)diameter=std::max(diameter,hypotl(R(vertex[parent[p].v[i]][0])-vertex[parent[p].v[j]][0],R(vertex[parent[p].v[i]][1])-vertex[parent[p].v[j]][1]));
            bool hit=minr<.5L-1e-10L?true:minr-diameter>.5L+1e-10L?false:core_exact(parent[p],vertex);
            core[f[0].cell]|=hit;
        }
        std::set<U> active_vertex;std::set<Key> active_edge;
        auto include=[&](const Triangle&t){for(int a=0;a<3;a++){active_vertex.insert(t.v[a]);active_edge.insert(std::minmax(t.v[a],t.v[(a+1)%3]));}};
        for(U s=0;s<nt;s++){
            const auto&t=fine[s];int owner=cell_owner[t.cell];
            for(int a=0;a<3;a++){
                fine_vertex_owner[t.v[a]]=std::min(fine_vertex_owner[t.v[a]],owner);
                auto&e=fine_edge[std::minmax(t.v[a],t.v[(a+1)%3])];e.first_triangle=std::min(e.first_triangle,s);e.owner=std::min(e.owner,owner);
            }
        }
        need(fine_edge.size()==ne,"saved fine edge count");
        for(auto&item:fine_edge){item.second.index=fine_edge_key.size();fine_edge_key.push_back(item.first);}
        for(U p=0;p<parent.size();p++){if(core[parent[p].cell])for(int k=0;k<4;k++)include(fine[4*p+k]);else include(parent[p]);}
        for(const auto&item:parent_edge){
            const auto&e=item.second;need(e.parents.size()==1||e.parents.size()==2,"parent manifold");
            if(e.parents.size()==2&&core[parent[e.parents[0]].cell]!=core[parent[e.parents[1]].cell]){
                interfaces++;slave_vertex[e.midpoint]=item.first;
                slave_edge[std::minmax(item.first.first,e.midpoint)]=item.first;
                slave_edge[std::minmax(item.first.second,e.midpoint)]=item.first;
            }
        }
        vertex_master.assign(nv,absent);fine_interior.assign(nt,absent);parent_interior.assign(parent.size(),absent);
        for(U id:active_vertex)if(!slave_vertex.count(id))vertex_master[id]=add(vertex[id]);
        for(Key key:active_edge)if(!slave_edge.count(key))edge_master[key]={add(edge_point(key,1)),add(edge_point(key,2))};
        for(U p=0;p<parent.size();p++){
            if(core[parent[p].cell])for(int k=0;k<4;k++)fine_interior[4*p+k]=add(center(fine[4*p+k]));
            else parent_interior[p]=add(center(parent[p]));
        }
        need(master_position.size()==1581181&&interfaces==240,"fixed candidate topology matches audited count");
        protected_cells=std::count(core.begin(),core.end(),uint8_t(1));need(protected_cells==5524,"fixed protected cell count");
        master_owner.assign(master_position.size()+2,INT_MAX);
        // Metadata only: choose each compact owner among actual source supports.
        for(U g=0;g<fine_poly;g++){
            int owner;
            if(g<nv)owner=fine_vertex_owner[g];
            else if(g<nv+2*ne)owner=fine_edge.at(fine_edge_key[(g-nv)/2]).owner;
            else owner=cell_owner[fine[g-(nv+2*ne)].cell];
            Row row=polynomial_row(g);row.check();
            for(int k=0;k<row.n;k++)master_owner[row.entry[k].column]=std::min(master_owner[row.entry[k].column],owner);
        }
        master_owner[master_position.size()]=master_owner[master_position.size()+1]=0;
        need(std::find(master_owner.begin(),master_owner.end(),INT_MAX)==master_owner.end(),"every compact master has actual owned support");
    }
    Row trace(Key key,int t6)const{
        Row out;const auto&e=edge_master.at(key);
        Index id[4]={vertex_master[key.first],e[0],e[1],vertex_master[key.second]};
        for(int j=0;j<4;j++){
            int numerator=1,denominator=1;for(int k=0;k<4;k++)if(k!=j){numerator*=t6-2*k;denominator*=2*(j-k);}
            double w=double(numerator)/denominator;out.add(id[j],w);
        }out.check();return out;
    }
    std::array<int,3> lambda2(U parent_id,U vertex_id)const{
        const auto&t=parent[parent_id];
        for(int a=0;a<3;a++)if(t.v[a]==vertex_id){std::array<int,3>x{};x[a]=2;return x;}
        for(int a=0;a<3;a++){
            int b=(a+1)%3;if(parent_edge.at(std::minmax(t.v[a],t.v[b])).midpoint==vertex_id){std::array<int,3>x{};x[a]=x[b]=1;return x;}
        }need(false,"fine vertex belongs to stored parent");return {};
    }
    Row parent_value(U p,const std::array<int,3>&lambda6)const{
        need(!core[parent[p].cell],"parent interpolation only outside protected fine core");Row out;
        const auto&t=parent[p];
        for(int a=0;a<10;a++){
            int m[3];gfs_material_p3_index(a,m);int numerator=1,denominator=1;
            for(int d=0;d<3;d++)for(int j=0;j<m[d];j++){numerator*=lambda6[d]-2*j;denominator*=2*(m[d]-j);}
            if(!numerator)continue;Index id=absent;
            if(a<3)id=vertex_master[t.v[a]];
            else if(a==9)id=parent_interior[p];
            else {int z=-1;for(int d=0;d<3;d++)if(m[d]==0)z=d;need(z>=0,"P3 edge node");int i=(z+1)%3,j=(z+2)%3;Key key=std::minmax(t.v[i],t.v[j]);int part=t.v[i]<t.v[j]?m[j]:m[i];need(part==1||part==2,"canonical P3 edge orientation");id=edge_master.at(key)[part-1];}
            double w=double(numerator)/denominator;out.add(id,w);
        }out.check();return out;
    }
    Row polynomial_row(U g)const{
        Row out;if(g>=fine_poly)return out;
        if(g<nv){
            if(vertex_master[g]!=absent){out.add(vertex_master[g],1);return out;}
            return trace(midpoint_parent.at(g),3);
        }
        if(g<nv+2*ne){
            U e=(g-nv)/2;int part=(g-nv)%2+1;Key key=fine_edge_key[e];
            auto master=edge_master.find(key);if(master!=edge_master.end()){out.add(master->second[part-1],1);return out;}
            auto slave=slave_edge.find(key);
            if(slave!=slave_edge.end()){
                Key whole=slave->second;U mid=parent_edge.at(whole).midpoint;
                auto param=[&](U a){need(a==whole.first||a==whole.second||a==mid,"interface endpoint");return a==whole.first?0:a==whole.second?6:3;};
                int t6=((3-part)*param(key.first)+part*param(key.second))/3;
                return trace(whole,t6);
            }
            U p=fine_edge.at(key).first_triangle/4;auto a=lambda2(p,key.first),b=lambda2(p,key.second);std::array<int,3>x;
            for(int d=0;d<3;d++)x[d]=(3-part)*a[d]+part*b[d];return parent_value(p,x);
        }
        U s=g-(nv+2*ne);if(fine_interior[s]!=absent){out.add(fine_interior[s],1);return out;}
        std::array<int,3>x{};for(U id:fine[s].v){auto l=lambda2(s/4,id);for(int d=0;d<3;d++)x[d]+=l[d];}
        return parent_value(s/4,x);
    }
};

static U hash_bytes(const void*p,size_t n,U h=14695981039346656037ULL){const auto*b=static_cast<const uint8_t*>(p);for(size_t k=0;k<n;k++){h^=b[k];h*=1099511628211ULL;}return h;}
struct Writer {
    FILE*f;U checksum=14695981039346656037ULL;
    explicit Writer(const std::string&path):f(fopen(path.c_str(),"wbx")){need(f,"new compact injection output");}
    template<class T>void put(const T*p,size_t n){need(fwrite(p,sizeof(T),n,f)==n,"compact injection write");checksum=hash_bytes(p,n*sizeof(T),checksum);}
    void finish(){need(fwrite(&checksum,sizeof(checksum),1,f)==1,"compact checksum write");material_full::durable(f);f=nullptr;}
    ~Writer(){if(f)fclose(f);}
};
}

int main(int argc,char**argv){
    MPI_Init(&argc,&argv);int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    try{
        need(argc==4&&ranks==4,"usage MPI4: builder reference_dir saved_initial_dir new_output_dir");
        double begin=MPI_Wtime();using namespace compact_injection;Topology topology(argv[1],ranks);
        auto fine=material_full::load_basis(std::string(argv[2])+"/full_conditioned_basis_r"+std::to_string(rank)+".bin");
        need(fine.global_vertices==topology.nv&&fine.global_edges==topology.ne&&fine.global_triangles==topology.nt&&fine.global_dofs==10816441,"saved conditioned field identity");
        std::vector<U>row(fine.global.size()+1),columns;
        std::vector<double>weight;
        std::array<std::vector<R>,2>residual;
        for(auto&r:residual)r.resize(fine.global.size());
        std::set<U>used;R local_position_error=0,local_velocity_error=0;U max_width=0;
        for(size_t i=0;i<fine.global.size();i++){
            U g=fine.global[i];Row r=topology.polynomial_row(g);R x[2]={},v[2]={};
            for(int k=0;k<r.n;k++){
                auto z=r.entry[k];columns.push_back(z.column);weight.push_back(z.weight);used.insert(z.column);
                for(int d=0;d<2;d++){x[d]+=R(z.weight)*topology.master_position[z.column][d];v[d]+=R(z.weight)*topology.master_velocity[z.column][d];}
            }
            row[i+1]=columns.size();max_width=std::max<U>(max_width,r.n);
            for(int d=0;d<2;d++){
                residual[d][i]=R(fine.f.coefficient[i].v[d])-v[d];
                local_position_error=std::max(local_position_error,fabsl(R(double(x[d]))-fine.f.coefficient[i].x[d]));
                local_velocity_error=std::max(local_velocity_error,fabsl(R(double(v[d]+residual[d][i]))-fine.f.coefficient[i].v[d]));
            }
        }
        U polynomial=topology.master_position.size();used.insert(polynomial);used.insert(polynomial+1);
        std::vector<U>global(used.begin(),used.end());std::vector<int>owner(global.size());std::vector<STCoefficient>initial(global.size());
        for(size_t i=0;i<global.size();i++){
            U g=global[i];owner[i]=topology.master_owner[g];
            if(g<polynomial)for(int d=0;d<2;d++){initial[i].x[d]=topology.master_position[g][d];initial[i].v[d]=topology.master_velocity[g][d];}
            else initial[i].v[g-polynomial]=1;
        }
        std::vector<Index>local_columns(columns.size());
        for(size_t k=0;k<columns.size();k++){auto at=std::lower_bound(global.begin(),global.end(),columns[k]);need(at!=global.end()&&*at==columns[k],"compact column remap");local_columns[k]=at-global.begin();}
        material_distributed::Ownership layout(global,owner,MPI_COMM_WORLD);
        need(layout.global_owned==polynomial+2,"actual unique compact MPI ownership");
        for(int d=0;d<4;d++){
            material_distributed::Vec v(initial.size());for(size_t i=0;i<v.size();i++)v[i]=d<2?initial[i].x[d]:initial[i].v[d-2];auto before=v;layout.pull(v);need(v==before,"compact initial owner/ghost bitmatch");
        }
        // Persist the actual injection, without evolving or mutating a field.
        Writer out(std::string(argv[3])+"/compact_injection_r"+std::to_string(rank)+".bin");
        U header[12]={0x4d43494e4a454331ULL,1,U(rank),U(ranks),fine.global.size(),global.size(),polynomial,polynomial+2,U(weight.size()),sizeof(R),sizeof(STCoefficient),fine.global_dofs};
        out.put(header,12);out.put(fine.global.data(),fine.global.size());out.put(row.data(),row.size());out.put(local_columns.data(),local_columns.size());out.put(weight.data(),weight.size());
        for(auto&r:residual)out.put(r.data(),r.size());out.put(global.data(),global.size());out.put(owner.data(),owner.size());out.put(initial.data(),initial.size());out.finish();
        U summary[5]={U(layout.owned.size()),layout.ghost_count,U(weight.size()),max_width,U(fine.global.size())};
        std::vector<U>all(rank==0?ranks*5:0);MPI_Gather(summary,5,MPI_UINT64_T,all.data(),5,MPI_UINT64_T,0,MPI_COMM_WORLD);
        R errors[2]={local_position_error,local_velocity_error};MPI_Allreduce(MPI_IN_PLACE,errors,2,MPI_LONG_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        need(errors[0]<1e-13L&&errors[1]<1e-13L,"reported initial coefficient roundoff only; no physical gate implied");
        if(rank==0){
            FILE*f=fopen((std::string(argv[3])+"/injection_result.json").c_str(),"wx");need(f,"new compact injection report");
            fprintf(f,"{\"scope\":\"actual shared polynomial injection and two dynamic initial-residual columns; no mass solve, quadrature or evolution\",\"MPI_ranks\":4,\"compact_scalar_dofs\":%llu,\"core_cells\":%llu,\"hanging_interfaces\":%llu,\"initial_position_coefficient_max_difference\":%.21Lg,\"initial_velocity_coefficient_max_difference\":%.21Lg,\"original_reference_rounding_kept_visible\":true,\"actual_MPI_owned_ghost_bitmatch\":true,\"rank_metadata\":[",(unsigned long long)(polynomial+2),(unsigned long long)topology.protected_cells,(unsigned long long)topology.interfaces,errors[0],errors[1]);
            for(int p=0;p<ranks;p++){auto*a=all.data()+5*p;fprintf(f,"%s{\"rank\":%d,\"owned\":%llu,\"ghost\":%llu,\"polynomial_nonzeros\":%llu,\"max_polynomial_row_width\":%llu,\"local_fine_rows\":%llu}",p?",":"",p,(unsigned long long)a[0],(unsigned long long)a[1],(unsigned long long)a[2],(unsigned long long)a[3],(unsigned long long)a[4]);}
            fprintf(f,"],\"seconds\":%.9g,\"initial_physical_gates_not_yet_evaluated\":true,\"mass_congruence_not_yet_assembled\":true,\"new_accepted_step\":false,\"long_Kepler_validated\":false}\n",MPI_Wtime()-begin);material_full::durable(f);
        }
    }catch(const std::exception&e){fprintf(stderr,"rank%d COMPACT_INJECTION_REJECTED %s\n",rank,e.what());MPI_Abort(MPI_COMM_WORLD,2);return 2;}
    MPI_Finalize();return 0;
}
