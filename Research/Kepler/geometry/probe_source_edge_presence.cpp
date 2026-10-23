#include "../material_compact_rezone_parallel.cpp"
#include <fstream>
#include <iostream>

int main(int argc,char**argv){
    int provided=0;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
    int rank=0,ranks=0;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    try{
        compact_p3::require(argc==8&&ranks==4&&provided>=MPI_THREAD_FUNNELED,
            "four-rank accepted-source edge incidence audit arguments");
        compact_source::Inputs input;
        input.initial=argv[1];input.immutable_manifest=argv[2];
        input.positive_directory=argv[3];input.whole_audit=argv[4];
        input.reference_receipt=argv[5];input.accepted_prefix=argv[6];
        input.accepted_steps=5;input.accepted_time=0.00011263377801805954;
        {
        auto source=compact_source::import_accepted_state(input,MPI_COMM_WORLD);
        const auto&owned=source.owned;
        compact_p3::require(owned.global_edges>0&&
            owned.triangle_edge.size()==3*owned.f.triangle.size()&&
            owned.external_edge.size()==owned.global_edges,
            "complete accepted original triangle-edge incidence registry");
        const size_t bytes=(owned.global_edges+7)/8;
        compact_p3::require(bytes<=INT_MAX,"finite MPI original edge bitmap");
        std::vector<unsigned char>present(bytes,0),global(bytes,0);
        for(uint64_t edge:owned.triangle_edge){
            compact_p3::require(edge<owned.global_edges,"original source edge in global registry");
            present[edge/8]|=static_cast<unsigned char>(1u<<(edge%8));
        }
        compact_p3::require(MPI_Allreduce(present.data(),global.data(),int(bytes),
            MPI_UNSIGNED_CHAR,MPI_BOR,MPI_COMM_WORLD)==MPI_SUCCESS,
            "all original source edge incidences assembled");
        size_t local_count=0,global_count=0;
        for(size_t edge=0;edge<owned.global_edges;edge++){
            local_count+=bool(present[edge/8]&(1u<<(edge%8)));
            global_count+=bool(global[edge/8]&(1u<<(edge%8)));
        }
        compact_p3::require(global_count==owned.global_edges,
            "every original global source edge has at least one triangle incidence");
        auto dir=std::filesystem::path(argv[7]);
        std::ofstream file(dir/("presence_rank"+std::to_string(rank)+".bin"),std::ios::binary);
        compact_p3::require(bool(file),"create source edge incidence bitmap");
        file.write(reinterpret_cast<const char*>(present.data()),present.size());
        file.close();compact_p3::require(bool(file),"complete source edge incidence bitmap");
        if(rank==0){
            for(uint8_t flag:owned.external_edge)
                compact_p3::require(flag<=1,"binary accepted source exterior edge flag");
            std::ofstream flags(dir/"external_edge.bin",std::ios::binary);
            compact_p3::require(bool(flags),"create original exterior edge flags");
            flags.write(reinterpret_cast<const char*>(owned.external_edge.data()),
                owned.external_edge.size());
            flags.close();compact_p3::require(bool(flags),
                "complete original exterior edge flags");
        }
        std::cout<<"SOURCE_EDGE_PRESENCE_PASS rank="<<rank
            <<" global_edges="<<owned.global_edges
            <<" local_unique="<<local_count
            <<" triangle_incidences="<<owned.triangle_edge.size()
            <<" global_covered="<<global_count
            <<" source_epoch="<<source.accepted.physical.epoch
            <<" source_hash="<<source.accepted.coefficient_hash
            <<std::endl;
        MPI_Barrier(MPI_COMM_WORLD);
        }
        MPI_Finalize();return 0;
    }catch(const std::exception&e){
        std::cerr<<"rank"<<rank<<" SOURCE_EDGE_PRESENCE_FAIL "<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD,2);return 2;
    }
}
