// Actual step5 compact-transaction setup through source context, target
// registry and immutable screened weight law. Read only; no remap or step.
#include "../material_compact_rezone.cpp"
#include <iomanip>
#include <iostream>
#include <unistd.h>

static void save_geometry(const compact_source::GeometryBound&geometry,
        int rank,const std::string&directory){
    static_assert(std::is_trivially_copyable<shared_function::RealBounds>::value,
        "fixed local source-bound binary payload");
    std::string path=directory+"/geometry_r"+std::to_string(rank)+".bin";
    FILE*file=fopen(path.c_str(),"wbx");
    compact_p3::require(file,"new actual source-bound output path");
    auto write=[&](const void*data,size_t size){
        compact_p3::require(fwrite(data,1,size,file)==size,
            "complete actual source geometry witness write");
    };
    uint64_t header[9]={0x4c46414354534731ULL,1,uint64_t(rank),4,
        uint64_t(geometry.triangle.size()),sizeof(shared_function::RealBounds),
        geometry.source_coefficient_hash,geometry.source_epoch,geometry.source_steps};
    compact_p3::require(geometry.canonical_global_injective&&
        geometry.triangle.size()==geometry.complete.size(),
        "persist only complete certified same-epoch source geometry");
    write(header,sizeof(header));write(&geometry.source_time,sizeof(geometry.source_time));
    write(&geometry.maximum_Dchi_minus_I,sizeof(geometry.maximum_Dchi_minus_I));
    write(&geometry.minimum_canonical_J,sizeof(geometry.minimum_canonical_J));
    write(geometry.triangle.data(),geometry.triangle.size()*sizeof(shared_function::RealBounds));
    write(geometry.complete.data(),geometry.complete.size());
    int flushed=fflush(file),synced=fsync(fileno(file)),closed=fclose(file);
    compact_p3::require(flushed==0&&synced==0&&closed==0,
        "durable actual source geometry witness");
}

int main(int argc,char**argv){
    int provided=0;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
    int rank=0,ranks=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    try{
        compact_p3::require(argc==8&&ranks==4&&provided>=MPI_THREAD_FUNNELED,
            "actual step5 setup needs six sealed inputs, new output and four MPI owners");
        compact_source::Inputs input;
        input.initial=argv[1];input.immutable_manifest=argv[2];
        input.positive_directory=argv[3];input.whole_audit=argv[4];
        input.reference_receipt=argv[5];input.accepted_prefix=argv[6];
        input.accepted_steps=5;input.accepted_time=0.00011263377801805954;
        const double begin=MPI_Wtime();
        {
            auto source=compact_source::import_accepted_state(input,MPI_COMM_WORLD);
            if(rank==0)std::cout<<"ACTUAL_STEP5_IMPORT_PASS seconds="
                <<MPI_Wtime()-begin<<std::endl;
            auto geometry=compact_source::certify_accepted_geometry(source);
            compact_p3::require(geometry.canonical_global_injective,
                "actual accepted step5 current source geometry prerequisite");
            if(rank==0)std::cout<<std::setprecision(19)
                <<"ACTUAL_STEP5_GEOMETRY_PASS max_Dchi_minus_I="
                <<geometry.maximum_Dchi_minus_I<<" min_J="
                <<geometry.minimum_canonical_J<<" seconds="
                <<MPI_Wtime()-begin<<std::endl;
            save_geometry(geometry,rank,argv[7]);
            auto context=compact_source::derive_accepted_context(source,input,geometry);
            if(rank==0)std::cout<<std::setprecision(19)
                <<"ACTUAL_STEP5_CONTEXT_PASS mass_relative="<<context.mass_difference
                <<" U_over_U="<<context.U_difference
                <<" Phi_over_U="<<context.Phi_difference
                <<" K_over_U="<<context.K_difference
                <<" L_over_2K="<<context.L_difference
                <<" corrected_J_lower="<<context.conditional_corrected_J_lower
                <<" seconds="<<MPI_Wtime()-begin<<std::endl;
            auto target=compact_source::extract_target_inputs(source,input,geometry,context);
            auto registry=compact_target::gather(target,MPI_COMM_WORLD);
            if(rank==0)std::cout<<std::setprecision(19)
                <<"ACTUAL_STEP5_TARGET_REGISTRY_PASS sites="<<registry.site.size()
                <<" exterior="<<registry.exterior.size()
                <<" max_position_error="<<target.maximum_position_error
                <<" seconds="<<MPI_Wtime()-begin<<std::endl;
            auto law=compact_mass_law::bind_actual(input,source,geometry,context,
                registry,MPI_COMM_WORLD);
            if(rank==0)std::cout<<std::setprecision(19)
                <<"ACTUAL_STEP5_FIXED_MASS_LAW_PASS entries="<<law.entries().size()
                <<" total_fixed_mass="<<law.total_fixed_mass()
                <<" total_mass_relative_difference="
                <<law.total_mass_relative_difference()
                <<" seconds="<<MPI_Wtime()-begin<<std::endl;
        }
        MPI_Finalize();return 0;
    }catch(const std::exception&e){
        std::cerr<<"rank"<<rank<<" ACTUAL_STEP5_SETUP_FAIL "<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD,2);return 2;
    }
}
