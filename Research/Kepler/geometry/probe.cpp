// Actual accepted step5 source gate for the compact transaction. Read only:
// no checkpoint write, target substitution, remap or time step.
#include "../material_compact_rezone.cpp"
#include <iomanip>
#include <iostream>

int main(int argc,char**argv){
    int provided=0;
    MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
    int rank=0,ranks=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    try{
        compact_p3::require(argc==7&&ranks==4&&provided>=MPI_THREAD_FUNNELED,
            "actual step5 source gate requires six sealed inputs and four MPI owners");
        compact_source::Inputs input;
        input.initial=argv[1];input.immutable_manifest=argv[2];
        input.positive_directory=argv[3];input.whole_audit=argv[4];
        input.reference_receipt=argv[5];input.accepted_prefix=argv[6];
        input.accepted_steps=5;input.accepted_time=0.00011263377801805954;
        const double begun=MPI_Wtime();
        int status=0;
        {
            auto source=compact_source::import_accepted_state(input,MPI_COMM_WORLD);
            const double imported=MPI_Wtime();
            if(rank==0)std::cout<<std::setprecision(17)
                <<"ACTUAL_STEP5_IMPORT_PASS time="<<source.accepted.physical.time
                <<" steps="<<source.accepted.physical.steps
                <<" seconds="<<imported-begun<<std::endl;
            auto geometry=compact_source::certify_accepted_geometry(source);
            const double proved=MPI_Wtime();
            if(rank==0)std::cout<<std::setprecision(19)
                <<"ACTUAL_STEP5_GEOMETRY_RESULT complete="
                <<int(geometry.canonical_global_injective)
                <<" unresolved="<<geometry.unresolved
                <<" first_unresolved="<<geometry.first_unresolved
                <<" max_Dchi_minus_I="<<geometry.maximum_Dchi_minus_I
                <<" min_J="<<geometry.minimum_canonical_J
                <<" seconds="<<proved-imported<<std::endl;
            if(!geometry.canonical_global_injective)status=3;
        }
        MPI_Finalize();return status;
    }catch(const std::exception&e){
        std::cerr<<"rank"<<rank<<" ACTUAL_STEP5_SOURCE_GATE_FAIL "<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD,2);return 2;
    }
}
