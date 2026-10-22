#include "../material_compact_rezone_parallel.cpp"
#include <iostream>
#include <iomanip>

int main(int argc,char**argv){
    int provided=0;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
    int rank=0,ranks=0;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    try{
        compact_p3::require(argc==8&&ranks==4&&provided>=MPI_THREAD_FUNNELED,
            "four-rank accepted-source edge-box witness arguments");
        compact_source::Inputs input;
        input.initial=argv[1];input.immutable_manifest=argv[2];
        input.positive_directory=argv[3];input.whole_audit=argv[4];
        input.reference_receipt=argv[5];input.accepted_prefix=argv[6];
        input.accepted_steps=5;input.accepted_time=0.00011263377801805954;
        const uint64_t wanted=std::stoull(argv[7]);
        {
            auto source=compact_source::import_accepted_state(input,MPI_COMM_WORLD);
            const auto&o=source.owned;
            size_t count=0;
            for(size_t s=0;s<o.f.triangle.size();s++)for(int k=0;k<3;k++)
                if(o.triangle_edge[3*s+k]==wanted){
                    auto curve=compact_power::source_edge_curve(o.f,s,k,wanted,
                        source.accepted.coefficient_hash,source.accepted.physical.epoch);
                    auto whole=compact_overlap_boundary::source_edge_box(curve);
                    std::cout<<std::setprecision(22)
                        <<"EDGE_BOX_WITNESS rank="<<rank<<" edge="<<wanted
                        <<" local_triangle="<<s<<" local_edge="<<k
                        <<" source_triangle="<<o.source_triangle[s]
                        <<" original_vertices="<<o.global[o.f.triangle[s].edge[k].poly[0]]
                        <<","<<o.global[o.f.triangle[s].edge[k].poly[3]]
                        <<" width0="<<whole.hi[0]-whole.lo[0]
                        <<" width1="<<whole.hi[1]-whole.lo[1]<<std::endl;
                    for(int parts: {2,4,8,16,32,64}){
                        compact_source_spatial::Box box{};
                        bool first=true;
                        for(int j=0;j<parts;j++){
                            auto interval=shared_function::range(
                                compact_power::Q(j)/parts,
                                compact_power::Q(j+1)/parts);
                            auto x=compact_power::position(curve,interval);
                            compact_source_spatial::Box piece;
                            for(int d=0;d<2;d++){
                                piece.lo[d]=nextafterl((long double)shared_function::lower(x[d].v),-INFINITY);
                                piece.hi[d]=nextafterl((long double)shared_function::upper(x[d].v),INFINITY);
                            }
                            box=first?piece:compact_source_spatial::hull(box,piece);
                            first=false;
                        }
                        std::cout<<"EDGE_BOX_SUBDIV rank="<<rank<<" edge="<<wanted
                            <<" parts="<<parts
                            <<" width0="<<box.hi[0]-box.lo[0]
                            <<" width1="<<box.hi[1]-box.lo[1]
                            <<" reduction_area="
                            <<(whole.hi[0]-whole.lo[0])*(whole.hi[1]-whole.lo[1])/
                                ((box.hi[0]-box.lo[0])*(box.hi[1]-box.lo[1]))
                            <<std::endl;
                    }
                    count++;
                }
            std::cout<<"EDGE_BOX_LOCAL_COUNT rank="<<rank<<" edge="<<wanted
                <<" incidences="<<count<<std::endl;
        }
        MPI_Finalize();return 0;
    }catch(const std::exception&e){
        std::cerr<<"rank"<<rank<<" EDGE_BOX_WITNESS_FAIL "<<e.what()<<std::endl;
        MPI_Abort(MPI_COMM_WORLD,2);return 2;
    }
}
