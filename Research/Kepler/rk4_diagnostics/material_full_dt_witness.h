// Optional provenance for the EXISTING material dt reduction. This does not
// choose a step, loosen a bound, inspect new points or claim a spectral CFL.
// Included inside full_evolution after Material/Evaluation declarations.
static void observe_dt_witness(const Material&m,Evaluation&e,const MFEvolvedBlock&b,size_t s){
    if(!m.save_dt_witness)return;
    R candidate[3]={R(m.eta)/std::max(R(b.max_grad_v),R(1e-30)),R(m.eta)/std::max(R(b.max_omega),R(1e-30)),.2L*b.min_length/(7*std::max(R(b.max_sound),R(1e-30)))};
    int branch=0;for(int j=1;j<3;j++)if(candidate[j]<candidate[branch])branch=j;
    uint64_t source=m.owned.source_triangle[s];
    if(candidate[branch]<e.dt_witness.bound||(candidate[branch]==e.dt_witness.bound&&(source<e.dt_witness.source||(source==e.dt_witness.source&&branch<e.dt_witness.branch)))){
        auto&w=e.dt_witness;w.bound=candidate[branch];w.source=source;w.cell=m.owned.source_cell[s];w.local=s;w.branch=branch;
        for(int j=0;j<3;j++)w.candidate[j]=candidate[j];w.max_grad=b.max_grad_v;w.max_omega=b.max_omega;w.min_length=b.min_length;w.max_sound=b.max_sound;
    }
}
static void report_dt_witness(const Material&m,const md::Ownership&layout,const Evaluation&e,const std::string&path){
    if(!m.save_dt_witness)return;
    // e.dt_bound has ALREADY been reduced by the unchanged controller code.
    // Resolve exact ties by canonical source ID, then unique owner rank.
    uint64_t source=e.dt_witness.bound==e.dt_bound?e.dt_witness.source:UINT64_MAX;
    MPI_Allreduce(MPI_IN_PLACE,&source,1,MPI_UINT64_T,MPI_MIN,layout.comm);
    int winner=source==e.dt_witness.source&&e.dt_witness.bound==e.dt_bound?layout.rank:INT_MAX;
    MPI_Allreduce(MPI_IN_PLACE,&winner,1,MPI_INT,MPI_MIN,layout.comm);
    need(winner>=0&&winner<layout.ranks&&source!=UINT64_MAX,"existing dt minimum has an actual recorded source");
    R data[9]={};uint64_t ids[3]={};
    if(layout.rank==winner){const auto&w=e.dt_witness;data[0]=w.bound;for(int j=0;j<3;j++)data[j+1]=w.candidate[j];data[4]=w.max_grad;data[5]=w.max_omega;data[6]=w.min_length;data[7]=w.max_sound;
        const auto&t=m.owned.f.triangle[w.local];R x=0,y=0;for(int j=0;j<3;j++){x+=R(t.vertex[j][0])/3;y+=R(t.vertex[j][1])/3;}data[8]=hypotl(x,y);ids[0]=w.cell;ids[1]=w.local;ids[2]=w.branch;
    }
    MPI_Bcast(data,9,MPI_LONG_DOUBLE,winner,layout.comm);MPI_Bcast(ids,3,MPI_UINT64_T,winner,layout.comm);
    need(data[0]==e.dt_bound&&ids[2]<3,"recorded dt provenance equals the existing global reduction");
    int io_bad=0;if(layout.rank==0){FILE*f=fopen((path+"/material_dt_witness.json").c_str(),"wx");if(!f)io_bad=1;else{
        const char*branch[3]={"full_velocity_gradient_Frobenius","central_gravity_Omega","acoustic_block_min_length_over_max_sound"};
        fprintf(f,"{\"scope\":\"argmin provenance of unchanged sampled material dt controller; no new point queries or spectral stability claim\",\"base_accepted_time\":%.17g,\"configured_stage\":%llu,\"coefficient_epoch\":%llu,\"h\":%.17g,\"source_triangle\":%llu,\"source_cell\":%llu,\"source_local_triangle\":%llu,\"owner_rank\":%d,\"reference_centroid_radius\":%.21Lg,\"limiting_branch\":\"%s\",\"existing_global_dt_bound\":%.21Lg,\"gradient_candidate\":%.21Lg,\"gravity_candidate\":%.21Lg,\"acoustic_candidate\":%.21Lg,\"block_max_grad_v\":%.21Lg,\"block_max_Omega\":%.21Lg,\"block_min_material_length\":%.21Lg,\"block_max_sound_speed\":%.21Lg,\"selected_positive_rule_active\":%s,\"rule_generation\":%llu,\"extremum_point_recorded\":false,\"same_single_point_for_acoustic_extrema_claimed\":false,\"dt_or_state_modified_by_diagnostic\":false,\"spectral_enriched_CFL_proven\":false}\n",m.ledger.time,(unsigned long long)m.ledger.stage,(unsigned long long)m.ledger.epoch,m.ledger.h,(unsigned long long)source,(unsigned long long)ids[0],(unsigned long long)ids[1],winner,data[8],branch[ids[2]],data[0],data[1],data[2],data[3],data[4],data[5],data[6],data[7],m.positive.active?"true":"false",(unsigned long long)m.positive.generation);
        if(fflush(f)!=0||fsync(fileno(f))!=0)io_bad=1;if(fclose(f)!=0)io_bad=1;
    }}MPI_Bcast(&io_bad,1,MPI_INT,0,layout.comm);need(!io_bad,"durable native dt minimum provenance report");
}
