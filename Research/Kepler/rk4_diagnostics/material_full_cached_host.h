// Stage-local cache tiles contain complete source cells. Every dependent pass
// keeps the accepted quadrature, reduction order, generating field and EOS.
struct PointTile {size_t first=0,count=0,cell_first=0,cell_end=0;std::vector<uint32_t>prefix{0};};
static std::vector<PointTile> point_tiles(const Material&m,const GfsEnrichConfig&cfg){
    const auto&o=m.owned;std::vector<PointTile>out;size_t at=0,last_cell=0;
    while(at<o.f.triangle.size()){
        PointTile tile;tile.first=at;tile.cell_first=o.f.triangle[at].cell;
        while(at<o.f.triangle.size()){
            size_t end=at,cell=o.f.triangle[at].cell;uint64_t points=0;
            need(cell>=last_cell,"source cells contiguous in retained triangle order");
            std::vector<uint32_t>counts;
            while(end<o.f.triangle.size()&&size_t(o.f.triangle[end].cell)==cell){
                const auto&t=o.f.triangle[end];const auto&r=o.rule[end];int active=0;
                for(int k=0;k<3;k++){int j=(k+1)%3;double x=t.vertex[k][0]-r.anchor[0],y=t.vertex[k][1]-r.anchor[1],u=t.vertex[j][0]-r.anchor[0],v=t.vertex[j][1]-r.anchor[1],jac=x*v-y*u;if(k==r.anchor_edge||(x==0&&y==0)||(u==0&&v==0))jac=0;need(jac>=0,"unchanged cached positive rule");if(jac>0)active++;}
                uint64_t n=uint64_t(active)*r.ns*cfg.na*cfg.nr;
                if(m.positive.active){need(m.positive.count.size()==o.f.triangle.size(),"complete loaded source point counts");n=m.positive.count[end];}
                need(n>0&&n<=UINT32_MAX,"bounded selected positive point prefix");counts.push_back(n);points+=n;end++;
            }
            need(end-at<=MF_BATCH_TRIANGLES&&points<=m.cache_capacity,"one whole source cell fits declared cache limit");
            if(tile.count&&(tile.count+end-at>MF_BATCH_TRIANGLES||uint64_t(tile.prefix.back())+points>m.cache_capacity))break;
            for(auto n:counts)tile.prefix.push_back(tile.prefix.back()+n);
            tile.count+=end-at;tile.cell_end=cell+1;last_cell=cell;at=end;
        }
        need(tile.count>0&&tile.prefix.size()==tile.count+1,"nonempty complete-source point tile");out.push_back(std::move(tile));
    }
    return out;
}
static void cached_passes(Material&m,const md::Ownership&layout,Device&gpu,const GfsEnrichConfig&cfg,
        Evaluation&e,bool forces,bool cpu_shadow,const std::vector<size_t>&witness,FILE*log,const std::string&path){
    auto&o=m.owned;auto tiles=point_tiles(m,cfg);uint64_t nt=tiles.size(),maxnt;
    MPI_Allreduce(&nt,&maxnt,1,MPI_UINT64_T,MPI_MAX,layout.comm);
    int bad=mf_evolved_stage(gpu.ptr,o.f.coefficient.data(),o.f.coefficient.size(),e.context.data(),e.context.size())!=0,any;
    MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,layout.comm);need(!any,"cached stage coefficient/context upload");
    std::vector<MFEvolvedBlock>batch(MF_BATCH_TRIANGLES),reference(MF_BATCH_TRIANGLES);
    uint64_t positive[3]={};R comparison=0;
    for(uint64_t tile_index=0;tile_index<maxnt;tile_index++)for(int pass=0;pass<(forces?3:2);pass++){
        bad=0;
        try{if(tile_index<tiles.size()){
            const auto&tile=tiles[tile_index];double seconds=0;
            int code=pass==0?mf_evolved_cache_build(gpu.ptr,tile.first,tile.count,tile.prefix.data(),batch.data(),&seconds):mf_evolved_cache_batch(gpu.ptr,pass,tile.first,tile.count,batch.data(),&seconds);
            need(!code,"stage-local native cached point pass");e.kernel[pass]+=seconds;
            // Only one tile of the first new RHS compares the original stream
            // against its cache at this SAME state. No accepted interval replay.
            if(!m.cache_checked&&tile_index==0){
                double ignored;need(!mf_evolved_batch(gpu.ptr,pass,tile.first,tile.count,reference.data(),&ignored),"same-stage original/cache witness");
                for(size_t j=0;j<tile.count;j++){
                    const auto&a=batch[j];const auto&b=reference[j];
                    need(a.bad_point==b.bad_point&&a.points==b.points,"same positive points in original/cache witness");
                    for(int k=0;k<MF_EVOL_METRICS;k++){R x=value(a.metric[k]),y=value(b.metric[k]);comparison=std::max(comparison,fabsl(x-y)/std::max({fabsl(x),fabsl(y),R(1e-28)}));}
                    if(pass==2){R scale=0,difference=0;for(int d=0;d<2;d++)for(int k=0;k<o.f.triangle[tile.first+j].n;k++){R x=value(a.force[d][k]),y=value(b.force[d][k]);scale=std::max(scale,fabsl(y));difference=std::max(difference,fabsl(x-y));}need(scale>0,"nonzero cached force witness scale");comparison=std::max(comparison,difference/scale);}
                    need(a.minJ==b.minJ&&(pass==0||(a.min_j==b.min_j&&a.min_U==b.min_U)),"original/cache positivity values unchanged");
                }
            }
            for(size_t c=tile.cell_first;c<tile.cell_end;c++)e.cell[c]={};
            for(size_t j=0;j<tile.count;j++){
                const auto&b=batch[j];size_t s=tile.first+j,c=o.f.triangle[s].cell;
                if(b.bad_point>=0){FILE*f=fopen((path+"/failed_point_r"+std::to_string(layout.rank)+".json").c_str(),"wx");if(f){fprintf(f,"{\"epoch\":%llu,\"stage\":%llu,\"pass\":%d,\"source_triangle\":%llu,\"cell\":%llu,\"compact_point\":%d,\"kind\":%d}\n",(unsigned long long)m.ledger.epoch,(unsigned long long)m.ledger.stage,pass,(unsigned long long)o.source_triangle[s],(unsigned long long)o.source_cell[s],b.bad_point,b.bad_kind);material_full::durable(f);}need(false,"cached actual changed-state failure witness retained");}
                positive[pass]+=b.points;
                for(int k=0;k<MF_EVOL_METRICS;k++){R v=value(b.metric[k]);need(std::isfinite(v),"finite cached source metric");e.cell[c][k]+=v;}
                e.minJ=std::min(e.minJ,R(b.minJ));e.maxJ=std::max(e.maxJ,R(b.max_J));e.maxGrad=std::max(e.maxGrad,R(b.max_grad_v));
                if(pass>0){e.minj=std::min(e.minj,R(b.min_j));e.minU=std::min(e.minU,R(b.min_U));need(b.min_length>0&&std::isfinite(b.max_grad_v)&&std::isfinite(b.max_sound),"cached material CFL bounds");e.dt_bound=std::min(e.dt_bound,std::min({R(m.eta)/std::max(R(b.max_grad_v),R(1e-30)),R(m.eta)/std::max(R(b.max_omega),R(1e-30)),.2L*b.min_length/(7*std::max(R(b.max_sound),R(1e-30)))}));observe_dt_witness(m,e,b,s);}
                if(pass==2){for(int d=0;d<2;d++)for(int k=0;k<o.f.triangle[s].n;k++){R v=value(b.force[d][k]);need(std::isfinite(v),"finite cached stress/gravity force");e.element[s][d*ST_MAX+k]=v;}if(cpu_shadow&&std::binary_search(witness.begin(),witness.end(),s)){auto diff=shadow(m,s,cfg,e.context[c],b);e.shadow_force=std::max(e.shadow_force,double(diff[0]));e.shadow_U=std::max(e.shadow_U,double(diff[1]));need(diff[0]<1e-11L&&diff[1]<1e-11L,"cached native CPU/CUDA force/heat witness");}}
            }
            for(size_t c=tile.cell_first;c<tile.cell_end;c++){auto&p=e.cell[c];need(p[0]>0&&p[3]>0,"positive cached source measure/mass");if(pass==0){e.context[c].delta=double((e.geometry.volume[c]-p[1])/p[0]);e.context[c].delta_dot=double((e.geometry.vdot[c]-p[2])/p[0]);}else if(pass==1)e.context[c].mean_deltaP=double(p[4]/p[0]);}
            if(pass<2)need(!mf_evolved_context_range(gpu.ptr,e.context.data(),tile.cell_first,tile.cell_end-tile.cell_first),"cache tile dependent source context upload");
            fprintf(log,"%d\t%llu\t%.9g\n",pass,(unsigned long long)(tile.first+tile.count),e.kernel[pass]);fflush(log);
        }}catch(const std::exception&x){fprintf(stderr,"rank%d CACHED_PASS_FAIL %s\n",layout.rank,x.what());bad=1;}
        MPI_Allreduce(&bad,&any,1,MPI_INT,MPI_MAX,layout.comm);need(!any,"collective cached current-state rejection");
    }
    MPI_Allreduce(MPI_IN_PLACE,positive,3,MPI_UINT64_T,MPI_SUM,layout.comm);
    uint64_t expected=m.positive.active?m.positive.global_points:2040256512ULL;
    for(int pass=0;pass<(forces?3:2);pass++)need(positive[pass]==expected,"every selected positive point accounted for, including full-rule fallbacks");
    MPI_Allreduce(MPI_IN_PLACE,&comparison,1,MPI_LONG_DOUBLE,MPI_MAX,layout.comm);need(comparison<1e-11L,"original/cache same-stage reduced operator agreement");
    if(!m.cache_checked){if(!layout.rank){FILE*f=fopen((path+"/cache_equivalence.json").c_str(),"wx");need(f,"actual new RHS cache equivalence report");fprintf(f,"{\"same_new_RHS_first_complete_source_tile_per_rank\":true,\"maximum_relative_reduced_difference\":%.21Lg,\"same_selected_rule_in_both_paths\":true,\"cache_point_capacity_per_rank\":%zu,\"full_point_array_stored\":false,\"pass\":true}\n",comparison,m.cache_capacity);material_full::durable(f);}m.cache_checked=true;}
    e.cache_difference=double(comparison);
}
