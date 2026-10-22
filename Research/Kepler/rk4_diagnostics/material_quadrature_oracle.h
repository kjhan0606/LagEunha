#ifndef MATERIAL_QUADRATURE_ORACLE_H
#define MATERIAL_QUADRATURE_ORACLE_H
// Same-stage reference comparison. Only the selected low-rule rate advances
// RK4. The high-rule k1 is restored; missing source summaries need no force/CG.
static Rate restore_high_k1(Material&m,const md::Ownership&layout,const std::string&directory){
    std::string report=textfile(directory+"/../next_RHS_result.json");need(report.find("\"actual_next_high_rule_RHS\":true")!=std::string::npos&&report.find("\"prepared_k1_reusable_for_same_high_rule\":true")!=std::string::npos,"accepted prepared high-rule k1 contract");
    FILE*f=fopen((directory+"/prepared_k1_r"+std::to_string(layout.rank)+".bin").c_str(),"rb");need(f,"actual retained high k1");uint64_t h[8];material_full::read(f,h,8);
    need(h[0]==0x4d464e4558544b31ULL&&h[1]==1&&h[2]==uint64_t(layout.rank)&&h[3]==uint64_t(layout.ranks)&&h[4]==m.owned.global.size()&&h[5]==sizeof(STCoefficient)&&h[6]==hash(m.owned.f.coefficient.data(),m.owned.f.coefficient.size()*sizeof(STCoefficient))&&h[7]==hash(m.immutable.data(),m.immutable.size()),"exact prepared state, owner and immutable thermal binding");
    Rate r;double time;material_full::read(f,&time);material_full::read(f,&r.dt_bound);need(time==m.ledger.time,"prepared rate time matches same-stage anchor");R rates[3];material_full::read(f,rates,3);r.work=rates[0];r.torque=rates[1];r.heat=rates[2];r.k.resize(h[4]);material_full::read(f,r.k.data(),r.k.size());need(fgetc(f)==EOF&&fclose(f)==0,"complete prepared k1 payload");
    for(size_t i=0;i<r.k.size();i++)for(int d=0;d<2;d++)need(memcmp(&r.k[i].x[d],&m.owned.f.coefficient[i].v[d],sizeof(double))==0&&std::isfinite(r.k[i].v[d]),"prepared xdot/current velocity bits and finite actual acceleration");return r;
}
static void attach_metadata(Rate&r,Evaluation&&e,const Endpoint&metric){
    r.endpoint=e.endpoint;r.endpoint.K=metric.K;r.endpoint.L=metric.L;r.minj=e.minj;r.maxJ=e.maxJ;r.maxGrad=e.maxGrad;r.dt_bound=std::min(r.dt_bound,e.dt_bound);
    r.source=std::move(e.cell);r.context=std::move(e.context);r.volume=std::move(e.geometry.volume);r.vdot=std::move(e.geometry.vdot);
}
struct QuadratureDifference {R pressure=0,density=0,pressure_rate=0,density_rate=0,velocity_rate=0,U=0,Phi=0;std::vector<std::array<R,20>>source;};
static QuadratureDifference compare_quadrature(const Material&m,const md::Mass&M,const Rate&low,const Rate&high,bool acceleration){
    const auto&layout=M.layout;const auto&o=m.owned;need(low.source.size()==o.f.cell.size()&&high.source.size()==low.source.size()&&low.context.size()==low.source.size()&&high.context.size()==high.source.size()&&low.volume.size()==low.source.size()&&high.volume.size()==high.source.size()&&high.vdot.size()==high.source.size()&&high.minj>0&&std::isfinite(high.maxJ)&&std::isfinite(high.maxGrad),"complete same-stage source/geometry summaries");QuadratureDifference q;q.source.resize(o.f.cell.size());
    for(size_t c=0;c<o.f.cell.size();c++){
        const auto&l=low.source[c];const auto&h=high.source[c];need(h[5]>0&&h[3]>0&&high.volume[c]>0,"positive source reference heat/mass/volume");
        R geometry_error=fabsl(low.volume[c]-high.volume[c])/high.volume[c];need(geometry_error<1e-14L,"same stage and canonical source boundary");
        R delta=R(low.context[c].delta)-high.context[c].delta,delta_dot=R(low.context[c].delta_dot)-high.context[c].delta_dot,jlow_min=high.minj-fabsl(delta);need(jlow_min>0,"low corrected map positive on every original high-rule node");
        R density_bound=fabsl(delta)/jlow_min,pressure_bound=expm1l(m.gamma*log1pl(density_bound));
        R pressure=std::max(fabsl(l[5]-h[5])/h[5],pressure_bound),density=std::max(fabsl(l[3]/h[3]-1),density_bound);
        R Vrate=high.vdot[c]/high.volume[c],pressure_rate=fabsl((l[6]-h[6])-(l[5]-h[5])*Vrate)/h[5];
        R nudot_bound=high.maxJ*sqrtl(2.L)*high.maxGrad+fabsl(high.context[c].delta_dot),log_density_rate=fabsl(delta_dot)/jlow_min+nudot_bound*fabsl(delta)/(high.minj*jlow_min);
        R high_log_rate=nudot_bound/high.minj,density_rate=density_bound*high_log_rate+(1+density_bound)*log_density_rate;
        pressure_rate=std::max(pressure_rate,R(m.gamma)*(pressure_bound*high_log_rate+(1+pressure_bound)*log_density_rate));
        q.pressure=std::max(q.pressure,pressure);q.density=std::max(q.density,density);q.pressure_rate=std::max(q.pressure_rate,pressure_rate);q.density_rate=std::max(q.density_rate,density_rate);q.U+=l[5]-h[5];q.Phi+=l[10]-h[10];
        q.source[c]={l[5],h[5],l[10],h[10],(m.gamma-1)*l[5]/high.volume[c],(m.gamma-1)*h[5]/high.volume[c],l[3]/high.volume[c],h[3]/high.volume[c],pressure,density,pressure_rate,density_rate,0,delta,delta_dot,high.volume[c],high.vdot[c],R(low.context[c].mean_deltaP),R(high.context[c].mean_deltaP),pressure_bound};
    }
    if(acceleration){
        need(low.k.size()==o.global.size()&&high.k.size()==low.k.size(),"same-stage native acceleration extents");std::vector<std::array<R,3>>element(o.f.triangle.size());int bad=0;
        #pragma omp parallel for schedule(static) reduction(+:bad)
        for(size_t s=0;s<o.f.triangle.size();s++){
            const auto&t=o.f.triangle[s];const R*mass=M.block.data()+s*ST_MAX*ST_MAX;R energy=0,absolute=0,K=0;
            for(int a=0;a<t.n;a++)for(int b=0;b<t.n;b++)for(int d=0;d<2;d++){R da=R(low.k[t.id[a]].v[d])-high.k[t.id[a]].v[d],db=R(low.k[t.id[b]].v[d])-high.k[t.id[b]].v[d],term=da*mass[a*ST_MAX+b]*db;energy+=term;absolute+=fabsl(term);K+=R(o.f.coefficient[t.id[a]].v[d])*mass[a*ST_MAX+b]*o.f.coefficient[t.id[b]].v[d];}
            R uncertainty=64*t.n*LDBL_EPSILON*absolute;if(energy< -uncertainty||!(K>0))bad++;element[s]={energy+uncertainty,uncertainty,K};
        }
        MPI_Allreduce(MPI_IN_PLACE,&bad,1,MPI_INT,MPI_SUM,layout.comm);need(!bad,"positive authoritative mass norm for every acceleration-difference element");
        std::vector<std::array<R,2>>cell(o.f.cell.size());for(size_t s=0;s<element.size();s++){auto&sum=cell[o.f.triangle[s].cell];sum[0]+=element[s][0];sum[1]+=element[s][2];}for(size_t c=0;c<cell.size();c++){R rate=sqrtl(cell[c][0]/cell[c][1]);q.source[c][12]=rate;q.velocity_rate=std::max(q.velocity_rate,rate);}
    }
    R maximum[5]={q.pressure,q.density,q.pressure_rate,q.density_rate,q.velocity_rate};MPI_Allreduce(MPI_IN_PLACE,maximum,5,MPI_LONG_DOUBLE,MPI_MAX,layout.comm);q.pressure=maximum[0];q.density=maximum[1];q.pressure_rate=maximum[2];q.density_rate=maximum[3];q.velocity_rate=maximum[4];R total[2]={q.U,q.Phi};MPI_Allreduce(MPI_IN_PLACE,total,2,MPI_LONG_DOUBLE,MPI_SUM,layout.comm);q.U=total[0];q.Phi=total[1];return q;
}
static void save_quadrature_cells(const Material&m,const md::Ownership&layout,const std::string&path,const QuadratureDifference&q,bool acceleration){
    FILE*out=fopen((path+"/quadrature_cells_r"+std::to_string(layout.rank)+".csv").c_str(),"wx");need(out,"retain every actual source quadrature discrepancy");
    fprintf(out,"cell,U_low,U_high,Phi_low,Phi_high,area_mean_P_low,area_mean_P_high,mean_rho_low,mean_rho_high,P_relative_bound,rho_relative_bound,P_relative_rate_difference_bound,rho_relative_rate_difference_bound,M_acceleration_difference_over_velocity,delta_difference,delta_dot_difference,volume,volume_dot,mean_deltaP_low,mean_deltaP_high,original_node_P_relative_bound,acceleration_compared\n");
    need(q.source.size()==m.owned.f.cell.size(),"complete ordered quadrature source output");for(size_t c=0;c<q.source.size();c++){fprintf(out,"%llu",(unsigned long long)m.owned.f.cell[c]);for(R value:q.source[c]){need(std::isfinite(value),"finite per-source quadrature discrepancy");fprintf(out,",%.21Lg",value);}fprintf(out,",%d\n",acceleration?1:0);}material_full::durable(out);
}
static Rate paired_rhs(Material&m,const md::Mass&M,Device&gpu,const GfsEnrichConfig&cfg,const std::string&path,Endpoint&endpoint,R&next_dt,full_gpu_mass::Operator*mass){
    const auto&layout=M.layout;need(m.quadrature_oracle&&m.positive.loaded,"explicit paired-rule comparison mode");
    std::string low_path=path+"/low",high_path=path+"/high";mkdir_collective(low_path,layout);mkdir_collective(high_path,layout);
    use_positive_rule(m,layout,gpu,true);R low_dt=1e99L;auto low=rhs(m,M,gpu,cfg,false,low_path,endpoint,low_dt,mass,false);
    double oracle_begin=MPI_Wtime();use_positive_rule(m,layout,gpu,false);Rate high;Endpoint reference;R high_dt=1e99L;
    const char*prepared=getenv("GFS_MATERIAL_PREPARED_NEXT_RHS");bool reuse=m.ledger.stage==1&&m.ledger.time==m.positive.anchor_time&&prepared&&*prepared;
    if(reuse){high=restore_high_k1(m,layout,prepared);auto summaries=evaluate(m,layout,gpu,cfg,false,false,high_path);attach_metadata(high,std::move(summaries),m.current);high_dt=high.dt_bound;}
    else high=rhs(m,M,gpu,cfg,false,high_path,reference,high_dt,mass,false);
    auto difference=compare_quadrature(m,M,low,high,true);save_quadrature_cells(m,layout,path,difference,true);
    // The map/velocity are identical. Only corrected j changes the sound
    // speed, by (j_high/j_low)^((gamma-1)/2). This conservative rescaling
    // keeps the ORIGINAL-node timestep bound for the low functional too.
    high_dt/=powl(1+difference.density,(m.gamma-1)/2);next_dt=std::min({next_dt,low_dt,high_dt});
    if(m.ledger.stage==1&&!m.numerical.reduced){
        m.numerical.U_bridge+=difference.U;m.numerical.Phi_bridge+=difference.Phi;m.numerical.mass_bridge+=low.endpoint.mass-high.endpoint.mass;
        m.numerical.pressure_error+=difference.pressure;m.numerical.density_error+=difference.density;
    }
    const R weights[4]={1.L/6,1.L/3,1.L/3,1.L/6};need(m.ledger.stage>=1&&m.ledger.stage<=4,"actual RK4 stage weight");R duration=R(m.ledger.h)*weights[m.ledger.stage-1];
    m.numerical.pressure_error+=duration*difference.pressure_rate;m.numerical.density_error+=duration*difference.density_rate;m.numerical.velocity_error+=duration*difference.velocity_rate;
    R allowance=1e-3L*std::min(10.L,R(m.ledger.time+m.ledger.h-m.positive.anchor_time))/10;
    bool pass=R(m.ledger.h)<=high_dt*(1+1e-12L)&&difference.pressure<1e-3L&&difference.density<1e-3L&&difference.pressure_rate<=1e-4L&&difference.density_rate<=1e-4L&&difference.velocity_rate<=1e-4L&&m.numerical.pressure_error<=allowance&&m.numerical.density_error<=allowance&&m.numerical.velocity_error<=allowance;
    use_positive_rule(m,layout,gpu,true);
    if(!layout.rank){FILE*out=fopen((path+"/quadrature_difference.json").c_str(),"wx");need(out,"actual same-stage full-domain quadrature comparison");fprintf(out,"{\"stage\":%llu,\"same_positions_velocities\":true,\"high_k1_force_and_solve_reused\":%s,\"high_missing_source_summaries_only\":%s,\"max_source_pressure_relative\":%.21Lg,\"max_point_density_relative_bound\":%.21Lg,\"max_pressure_relative_rate_difference\":%.21Lg,\"max_density_relative_rate_difference_bound\":%.21Lg,\"max_source_mass_metric_acceleration_difference_over_velocity\":%.21Lg,\"delta_U\":%.21Lg,\"delta_Phi\":%.21Lg,\"error_allowance_at_step_end\":%.21Lg,\"cumulative_pressure_error\":%.21Lg,\"cumulative_density_error\":%.21Lg,\"cumulative_velocity_error\":%.21Lg,\"same_stage_quadrature_gate_pass\":%s,\"no_physical_state_correction\":true}\n",(unsigned long long)m.ledger.stage,reuse?"true":"false",reuse?"true":"false",difference.pressure,difference.density,difference.pressure_rate,difference.density_rate,difference.velocity_rate,difference.U,difference.Phi,allowance,m.numerical.pressure_error,m.numerical.density_error,m.numerical.velocity_error,pass?"true":"false");material_full::durable(out);}
    m.oracle_wall+=MPI_Wtime()-oracle_begin;need(pass,"same-stage pressure/density/orbit quadrature budget and original timestep gate");endpoint=low.endpoint;return low;
}
static void paired_endpoint(Material&m,const md::Mass&M,Device&gpu,const GfsEnrichConfig&cfg,const std::string&path,Evaluation&state,R&next_dt){
    double oracle_begin=MPI_Wtime();const auto&layout=M.layout;need(m.positive.active&&m.quadrature_oracle,"reduced endpoint only after four actual paired stages");
    Rate low;low.endpoint=state.endpoint;low.minj=state.minj;low.maxJ=state.maxJ;low.maxGrad=state.maxGrad;low.dt_bound=state.dt_bound;
    low.source=state.cell;low.context=state.context;low.volume=state.geometry.volume;low.vdot=state.geometry.vdot;
    std::string high_path=path+"/high";mkdir_collective(high_path,layout);use_positive_rule(m,layout,gpu,false);auto high_state=evaluate(m,layout,gpu,cfg,false,false,high_path);R high_minJ=high_state.minJ;Rate high;attach_metadata(high,std::move(high_state),state.endpoint);
    auto difference=compare_quadrature(m,M,low,high,false);save_quadrature_cells(m,layout,path,difference,false);R original_bound=high.dt_bound/powl(1+difference.density,(m.gamma-1)/2);next_dt=std::min(next_dt,original_bound);state.dt_bound=std::min(state.dt_bound,original_bound);
    R max_delta=0;for(size_t c=0;c<low.context.size();c++)max_delta=std::max(max_delta,fabsl(R(low.context[c].delta)-high.context[c].delta));state.minJ=std::min(state.minJ,high_minJ);state.minj=std::min(state.minj,high.minj-max_delta);
    R allowance=1e-3L*std::min(10.L,R(m.ledger.time+m.ledger.h-m.positive.anchor_time))/10;
    m.numerical.pressure_error=std::max(m.numerical.pressure_error,difference.pressure);m.numerical.density_error=std::max(m.numerical.density_error,difference.density);
    bool pass=state.minj>0&&R(m.ledger.h)<=original_bound*(1+1e-12L)&&difference.pressure<1e-3L&&difference.density<1e-3L&&m.numerical.pressure_error<=allowance&&m.numerical.density_error<=allowance&&m.numerical.velocity_error<=allowance;
    use_positive_rule(m,layout,gpu,true);
    if(!layout.rank){FILE*out=fopen((path+"/quadrature_difference.json").c_str(),"wx");need(out,"same-state endpoint and explicit functional transition ledger");fprintf(out,"{\"same_positions_velocities\":true,\"high_force_or_solve_repeated\":false,\"max_source_pressure_relative\":%.21Lg,\"max_point_density_relative_bound\":%.21Lg,\"delta_U\":%.21Lg,\"delta_Phi\":%.21Lg,\"cumulative_transition_U_difference\":%.21Lg,\"cumulative_transition_Phi_difference\":%.21Lg,\"cumulative_transition_mass_difference\":%.21Lg,\"cumulative_pressure_error\":%.21Lg,\"cumulative_density_error\":%.21Lg,\"cumulative_velocity_error\":%.21Lg,\"error_allowance_at_step_end\":%.21Lg,\"original_node_timestep_bound_with_sound_rescaling\":%.21Lg,\"original_node_corrected_J_lower_bound\":%.21Lg,\"endpoint_quadrature_gate_pass\":%s,\"transition_difference_is_physical_heat\":false,\"no_physical_state_correction\":true,\"global_map_injectivity_proven\":false}\n",difference.pressure,difference.density,difference.U,difference.Phi,m.numerical.U_bridge,m.numerical.Phi_bridge,m.numerical.mass_bridge,m.numerical.pressure_error,m.numerical.density_error,m.numerical.velocity_error,allowance,original_bound,state.minj,pass?"true":"false");material_full::durable(out);}
    m.oracle_wall+=MPI_Wtime()-oracle_begin;need(pass,"endpoint all-source pressure/density budget and original-node positivity/timestep bounds");
}
#endif
