#define GFS_SHARED_CLIPPED_LIBRARY
#include "test_shared_clipped_target_v1.cpp"
#include <atomic>

namespace {
struct FrozenPanel {uint64_t source,index;int32_t policy,dependent,lower,upper;target::ClipPanel geometry;};
struct Contribution {bool seen=false,complete=false;uint64_t points=0;std::array<bool,12>present={};S moment[12];Q remainder[12]={};};
using NativeBall=material_forward::B<long double>;
static void merge_panel(target::Reference80Accumulator&out,const target::Reference80Accumulator&panel,target::Reference80Positive&last){
    need(panel.points==256&&panel.occupied==(uint64_t(1)<<8),"one complete power-of-two panel reduction");
    NativeBall carry[12];for(int k=0;k<12;k++)carry[k]=panel.level[8][k];int level=8;
    while(level<64&&(out.occupied&(uint64_t(1)<<level))){
        for(int k=0;k<12;k++){
            carry[k]=out.level[level][k]+carry[k];
            last.sum_local_rounding[k]=sf::add(last.sum_local_rounding[k],Q(material_forward::round_error(carry[k].n)));
        }
        out.occupied&=~(uint64_t(1)<<level);level++;
    }
    need(level<64,"power-of-two source reduction capacity");
    for(int k=0;k<12;k++)out.level[level][k]=carry[k];out.occupied|=uint64_t(1)<<level;out.points+=256;
}
static void replay_moments_only(target::Reference80Accumulator&out,const target::Reference80Positive&row){
    NativeBall carry[12];for(int k=0;k<12;k++){
        long double n=static_cast<long double>(row.moment[k]),e=static_cast<long double>(row.radius[k]);
        need(Q(n)==row.moment[k]&&Q(e)==row.radius[k],"saved binary80 point value/radius recovered exactly");carry[k]=NativeBall(n,e);
    }
    int level=0;while(level<64&&(out.occupied&(uint64_t(1)<<level))){
        for(int k=0;k<12;k++)carry[k]=out.level[level][k]+carry[k];
        out.occupied&=~(uint64_t(1)<<level);level++;
    }
    need(level<64,"saved-moment reduction capacity");
    for(int k=0;k<12;k++)out.level[level][k]=carry[k];out.occupied|=uint64_t(1)<<level;out.points++;
}
static std::string csv_string(const std::string&text){std::string out="\"";for(char c:text){if(c=='\"')out+='\"';out+=c;}return out+'\"';}
}

int main(int argc,char**argv){
    int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);int code=2;
    try{
        need(argc==12&&ranks==1,"packet saved sites weights originalcall realbounds originating whole4 ORIGINAL_tau P2_proofs NEW_output; MPI1 with OpenMP panels");
        auto begin=std::chrono::steady_clock::now();auto p=curved::read(argv[1]);
        const auto f=load_saved(std::string(argv[2])+"/shared_state.bin",p);const sf::FrozenThermal heat(f,p.parameters[2],p.parameters[1],1);
        auto context=target::retained_context(argv[5],f,heat);auto source=target::retained_sources(argv[5],argv[6],f);
        auto group=target::retained_whole_group(argv[8],source);auto tau=target::original_tolerances(argv[9]);
        const std::string originating=argv[7],proofdir=argv[10],dest=argv[11];
        std::ifstream terminal(originating+"/native_exit_code.txt");int previous=-1;need(bool(terminal)&&bool(terminal>>previous)&&previous==2,"originating target terminal and rejected");
        std::ifstream proofexit(proofdir+"/native_exit_code.txt");int ready=-1;need(bool(proofexit)&&bool(proofexit>>ready)&&ready==0,"fixed P2 known budget compatible before missing quadrature");
        need(target::publication_unit_contract(),"original final R80 publication unit gates");
        std::vector<std::vector<std::string>>coverage(384);std::vector<bool>seen(384,false);std::string line;
        std::ifstream complete(originating+"/all384_source_completion.csv");need(bool(complete),"all384 saved coverage");std::getline(complete,line);std::string source_header=line;
        while(std::getline(complete,line)){
            auto c=target::columns(line);int s=std::stoi(c[0]);need(s>=0&&s<384&&!seen[s]&&c[3]=="1"&&c[4]=="1"&&c[5]=="1","exact original source coverage/atlas/analytic complete");
            need(c[6]==(s==301?"0":"1"),"only original actual source301 incomplete");seen[s]=true;coverage[s]=c;
        }
        need(std::all_of(seen.begin(),seen.end(),[](bool b){return b;}),"all384 saved source identities");
        std::vector<Contribution>value(384);std::ifstream values(originating+"/actual_clipped_source_integrals.csv");need(bool(values),"old25 actual values preserved");std::getline(values,line);size_t value_rows=0;uint64_t old_points=0;
        while(std::getline(values,line)){
            auto c=target::columns(line);int s=std::stoi(c[0]),k=-1;for(int j=0;j<12;j++)if(c[1]==target::moment_name[j])k=j;
            need(s>=0&&s<384&&k>=0,"saved actual source/moment identity");auto&r=value[s];bool done=c[2]=="1";
            need(!r.present[k]&&(!r.seen||r.complete==done),"one consistent saved all12 contribution");
            if(!r.seen){r.points=std::stoull(c[3]);old_points+=r.points;}else need(r.points==std::stoull(c[3]),"consistent saved source point count");
            r.seen=true;r.complete=done;r.present[k]=true;if(done)r.moment[k]=S(target::number(c[4]),target::number(c[5]));r.remainder[k]=target::number(c[6]);value_rows++;
        }
        need(value_rows==312&&old_points==220160&&value[301].seen&&!value[301].complete&&value[301].points==0,"retained exactly old25, no invented missing values");
        std::map<std::pair<int,int>,Q>old_R,chosen;std::array<Q,384>potential={};
        for(int owner=0;owner<4;owner++){
            std::ifstream input(originating+"/actual_clipped_panel_certificates_rank"+std::to_string(owner)+".csv");need(bool(input),"original valid panel bounds");std::getline(input,line);auto names=target::columns(line);
            auto it=std::find(names.begin(),names.end(),"potential_remainder");need(it!=names.end(),"original potential column");size_t column=it-names.begin();
            while(std::getline(input,line)){auto c=target::columns(line);need(c[8]=="1"&&old_R.emplace(std::make_pair(std::stoi(c[0]),std::stoi(c[1])),target::number(c[column])).second,"unique original valid bound");}
        }
        for(int owner=0;owner<4;owner++){
            std::ifstream input(proofdir+"/potential_P2_certificates_rank"+std::to_string(owner)+".csv");need(bool(input),"fixed P2 proof input");std::getline(input,line);auto names=target::columns(line);
            auto column=[&](const char*name){auto it=std::find(names.begin(),names.end(),name);need(it!=names.end(),"P2 schema");return size_t(it-names.begin());};
            size_t c_old=column("old_R"),c_new=column("new_R"),c_chosen=column("chosen_R"),c_valid=column("new_proof");
            while(std::getline(input,line)){
                auto c=target::columns(line);int s=std::stoi(c[0]),j=std::stoi(c[1]);Q old=old_R.at({s,j}),selected=target::number(c[c_chosen]);
                need(target::number(c[c_old])==old&&selected==(c[c_valid]=="1"?std::min(old,target::number(c[c_new])):old),"only valid minimum of original and fixed P2 proofs");
                need(chosen.emplace(std::make_pair(s,j),selected).second,"one chosen panel proof");potential[s]=sf::add(potential[s],selected);
            }
        }
        need(old_R.size()==871&&chosen.size()==871,"all871 potential proof patches retained");
        for(int s=0;s<384;s++)if(value[s].seen){need(potential[s]<=value[s].remainder[target::POTENTIAL],"only potential remainder tightened");value[s].remainder[target::POTENTIAL]=potential[s];}
        std::vector<FrozenPanel>panels;FILE*packet=fopen((originating+"/actual_reference_panel_packets_rank1.bin").c_str(),"rb");need(packet,"actual missing source geometry");uint64_t header[8];
        need(fread(header,sizeof(header),1,packet)==1&&header[0]==0x47514650414e4c31ULL&&header[1]==1&&header[2]==1&&header[3]==sizeof(FrozenPanel)&&header[5]==1,"original missing-panel ABI/epoch/owner");
        for(uint64_t i=0;i<header[4];i++){FrozenPanel row;need(fread(&row,sizeof(row),1,packet)==1,"actual panel record");if(row.source==301){need(row.index==panels.size(),"original missing source panel order");panels.push_back(row);}}
        need(fgetc(packet)==EOF,"exact typed packet length");fclose(packet);
        need(panels.size()==11&&old_points+group.integral.points+panels.size()*256<=817152,"ONLY original missing2816 within unchanged positive-node budget");
        auto m=sd::compile(f,301);auto planes=target::planes(sites(argv[3],argv[4]),68526);std::vector<target::Face>faces;for(const auto&plane:planes)faces.push_back(target::pullback(m,plane));
        std::array<target::Axis,6>axis;std::array<bool,6>axis_seen={};
        for(const auto&row:panels)if(!row.geometry.whole){int ci=row.geometry.strip.chart;if(!axis_seen[ci]){axis[ci]=target::certified_axis(m,row.geometry.reference,faces,planes,source[301].real);axis_seen[ci]=true;}need(axis[ci].dependent==row.dependent,"same frozen missing-source chart axis");}
        auto gauss=target::gauss16();std::vector<target::Reference80Accumulator>local(panels.size());std::vector<std::string>failure(panels.size());std::vector<uint64_t>fallback(panels.size(),0);std::atomic<bool>cancel(false);
        #pragma omp parallel for schedule(dynamic,1)
        for(size_t j=0;j<panels.size();j++){
            if(cancel.load()){failure[j]="not_attempted_after_actual_query_failure";continue;}
            FILE*out=fopen((dest+"/source301_panel"+std::to_string(j)+".bin").c_str(),"wb");
            uint64_t before=shared_reference80_ball::division_free_tail_calls;
            try{
                need(out,"missing-source panel query file");const auto&panel=panels[j].geometry;
                for(int iu=0;iu<16;iu++)for(int iv=0;iv<16;iv++){
                    sd::Point point;S weight;
                    if(panel.whole){S ds=panel.ordinary.shi-panel.ordinary.slo,dt=panel.ordinary.thi-panel.ordinary.tlo;point=sd::evaluate(m,panel.reference,panel.ordinary.slo+ds*gauss.x[iu],panel.ordinary.tlo+dt*gauss.x[iv]);weight=point.reference_jacobian*ds*dt*gauss.w[iu]*gauss.w[iv];}
                    else{S free=panel.strip.lo+(panel.strip.hi-panel.strip.lo)*gauss.x[iu],eta=panel.eta_lo+(panel.eta_hi-panel.eta_lo)*gauss.x[iv];auto q=target::strip_point(m,panel.reference,faces,axis[panel.strip.chart],panel.strip,free,eta,false);point=q.value;weight=q.jacobian*(panel.eta_hi-panel.eta_lo)*gauss.w[iu]*gauss.w[iv];}
                    for(const auto&plane:planes)need(sf::upper(target::value(plane,point.field.position[0],point.field.position[1]))<0,"every missing query inside every original target plane");
                    target::point_reference80(p,m,point,weight,context[f.triangle[301].cell],local[j],out);
                }
                target::finish_reference80(local[j]);
            }catch(const std::exception&e){failure[j]=e.what();cancel.store(true);}
            fallback[j]=shared_reference80_ball::division_free_tail_calls-before;if(out)fclose(out);
        }
        bool all=true;uint64_t missing_points=0,fallback_total=0;FILE*log=fopen((dest+"/missing_source_panels.csv").c_str(),"w");need(log,"actual missing panel outcomes");fprintf(log,"source,panel,points,complete,division_free_tail_calls,reason\n");
        for(size_t j=0;j<panels.size();j++){bool ok=failure[j].empty()&&local[j].points==256;all=all&&ok;missing_points+=local[j].points;fallback_total+=fallback[j];fprintf(log,"301,%zu,%llu,%d,%llu,%s\n",j,(unsigned long long)local[j].points,ok,(unsigned long long)fallback[j],csv_string(ok?"actual_panel_all12_complete":failure[j]).c_str());}fclose(log);
        value[301].points=missing_points;
        bool reduction_verified=false;
        if(all){
            target::Reference80Accumulator result,serial;FILE*combined=fopen((dest+"/NEW_missing_source301_queries.bin").c_str(),"wb");need(combined,"ordered missing query packet");
            for(size_t j=0;j<panels.size();j++){
                FILE*in=fopen((dest+"/source301_panel"+std::to_string(j)+".bin").c_str(),"r+b");need(in,"complete panel packet");target::Reference80Positive last;
                need(fseek(in,255*sizeof(last),SEEK_SET)==0&&fread(&last,sizeof(last),1,in)==1,"actual final panel point");
                merge_panel(result,local[j],last);
                need(fseek(in,255*sizeof(last),SEEK_SET)==0&&fwrite(&last,sizeof(last),1,in)==1&&fflush(in)==0&&fseek(in,0,SEEK_SET)==0,"cross-panel carry rounding recorded on original last point");
                for(int point=0;point<256;point++){target::Reference80Positive row;need(fread(&row,sizeof(row),1,in)==1&&row.source==301,"complete ordered native missing-source record");replay_moments_only(serial,row);need(fwrite(&row,sizeof(row),1,combined)==1,"ordered actual point packet");}
                need(fgetc(in)==EOF,"exact panel query count");fclose(in);
            }
            fclose(combined);target::finish_reference80(result);target::finish_reference80(serial);
            need(result.points==2816&&serial.points==2816,"exact missing-only query count");
            for(int k=0;k<12;k++)need(memcmp(&result.native_sum[k].n,&serial.native_sum[k].n,10)==0&&memcmp(&result.native_sum[k].e,&serial.native_sum[k].e,10)==0,"parallel power-of-two reduction is payload-bitwise identical to serial saved-moment-only tree");
            auto&r=value[301];r.complete=true;r.points=result.points;for(int k=0;k<12;k++)r.moment[k]=result.moment[k];reduction_verified=true;
        }
        FILE*integrals=fopen((dest+"/actual_clipped_source_integrals.csv").c_str(),"w");need(integrals,"reused25 plus actually completed missing source");fprintf(integrals,"source,moment,all12_complete,positive_points,value,native_geometry_pairwise_sum_error,analytic_remainder,arithmetic\n");
        for(int s=0;s<384;s++)if(value[s].seen)for(int k=0;k<12;k++){const auto&r=value[s];fprintf(integrals,"%d,%s,%d,%llu,%s,%s,%s,%s\n",s,target::moment_name[k],r.complete,(unsigned long long)r.points,r.complete?qs(r.moment[k].n).c_str():"unproved",r.complete?qs(r.moment[k].e).c_str():"unproved",qs(r.remainder[k]).c_str(),s==301?"conditional_division_free_binary80":"retained_binary80_values_unchanged");}fclose(integrals);
        S sum[12];Q remainder[12];for(int k=0;k<12;k++){sum[k]=group.integral.moment[k];remainder[k]=group.integral.remainder[k];}
        for(int s=0;s<384;s++)if(value[s].seen){all=all&&value[s].complete;for(int k=0;k<12;k++){if(value[s].complete)sum[k]=sum[k]+value[s].moment[k];remainder[k]=sf::add(remainder[k],value[s].remainder[k]);}}
        FILE*completion=fopen((dest+"/all384_source_completion.csv").c_str(),"w");need(completion,"explicit all384 completion");fprintf(completion,"%s\n",source_header.c_str());
        for(int s=0;s<384;s++){auto c=coverage[s];if(s==301){c[6]=all?"1":"0";c[10]=std::to_string(missing_points);c[12]=all?"actual_missing_source301_all12_complete":"actual_missing_source301_still_unproved";}for(size_t j=0;j<c.size();j++)fprintf(completion,"%s%s",j?",":"",c[j].c_str());fputc('\n',completion);}fclose(completion);
        target::SerializedPacket caller;caller.epoch=9;caller.target=1234;const auto original=caller;
        auto publication=target::publish(1,68526,sum,remainder,tau.data(),all,caller);if(!publication.published)need(!memcmp(&caller,&original,sizeof(caller)),"failed repaired target preserves whole caller packet");
        FILE*caller_file=fopen((dest+"/caller_target_packet.bin").c_str(),"wb");need(caller_file&&fwrite(&caller,sizeof(caller),1,caller_file)==1,"transactional caller packet");fclose(caller_file);
        FILE*budget=fopen((dest+"/all12_ORIGINAL_whole_target_budget.csv").c_str(),"w");need(budget,"original complete post-serialization physical budget");fprintf(budget,"moment,all384_complete,candidate,completed_native_geometry_sum_error,completed_remainder,final_R80_conversion_error,total_reported_error,ORIGINAL_tau,published\n");
        for(int k=0;k<12;k++)fprintf(budget,"%s,%d,%s,%s,%s,%s,%s,%s,%d\n",target::moment_name[k],all,qs(sum[k].n).c_str(),qs(sum[k].e).c_str(),qs(remainder[k]).c_str(),qs(publication.conversion[k]).c_str(),all?qs(publication.reported_error[k]).c_str():"unproved",qs(tau[k]).c_str(),publication.published);fclose(budget);
        checkpoint(dest+"/caller_unchanged_state.bin",f,1);checkpoint(dest+"/immutable_heat_generating_state.bin",heat.generating,1);
        printf("{\"scope\":\"ONLY missing source301 actual2816 queries, reuse old25 native values and retained4 whole group once, saved871 P2 potential proof patch, final original12tau transactional publication, no evolution\",\"MPI_ranks\":1,\"OMP_threads\":%d,\"new_query_source\":301,\"new_actual_queries\":%llu,\"retained_clipped_queries\":220160,\"retained_whole_queries\":45056,\"total_actual_queries\":%llu,\"division_free_tail_calls\":%llu,\"parallel_panel_reduction_matches_serial_saved_moment_tree\":%s,\"all384_complete\":%s,\"all_original12_tau_include_final_R80_conversion\":true,\"target_published\":%s,\"caller_preserved_on_rejection\":%s,\"accepted_RK4_evolution\":false,\"production_CUDA_transfer_certified\":false,\"long_Kepler_validated\":false,\"seconds\":%.17g}\n",omp_get_max_threads(),(unsigned long long)missing_points,(unsigned long long)(old_points+group.integral.points+missing_points),(unsigned long long)fallback_total,reduction_verified?"true":"false",all?"true":"false",publication.published?"true":"false",publication.published?"false":"true",std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count());
        code=publication.published?0:2;
    }catch(const std::exception&e){fprintf(stderr,"rank%d FAIL %s\n",rank,e.what());MPI_Abort(MPI_COMM_WORLD,2);}
    MPI_Finalize();return code;
}
