#define GFS_SHARED_CLIPPED_LIBRARY
#include "test_shared_clipped_target_v1.cpp"
#include "material_shared_potential_p2.h"

namespace p2=target::potential_p2;
namespace {
struct FrozenPanel {uint64_t source,index;int32_t policy,dependent,lower,upper;target::ClipPanel geometry;};
struct PanelProof {FrozenPanel panel;uint64_t record=0;Q old_R=0,chosen_R=0;p2::Certificate proof;};
struct SourceProof {uint64_t panels=0,proved=0,improved=0,evaluations=0;Q chosen_R=0;int complete=0;};
struct NativeSource {bool seen=false,complete=false;std::array<bool,12>moment_seen={};S moment[12];Q remainder[12]={};};
}

int main(int argc,char**argv){
    int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
    int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);int code=2;
    try{
        need(argc==11,"packet saved sites weights originalcall realbounds originating whole4 ORIGINAL_tau NEW_output");
        auto begin=std::chrono::steady_clock::now();auto p=curved::read(argv[1]);
        const auto f=load_saved(std::string(argv[2])+"/shared_state.bin",p);
        const sf::FrozenThermal heat(f,p.parameters[2],p.parameters[1],1);
        auto context=target::retained_context(argv[5],f,heat);
        auto sources=target::retained_sources(argv[5],argv[6],f);
        auto group=target::retained_whole_group(argv[8],sources);
        auto tau=target::original_tolerances(argv[9]);
        const std::string originating=argv[7],dest=argv[10];
        std::ifstream terminal(originating+"/native_exit_code.txt");int prior_exit=-1;
        need(bool(terminal)&&bool(terminal>>prior_exit)&&prior_exit==2,"one terminal rejected originating call");
        need(p2::unit_contract(),"normalized P2 manufactured unit contract");
        auto planes=target::planes(sites(argv[3],argv[4]),68526);
        std::vector<std::vector<PanelProof>>proofs(384);
        uint64_t loaded=0;
        for(int owner=0;owner<4;owner++){
            std::map<std::pair<int,int>,Q>old;
            std::ifstream certificate(originating+"/actual_clipped_panel_certificates_rank"+std::to_string(owner)+".csv");
            need(bool(certificate),"saved original valid panel certificates");
            std::string line;std::getline(certificate,line);auto names=target::columns(line);
            auto field=std::find(names.begin(),names.end(),"potential_remainder");
            need(field!=names.end(),"saved potential remainder column");size_t column=field-names.begin();
            while(std::getline(certificate,line)){
                auto c=target::columns(line);need(c[8]=="1","every old analytic certificate valid");
                need(old.emplace(std::make_pair(std::stoi(c[0]),std::stoi(c[1])),target::number(c[column])).second,"unique saved source/panel bound");
            }
            FILE*input=fopen((originating+"/actual_reference_panel_packets_rank"+std::to_string(owner)+".bin").c_str(),"rb");
            need(input,"immutable typed panel packet");uint64_t h[8];
            need(fread(h,sizeof(h),1,input)==1,"typed panel header");
            need(h[0]==0x47514650414e4c31ULL&&h[1]==1&&h[2]==1&&h[3]==sizeof(FrozenPanel)&&h[5]==uint64_t(owner)&&h[6]==0&&h[7]==0,"actual original epoch/ABI/owner");
            need(h[4]==old.size(),"certificate and actual packet panel count match");
            for(uint64_t i=0;i<h[4];i++){
                PanelProof row;need(fread(&row.panel,sizeof(row.panel),1,input)==1,"complete original typed panel");
                int s=row.panel.source,j=row.panel.index;need(s>=0&&s<384&&s%4==owner&&j==int(proofs[s].size()),"canonical original source/panel ordering");
                row.old_R=old.at({s,j});row.chosen_R=row.old_R;row.record=i;proofs[s].push_back(row);loaded++;
            }
            need(fgetc(input)==EOF,"typed packet length exact");fclose(input);
        }
        need(loaded==871,"exactly original871 potential certificates, no new panel topology");
        std::vector<NativeSource>native(384);std::ifstream values(originating+"/actual_clipped_source_integrals.csv");
        need(bool(values),"retained native source values/errors");std::string line;std::getline(values,line);size_t rows=0;
        while(std::getline(values,line)){
            auto c=target::columns(line);int s=std::stoi(c[0]),k=-1;
            for(int j=0;j<12;j++)if(c[1]==target::moment_name[j])k=j;
            need(s>=0&&s<384&&!proofs[s].empty()&&k>=0,"retained exact source/moment identity");
            auto&r=native[s];bool complete=c[2]=="1";
            need(!r.moment_seen[k]&&(!r.seen||r.complete==complete),"unique consistently complete saved source moments");
            r.seen=true;r.complete=complete;r.moment_seen[k]=true;
            if(complete)r.moment[k]=S(target::number(c[4]),target::number(c[5]));
            r.remainder[k]=target::number(c[6]);rows++;
        }
        need(rows==312,"all26x12 original source contributions retained");
        uint64_t source_count=0,missing=0;
        for(int s=0;s<384;s++)if(!proofs[s].empty()){
            need(native[s].seen&&std::all_of(native[s].moment_seen.begin(),native[s].moment_seen.end(),[](bool b){return b;}),"saved complete all12 schema");
            source_count++;missing+=!native[s].complete;
        }
        need(source_count==26&&missing==1&&!native[301].complete,"same actual missing source301, not a new test input");
        std::vector<SourceProof>local(384);
        #pragma omp parallel for schedule(dynamic,1)
        for(int s=0;s<384;s++)if(s%ranks==rank&&!proofs[s].empty()){
            auto m=sd::compile(f,s);std::vector<target::Face>faces;for(const auto&plane:planes)faces.push_back(target::pullback(m,plane));
            std::array<target::Axis,6>axis;std::array<bool,6>seen={};
            for(auto&row:proofs[s]){
                auto&panel=row.panel.geometry;int ci=panel.strip.chart;
                need(ci>=0&&ci<6,"original chart identity");
                if(!panel.whole&&!seen[ci]){axis[ci]=target::certified_axis(m,panel.reference,faces,planes,sources[s].real);seen[ci]=true;}
                if(!panel.whole)need(!axis[ci].empty&&!axis[ci].whole&&axis[ci].dependent==row.panel.dependent,"same original actual chart axis and branch");
                row.proof=p2::certify(p,m,faces,axis[ci],panel,row.panel.policy==1,context[f.triangle[s].cell]);
                if(row.proof.pass)row.chosen_R=std::min(row.old_R,row.proof.remainder);
                auto&r=local[s];r.panels++;r.proved+=row.proof.pass;r.improved+=row.chosen_R<row.old_R;r.evaluations+=row.proof.evaluations;
                r.chosen_R=sf::add(r.chosen_R,row.chosen_R);
            }
            need(local[s].chosen_R<=native[s].remainder[target::POTENTIAL],"valid-min chosen source remainder cannot exceed unchanged original bound");
            local[s].complete=1;
        }
        FILE*csv=fopen((dest+"/potential_P2_certificates_rank"+std::to_string(rank)+".csv").c_str(),"w");need(csv,"per-panel fixed P2 proof ledger");
        fprintf(csv,"source,panel,origin_rank,origin_record,coordinate_policy,dependent,lower,upper,new_proof,center_complete,axis0_complete,axis1_complete,old_R,new_R,chosen_R,a0,a1,a2,A00,A01,A11,center_value_error,center_grad0_error,center_grad1_error,axis0_majorant,axis1_majorant,Hdiff0_00,Hdiff0_01,Hdiff0_10,Hdiff0_11,Hdiff1_00,Hdiff1_01,Hdiff1_10,Hdiff1_11,AD2_evaluations,reason\n");
        for(int s=rank;s<384;s+=ranks)for(const auto&row:proofs[s]){
            const auto&r=row.proof;
            fprintf(csv,"%d,%llu,%d,%llu,%d,%d,%d,%d,%d,%d,%d,%d,%s,%s,%s",s,(unsigned long long)row.panel.index,s%4,(unsigned long long)row.record,row.panel.policy,row.panel.dependent,row.panel.lower,row.panel.upper,r.pass,r.center_complete,r.axis_complete[0],r.axis_complete[1],qs(row.old_R).c_str(),r.pass?qs(r.remainder).c_str():"unproved",qs(row.chosen_R).c_str());
            for(Q q:r.coefficient)fprintf(csv,",%s",qs(q).c_str());
            for(Q q:r.center_error)fprintf(csv,",%s",qs(q).c_str());
            for(Q q:r.axis_majorant)fprintf(csv,",%s",qs(q).c_str());
            for(const auto&axis:r.hessian_difference)for(Q q:axis)fprintf(csv,",%s",qs(q).c_str());
            fprintf(csv,",%llu,%s\n",(unsigned long long)r.evaluations,r.reason);
        }
        fclose(csv);std::vector<SourceProof>gathered(rank==0?384*ranks:0);
        MPI_Gather(local.data(),384*sizeof(SourceProof),MPI_BYTE,gathered.data(),384*sizeof(SourceProof),MPI_BYTE,0,MPI_COMM_WORLD);
        if(rank==0){
            S partial[12];Q remainder[12];for(int k=0;k<12;k++){partial[k]=group.integral.moment[k];remainder[k]=group.integral.remainder[k];}
            uint64_t proven=0,improved=0,evaluations=0,panels=0;
            for(int s=0;s<384;s++)if(!proofs[s].empty()){
                const auto&r=gathered[(s%ranks)*384+s];need(r.complete,"every source fixed proof pass attempted exactly once");
                proven+=r.proved;improved+=r.improved;evaluations+=r.evaluations;panels+=r.panels;
                for(int k=0;k<12;k++){
                    if(native[s].complete)partial[k]=partial[k]+native[s].moment[k];
                    remainder[k]=sf::add(remainder[k],k==target::POTENTIAL?r.chosen_R:native[s].remainder[k]);
                }
            }
            need(panels==871,"all original potential proof attempts complete");
            bool possible=true;FILE*budget=fopen((dest+"/known_all12_budget.csv").c_str(),"w");need(budget,"known subbudget, never a complete target verdict");
            fprintf(budget,"moment,completed_source_forward_and_Qsum_error,all_source_analytic_remainder,known_total_without_missing_source_or_final_conversion,ORIGINAL_tau,known_terms_strictly_below_tau\n");
            for(int k=0;k<12;k++){
                Q known=sf::add(partial[k].e,remainder[k]);bool ok=known<tau[k];possible=possible&&ok;
                fprintf(budget,"%s,%s,%s,%s,%s,%d\n",target::moment_name[k],qs(partial[k].e).c_str(),qs(remainder[k]).c_str(),qs(known).c_str(),qs(tau[k]).c_str(),ok);
            }
            fclose(budget);checkpoint(dest+"/caller_unchanged_state.bin",f,1);checkpoint(dest+"/immutable_heat_generating_state.bin",heat.generating,1);
            printf("{\"scope\":\"ONLY fixed P2 proof of normalized potential integrand on original871 stored panels; no quadrature or field-profile replay of original220160 points; no target publication or evolution\",\"MPI_ranks\":%d,\"OMP_threads\":%d,\"original_panels\":%llu,\"new_valid_P2_proofs\":%llu,\"strictly_improved_potential_remainders\":%llu,\"full_potential_AD2_evaluations_attempted\":%llu,\"normalized_polynomial_unit_pass\":true,\"all_stored_native_values_and_other11_remainders_unchanged\":true,\"original_missing_source\":301,\"known_budget_compatible\":%s,\"missing_source_forward_and_final_serialization_still_unproved\":true,\"target_published\":false,\"seconds\":%.17g}\n",ranks,omp_get_max_threads(),(unsigned long long)panels,(unsigned long long)proven,(unsigned long long)improved,(unsigned long long)evaluations,possible?"true":"false",std::chrono::duration<double>(std::chrono::steady_clock::now()-begin).count());
            code=possible?0:2;
        }
        MPI_Bcast(&code,1,MPI_INT,0,MPI_COMM_WORLD);
    }catch(const std::exception&e){fprintf(stderr,"rank%d FAIL %s\n",rank,e.what());MPI_Abort(MPI_COMM_WORLD,2);}
    MPI_Finalize();return code;
}
