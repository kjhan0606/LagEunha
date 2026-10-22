// Diagnostic only: one first query of an actually incomplete saved source.
// Build against the frozen originating run, never reconstruct its atlas.
#include <cstdio>
#include <quadmath.h>
#include <execinfo.h>
static void gfs_forward_quotient_failure(__float128 numerator,__float128 numerator_radius,__float128 denominator,__float128 denominator_radius){
    char n[128],ne[128],d[128],de[128];
    quadmath_snprintf(n,sizeof(n),"%.36Qg",numerator);
    quadmath_snprintf(ne,sizeof(ne),"%.36Qg",numerator_radius);
    quadmath_snprintf(d,sizeof(d),"%.36Qg",denominator);
    quadmath_snprintf(de,sizeof(de),"%.36Qg",denominator_radius);
    fprintf(stderr,"ACTUAL_QDIV numerator=%s numerator_radius=%s denominator=%s denominator_radius=%s\n",n,ne,d,de);
    void*frames[32];int count=backtrace(frames,32);backtrace_symbols_fd(frames,count,2);
}
#define GFS_FORWARD_QUOTIENT_FAILURE(a,ae,b,be) gfs_forward_quotient_failure(__float128(a),__float128(ae),__float128(b),__float128(be))
#define GFS_SHARED_CLIPPED_LIBRARY
#include "test_shared_clipped_target_v1.cpp"

namespace {
struct FrozenPanel {
    uint64_t source,index;
    int32_t policy,dependent,lower,upper;
    target::ClipPanel geometry;
};
std::string json_string(const std::string&input){
    const char*hex="0123456789abcdef";std::string out="\"";
    for(unsigned char c:input){
        if(c=='\"'||c=='\\'){out+='\\';out+=char(c);}
        else if(c<32){out+="\\u00";out+=hex[c>>4];out+=hex[c&15];}
        else out+=char(c);
    }
    return out+'\"';
}
}

int main(int argc,char**argv){
    int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);
    int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    int code=2;
    try{
        need(argc==10&&ranks==1,"packet saved sites weights originalcall realbounds TERMINAL originating call source NEW diagnostic output; MPI1 only");
        auto p=curved::read(argv[1]);
        const auto f=load_saved(std::string(argv[2])+"/shared_state.bin",p);
        const sf::FrozenThermal heat(f,p.parameters[2],p.parameters[1],1);
        auto context=target::retained_context(argv[5],f,heat);
        auto source=target::retained_sources(argv[5],argv[6],f);
        const std::string originating=argv[7],dest=argv[9];const int selected=std::stoi(argv[8]);
        need(selected>=0&&selected<384,"actual selected source identity");
        std::ifstream terminal(originating+"/native_exit_code.txt");int exit_code=-1;
        need(bool(terminal)&&bool(terminal>>exit_code)&&exit_code==2,"originating call must already be terminal and rejected before one-query diagnosis");
        std::ifstream ledger(originating+"/all384_source_completion.csv");
        need(bool(ledger),"actual originating source completion ledger");
        std::string line;std::getline(ledger,line);bool found=false;std::string actual_reason;
        while(std::getline(ledger,line)){
            auto c=target::columns(line);if(std::stoi(c[0])!=selected)continue;
            need(!found&&c[6]=="0"&&c[10]=="0"&&c[5]=="1","only an actually incomplete zero-record source after analytic success; first query identity is then exact");
            found=true;actual_reason=c[12];
        }
        need(found,"selected actual source row exists");
        FILE*packet=fopen((originating+"/actual_reference_panel_packets_rank"+std::to_string(selected%4)+".bin").c_str(),"rb");
        need(packet,"frozen actual full panel geometry");
        uint64_t header[8];need(fread(header,sizeof(header),1,packet)==1,"typed geometry header");
        need(header[0]==0x47514650414e4c31ULL&&header[1]==1&&header[2]==1&&header[3]==sizeof(FrozenPanel)&&header[5]==uint64_t(selected%4),"originating panel ABI, epoch and owner");
        FrozenPanel first={};bool panel_found=false;
        for(uint64_t i=0;i<header[4];i++){
            FrozenPanel row;need(fread(&row,sizeof(row),1,packet)==1,"complete typed panel record");
            if(row.source==uint64_t(selected)&&row.index==0){need(!panel_found,"unique first panel");first=row;panel_found=true;}
        }
        need(fgetc(packet)==EOF,"typed panel packet has exact terminal length");fclose(packet);
        need(panel_found,"actual first panel retained");
        auto m=sd::compile(f,selected);auto sites0=sites(argv[3],argv[4]);
        auto planes=target::planes(sites0,68526);std::vector<target::Face>faces;
        for(const auto&plane:planes)faces.push_back(target::pullback(m,plane));
        const auto&panel=first.geometry;auto gauss=target::gauss16();
        FILE*trace=fopen((dest+"/first_actual_query.csv").c_str(),"w");
        FILE*positive=fopen((dest+"/diagnostic_positive_query.bin").c_str(),"wb");
        need(trace&&positive,"one actual query diagnostic outputs");
        fprintf(trace,"quantity,value,radius\n");
        auto record=[&](const std::string&name,S value){fprintf(trace,"%s,%s,%s\n",name.c_str(),qs(value.n).c_str(),qs(value.e).c_str());};
        record("source",S(selected));record("panel",S(0));record("iu",S(0));record("iv",S(0));
        record("coordinate_policy",S(first.policy));record("chart",S(panel.strip.chart));
        record("dependent_axis",S(first.dependent));record("chart_det",panel.reference.det);
        std::string phase="real_quadrature_map",failure;bool passed=false;
        sd::Point point;S weight;
        try{
            if(panel.whole){
                const auto&b=panel.ordinary;S ds=b.shi-b.slo,dt=b.thi-b.tlo;
                S s=b.slo+ds*gauss.x[0],t=b.tlo+dt*gauss.x[0];
                record("s",s);record("t",t);
                point=sd::evaluate(m,panel.reference,s,t);
                weight=point.reference_jacobian*ds*dt*gauss.w[0]*gauss.w[0];
            }else{
                auto a=target::certified_axis(m,panel.reference,faces,planes,source[selected].real);
                need(!a.empty&&!a.whole&&a.dependent==first.dependent,"same frozen actual chart axis; no new atlas");
                S free=panel.strip.lo+(panel.strip.hi-panel.strip.lo)*gauss.x[0];
                S eta=panel.eta_lo+(panel.eta_hi-panel.eta_lo)*gauss.x[0];
                record("free",free);record("eta",eta);
                auto q=target::strip_point(m,panel.reference,faces,a,panel.strip,free,eta,false);
                record("s",q.s);record("t",q.t);record("strip_jacobian",q.jacobian);
                point=q.value;weight=q.jacobian*(panel.eta_hi-panel.eta_lo)*gauss.w[0]*gauss.w[0];
            }
            record("X",point.X[0]);record("Y",point.X[1]);record("positive_weight",weight);
            for(int k=0;k<2;k++){record("position"+std::to_string(k),point.field.position[k]);record("velocity"+std::to_string(k),point.field.velocity[k]);}
            for(int k=0;k<4;k++){record("F"+std::to_string(k),point.field.F[k]);record("V"+std::to_string(k),point.field.V[k]);record("grad_v"+std::to_string(k),point.field.grad_v[k]);}
            record("J",point.field.J);fflush(trace);
            phase="original_target_plane_ownership";int failed_plane=-1;
            for(size_t i=0;i<planes.size();i++){
                S value=target::value(planes[i],point.field.position[0],point.field.position[1]);
                record("target_plane"+std::to_string(i),value);
                record("plane_competitor"+std::to_string(i),S(Q(planes[i].competitor)));
                if(failed_plane<0&&!(sf::upper(value)<0))failed_plane=int(i);
            }
            fflush(trace);need(failed_plane<0,"every new common quadrature point inside all actual target planes");
            phase="original_binary80_field_profile_all12";
            target::Reference80Accumulator result;
            target::point_reference80(p,m,point,weight,context[f.triangle[selected].cell],result,positive);
            passed=true;
        }catch(const std::exception&e){failure=e.what();}
        fclose(trace);fclose(positive);
        checkpoint(dest+"/caller_unchanged_state.bin",f,1);
        checkpoint(dest+"/immutable_heat_generating_state.bin",heat.generating,1);
        uint64_t fallback_calls=0;
#ifdef GFS_REFERENCE80_DIVISION_FREE_FALLBACK
        fallback_calls=shared_reference80_ball::division_free_tail_calls;
#endif
        printf("{\"scope\":\"ONLY first actual query source/panel0/iu0/iv0, originating frozen geometry, no atlas or other-source quadrature replay, no evolution\",\"source\":%d,\"originating_source_reason\":%s,\"query_passed\":%s,\"phase\":%s,\"failure\":%s,\"failure_message_reproduced\":%s,\"division_free_tail_calls\":%llu,\"target_published\":false}\n",selected,json_string(actual_reason).c_str(),passed?"true":"false",json_string(phase).c_str(),json_string(failure).c_str(),failure==actual_reason?"true":"false",(unsigned long long)fallback_calls);
        code=0; // A recorded diagnosis is not an accepted target query.
    }catch(const std::exception&e){fprintf(stderr,"rank%d FAIL %s\n",rank,e.what());MPI_Abort(MPI_COMM_WORLD,2);}
    MPI_Finalize();return code;
}
