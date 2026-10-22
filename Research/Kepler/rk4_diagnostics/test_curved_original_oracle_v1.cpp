#define GFS_CURVED_BRANCH_LIBRARY
#include "test_curved_branch_v1.cpp"
#include "material_original_oracle.h"
#include "material_native_forward.h"
#include <atomic>
namespace mf=material_forward;
using Q=mf::Q;
template<class T> static std::array<mf::B<T>,13> fields(const mf::Point<T>&p){return {p.position[0],p.position[1],p.velocity[0],p.velocity[1],p.F[0],p.F[1],p.F[2],p.F[3],p.J,p.grad_v[0],p.grad_v[1],p.grad_v[2],p.grad_v[3]};}
static const char *field_names[13]={"x","y","vx","vy","F00","F01","F10","F11","J","grad_v00","grad_v01","grad_v10","grad_v11"};
static mf::B<Q> original_density(const curved::Packet&p,const QuadPoint&q){
    using S=mf::B<Q>;S X(Q(q.X),Q(mf::cast_radius(q.X))),Y(Q(q.Y),Q(mf::cast_radius(q.Y))),r=mf::root(X*X+Y*Y);
    return S(Q(p.parameters[6]))+S(Q(p.parameters[9]))/(r+S(Q(p.parameters[7])))*S(.5L)*(S(1)+mf::tanh_((r-S(Q(p.parameters[4])))/S(Q(p.parameters[8]))))*S(.5L)*(S(1)+mf::tanh_((S(Q(p.parameters[5]))-r)/S(Q(p.parameters[8]))));
}
template<class T> static std::array<mf::B<T>,NMOMENT> moment(const curved::Packet&p,const mf::Point<T>&v,mf::B<T>rho,mf::B<T>w,mf::B<T>div,mf::B<T>eps2){
    using S=mf::B<T>;std::array<S,NMOMENT> z;S dm=rho*w,x=v.position[0],y=v.position[1],vx=v.velocity[0],vy=v.velocity[1],P=S(T(p.parameters[2]))*mf::power(v.J,T(-p.parameters[1]));
    z[MASS]=dm;z[PX]=dm*vx;z[PY]=dm*vy;z[ANGULAR]=dm*(x*vy-y*vx);z[KINETIC]=S(.5L)*dm*(vx*vx+vy*vy);z[THERMAL]=P*v.J*w/(S(T(p.parameters[1]))-S(1));z[VOLUME]=v.J*w;z[QX]=dm*x;z[QY]=dm*y;z[POLAR]=dm*(x*x+y*y);z[THERMAL_RATE]=-P*v.J*w*div;z[POTENTIAL]=-dm/mf::root(x*x+y*y+eps2);return z;
}
static mf::Point<R> lift(const mf::Point<double>&p){mf::Point<R> out;auto c=[](mf::B<double>a){return mf::B<R>(a.n,(R)a.e);};for(int k=0;k<2;k++){out.position[k]=c(p.position[k]);out.velocity[k]=c(p.velocity[k]);}for(int k=0;k<4;k++){out.F[k]=c(p.F[k]);out.V[k]=c(p.V[k]);out.grad_v[k]=c(p.grad_v[k]);}out.J=c(p.J);return out;}
struct Stats {R error=0,radius=0,ratio=0,oracle_radius=0;uint64_t max_index=0;};
static void observe(Stats&s,uint64_t index,Q error,R radius,Q oracle){s.error=std::max(s.error,(R)error);s.radius=std::max(s.radius,radius);s.oracle_radius=std::max(s.oracle_radius,(R)oracle);R ratio=radius>0?(R)std::max(Q(0),error-oracle)/radius:(error>oracle?INFINITY:0);if(ratio>s.ratio){s.ratio=ratio;s.max_index=index;}}
#ifndef GFS_CURVED_ORIGINAL_ORACLE_LIBRARY
int main(int argc,char **argv){
    int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);active_mpi=true;need(argc==7,"original packet sites guard fixed-positive-packet new-output case");
    auto p=curved::read(argv[1]);auto sites=curved_partition::read_sites(argv[2]);auto weights=read_weights(argv[3],sites);std::string name=argv[6],dest=argv[5];bool lens=name=="lens";size_t it=lens?0:atoi(argv[6]);need(lens||it==225||it==137,"same three accepted controls only");if(lens)p=lens_setup(sites,weights);const auto&r=p.record[it];
    FILE*f=fopen(argv[4],"rb");uint64_t h[6];need(f&&fread(h,sizeof(h),1,f)==1&&h[0]==0x4746503346435131ULL&&h[1]==1&&h[2]==it&&h[4]==56,"unchanged v5 positive packet");std::vector<QuadPoint> points(h[3]);need(fread(points.data(),sizeof(QuadPoint),points.size(),f)==points.size()&&fgetc(f)==EOF,"all positive points preserved");fclose(f);
    auto map=curved_compiled::compile(r,p.parameters[3]);Stats native_stats[13],representative_stats[13],conversion_stats[13],moment_stats[NMOMENT];uint64_t used=0,native_graph_mismatch=0,moment_graph_mismatch=0,model_excess=0;std::atomic<bool> stop(false);uint64_t witness_index=0;std::string witness_reason;
    mf::B<Q> reference[24][NMOMENT],canonical[24][NMOMENT];R forward[24][NMOMENT]={},absolute[24][NMOMENT]={};uint64_t counts[24]={};auto started=std::chrono::steady_clock::now();
    #pragma omp parallel
    {
        Stats ns[13],rs[13],cs[13],ms[NMOMENT];uint64_t local_used=0,ng=0,mg=0,ex=0;
        mf::B<Q> ref[24][NMOMENT],nom[24][NMOMENT];R bound[24][NMOMENT]={},ab[24][NMOMENT]={};uint64_t count[24]={};
        #pragma omp for schedule(static)
        for(size_t ip=0;ip<points.size();ip++){
            if(stop.load(std::memory_order_relaxed))continue;const auto&q=points[ip];need(q.triangle==it&&q.target<sites.size()&&q.reference_weight>0&&q.rho0>0,"saved positive canonical point");
            auto cpu=curved::native_point(r,q.X,q.Y,p.parameters[3]);auto native=mf::native(r,q.X,q.Y,p.parameters[3],true);auto oracle=mf::original(r,q.X,q.Y,p.parameters[3]),rep=mf::representative(map,q.X,q.Y,p.parameters[3]);auto stable=mf::stable(map,q.X,q.Y,p.parameters[3],false),converted=mf::stable(map,q.X,q.Y,p.parameters[3],true);
            auto n=fields(native);auto o=fields(oracle),rr=fields(rep);auto ss=fields(stable),cc=fields(converted);double actual[13]={cpu.position[0],cpu.position[1],cpu.velocity[0],cpu.velocity[1],cpu.F[0],cpu.F[1],cpu.F[2],cpu.F[3],cpu.J,cpu.grad_v[0],cpu.grad_v[1],cpu.grad_v[2],cpu.grad_v[3]};bool bad=false;
            for(int k=0;k<13;k++){
                if(n[k].n!=actual[k]){ng++;bad=true;}
                Q en=fabsq(Q(n[k].n)-o[k].n),er=fabsq(Q(ss[k].n)-rr[k].n),ec=fabsq(Q(cc[k].n)-o[k].n);observe(ns[k],ip,en,n[k].e,o[k].e);observe(rs[k],ip,er,ss[k].e,rr[k].e);observe(cs[k],ip,ec,cc[k].e,o[k].e);
                if(en>mf::FP<Q>::up(Q(n[k].e)+o[k].e)||er>mf::FP<Q>::up(Q(ss[k].e)+rr[k].e)||ec>mf::FP<Q>::up(Q(cc[k].e)+o[k].e)){ex++;bad=true;}
            }
            auto rho=original_density(p,q);R density_error=mf::FP<R>::up((R)fabsq(Q(q.rho0)-rho.n)+(R)rho.e);mf::B<double> eps2native=mf::B<double>(p.parameters[3])*mf::B<double>(p.parameters[3]),divnative=native.grad_v[0]+native.grad_v[3];
            auto nn=moment(p,lift(native),mf::B<R>(q.rho0,density_error),mf::B<R>(q.reference_weight),mf::B<R>(divnative.n,divnative.e),mf::B<R>(eps2native.n,eps2native.e));
            auto oo=moment(p,oracle,rho,mf::B<Q>(Q(q.reference_weight)),oracle.grad_v[0]+oracle.grad_v[3],mf::B<Q>(Q(p.parameters[3]))*mf::B<Q>(Q(p.parameters[3])));auto ordinary=integrand(p,cpu,q.rho0,q.reference_weight);
            for(int k=0;k<NMOMENT;k++){
                if(nn[k].n!=ordinary.a[k]){mg++;bad=true;}Q e=fabsq(Q(nn[k].n)-oo[k].n);observe(ms[k],ip,e,nn[k].e,oo[k].e);if(e>mf::FP<Q>::up(Q(nn[k].e)+oo[k].e)){ex++;bad=true;}
                ref[q.target][k]=ref[q.target][k]+oo[k];nom[q.target][k]=nom[q.target][k]+mf::B<Q>(Q(nn[k].n));bound[q.target][k]=mf::plus_up(bound[q.target][k],nn[k].e);ab[q.target][k]=mf::plus_up(ab[q.target][k],fabsl(nn[k].n));
            }count[q.target]++;local_used++;
            if(bad){bool expected=false;if(stop.compare_exchange_strong(expected,true)){
                #pragma omp critical(oracle_witness)
                {witness_index=ip;witness_reason=ng?"native scalar graph mismatch":(mg?"moment scalar graph mismatch":"independent oracle outside proposed forward model");}
            }}
        }
        #pragma omp critical(oracle_merge)
        {
            auto merge=[](Stats&a,const Stats&b){a.error=std::max(a.error,b.error);a.radius=std::max(a.radius,b.radius);a.oracle_radius=std::max(a.oracle_radius,b.oracle_radius);if(b.ratio>a.ratio){a.ratio=b.ratio;a.max_index=b.max_index;}};
            for(int k=0;k<13;k++){merge(native_stats[k],ns[k]);merge(representative_stats[k],rs[k]);merge(conversion_stats[k],cs[k]);}for(int k=0;k<NMOMENT;k++)merge(moment_stats[k],ms[k]);used+=local_used;native_graph_mismatch+=ng;moment_graph_mismatch+=mg;model_excess+=ex;
            for(size_t j=0;j<sites.size();j++){counts[j]+=count[j];for(int k=0;k<NMOMENT;k++){reference[j][k]=reference[j][k]+ref[j][k];canonical[j][k]=canonical[j][k]+nom[j][k];forward[j][k]=mf::plus_up(forward[j][k],bound[j][k]);absolute[j][k]=mf::plus_up(absolute[j][k],ab[j][k]);}}
        }
    }
    FILE*out=fopen((dest+"/field_forward.csv").c_str(),"w");need(out,"forward field error ledger");fprintf(out,"path,field,max_observed_error,max_forward_radius,max_oracle_radius,max_error_over_radius,witness_point\n");for(int mode=0;mode<3;mode++)for(int k=0;k<13;k++){auto s=(mode==0?native_stats:(mode==1?representative_stats:conversion_stats))[k];fprintf(out,"%s,%s,%.21Lg,%.21Lg,%.21Lg,%.21Lg,%llu\n",mode==0?"actual_native":"stable" ,field_names[k],s.error,s.radius,s.oracle_radius,s.ratio,(unsigned long long)s.max_index);}for(int k=0;k<NMOMENT;k++){auto s=moment_stats[k];fprintf(out,"moment,%s,%.21Lg,%.21Lg,%.21Lg,%.21Lg,%llu\n",moment_names[k],s.error,s.radius,s.oracle_radius,s.ratio,(unsigned long long)s.max_index);}fclose(out);
    out=fopen((dest+"/canonical_target_fields.csv").c_str(),"w");need(out,"per-target original reference field ledger");fprintf(out,"source,target,moment,points,canonical_native,original_reference,reference_radius,observed_signed_difference,field_forward_bound,absolute_term_sum\n");for(size_t j=0;j<sites.size();j++)if(counts[j])for(int k=0;k<NMOMENT;k++)fprintf(out,"%llu,%llu,%s,%llu,%.21Lg,%.21Lg,%.21Lg,%.21Lg,%.21Lg,%.21Lg\n",(unsigned long long)r.cell,(unsigned long long)sites[j].id,moment_names[k],(unsigned long long)counts[j],(R)canonical[j][k].n,(R)reference[j][k].n,(R)reference[j][k].e,(R)(canonical[j][k].n-reference[j][k].n),forward[j][k],absolute[j][k]);fclose(out);
    bool pass=!stop&&used==points.size();printf("{\"scope\":\"same saved positive points; original stored-coefficient mathematical basis; actual native scalar graph and stable FMA Horner; no state or nominal packet update\",\"case\":\"%s\",\"points_available\":%zu,\"points_checked\":%llu,\"native_graph_mismatches\":%llu,\"moment_graph_mismatches\":%llu,\"oracle_model_excess\":%llu,\"witness_point\":%llu,\"reason\":\"%s\",\"fixed_binary128_oracle\":true,\"libm_model\":\"normal IEEE nearest arithmetic; elementary functions <=8ulp engineering assumption; CPU only\",\"point_forward_model_pass\":%s,\"source_quadrature_geometry_gluing_serialization_sum_and_target_budget_complete\":false,\"full_cold_transfer_proved\":false,\"seconds\":%.17g}\n",name.c_str(),points.size(),(unsigned long long)used,(unsigned long long)native_graph_mismatch,(unsigned long long)moment_graph_mismatch,(unsigned long long)model_excess,(unsigned long long)witness_index,witness_reason.c_str(),pass?"true":"false",std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count());active_mpi=false;MPI_Finalize();return pass?0:2;
}
#endif
