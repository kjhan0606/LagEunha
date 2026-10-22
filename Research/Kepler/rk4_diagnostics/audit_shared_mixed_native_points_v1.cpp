// Integer audit of stored point packets and replay of the summation tree.
// No field evaluator, geometry, physical query, CUDA or MPI is linked here.
#include <algorithm>
#include <array>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include <omp.h>
static void need(bool ok,const char*message){if(!ok)throw std::runtime_error(message);}
#include "material_forward_ball.h"
using Q=__float128;using S=material_forward::B<Q>;using B=material_forward::B<long double>;
using Wide=__int128;
struct SourcePacket {uint64_t slot,source,target,planned_points,points,complete,oracles;S moment[12];Q remainder[12];char reason[192];};
struct Row {uint64_t source;Q X[2],weight,input_radius[3],cache_radius_max[3],field_radius[17],moment[12],radius[12],sum_local_rounding[12];};
struct Commit {uint64_t slot,panel,next_query,points,occupied,record_hash,successful,version;S level[64][12];};
static_assert(sizeof(SourcePacket)==832&&sizeof(Row)==1008&&sizeof(Commit)==24640,"stored packet ABI");
// Exact dyadic superaccumulator. At most13074432 64-bit significands
// contribute per exponent, so signed128 bins have over35 bits of headroom.
// Sign uses mathematical floor carries: every discarded remainder is
// nonnegative, and all lower remainders together are smaller than one unit
// of the next bin. No floating addition or external bignum package is used.
struct Bank {std::map<int,Wide>values;void add(Q x,bool nonnegative=false,int sign=1){uint64_t bits[2];memcpy(bits,&x,16);int exponent=(bits[1]>>48)&32767;need(exponent<32767&&(!nonnegative||!(bits[1]>>63)),"finite native value/nonnegative radius");unsigned __int128 mantissa=(static_cast<unsigned __int128>(bits[1]&((1ULL<<48)-1))<<64)|bits[0];if(!exponent){need(!mantissa,"normal or zero native point value");return;}mantissa|=static_cast<unsigned __int128>(1)<<112;need(!(mantissa&((static_cast<unsigned __int128>(1)<<49)-1)),"exact binary80 point serialization");Wide n=uint64_t(mantissa>>49);if(bits[1]>>63)n=-n;values[exponent-16446]+=sign*n;}
 void merge(const Bank&b,int sign){for(const auto&v:b.values)values[v.first]+=sign*v.second;}
 int sign()const{if(values.empty())return 0;auto it=values.begin();Wide n=it->second;int exponent=it->first;bool lower_nonzero=false;while(++it!=values.end()){int gap=it->first-exponent;need(gap>0,"strict superaccumulator exponent order");Wide carry;if(gap>=127){carry=n<0?-1:0;lower_nonzero=lower_nonzero||n!=0;}else{Wide base=Wide(1)<<gap;carry=n/base;Wide remainder=n%base;if(remainder<0){carry--;remainder+=base;}lower_nonzero=lower_nonzero||remainder!=0;}n=it->second+carry;exponent=it->first;}return n>0?1:n<0?-1:lower_nonzero?1:0;}};
template<class T>static std::vector<T>read_packet(const std::string&path){FILE*f=fopen(path.c_str(),"rb");need(f,"existing source packet");uint64_t h[8];need(fread(h,sizeof(h),1,f)==1&&h[1]==1&&h[2]==1&&h[3]==sizeof(T)&&h[6]==24&&h[7]==384,"source packet header");std::vector<T>v(h[4]);need(fread(v.data(),sizeof(T),v.size(),f)==v.size()&&fgetc(f)==EOF,"source packet length");fclose(f);return v;}
int main(int argc,char**argv){try{need(argc==3,"run NEW audit directory");const double start=omp_get_wtime();std::string run=argv[1],dest=argv[2];auto packets=read_packet<SourcePacket>(run+"/native_source_packets_rank0.bin");need(packets.size()==465,"all465source overlap summaries");uint64_t all_points=0,complete=0,all_commits=0,failed=0;
 std::ofstream csv(dest+"/exact_native_prefix_accounting.csv");csv<<"slot_id,source,target,points,complete,checked_commits,final_FNV64,exact_source_enclosure_pass,reason\n";
 #pragma omp parallel for schedule(dynamic,1) reduction(+:all_points,complete,all_commits,failed)
 for(size_t i=0;i<packets.size();i++){const auto&p=packets[i];FILE*points=nullptr,*commits=nullptr;uint64_t count=0,checks=0,hash=14695981039346656037ULL;bool pass=false;std::string reason="PASS";try{
  need(p.planned_points<13074432&&p.points<=p.planned_points,"fixed count bound below int128 overflow limit");points=fopen((run+"/positive_slot"+std::to_string(p.slot)+".bin").c_str(),"rb");commits=fopen((run+"/tree_commits_slot"+std::to_string(p.slot)+".bin").c_str(),"rb");need(points&&commits,"retained positive and prefix commit files");need(fseek(commits,0,SEEK_END)==0&&ftell(commits)%sizeof(Commit)==0&&fseek(commits,0,SEEK_SET)==0,"complete commit record boundaries");setvbuf(points,nullptr,_IOFBF,1<<20);setvbuf(commits,nullptr,_IOFBF,1<<20);std::array<Bank,12>values,radii;B level[64][12];uint64_t occupied=0;Commit c;bool rejected=false;
  while(fread(&c,sizeof(c),1,commits)==1){need(c.version==1&&c.slot==p.slot&&c.next_query==c.points&&c.points>=count&&c.points<=p.points&&c.occupied==c.points&&!rejected,"committed exact query range and occupancy");if(c.successful)need(c.points==(c.panel+1)*256&&c.points-count==256,"complete original panel commit");else rejected=true;
   while(count<c.points){Row row;need(fread(&row,sizeof(row),1,points)==1&&row.source==p.source,"actual point in exact source prefix");const auto*bytes=reinterpret_cast<const unsigned char*>(&row);for(size_t at=0;at<sizeof(row);at++){hash^=bytes[at];hash*=1099511628211ULL;}B carry[12];Q local[12]={};for(int k=0;k<12;k++){values[k].add(row.moment[k]);radii[k].add(row.radius[k],true);carry[k]=B(static_cast<long double>(row.moment[k]),static_cast<long double>(row.radius[k]));}for(int k:{0,5})need(row.moment[k]>row.radius[k]&&row.radius[k]>=0,"point mass and physical U positive");
    int l=0;while(occupied&(uint64_t(1)<<l)){for(int k=0;k<12;k++){carry[k]=level[l][k]+carry[k];Q round=Q(material_forward::round_error(carry[k].n));local[k]=material_forward::plus_up(local[k],round);}occupied&=~(uint64_t(1)<<l);l++;need(l<64,"tree audit capacity");}for(int k=0;k<12;k++){level[l][k]=carry[k];need(local[k]==row.sum_local_rounding[k],"saved actual local tree-rounding diagnostic");}occupied|=uint64_t(1)<<l;count++;
   }
   need(hash==c.record_hash&&occupied==c.occupied,"positive prefix rolling checksum and tree occupancy");for(int l=0;l<64;l++)if(occupied&(uint64_t(1)<<l))for(int k=0;k<12;k++)need(Q(level[l][k].n)==c.level[l][k].n&&Q(level[l][k].e)==c.level[l][k].e,"committed full actual tree state");checks++;
  }need(feof(commits)&&fgetc(points)==EOF&&count==p.points&&checks>0,"exact committed point prefix, no trailing point/commit data");need(!p.complete||(!rejected&&p.points==p.planned_points),"complete overlap has no rejected final commit");
  for(int k=0;k<12;k++){B sum;for(int l=0;l<64;l++)if(occupied&(uint64_t(1)<<l))sum=sum+level[l][k];Bank lower=values[k];lower.merge(radii[k],-1);lower.add(Q(sum.e),true);lower.add(Q(sum.n),false,-1);Bank upper;upper.merge(values[k],-1);upper.merge(radii[k],-1);upper.add(Q(sum.e),true);upper.add(Q(sum.n));need(lower.sign()>=0&&upper.sign()>=0,"independent integer source value/radius enclosure");if(p.complete)need(Q(sum.n)==p.moment[k].n&&Q(sum.e)==p.moment[k].e,"published overlap matches exact committed final tree");}pass=true;
 }catch(const std::exception&e){reason=e.what();}if(points)fclose(points);if(commits)fclose(commits);all_points+=count;all_commits+=checks;complete+=p.complete;failed+=!pass;
  #pragma omp critical(gfs_integer_prefix_audit)
  {csv<<p.slot<<','<<p.source<<','<<p.target<<','<<count<<','<<p.complete<<','<<checks<<','<<hash<<','<<pass<<','<<reason<<'\n';fprintf(stderr,"native_integer_prefix_audit slot=%llu points=%llu complete=%llu pass=%d reason=%s\n",(unsigned long long)p.slot,(unsigned long long)count,(unsigned long long)p.complete,pass,reason.c_str());fflush(stderr);}
 }
 csv.close();printf("{\"scope\":\"all saved native point centers/radii exact integer source enclosure checks and every committed query range/FNV/tree state; native field and libm model retained, no field replay\",\"audited_overlap_slots\":465,\"completed_overlap_slots\":%llu,\"accounted_actual_points\":%llu,\"checked_tree_commits\":%llu,\"failed_accounting_slots\":%llu,\"independent_integer_source_enclosures_pass\":%s,\"new_field_or_physical_queries\":0,\"accepted_RK4_evolution\":false,\"long_Kepler_validated\":false,\"seconds\":%.17g}\n",(unsigned long long)complete,(unsigned long long)all_points,(unsigned long long)all_commits,(unsigned long long)failed,failed?"false":"true",omp_get_wtime()-start);return failed?1:0;
 }catch(const std::exception&e){fprintf(stderr,"native integer prefix audit FAIL %s\n",e.what());return 2;}}
