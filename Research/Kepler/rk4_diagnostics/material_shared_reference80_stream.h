#ifndef MATERIAL_SHARED_REFERENCE80_STREAM_H
#define MATERIAL_SHARED_REFERENCE80_STREAM_H
#include "material_shared_reference80_integral.h"
namespace shared_target {
// This is the existing binary80 field/query/tree graph with geometry inputs
// supplied directly. The expensive Q128 Hessian oracle is optional; the
// native forward error graph and all positivity/denominator guards are not.
struct StreamOracle {uint64_t slot,source,target,panel,iu,iv,point_index,passed;S native_field[17],oracle_field[17],native_moment[12],oracle_moment[12];};
static_assert(sizeof(StreamOracle)==1920,"native streaming oracle ABI");
template<class T>static std::array<mf::B<T>,17>stream_fields(const mf::Point<T>&p){std::array<mf::B<T>,17>a;int j=0;for(auto x:p.position)a[j++]=x;for(auto x:p.velocity)a[j++]=x;for(auto x:p.F)a[j++]=x;for(auto x:p.V)a[j++]=x;a[j++]=p.J;for(auto x:p.grad_v)a[j++]=x;return a;}
static bool stream_reference80(const curved::Packet&p,const sf::Map&m,S original_X,S original_Y,S original_weight,const Thermal&thermal,const Q*cache_radius,Reference80Accumulator&out,FILE*positive,const sd::Point*oracle=nullptr,StreamOracle*proof=nullptr,const S*reference_lambda=nullptr,uint64_t*record_hash=nullptr){
 namespace r80=shared_reference80_ball;using T=long double;
 same_control(thermal,m.state);need(sf::lower(original_weight)>0,"positive whole-panel reference weight");r80::PointInput q;q.triangle=m.source;q.X=T(original_X.n);q.Y=T(original_Y.n);const auto&f=m.state->triangle[m.source];
 if(reference_lambda){for(int j=0;j<3;j++)q.l[j]=T(reference_lambda[j].n);}else{T dx=q.X-T(f.vertex[0][0]),dy=q.Y-T(f.vertex[0][1]),ux=T(f.vertex[1][0])-T(f.vertex[0][0]),uy=T(f.vertex[1][1])-T(f.vertex[0][1]),vx=T(f.vertex[2][0])-T(f.vertex[0][0]),vy=T(f.vertex[2][1])-T(f.vertex[0][1]),det=ux*vy-uy*vx;
 q.l[1]=(dx*vy-dy*vx)/det;q.l[2]=(ux*dy-uy*dx)/det;q.l[0]=1-q.l[1]-q.l[2];}
 auto native=r80::evaluate(m,q,original_X,original_Y,reference_lambda);auto X=r80::cache(q.X,original_X),Y=r80::cache(q.Y,original_Y),weight=r80::cache(T(original_weight.n),original_weight);
 auto actual=physical_query(p,native,X,Y,weight,thermal);Reference80Positive row={};row.source=m.source;row.X[0]=original_X.n;row.X[1]=original_Y.n;row.weight=original_weight.n;row.input_radius[0]=Q(X.e);row.input_radius[1]=Q(Y.e);row.input_radius[2]=Q(weight.e);for(int j=0;j<3;j++)row.cache_radius_max[j]=cache_radius[j];auto nf=stream_fields(native);for(int j=0;j<17;j++)row.field_radius[j]=Q(nf[j].e);
 bool valid=true;if(oracle){need(proof,"oracle proof destination present");auto of=stream_fields(oracle->field);auto om=physical_query(p,oracle->field,original_X,original_Y,original_weight,thermal);for(int j=0;j<17;j++){proof->native_field[j]=S(Q(nf[j].n),Q(nf[j].e));proof->oracle_field[j]=of[j];valid=valid&&fabsq(Q(nf[j].n)-of[j].n)<=sf::add(Q(nf[j].e),of[j].e);}for(int k=0;k<12;k++){proof->native_moment[k]=S(Q(actual[k].n),Q(actual[k].e));proof->oracle_moment[k]=om[k];valid=valid&&fabsq(Q(actual[k].n)-om[k].n)<=sf::add(Q(actual[k].e),om[k].e);}proof->passed=valid;}
 mf::B<T>carry[12];for(int k=0;k<12;k++){row.moment[k]=Q(actual[k].n);row.radius[k]=Q(actual[k].e);carry[k]=actual[k];}int level=0;
 while(level<64&&(out.occupied&(uint64_t(1)<<level))){for(int k=0;k<12;k++){carry[k]=out.level[level][k]+carry[k];row.sum_local_rounding[k]=sf::add(row.sum_local_rounding[k],Q(mf::round_error(carry[k].n)));}out.occupied&=~(uint64_t(1)<<level);level++;}
 need(level<64,"native streaming binary80 tree capacity");for(int k=0;k<12;k++)out.level[level][k]=carry[k];out.occupied|=uint64_t(1)<<level;
 need(fwrite(&row,sizeof(row),1,positive)==1,"native streaming positive record retained");if(record_hash){const auto*bytes=reinterpret_cast<const unsigned char*>(&row);for(size_t i=0;i<sizeof(row);i++){*record_hash^=bytes[i];*record_hash*=1099511628211ULL;}}out.points++;return valid;
}
}
#endif
