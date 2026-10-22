#ifndef MATERIAL_SHARED_REFERENCE80_INTEGRAL_H
#define MATERIAL_SHARED_REFERENCE80_INTEGRAL_H
#include "material_shared_target_integral.h"
#include "material_shared_reference80_ball.h"
namespace shared_target {
struct Reference80Accumulator:Accumulator {mf::B<long double> native_sum[12],level[64][12];uint64_t occupied=0;};
static void finish_reference80(Reference80Accumulator&out){for(int k=0;k<12;k++){mf::B<long double>sum;for(int level=0;level<64;level++)if(out.occupied&(uint64_t(1)<<level))sum=sum+out.level[level][k];out.native_sum[k]=sum;out.moment[k]=S(Q(sum.n),Q(sum.e));}}
// These are separate cumulative diagnostics with different units. They
// are not an additive decomposition of the final all12 error envelope.
struct Reference80Positive {uint64_t source;Q X[2],weight,input_radius[3],cache_radius_max[3],field_radius[17],moment[12],radius[12],sum_local_rounding[12];};
static void point_reference80(const curved::Packet&p,const sd::Compiled&m,const sd::Point&field,S w,const Thermal&thermal,Reference80Accumulator&out,FILE*positive){
    using R80=long double;namespace r80=shared_reference80_ball;
    static_assert(LDBL_MANT_DIG==64,"CPU reference requires binary80 mantissa");
    same_control(thermal,m.base.state);need(sf::lower(w)>0,"one positive CPU binary80 all12 reference packet");
    const auto&f=m.base.state->triangle[m.base.source];r80::PointInput q;q.triangle=m.base.source;q.X=R80(field.X[0].n);q.Y=R80(field.X[1].n);
    R80 dx=q.X-R80(f.vertex[0][0]),dy=q.Y-R80(f.vertex[0][1]),ux=R80(f.vertex[1][0])-R80(f.vertex[0][0]),uy=R80(f.vertex[1][1])-R80(f.vertex[0][1]),vx=R80(f.vertex[2][0])-R80(f.vertex[0][0]),vy=R80(f.vertex[2][1])-R80(f.vertex[0][1]),det=ux*vy-uy*vx;
    q.l[1]=(dx*vy-dy*vx)/det;q.l[2]=(ux*dy-uy*dx)/det;q.l[0]=1-q.l[1]-q.l[2];
    auto native=r80::evaluate(m.base,q,field.X[0],field.X[1]);
    auto X=r80::cache(q.X,field.X[0]),Y=r80::cache(q.Y,field.X[1]),weight=r80::cache(R80(w.n),w);
    auto actual=physical_query(p,native,X,Y,weight,thermal);auto oracle=physical_query(p,field.field,field.X[0],field.X[1],w,thermal);
    Reference80Positive row={};row.source=m.base.source;row.X[0]=field.X[0].n;row.X[1]=field.X[1].n;row.weight=w.n;row.input_radius[0]=Q(X.e);row.input_radius[1]=Q(Y.e);row.input_radius[2]=Q(weight.e);
    for(int a=0;a<3;a++)for(int k=0;k<2;k++)row.cache_radius_max[0]=std::max(row.cache_radius_max[0],Q(m.base.gradcache[a][k].e));
    for(int a=0;a<10;a++)for(int k=0;k<2;k++)row.cache_radius_max[1]=std::max(row.cache_radius_max[1],Q(m.base.qcache[a][k].e));
    for(int e=0;e<3;e++){for(int k=0;k<2;k++)row.cache_radius_max[2]=std::max(row.cache_radius_max[2],Q(m.base.edge[e].Dcache[k].e));for(int a=0;a<4;a++)for(int k=0;k<2;k++)row.cache_radius_max[2]=std::max(row.cache_radius_max[2],Q(m.base.edge[e].qcache[a][k].e));}
    auto compare_field=[&](mf::B<R80> a,S b){need(fabsq(Q(a.n)-b.n)<=sf::add(Q(a.e),b.e),"binary80 field graph encloses independent canonical value/first derivative");};for(int k=0;k<2;k++){compare_field(native.position[k],field.field.position[k]);compare_field(native.velocity[k],field.field.velocity[k]);}for(int k=0;k<4;k++){compare_field(native.F[k],field.field.F[k]);compare_field(native.V[k],field.field.V[k]);compare_field(native.grad_v[k],field.field.grad_v[k]);}compare_field(native.J,field.field.J);
    int j=0;for(const auto&v:native.position)row.field_radius[j++]=Q(v.e);for(const auto&v:native.velocity)row.field_radius[j++]=Q(v.e);for(const auto&v:native.F)row.field_radius[j++]=Q(v.e);for(const auto&v:native.V)row.field_radius[j++]=Q(v.e);row.field_radius[j++]=Q(native.J.e);for(const auto&v:native.grad_v)row.field_radius[j++]=Q(v.e);need(j==17,"reference field diagnostic schema");
    for(int k=0;k<12;k++){
        Q error=fabsq(Q(actual[k].n)-oracle[k].n),radius=sf::add(Q(actual[k].e),oracle[k].e);need(error<=radius,"actual binary80 shared graph/cache/all12 encloses Q128 same-function oracle");
        row.moment[k]=Q(actual[k].n);row.radius[k]=Q(actual[k].e);
    }
    // Streaming binary tree: actual binary80 additions with their forward
    // radii, O(log N) depth. No change to any point, weight or physics.
    mf::B<R80>carry[12];for(int k=0;k<12;k++)carry[k]=actual[k];int level=0;
    while(level<64&&(out.occupied&(uint64_t(1)<<level))){for(int k=0;k<12;k++){carry[k]=out.level[level][k]+carry[k];row.sum_local_rounding[k]=sf::add(row.sum_local_rounding[k],Q(mf::round_error(carry[k].n)));}out.occupied&=~(uint64_t(1)<<level);level++;}
    need(level<64,"reference summation tree capacity");for(int k=0;k<12;k++)out.level[level][k]=carry[k];out.occupied|=uint64_t(1)<<level;
    size_t written=0;
    #pragma omp critical(gfs_reference80_positive_packet)
    {written=fwrite(&row,sizeof(row),1,positive);}
    need(written==1,"binary80 positive packet with input/cache/field/moment/sum diagnostics");out.points++;
}
}
#endif
