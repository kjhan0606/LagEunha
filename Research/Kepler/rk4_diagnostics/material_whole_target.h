#ifndef MATERIAL_WHOLE_TARGET_H
#define MATERIAL_WHOLE_TARGET_H
#include "material_whole_target_geometry.h"
#include <cstring>

namespace whole_target {
namespace wg=whole_geometry;
struct Input {
    const curved::Packet& packet;
    const std::vector<curved_partition::Site>& sites;
    const std::vector<R>& weights;
    uint64_t epoch;
};
// CPU binary80 integration/Q128 geometry contract. Device code will have to
// supply the SAME positive node packet, all12 values AND absolute bounds.
// No device backend is registered by this reference implementation.
struct TargetPacket {uint64_t epoch=0,target=0;R moment[12]={},absolute_error[12]={};};
struct Contribution {
    uint64_t source=0,cell=0,points=0,outgoing=0;
    int classification=-1,complete=0;
    R moment[12]={},native_error[12]={};Q integral_error[12]={};
    char reason[192]={};
};
struct Value {mf::B<R> z[12];Q remainder[12]={};uint64_t points=0;};
struct Logs {std::ofstream&sectors;FILE*positive;};
static void point(const curved::Packet&p,const curved::Record&r,size_t source,uint64_t target,
                  const tg::Map&m,S X,S Y,S x,S y,S w,const tg::Inverse*inverse,
                  uint64_t loop,uint64_t boundary,Value&out,Logs&logs) {
    need(ti::lower(w)>0,"positive common all12 physical weight");
    auto native=tg::evaluate<R>(m,X,Y);auto oracle=mf::original_at(r,X,Y,p.parameters[3]);
    auto nn=physical_moment(p,native,profile(p,ti::cast(X),ti::cast(Y),true),ti::cast(x),ti::cast(y),ti::cast(w));
    auto oo=physical_moment(p,oracle,profile(p,X,Y,false),x,y,w);
    Positive record={};record.triangle=source;record.target=target;record.loop=loop;record.boundary=boundary;
    record.x[0]=x.n;record.x[1]=y.n;record.X[0]=X.n;record.X[1]=Y.n;record.weight=w.n;
    record.physical_radius[0]=x.e;record.physical_radius[1]=y.e;record.weight_radius=w.e;
    record.inverse_radius=inverse?inverse->radius:std::max(X.e,Y.e);
    for(int k=0;k<12;k++) {
        need(fabsq(Q(nn[k].n)-oo[k].n)<=ti::add(Q(nn[k].e),oo[k].e),"original basis/profile oracle contained by native forward bound");
        out.z[k]=out.z[k]+nn[k];record.moment[k]=nn[k].n;record.error[k]=nn[k].e;
    }
    need(fwrite(&record,sizeof(record),1,logs.positive)==1,"immutable positive-point record write");out.points++;
}
static Value intersection(const Input&in,size_t source,int target,const wg::Source&geometry,Logs&logs) {
    const auto&p=in.packet;const auto&r=p.record[source];const auto&m=geometry.map;
    auto d=td::build(m,r,in.sites,in.weights,target);Value out;auto gauss=ti::gauss16();
    for(size_t il=0;il<d.loop.size();il++)for(int ib:d.loop[il].boundary)for(int sector=0;sector<4;sector++) {
        const auto&loop=d.loop[il];const auto&boundary=d.boundary[ib];
        auto c=ta::sector(p,r,m,d,loop,boundary,sector);
        logs.sectors<<source<<','<<in.sites[target].id<<','<<il<<','<<ib<<','<<sector<<','<<c.pass<<','<<c.reason<<','<<qs(c.contraction)<<','<<qs(c.invariance_ratio)<<'\n';
        need(c.pass,c.reason);for(int k=0;k<12;k++)out.remainder[k]=ti::add(out.remainder[k],c.remainder[k]);
        for(int iu=0;iu<16;iu++)for(int iv=0;iv<16;iv++) {
            S u=(S(sector)+gauss.x[iu])/S(4),s=gauss.x[iv];auto curve=td::curve(d,m,r,boundary,u);
            S delta[2]={curve.x[0]-S(loop.anchor[0]),curve.x[1]-S(loop.anchor[1])},J=td::cross(delta,curve.d);
            S x=S(loop.anchor[0])+s*delta[0],y=S(loop.anchor[1])+s*delta[1],w=s*J*gauss.w[iu]*gauss.w[iv]/S(4);
            auto inverse=tg::inverse(m,r,x,y);need(inverse.inside,"every positive physical target node inside original source");
            point(p,r,source,in.sites[target].id,m,inverse.X,inverse.Y,x,y,w,&inverse,il,ib,out,logs);
        }
    }
    return out;
}
struct Piece {std::array<S,2> vertex[3];};
static std::array<Q,12> reference_remainder(const Input&in,size_t source,const Piece&piece) {
    const auto&p=in.packet;const auto&r=p.record[source];const Q radius=Q(.875L),h=Q(.5L)+radius;
    S D[2][2];Q anchor[2],rr[2];
    for(int k=0;k<2;k++) {
        D[0][k]=piece.vertex[1][k]-piece.vertex[0][k];D[1][k]=piece.vertex[2][k]-piece.vertex[0][k];
        S center=piece.vertex[0][k]+S(Q(.5L))*D[0][k]+S(Q(.25L))*D[1][k];anchor[k]=center.n;
        Q a=ti::mag(D[0][k]-S(Q(.5L))*D[1][k]),b=ti::mag(D[1][k]);
        rr[k]=ti::add(center.e,std::max(ti::add(ti::mul(h,a),ti::mul(ti::add(Q(.25L),ti::mul(h,Q(.5L))),b)),ti::add(ti::mul(Q(.5L),a),ti::mul(h,b))));
    }
    auto field=ta::bounds(tg::compile(r,anchor,p.parameters[3]),p,rr);need(field.pass,field.reason);
    S jac=td::cross(D[0],D[1]),f16(1),f32(1);need(ti::lower(jac)>0,"positive reference red subtriangle");
    for(int k=2;k<=32;k++){f32=f32*S(k);if(k<=16)f16=f16*S(k);}S ratio=f16*f16/f32,C=ratio*ratio/S(33)*mf::power(S(radius),Q(-32));
    std::array<Q,12> out;Q coefficient=ti::mul(Q(2),ti::mag(C));
    for(int k=0;k<12;k++)out[k]=ti::mul(coefficient,ti::mul(ti::mul(Q(1)+radius,ti::mag(jac)),ti::mul(field.moment[k],field.Jhi)));
    return out;
}
static Value baseline(const Input&in,size_t source,Logs&logs) {
    const auto&p=in.packet;const auto&r=p.record[source];std::array<S,2> v[3],mid[3];
    for(int a=0;a<3;a++)for(int k=0;k<2;k++)v[a][k]=S(Q(r.triangle.geometry.x[a][k]));
    for(int a=0;a<3;a++)for(int k=0;k<2;k++)mid[a][k]=(v[a][k]+v[(a+1)%3][k])/S(2);
    std::vector<Piece> pieces={{{v[0],mid[0],mid[2]}},{{v[1],mid[1],mid[0]}},{{v[2],mid[2],mid[1]}},{{mid[0],mid[1],mid[2]}}};
    auto gauss=ti::gauss16();Value out;
    for(size_t ip=0;ip<pieces.size();ip++) {
        const auto&piece=pieces[ip];auto error=reference_remainder(in,source,piece);for(int k=0;k<12;k++)out.remainder[k]=ti::add(out.remainder[k],error[k]);
        logs.sectors<<source<<','<<r.cell<<",baseline,"<<ip<<",0,true,pass,0,0\n";
        Q anchor[2];S D[2][2];for(int k=0;k<2;k++){S c;for(int a=0;a<3;a++)c=c+piece.vertex[a][k]/S(3);anchor[k]=c.n;D[0][k]=piece.vertex[1][k]-piece.vertex[0][k];D[1][k]=piece.vertex[2][k]-piece.vertex[0][k];}
        auto m=tg::compile(r,anchor,p.parameters[3]);S jac=td::cross(D[0],D[1]);
        for(int iu=0;iu<16;iu++)for(int iv=0;iv<16;iv++) {
            S u=gauss.x[iu],v=gauss.x[iv],X=piece.vertex[0][0]+u*D[0][0]+(S(1)-u)*v*D[1][0],Y=piece.vertex[0][1]+u*D[0][1]+(S(1)-u)*v*D[1][1],w=jac*(S(1)-u)*gauss.w[iu]*gauss.w[iv];
            auto q=tg::evaluate<Q>(m,X,Y);
            point(p,r,source,r.cell,m,X,Y,q.position[0],q.position[1],w*q.J,nullptr,ip,UINT64_MAX,out,logs);
        }
    }
    return out;
}
static Contribution integrate(const Input&in,size_t source,int target,const wg::Source&g,Logs&logs) {
    Contribution row={};row.source=source;row.cell=in.packet.record[source].cell;
    try {
        row.classification=wg::classify(g,in.sites,in.weights,target);
        if(row.classification==0){row.complete=1;std::strcpy(row.reason,"whole_image_plane_exclusion");return row;}
        Value v;
        if(row.classification==1)v=baseline(in,source,logs);
        else if(row.cell==in.sites[target].id) {
            // Large own component: complete original source minus ALL other
            // target intersections. Certified zero terms use whole-image
            // half-plane exclusion, never a retained flow table.
            v=baseline(in,source,logs);
            for(int j=0;j<int(in.sites.size());j++)if(j!=target&&wg::classify(g,in.sites,in.weights,j)!=0) {
                Value outgoing=intersection(in,source,j,g,logs);row.outgoing++;
                for(int k=0;k<12;k++){v.z[k]=v.z[k]-outgoing.z[k];v.remainder[k]=ti::add(v.remainder[k],outgoing.remainder[k]);}v.points+=outgoing.points;
            }
            need(Q(v.z[MASS].n)-Q(v.z[MASS].e)-v.remainder[MASS]>0,"complete own complement positive mass");
            need(Q(v.z[THERMAL].n)-Q(v.z[THERMAL].e)-v.remainder[THERMAL]>0,"complete own complement positive physical U");
        }else v=intersection(in,source,target,g,logs);
        for(int k=0;k<12;k++){row.moment[k]=v.z[k].n;row.native_error[k]=v.z[k].e;row.integral_error[k]=v.remainder[k];}
        row.points=v.points;row.complete=1;std::strcpy(row.reason,v.points?"all12_original_pushforward_integrated":"complete_empty_boundary_arrangement");
    }catch(const std::exception&e){std::snprintf(row.reason,sizeof(row.reason),"%s",e.what());}
    return row;
}
struct Report {
    bool published=false,epoch_valid=false,sources_complete=false,geometry_complete=false,budgets_pass=false,positive=false;
    uint64_t failed_sources=0,points=0,excluded=0,inside=0,intersected=0;
    Q error[12]={},tolerance[12]={},integral[12]={},native[12]={},sum[12]={};
    TargetPacket candidate;
    wg::Audit geometry;
    std::string witness;
};
static Report assemble(const Input&in,uint64_t expected_epoch,int target,const std::vector<Contribution>&rows,
                       const std::vector<wg::Source>&geometry,const R*normalizers,
                       const std::string&dest,TargetPacket&caller) {
    Report report;
    if(in.epoch!=expected_epoch){report.witness="immutable_input_epoch_mismatch";return report;}
    report.epoch_valid=true;report.candidate.epoch=in.epoch;report.candidate.target=in.sites[target].id;
    if(rows.size()!=in.packet.record.size()){report.witness="missing_source_rows";return report;}
    report.sources_complete=true;mf::B<R> sum[12];Q vscale=__builtin_huge_valq(),rscale=vscale;
    for(size_t s=0;s<rows.size();s++) {
        const auto&r=rows[s];if(r.source!=s||!r.complete){report.sources_complete=false;report.failed_sources++;if(report.witness.empty())report.witness="source_"+std::to_string(s)+":"+r.reason;continue;}
        report.excluded+=r.classification==0;report.inside+=r.classification==1;report.intersected+=r.classification==2;report.points+=r.points;
        if(r.classification!=0){auto q=tg::evaluate<Q>(geometry[s].map,geometry[s].reference[0],geometry[s].reference[1]);Q vv[2]={ti::mag(q.velocity[0]),ti::mag(q.velocity[1])};vscale=std::min(vscale,ti::mag(mf::root(S(vv[0])*S(vv[0])+S(vv[1])*S(vv[1]))));Q radius=0;for(int a=0;a<3;a++){S x(Q(in.packet.record[s].triangle.geometry.x[a][0])),y(Q(in.packet.record[s].triangle.geometry.x[a][1]));radius=std::max(radius,ti::mag(mf::root(x*x+y*y)));}rscale=std::min(rscale,radius);}
        for(int k=0;k<12;k++){sum[k]=sum[k]+mf::B<R>(r.moment[k]);report.native[k]=ti::add(report.native[k],Q(r.native_error[k]));report.integral[k]=ti::add(report.integral[k],r.integral_error[k]);}
    }
    try{report.geometry=wg::audit(in.packet,geometry,in.sites,in.weights,target,dest);report.geometry_complete=true;}
    catch(const std::exception&e){if(report.witness.empty())report.witness=std::string("geometry:")+e.what();}
    // Same 1e-12 moment / 1e-10 physical energy / 1e-9 Udot
    // coefficients. Energy U lower is read from the UNCHANGED old proof.
    // For the new complete nonenergy budget use the MINIMUM of relevant
    // old source scale formulas, never a generous packet-wide maximum.
    for(int k=0;k<12;k++)report.tolerance[k]=ti::lower(S(Q(normalizers[1]))/S(Q(1000000000000LL)));
    report.tolerance[PX]=report.tolerance[PY]=ti::lower(S(report.tolerance[MASS])*S(vscale));
    report.tolerance[ANGULAR]=ti::lower(S(report.tolerance[MASS])*S(rscale)*S(vscale));
    for(int k:{KINETIC,THERMAL,POTENTIAL})report.tolerance[k]=ti::lower(S(Q(normalizers[0]))/S(Q(10000000000LL)));
    report.tolerance[VOLUME]=ti::lower(S(Q(normalizers[2]))/S(Q(1000000000000LL)));
    report.tolerance[QX]=report.tolerance[QY]=ti::lower(S(report.tolerance[MASS])*S(rscale));
    report.tolerance[POLAR]=ti::lower(S(report.tolerance[MASS])*S(rscale)*S(rscale));
    report.tolerance[THERMAL_RATE]=ti::lower(S(Q(normalizers[0]))/S(Q(1000000000)));
    report.budgets_pass=report.sources_complete&&report.geometry_complete;
    for(int k=0;k<12;k++) {
        report.candidate.moment[k]=sum[k].n;report.sum[k]=Q(sum[k].e);
        report.error[k]=ti::add(ti::add(report.native[k],report.integral[k]),ti::add(report.sum[k],ti::add(report.geometry.seam[k],report.geometry.cap[k])));
        auto e=ti::cast(S(report.error[k]));report.candidate.absolute_error[k]=mf::FP<R>::up(e.n+e.e);
        if(report.error[k]>report.tolerance[k]){report.budgets_pass=false;if(report.witness.empty())report.witness=std::string("whole_target_budget:")+moment_names[k];}
    }
    report.positive=Q(report.candidate.moment[MASS])-report.error[MASS]>0&&Q(report.candidate.moment[THERMAL])-report.error[THERMAL]>0;
    if(!report.positive&&report.witness.empty())report.witness="target_mass_or_physical_U_not_certified_positive";
    if(report.sources_complete&&report.geometry_complete&&report.budgets_pass&&report.positive){caller=report.candidate;report.published=true;}
    return report;
}
}
#endif
