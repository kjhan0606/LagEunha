#ifndef MATERIAL_WHOLE_TARGET_GEOMETRY_H
#define MATERIAL_WHOLE_TARGET_GEOMETRY_H
#include "material_target_gluing.h"

// CPU reference certificates for ONE target of the unchanged finite packet.
// The reference is the sum of original source push-forwards. No gap is filled
// and no overlapping original source measure is removed.
namespace whole_geometry {
namespace gl=material_target_gluing; namespace ta=material_target_analytic;
struct Source {
    tg::Map map;
    S reference[2],image[2],displacement[2];
    Q L=0,L2=0;
};
static Source source(const curved::Record&r,double eps) {
    Source z; Q anchor[2];
    for(int k=0;k<2;k++) {
        Q lo=r.triangle.geometry.x[0][k],hi=lo;
        for(int a=1;a<3;a++){lo=std::min(lo,Q(r.triangle.geometry.x[a][k]));hi=std::max(hi,Q(r.triangle.geometry.x[a][k]));}
        S c=(S(lo)+S(hi))/S(2),h=(S(hi)-S(lo))/S(2);
        z.reference[k]=S(c.n,ti::add(c.e,ti::mag(h)));anchor[k]=c.n;
    }
    z.map=tg::compile(r,anchor,eps);
    auto p=tg::evaluate<Q>(z.map,z.reference[0],z.reference[1]);
    S dx=z.reference[0]-S(anchor[0]),dy=z.reference[1]-S(anchor[1]);
    S omega=mf::power(z.reference[0]*z.reference[0]+z.reference[1]*z.reference[1]+S(Q(eps))*S(Q(eps)),Q(-.75L));
    S initial[2]={-z.reference[1]*omega,z.reference[0]*omega};
    for(int k=0;k<2;k++) {
        z.image[k]=p.position[k]; auto d=z.map.x[k];
        d.c[0][0]=d.c[0][0]-S(anchor[k]);d.c[k==0?1:0][k==1?1:0]=d.c[k==0?1:0][k==1?1:0]-S(1);
        z.displacement[k]=ti::value(d,dx,dy);
        for(int j=0;j<2;j++)z.displacement[k]=z.displacement[k]+z.map.Bx[k][j]*initial[j];
        Q row=0;for(int j=0;j<2;j++)row=ti::add(row,ti::mag(p.F[2*k+j]-S(k==j?1:0)));
        z.L=std::max(z.L,row);
    }
    z.L2=ti::mul(ti::mag(mf::root(S(2))),z.L);
    need(z.L2<1,"original whole source injectivity (Euclidean norm)");return z;
}
static int classify(const Source&s,const std::vector<curved_partition::Site>&site,const std::vector<R>&w,int target) {
    bool inside=true;Q origin[2]={s.image[0].n,s.image[1].n};
    for(int j=0;j<int(site.size());j++)if(j!=target){S f=tg::value(tg::plane(site,w,target,j,origin),s.image[0],s.image[1],origin);if(ti::lower(f)>0)return 0;inside=inside&&ti::upper(f)<0;}
    return inside?1:2; // Excluded / whole inclusion / needs exact arrangement.
}
static bool excluded(const S*box,const std::vector<curved_partition::Site>&site,const std::vector<R>&w,int target) {
    Q origin[2]={box[0].n,box[1].n};
    for(int j=0;j<int(site.size());j++)if(j!=target&&ti::lower(tg::value(tg::plane(site,w,target,j,origin),box[0],box[1],origin))>0)return true;
    return false;
}
static S interval(Q lo,Q hi){S c=(S(lo)+S(hi))/S(2),h=(S(hi)-S(lo))/S(2);return S(c.n,ti::add(c.e,ti::mag(h)));}
static Pt vertex(const curved::Record&r,int a){return {r.triangle.geometry.x[a][0],r.triangle.geometry.x[a][1]};}
static Q norm(const S*x){return gl::norm2(x[0],x[1]);}
static std::map<Key,std::vector<Edge>> edges(const curved::Packet&p){std::map<Key,std::vector<Edge>> out;for(size_t s=0;s<p.record.size();s++)for(int a=0;a<3;a++){Pt x=vertex(p.record[s],a),y=vertex(p.record[s],(a+1)%3);out[x<y?Key(x,y):Key(y,x)].push_back({s,p.record[s].component,x,y});}return out;}
struct Audit {
    Q seam[12]={},cap[12]={},duplicated_seam[12]={};
    uint64_t paired=0,seam_parts=0,external=0,relevant_external=0,pairs=0,bbox_excluded=0,target_excluded=0,edge_pairs=0,vertex_pairs=0,separated_pairs=0;
};
static bool convex_union(const curved::Record&a,const curved::Record&b,const Key&e) {
    Pt u={},v={};for(int i=0;i<3;i++){Pt p=vertex(a,i);if(p!=e.first&&p!=e.second)u=p;p=vertex(b,i);if(p!=e.first&&p!=e.second)v=p;}
    S D[2],U[2],V[2],W[2],DU[2];for(int k=0;k<2;k++){D[k]=tg::exact_difference(Q(e.second[k]),Q(e.first[k]));U[k]=tg::exact_difference(Q(u[k]),Q(e.first[k]));V[k]=tg::exact_difference(Q(v[k]),Q(e.first[k]));W[k]=tg::exact_difference(Q(u[k]),Q(v[k]));DU[k]=tg::exact_difference(Q(e.second[k]),Q(u[k]));}
    auto cross_exact=[](const S*a,const S*b){return td::exact_det(a[0],a[1],b[0],b[1]);};
    S ca=cross_exact(D,U),cb=cross_exact(D,V);if(!(ti::lower(ca)>0&&ti::upper(cb)<0||ti::upper(ca)<0&&ti::lower(cb)>0))return false;
    // Opposite diagonal intersects the CLOSED common segment. A degenerate
    // quadrilateral may be a convex triangle; exact zero is not uncertainty.
    S den=cross_exact(D,W),num=cross_exact(U,W),remaining=cross_exact(DU,W);
    auto nonnegative=[](S x){return x.n==0&&x.e==0||ti::lower(x)>=0;};
    return ti::lower(den)>0&&nonnegative(num)&&nonnegative(remaining)||ti::upper(den)<0&&nonnegative(-num)&&nonnegative(-remaining);
}
static Q local_L(const tg::Map&m,const Q*anchor,const Q*rr){auto p=tg::evaluate<Q>(m,S(anchor[0],rr[0]),S(anchor[1],rr[1]));Q L=0;for(int i=0;i<2;i++){Q row=0;for(int j=0;j<2;j++)row=ti::add(row,ti::mag(p.F[2*i+j]-S(i==j?1:0)));L=std::max(L,row);}return L;}
static void add_fields(Audit&out,Q area,const curved::Packet&p,const curved::Record&a,const curved::Record&b,const Q*anchor,const Q*rr,bool seam) {
    auto ma=tg::compile(a,anchor,p.parameters[3]),mb=tg::compile(b,anchor,p.parameters[3]);auto fa=ta::bounds(ma,p,rr),fb=ta::bounds(mb,p,rr);
    need(fa.pass&&fb.pass,"canonical ambiguity field envelopes complete");
    Q*dest=seam?out.seam:out.cap;
    for(int k=0;k<12;k++){Q error=ti::mul(area,ti::add(fa.moment[k],fb.moment[k]));dest[k]=ti::add(dest[k],error);if(seam)out.duplicated_seam[k]=ti::add(out.duplicated_seam[k],ti::mul(Q(2),error));}
}
// A separated reference triangle pair has separated image projections if
// its SAT gap exceeds the ORIGINAL displacement projection difference.
static bool separated(const curved::Record&a,const curved::Record&b,const Source&sa,const Source&sb) {
    for(const auto*triangle:{&a,&b})for(int edge=0;edge<3;edge++)for(int sign:{-1,1}){
        Pt x=vertex(*triangle,edge),y=vertex(*triangle,(edge+1)%3);S n[2]={S(sign)*(S(Q(y[1]))-S(Q(x[1]))),S(-sign)*(S(Q(y[0]))-S(Q(x[0])))};
        Q amax=-__builtin_huge_valq(),bmin=__builtin_huge_valq();
        for(int v=0;v<3;v++){S va=n[0]*S(Q(a.triangle.geometry.x[v][0]))+n[1]*S(Q(a.triangle.geometry.x[v][1])),vb=n[0]*S(Q(b.triangle.geometry.x[v][0]))+n[1]*S(Q(b.triangle.geometry.x[v][1]));amax=std::max(amax,ti::upper(va));bmin=std::min(bmin,ti::lower(vb));}
        S da=n[0]*sa.displacement[0]+n[1]*sa.displacement[1],db=n[0]*sb.displacement[0]+n[1]*sb.displacement[1];
        if(ti::lower(S(bmin)-S(amax)+db-da)>0)return true;
    }return false;
}
static Audit audit(const curved::Packet&p,const std::vector<Source>&source,const std::vector<curved_partition::Site>&sites,const std::vector<R>&weights,int target,const std::string&dest) {
    Audit out;auto paired=edges(p);
    std::ofstream rows(dest+"/canonical_seams.csv"),pairs(dest+"/nonadjacent_pairs.csv"),external(dest+"/external_edges.csv");
    rows<<"left,right,part,position_gap,tube_radius,tube_area,possible_target,local_extension_complete\n";
    pairs<<"left,right,shared_vertices,certificate,cap_radius,cap_area\n";external<<"source,edge_x0,edge_y0,edge_x1,edge_y1,possibly_touches_target,policy\n";
    for(auto&entry:paired){const auto&e=entry.second;need(e.size()==1||e.size()==2,"manifold reference edge incidence");if(e.size()==1){out.external++;size_t s=e[0].triangle;bool possible=classify(source[s],sites,weights,target)!=0;out.relevant_external+=possible;external<<s<<','<<entry.first.first[0]<<','<<entry.first.first[1]<<','<<entry.first.second[0]<<','<<entry.first.second[1]<<','<<possible<<",original_pushforward_boundary_no_fill\n";continue;}
        out.paired++;size_t ia=e[0].triangle,ib=e[1].triangle;const auto&a=p.record[ia];const auto&b=p.record[ib];
        S A[2]={S(Q(entry.first.first[0])),S(Q(entry.first.first[1]))},D[2]={S(Q(entry.first.second[0]))-A[0],S(Q(entry.first.second[1]))-A[1]};Q anchor[2]={(A[0]+D[0]/S(2)).n,(A[1]+D[1]/S(2)).n};
        auto ma=tg::compile(a,anchor,p.parameters[3]),mb=tg::compile(b,anchor,p.parameters[3]);auto delta=gl::trace_difference(ma,mb,D);Q L=std::max(source[ia].L2,source[ib].L2);
        for(int part=0;part<16;part++){
            S t((Q(part)+Q(.5L))/16-Q(.5L),Q(1)/32),X=S(anchor[0])+D[0]*t,Y=S(anchor[1])+D[1]*t;
            auto left=tg::evaluate<Q>(ma,X,Y),right=tg::evaluate<Q>(mb,X,Y);S z=X*X+Y*Y+S(Q(p.parameters[3]))*S(Q(p.parameters[3])),omega=mf::power(z,Q(-.75L)),v0[2]={-Y*omega,X*omega},gap[2];
            for(int k=0;k<2;k++){gap[k]=gl::value(delta.p[k],t);for(int j=0;j<2;j++)gap[k]=gap[k]+delta.B[k][j]*v0[j];}
            Q dx=norm(gap),radius=ta::quotient_up(ti::mul(ti::add(1,L),dx),ta::difference_down(1,L));S physical[2];
            for(int k=0;k<2;k++)physical[k]=interval(mf::FP<Q>::down(std::min(ti::lower(left.position[k]),ti::lower(right.position[k]))-radius),mf::FP<Q>::up(std::max(ti::upper(left.position[k]),ti::upper(right.position[k]))+radius));
            bool possible=!excluded(physical,sites,weights,target);if(!possible)continue;
            need(convex_union(a,b,entry.first),"paired triangle union convex for overlap collar coverage");
            S tangent[2];for(int k=0;k<2;k++)tangent[k]=left.F[2*k]*D[0]+left.F[2*k+1]*D[1];
            // Tubular neighbourhood of the canonical original edge. Endpoint
            // caps intentionally overlap; no cancellation or min-of-audits.
            Q area=ti::add(ti::mul(ti::mul(Q(2),radius),ti::mul(norm(tangent),Q(1)/16)),ti::mul(Q(8),ti::mul(radius,radius)));
            Q pad=ta::quotient_up(ti::mul(Q(2),ti::add(radius,dx)),ta::difference_down(1,std::max(source[ia].L,source[ib].L))),local[2]={X.n,Y.n},rr[2]={ti::add(X.e,pad),ti::add(Y.e,pad)};
            Q extension_limit=ti::lower((S(1)+S(std::max(source[ia].L,source[ib].L)))/S(2));
            need(std::max(local_L(ma,local,rr),local_L(mb,local,rr))<extension_limit,"canonical tube source extensions cover both fields");
            add_fields(out,area,p,a,b,local,rr,true);out.seam_parts++;rows<<ia<<','<<ib<<','<<part<<','<<qs(dx)<<','<<qs(radius)<<','<<qs(area)<<",true,true\n";
        }
    }
    for(size_t ia=0;ia<p.record.size();ia++)for(size_t ib=ia+1;ib<p.record.size();ib++){
        out.pairs++;const auto&a=p.record[ia];const auto&b=p.record[ib];std::vector<Pt> shared;for(int i=0;i<3;i++)for(int j=0;j<3;j++)if(vertex(a,i)==vertex(b,j))shared.push_back(vertex(a,i));
        S box[2];bool disjoint=false;
        for(int k=0;k<2;k++){
            // Shared reference vertices can leave a central gap even when
            // their images are disjoint. Use the UNION box for those pairs;
            // intersection-box exclusion applies only to nonadjacent images.
            Q lo=shared.empty()?std::max(ti::lower(source[ia].image[k]),ti::lower(source[ib].image[k])):std::min(ti::lower(source[ia].image[k]),ti::lower(source[ib].image[k]));
            Q hi=shared.empty()?std::min(ti::upper(source[ia].image[k]),ti::upper(source[ib].image[k])):std::max(ti::upper(source[ia].image[k]),ti::upper(source[ib].image[k]));
            if(lo>hi){disjoint=true;break;}box[k]=interval(lo,hi);
        }
        if(disjoint){out.bbox_excluded++;continue;}if(excluded(box,sites,weights,target)){out.target_excluded++;continue;}
        if(shared.size()==2){out.edge_pairs++;need(convex_union(a,b,shared[0]<shared[1]?Key(shared[0],shared[1]):Key(shared[1],shared[0])),"all relevant adjacent overlaps covered by canonical collars");pairs<<ia<<','<<ib<<",2,canonical_edge_tube,0,0\n";continue;}
        if(shared.empty()){need(separated(a,b,source[ia],source[ib]),"nonadjacent original source images unresolved: no displacement SAT certificate");out.separated_pairs++;pairs<<ia<<','<<ib<<",0,displacement_SAT,0,0\n";continue;}
        need(shared.size()==1,"distinct source triangles");Pt v=shared[0];S rays[2][2][2];int counts[2]={};for(int side=0;side<2;side++){const auto&r=side?b:a;for(int j=0;j<3;j++)if(vertex(r,j)!=v){int q=counts[side]++;for(int k=0;k<2;k++)rays[side][q][k]=S(Q(r.triangle.geometry.x[j][k]))-S(Q(v[k]));}}
        Q alpha=0;for(int i=0;i<2&&alpha==0;i++)for(int j=0;j<2&&alpha==0;j++)for(int sign:{-1,1}){
            S la=mf::root(rays[0][i][0]*rays[0][i][0]+rays[0][i][1]*rays[0][i][1]),lb=mf::root(rays[1][j][0]*rays[1][j][0]+rays[1][j][1]*rays[1][j][1]);S h[2]={rays[0][i][0]/la+rays[1][j][0]/lb,rays[0][i][1]/la+rays[1][j][1]/lb};
            if(h[0].n==0&&h[1].n==0)continue;Q len=hypotq(h[0].n,h[1].n);S n[2]={S(sign*(-h[1].n/len)),S(sign*h[0].n/len)};Q ca=__builtin_huge_valq(),cb=ca;
            for(int q=0;q<2;q++){S norm_a=mf::root(rays[0][q][0]*rays[0][q][0]+rays[0][q][1]*rays[0][q][1]),norm_b=mf::root(rays[1][q][0]*rays[1][q][0]+rays[1][q][1]*rays[1][q][1]);ca=std::min(ca,ti::lower(-(n[0]*rays[0][q][0]+n[1]*rays[0][q][1])/norm_a));cb=std::min(cb,ti::lower((n[0]*rays[1][q][0]+n[1]*rays[1][q][1])/norm_b));}
            Q norm_n=norm(n),value=std::min(ta::difference_down(ca,ti::mul(norm_n,source[ia].L2)),ta::difference_down(cb,ti::mul(norm_n,source[ib].L2)));if(value>0){alpha=mf::FP<Q>::down(value/norm_n);break;}
        }
        need(alpha>0,"vertex-only reference wedges separate faster than original displacement");
        auto pa=tg::evaluate<Q>(source[ia].map,S(Q(v[0])),S(Q(v[1]))),pb=tg::evaluate<Q>(source[ib].map,S(Q(v[0])),S(Q(v[1])));S delta[2]={pa.position[0]-pb.position[0],pa.position[1]-pb.position[1]};Q dx=norm(delta),L=std::max(source[ia].L2,source[ib].L2),radius=ti::add(ta::quotient_up(ti::mul(ti::add(1,L),dx),alpha),dx);S capbox[2];for(int k=0;k<2;k++)capbox[k]=S(pa.position[k].n,ti::add(pa.position[k].e,radius));
        if(excluded(capbox,sites,weights,target)){out.vertex_pairs++;pairs<<ia<<','<<ib<<",1,vertex_cap_excluded,"<<qs(radius)<<",0\n";continue;}
        Q area=ti::mul(Q(8),ti::mul(radius,radius)),anchor[2]={Q(v[0]),Q(v[1])},pad=ta::quotient_up(ti::mul(Q(2),ti::add(radius,dx)),ta::difference_down(1,std::max(source[ia].L,source[ib].L))),rr[2]={pad,pad};
        need(std::max(local_L(source[ia].map,anchor,rr),local_L(source[ib].map,anchor,rr))<ti::lower((S(1)+S(std::max(source[ia].L,source[ib].L)))/S(2)),"vertex cap source extensions complete");
        add_fields(out,area,p,a,b,anchor,rr,false);out.vertex_pairs++;pairs<<ia<<','<<ib<<",1,vertex_cap,"<<qs(radius)<<','<<qs(area)<<'\n';
    }
    return out;
}
}
#endif
