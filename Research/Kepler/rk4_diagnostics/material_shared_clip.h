#ifndef MATERIAL_SHARED_CLIP_H
#define MATERIAL_SHARED_CLIP_H
#include "material_shared_atlas.h"
namespace shared_target {
struct Axis {int dependent=-1;bool empty=false,whole=false;std::vector<int>active;std::vector<S>derivative;char reason[192]={};};
static Axis axis(const sd::Compiled&m,const sd::Chart&chart,const std::vector<Plane>&planes,const sf::RealBounds&bound){Axis out;out.derivative.resize(planes.size());S reference[3][2];for(int a=0;a<3;a++)for(int j=0;j<3;j++)for(int k=0;k<2;k++)reference[a][k]=reference[a][k]+chart.vertex[a][j]*m.base.vertex[j][k];for(size_t i=0;i<planes.size();i++){Q lo=sf::lower(value(planes[i],reference[0][0],reference[0][1])),hi=sf::upper(value(planes[i],reference[0][0],reference[0][1])),pad=0;for(int a=1;a<3;a++){S f=value(planes[i],reference[a][0],reference[a][1]);lo=std::min(lo,sf::lower(f));hi=std::max(hi,sf::upper(f));}for(int k=0;k<2;k++)pad=sf::add(pad,sf::mul(sf::mag(planes[i].n[k]),bound.displacement[k]));if(lo>pad){out.empty=true;snprintf(out.reason,sizeof(out.reason),"whole_chart_excluded_by_original_target_plane");return out;}if(hi+pad>=0)out.active.push_back(i);}if(out.active.empty()){out.whole=true;snprintf(out.reason,sizeof(out.reason),"whole_chart_inside_all_original_target_planes");return out;}
    auto f=frame(m,chart,sf::range(0,1));for(int dependent:{1,0}){bool valid=true;for(int i:out.active){const auto&p=planes[i];S coefficient;Q length=0;for(int k=0;k<2;k++){S direction=dependent?f.E[k]:f.D[k];coefficient=coefficient+p.n[k]*direction;length=std::max(length,sf::mag(direction));}Q radius=sf::mul(sf::add(sf::mag(p.n[0]),sf::mag(p.n[1])),sf::mul(bound.L,length));out.derivative[i]=S(coefficient.n,sf::add(coefficient.e,radius));valid=valid&&(sf::lower(out.derivative[i])>0||sf::upper(out.derivative[i])<0);}if(valid){out.dependent=dependent;snprintf(out.reason,sizeof(out.reason),"whole_chart_target_faces_transverse_in_selected_reference_axis");return out;}}snprintf(out.reason,sizeof(out.reason),"no_single_transverse_reference_axis_for_all_active_faces_event_chart_required");return out;}
static Axis certified_axis(const sd::Compiled&m,const sd::Chart&chart,const std::vector<Face>&faces,const std::vector<Plane>&planes,const sf::RealBounds&bound){auto out=axis(m,chart,planes,bound);if(out.empty||out.whole)return out;out.dependent=-1;auto shape=shape_bounds(m);std::vector<GradientCone>cone(planes.size());std::vector<int>active;for(int index:out.active){cone[index]=gradient_cone(faces[index],shape);S actual=actual_face_range(faces[index],chart,shape);if(sf::lower(actual)>0){out.empty=true;out.active.clear();snprintf(out.reason,sizeof(out.reason),"actual_restricted_cubic_and_analytic_residual_proves_whole_chart_excluded");return out;}if(sf::upper(actual)>=0)active.push_back(index);}out.active=active;if(active.empty()){out.whole=true;snprintf(out.reason,sizeof(out.reason),"actual_restricted_cubic_and_analytic_residual_proves_whole_chart_inside");return out;}auto f=frame(m,chart,sf::range(0,1));Q best=-1;std::vector<S>chosen;for(int dependent:{1,0}){bool valid=true;Q margin=1;std::vector<S>derivatives(planes.size());for(int index:active){S derivative;Q scale=0;for(int k=0;k<2;k++){S direction=dependent?f.E[k]:f.D[k];derivative=derivative+cone[index].gradient[k]*direction;scale=sf::add(scale,sf::mul(sf::mag(cone[index].gradient[k]),sf::mag(direction)));}derivatives[index]=derivative;Q separation=sf::lower(derivative)>0?sf::lower(derivative):-sf::upper(derivative);valid=valid&&separation>0;if(scale>0)margin=std::min(margin,separation/scale);}if(valid&&margin>best){best=margin;chosen=derivatives;out.dependent=dependent;}}if(out.dependent>=0){out.derivative=chosen;snprintf(out.reason,sizeof(out.reason),"actual_gradient_cone_maximum_normalized_transverse_margin");return out;}snprintf(out.reason,sizeof(out.reason),"actual_gradient_cone_needs_oblique_event_projection");return out;}
struct ScalarRoots {std::vector<S>root;bool complete=false;uint64_t windows=0;char reason[192]={};};
template<class Function>static ScalarRoots scalar_roots(Function evaluate){ScalarRoots out;try{struct W {Q lo,hi;int depth;bool exists;};std::vector<W>stack={{0,1,0,false}};while(!stack.empty()){auto w=stack.back();stack.pop_back();need(out.windows++<256,"fixed256 reference chart edge event budget");S parameter=sf::range(w.lo,w.hi);auto f=evaluate(parameter);if(sf::lower(f.v)>0||sf::upper(f.v)<0)continue;auto left=evaluate(S(w.lo)),right=evaluate(S(w.hi));bool transverse=sf::lower(f.g[0])>0||sf::upper(f.g[0])<0,opposite=sf::upper(left.v)<0&&sf::lower(right.v)>0||sf::lower(left.v)>0&&sf::upper(right.v)<0;if(transverse){if(!opposite&&((sf::lower(left.v)>0&&sf::lower(right.v)>0)||(sf::upper(left.v)<0&&sf::upper(right.v)<0)))continue;auto center=evaluate(S(parameter.n));S N=S(parameter.n)-center.v/f.g[0];Q lo=std::max(w.lo,sf::lower(N)),hi=std::min(w.hi,sf::upper(N));if(lo>hi)continue;if((opposite||w.exists)&&hi-lo<1e-27Q){out.root.push_back(sf::range(lo,hi));continue;}if(hi-lo<Q(.75L)*(w.hi-w.lo)){stack.push_back({lo,hi,w.depth+1,opposite||w.exists});continue;}}need(w.depth<64&&parameter.n>w.lo&&parameter.n<w.hi,"reference edge endpoint/tangency event not separated at fixed budget");stack.push_back({parameter.n,w.hi,w.depth+1,false});stack.push_back({w.lo,parameter.n,w.depth+1,false});}std::sort(out.root.begin(),out.root.end(),[](S a,S b){return a.n<b.n;});for(size_t i=1;i<out.root.size();i++)need(sf::upper(out.root[i-1])<sf::lower(out.root[i]),"distinct reference edge events");out.complete=true;snprintf(out.reason,sizeof(out.reason),"complete_same_function_reference_edge_roots");}catch(const std::exception&e){snprintf(out.reason,sizeof(out.reason),"%s",e.what());}return out;}
static sd::Jet in_axis(const Face&f,const sd::Chart&chart,int dependent,S free,S bound){return dependent==1?face(f,chart,free,bound):face(f,chart,bound,free);}
// A monotone implicit graph at one real free-coordinate box. Its retained
// box covers every root of that input family; no query uncertainty is
// discarded merely to meet a parameter-width stopping criterion.
static S graph(const Face&f,const sd::Chart&chart,const Axis&axis,int index,S free){
    S derivative=axis.derivative[index];if(axis.dependent==1)derivative=free*derivative;
    need(sf::lower(derivative)>0||sf::upper(derivative)<0,"implicit graph real derivative separated, including radial factor");
    auto left=in_axis(f,chart,axis.dependent,free,S(0)),right=in_axis(f,chart,axis.dependent,free,S(1));
    bool exists=sf::upper(left.v)<0&&sf::lower(right.v)>0||sf::lower(left.v)>0&&sf::upper(right.v)<0;
    need(exists,"implicit graph has a root for the entire real query box");Q lo=0,hi=1;
    for(int it=0;it<256;it++){
        S box=sf::range(lo,hi);auto center=in_axis(f,chart,axis.dependent,free,S(box.n));
        S N=S(box.n)-center.v/derivative;Q a=std::max(lo,sf::lower(N)),b=std::min(hi,sf::upper(N));
        need(a<=b,"real implicit graph inclusion");if(b-a<1e-27Q)return sf::range(a,b);
        if(b-a<Q(.75L)*(hi-lo)){lo=a;hi=b;continue;}
        // A wide whole-real derivative cone can make the first Newton
        // inclusion unchanged even for one exact free coordinate. Bisect
        // by certified monotone signs before treating it as a family floor.
        bool increasing=sf::lower(derivative)>0;
        if(sf::lower(center.v)>0){if(increasing)hi=box.n;else lo=box.n;continue;}
        if(sf::upper(center.v)<0){if(increasing)lo=box.n;else hi=box.n;continue;}
        return sf::range(a,b); // actual unresolved input-family enclosure
    }
    throw std::runtime_error("fixed256 implicit graph work budget");
}

struct Strip {int chart=0,dependent=1,lower=-1,upper=-1;S lo,hi;};
struct ClippedChart {sd::Chart reference;Q radial_cut=1,radial_derivative_lower=0;S apex_face;int trim_face=-1;S unresolved_free,candidate,selected;int unresolved_face=-1,selected_face=-1;Axis coordinate;std::vector<Strip>strip;bool complete=false;uint64_t edge_windows=0;char reason[192]={};};
static ClippedChart clipped_chart(const sd::Compiled&m,int chart_index,const std::vector<Face>&faces,const std::vector<Plane>&planes,const sf::RealBounds&bound,const Domain&domain){ClippedChart out;try{auto chart=sd::charts(m)[chart_index];out.reference=chart;out.coordinate=certified_axis(m,chart,faces,planes,bound);
        if(!out.coordinate.empty&&!out.coordinate.whole&&out.coordinate.dependent<0){
            // An outward face can prove that the actual clipped portion is
            // confined to a small original-apex neighborhood. This is a
            // geometric event chart, not accuracy-driven subdivision: the
            // discarded radial tail violates that one original target face.
            auto shape=shape_bounds(m);auto r=frame(m,chart,sf::range(0,1));Q cut=1;
            for(int index:out.coordinate.active){auto cone=gradient_cone(faces[index],shape);S derivative;for(int k=0;k<2;k++)derivative=derivative+cone.gradient[k]*r.D[k];S A=faces[index].vertex[chart.apex];A.e=sf::add(A.e,planes[index].input_offset_radius);
                if(sf::lower(derivative)>0&&sf::upper(A)<0){Q candidate=sf::quotient(-sf::lower(A),sf::lower(derivative));if(candidate<cut){cut=candidate;out.trim_face=index;out.apex_face=A;out.radial_derivative_lower=sf::lower(derivative);}}
            }
            if(cut<1){need(cut>0,"positive original-apex target event radius");out.radial_cut=cut;for(int a=1;a<3;a++)for(int j=0;j<3;j++)chart.vertex[a][j]=chart.vertex[0][j]+S(cut)*(chart.vertex[a][j]-chart.vertex[0][j]);chart.det=chart.det*S(cut)*S(cut);out.reference=chart;out.coordinate=certified_axis(m,chart,faces,planes,bound);}
        }
        const auto&a=out.coordinate;if(a.empty||a.whole){out.complete=true;snprintf(out.reason,sizeof(out.reason),"%s",a.reason);return out;}need(a.dependent>=0,a.reason);std::vector<S>events={S(0),S(1)};for(int index:a.active)for(int boundary=0;boundary<2;boundary++){auto roots=scalar_roots([&](S u){auto j=in_axis(faces[index],chart,a.dependent,u,S(boundary));sd::Jet f(j.v,j.g[1-a.dependent]);f.H[0][0]=j.H[1-a.dependent][1-a.dependent];return f;});out.edge_windows+=roots.windows;need(roots.complete,roots.reason);events.insert(events.end(),roots.root.begin(),roots.root.end());}
        for(const auto&node:domain.node)if(node.active.size()>=2&&node.vertex<0){auto inverse=inverse_reference(m,bound,node.p[0],node.p[1]);auto coordinate=coordinates(m,chart,inverse.X,inverse.Y);if(sf::upper(coordinate.s)<0||sf::lower(coordinate.s)>1||sf::upper(coordinate.t)<0||sf::lower(coordinate.t)>1)continue;S u=a.dependent==1?coordinate.s:coordinate.t;need(sf::lower(u)>0&&sf::upper(u)<1,"target vertex on reference chart boundary needs a shared event identity");events.push_back(u);}std::sort(events.begin(),events.end(),[](S a,S b){return a.n<b.n;});for(size_t i=0;i+1<events.size();i++){need(sf::upper(events[i])<sf::lower(events[i+1]),"reference edge/target vertex event ordering must be proved without coalescing");S free=(events[i]+events[i+1])/S(2),lo(0),hi(1);int lower=-1,upper=-1;bool empty=false;for(int index:a.active){auto left=in_axis(faces[index],chart,a.dependent,S(free.n),S(0)),right=in_axis(faces[index],chart,a.dependent,S(free.n),S(1));if(sf::upper(left.v)<0&&sf::upper(right.v)<0)continue;if(sf::lower(left.v)>0&&sf::lower(right.v)>0){empty=true;break;}S root=graph(faces[index],chart,a,index,S(free.n));out.unresolved_free=free;out.candidate=root;out.unresolved_face=index;out.selected=sf::lower(a.derivative[index])>0?hi:lo;out.selected_face=sf::lower(a.derivative[index])>0?upper:lower;if(sf::lower(a.derivative[index])>0){if(sf::upper(root)<sf::lower(hi)){hi=root;upper=index;}else need(sf::lower(root)>sf::upper(hi),"distinct selected upper implicit target graphs");}else{if(sf::lower(root)>sf::upper(lo)){lo=root;lower=index;}else need(sf::upper(root)<sf::lower(lo),"distinct selected lower implicit target graphs");}}if(empty||sf::lower(lo)>sf::upper(hi))continue;need(sf::upper(lo)<sf::lower(hi),"positive actual clipped interval in dependent coordinate");out.strip.push_back({chart_index,a.dependent,lower,upper,events[i],events[i+1]});}out.complete=true;snprintf(out.reason,sizeof(out.reason),"complete_actual_clipped_reference_graph_arrangement");
    }catch(const std::exception&e){snprintf(out.reason,sizeof(out.reason),"%s",e.what());}return out;}
struct Implicit {S value;Q contraction=0,invariance=0,forcing=0,correction_radius=0,prior_radius=0,real_root_radius=0;};
struct GraphFailure:std::runtime_error {Q contraction,invariance,forcing,radius,prior=0,real_root_radius=0;GraphFailure(const char*reason,Q q,Q inv,Q residual,Q R,Q prior_radius=0,Q real_radius=0):std::runtime_error(reason),contraction(q),invariance(inv),forcing(residual),radius(R),prior(prior_radius),real_root_radius(real_radius){}};
static Implicit graph_disk(const Face&f,const sd::Chart&chart,const Axis&a,int index,S free,S center){
    auto at=[&]{sd::TracePhase phase(41);return in_axis(f,chart,a.dependent,S(free.n),S(center.n));}();
    Q A=Q(1)/at.g[a.dependent].n,slope=(-at.g[1-a.dependent]/at.g[a.dependent]).n;
    S delta=free-S(free.n),predictor=S(center.n)+S(slope)*delta;
    auto constant_forcing=[&]{sd::TracePhase phase(42);return in_axis(f,chart,a.dependent,free,S(center.n));}();
    Q prior=sf::add(sf::mul(Q(1.5L),sf::mag(S(A)*constant_forcing.v)),center.e);
    need(prior>0,"positive geometric prior for implicit correction");
    auto enclosure=[&]{sd::TracePhase phase(43);return in_axis(f,chart,a.dependent,free,S(predictor.n,sf::add(predictor.e,prior)));}();
    int d=a.dependent,u=1-d;
    // Preserve the affine predictor's first-order cancellation BEFORE
    // magnitudes. The remaining forcing is a directional Hessian term.
    S first=at.g[u]+S(slope)*at.g[d];
    S second=enclosure.H[u][u]+S(2)*S(slope)*enclosure.H[u][d]+S(slope)*S(slope)*enclosure.H[d][d];
    Q D=sf::mag(delta),Aabs=fabsq(A);
    Q residual=sf::mul(Aabs,sf::add(sf::mag(at.v),sf::add(sf::mul(sf::mag(first),D),sf::mul(Q(.5L),sf::mul(sf::mag(second),sf::mul(D,D))))));
    S mixed=enclosure.H[d][u]+S(slope)*enclosure.H[d][d];
    Q q0=sf::mag(S(1)-S(A)*at.g[d]),g=sf::add(q0,sf::mul(Aabs,sf::mul(sf::mag(mixed),D))),h=sf::mul(Aabs,sf::mag(enclosure.H[d][d]));
    Q gap=mf::FP<Q>::down(Q(1)-g);
    if(!(gap>0))throw GraphFailure("implicit affine-predictor derivative cone unproved",g,0,residual,prior,prior,center.e);
    Q R=std::max(center.e,sf::quotient(sf::mul(Q(2),residual),gap));
    // The 2*residual/gap choice is only a sufficient construction. If it
    // exceeds the ALREADY certified Hessian neighborhood, test that prior
    // itself as an invariant contraction domain. Never extrapolate H.
    if(R>prior){Q qp=sf::add(g,sf::mul(h,prior)),ip=sf::quotient(sf::add(residual,sf::mul(qp,prior)),prior);if(qp<1&&ip<1)R=prior;}
    Q q=sf::add(g,sf::mul(h,R));
    Q invariant=R>0?sf::quotient(sf::add(residual,sf::mul(q,R)),R):__builtin_huge_valq();
    if(!(R>0&&R<=prior))throw GraphFailure("implicit correction radius exceeds certified Hessian neighborhood",q,invariant,residual,R,prior,center.e);
    if(!(q<1&&invariant<1))throw GraphFailure("implicit affine-predictor correction contraction/invariance unproved",q,invariant,residual,R,prior,center.e);
    return {S(predictor.n,sf::add(predictor.e,R)),q,invariant,residual,R,prior,center.e};
}
struct StripPoint {sd::Point value;S s,t,jacobian;Q contraction=0,invariance=0,forcing=0,correction_radius=0,prior_radius=0,real_root_radius=0;};
struct RegularGraph {sd::AffineAngular coordinate;Implicit certificate;};
// A finite correction neighborhood from distances to the actual
// canonical denominator and fixed binomial-domain boundaries. This is
// chosen before a Hessian evaluation, independently of quadrature error.
static Q binomial_clearance(const double*anchor,S X,S Y,S Dx,S Dy,double eps){
    S ax=sf::exact(anchor[0]),ay=sf::exact(anchor[1]),dx=X-ax,dy=Y-ay,z=ax*ax+ay*ay+sf::exact(eps)*sf::exact(eps),argument=(S(2)*(ax*dx+ay*dy)+dx*dx+dy*dy)/z;
    Q gap=mf::FP<Q>::down(Q(.75L)-sf::mag(argument));if(!(gap>0))throw sd::DomainFailure("geometric regular-q base binomial domain",argument,Q(.75L),X,Y,anchor);
    S first=S(2)*(X*Dx+Y*Dy)/z,second=(Dx*Dx+Dy*Dy)/z;Q derivative=sf::add(sf::mag(first),sf::mag(second));return derivative>0?std::min(Q(1),mf::FP<Q>::down(gap/sf::mul(Q(2),derivative))):Q(1);
}
static Q geometric_q_prior(const Face&face,const sd::Chart&chart,S free,const sd::AffineAngular&prediction){
    const auto&m=*face.map;const auto&t=m.base.state->triangle[m.base.source];auto frame0=frame(m,chart,S(0));S X[2];for(int k=0;k<2;k++)X[k]=sd::affine_coordinate(frame0.A[k],frame0.A[k]+frame0.D[k],frame0.A[k]+frame0.D[k]+frame0.E[k],free,prediction);
    Q radius=binomial_clearance(t.vertex[chart.apex],X[0],X[1],frame0.E[0],frame0.E[1],m.base.eps);
    for(int e=0;e<3;e++){
        int a=t.endpoint[e][0],b=t.endpoint[e][1];S D;S parameter=sd::regular_edge_parameter(chart,free,prediction,a,b,D),Ea=chart.vertex[2][a]-chart.vertex[1][a],Eb=chart.vertex[2][b]-chart.vertex[1][b],Ed=Ea+Eb;
        Q floor=sf::lower(D),change=sf::mag(Ed);if(change>0)radius=std::min(radius,mf::FP<Q>::down(floor/sf::mul(Q(2),change)));
        const auto&edge=m.base.edge[e];S ex=edge.A[0]+parameter*edge.Dref[0],ey=edge.A[1]+parameter*edge.Dref[1];const double*anchor=parameter.n>Q(.5L)?t.vertex[t.endpoint[e][1]]:edge.native.A;Q p_radius;
        try{p_radius=binomial_clearance(anchor,ex,ey,edge.Dref[0],edge.Dref[1],m.base.eps);}catch(sd::DomainFailure&failure){failure.edge=e;failure.edge_parameter=parameter;throw;}
        S determinant_weights[3];determinant_weights[a]=Eb;determinant_weights[b]=Eb-Ed;S determinant=sd::barycentric_combination(chart,free,prediction,determinant_weights);
        Q denominator=mf::FP<Q>::down(floor*floor),parameter_lipschitz=sf::quotient(sf::mul(Q(2),sf::mag(determinant)),denominator);
        if(parameter_lipschitz>0)radius=std::min(radius,mf::FP<Q>::down(p_radius/parameter_lipschitz));
    }
    need(radius>0,"positive geometric regular-q correction neighborhood");return radius;
}
static RegularGraph graph_q_disk(const Face&f,const sd::Chart&chart,const Axis&a,int index,S free,S center){
    need(a.dependent==1,"regular q graph is for an angular dependent boundary");S s0(free.n),q0=s0*S(center.n);
    sd::AffineAngular constant;constant.intercept=q0;
    auto at=[&]{sd::TracePhase phase(51);return face_q(f,chart,s0,constant);}();
    Q A=Q(1)/at.g[1].n,beta=(-at.g[0]/at.g[1]).n;
    RegularGraph out;out.coordinate.slope=S(beta);out.coordinate.intercept=q0-S(beta)*s0;
    S delta=free-s0;Q real_radius=sf::mag(s0*center-S(q0.n));
    Q prior=[&]{sd::TracePhase phase(55);return geometric_q_prior(f,chart,free,out.coordinate);}();
    need(prior>0,"positive geometric regular-q correction prior");auto enclosure_coordinate=out.coordinate;enclosure_coordinate.correction=S(0,prior);
    auto enclosure=[&]{sd::TracePhase phase(53);return face_q(f,chart,free,enclosure_coordinate);}();
    S first=at.g[0]+S(beta)*at.g[1],second=enclosure.H[0][0]+S(2)*S(beta)*enclosure.H[0][1]+S(beta)*S(beta)*enclosure.H[1][1];
    Q D=sf::mag(delta),Aabs=fabsq(A),residual=sf::mul(Aabs,sf::add(sf::mag(at.v),sf::add(sf::mul(sf::mag(first),D),sf::mul(Q(.5L),sf::mul(sf::mag(second),sf::mul(D,D))))));
    S mixed=enclosure.H[1][0]+S(beta)*enclosure.H[1][1];Q g=sf::add(sf::mag(S(1)-S(A)*at.g[1]),sf::mul(Aabs,sf::mul(sf::mag(mixed),D))),h=sf::mul(Aabs,sf::mag(enclosure.H[1][1])),gap=mf::FP<Q>::down(Q(1)-g);
    if(!(gap>0))throw GraphFailure("regular q implicit derivative gap unproved",g,0,residual,prior,prior,real_radius);
    Q R=std::max(real_radius,sf::quotient(sf::mul(Q(2),residual),gap)),q=sf::add(g,sf::mul(h,R)),invariance=R>0?sf::quotient(sf::add(residual,sf::mul(q,R)),R):__builtin_huge_valq();
    if(!(R>0&&R<=prior))throw GraphFailure("regular q correction outside certified Hessian neighborhood",q,invariance,residual,R,prior,real_radius);
    if(!(q<1&&invariance<1))throw GraphFailure("regular q correction contraction/invariance unproved",q,invariance,residual,R,prior,real_radius);
    out.coordinate.correction=S(0,R);S predicted=out.coordinate.slope*free+out.coordinate.intercept;
    out.certificate={S(predicted.n,sf::add(predicted.e,R)),q,invariance,residual,R,prior,real_radius};return out;
}
static StripPoint regular_angular_strip_point(const sd::Compiled&m,const sd::Chart&chart,const std::vector<Face>&faces,const Axis&a,const Strip&strip,S free,S eta){
    need(sf::lower(free)>0,"regular angular strip keeps its genuine s=0 pole excluded");RegularGraph lo,hi;hi.coordinate.slope=S(1);StripPoint out;
    auto branch=[&](int index){auto center=graph(faces[index],chart,a,index,S(free.n));auto result=graph_q_disk(faces[index],chart,a,index,free,center);const auto&b=result.certificate;out.contraction=std::max(out.contraction,b.contraction);out.invariance=std::max(out.invariance,b.invariance);out.forcing=std::max(out.forcing,b.forcing);out.correction_radius=std::max(out.correction_radius,b.correction_radius);out.prior_radius=std::max(out.prior_radius,b.prior_radius);out.real_root_radius=std::max(out.real_root_radius,b.real_root_radius);return result;};
    if(strip.lower>=0)lo=branch(strip.lower);if(strip.upper>=0)hi=branch(strip.upper);
    sd::AffineAngular mixed=lo.coordinate;mixed.interpolated=true;mixed.eta=eta;mixed.width_slope=hi.coordinate.slope-lo.coordinate.slope;mixed.width_intercept=hi.coordinate.intercept-lo.coordinate.intercept;mixed.width_correction=hi.coordinate.correction-lo.coordinate.correction;
    S width=free*(hi.coordinate.slope-lo.coordinate.slope)+(hi.coordinate.intercept-lo.coordinate.intercept)+(hi.coordinate.correction-lo.coordinate.correction);
    out.s=free;out.t=mixed.slope+(mixed.intercept+mixed.correction)/free+eta*width/free;
    out.value=[&]{sd::TracePhase phase(54);return sd::evaluate(m,chart,free,out.t,&mixed);}();
    // det(dX/d(s,q)) is chart.det.  s*(t_hi-t_lo)=q_hi-q_lo
    // is cancelled analytically before the Jacobian is bounded.
    out.jacobian=chart.det*width*(strip.hi-strip.lo);return out;
}
static StripPoint strip_point(const sd::Compiled&m,const sd::Chart&chart,const std::vector<Face>&faces,const Axis&a,const Strip&strip,S free,S eta,bool complex,bool regular_q=true){if(complex&&regular_q&&a.dependent==1&&(strip.lower>=0||strip.upper>=0))return regular_angular_strip_point(m,chart,faces,a,strip,free,eta);S lo(0),hi(1);Q qmax=0,imax=0,forcing=0,radius=0,prior=0,real_radius=0;auto extend=[&](const Face&face,int index,S center){auto bound=graph_disk(face,chart,a,index,free,center);qmax=std::max(qmax,bound.contraction);imax=std::max(imax,bound.invariance);forcing=std::max(forcing,bound.forcing);radius=std::max(radius,bound.correction_radius);prior=std::max(prior,bound.prior_radius);real_radius=std::max(real_radius,bound.real_root_radius);return bound.value;};if(strip.lower>=0){S center=graph(faces[strip.lower],chart,a,strip.lower,S(free.n));lo=complex?extend(faces[strip.lower],strip.lower,center):graph(faces[strip.lower],chart,a,strip.lower,free);}if(strip.upper>=0){S center=graph(faces[strip.upper],chart,a,strip.upper,S(free.n));hi=complex?extend(faces[strip.upper],strip.upper,center):graph(faces[strip.upper],chart,a,strip.upper,free);}S width=hi-lo;if(!complex)need(sf::lower(width)>0,"positive actual clipped-dependent-coordinate width");S dependent=lo+eta*width;StripPoint out;out.contraction=qmax;out.invariance=imax;out.forcing=forcing;out.correction_radius=radius;out.prior_radius=prior;out.real_root_radius=real_radius;out.s=a.dependent==1?free:dependent;out.t=a.dependent==1?dependent:free;out.value=[&]{sd::TracePhase phase(44);return sd::evaluate(m,chart,out.s,out.t);}();out.jacobian=out.value.reference_jacobian*width*(strip.hi-strip.lo);if(!complex)need(sf::lower(out.s)>0&&sf::lower(out.t)>0&&sf::upper(out.s)<1&&sf::upper(out.t)<1&&sf::lower(out.jacobian)>0,"positive interior reference strip point with tracked orientation");return out;}
static Certificate strip_certificate(const curved::Packet&p,const sd::Compiled&m,const std::vector<Face>&faces,const Axis&a,const Strip&strip,const Thermal&context){Certificate out;try{auto chart=sd::charts(m)[strip.chart];S width=strip.hi-strip.lo,center_free=(strip.lo+strip.hi)/S(2);auto center=strip_point(m,chart,faces,a,strip,center_free,S(Q(.5L)),false);auto zc=reference_moments(p,center.value,reference_radius(m,chart,center.s,center.t,center.value),context);Q constant[12];for(int k=0;k<12;k++)constant[k]=(zc[k]*center.jacobian).n;Q maxima[12]={};for(int axis=0;axis<2;axis++){S u=axis?sf::range(0,1):S(Q(.5L),Q(1.25L)),eta=axis?S(Q(.5L),Q(1.25L)):sf::range(0,1),free=strip.lo+width*u;auto q=strip_point(m,chart,faces,a,strip,free,eta,true);auto z=reference_moments(p,q.value,reference_radius(m,chart,q.s,q.t,q.value),context);for(int k=0;k<12;k++){S g=z[k]*q.jacobian;maxima[k]=std::max(maxima[k],sf::add(fabsq(g.n-constant[k]),g.e));}}Q C=sf::mag(cauchy_coefficient());for(int k=0;k<12;k++)out.remainder[k]=sf::mul(C,maxima[k]);out.pass=true;snprintf(out.reason,sizeof(out.reason),"actual_clipped_reference_strip_implicit_graph_all12_remainder_pass");}catch(const std::exception&e){snprintf(out.reason,sizeof(out.reason),"%s",e.what());}return out;}
struct ClipPanel {Strip strip;S eta_lo,eta_hi;bool whole=false;Panel ordinary;sd::Chart reference;};
struct ClipAtlas {std::vector<ClipPanel>panel;bool complete=false;Q reference_radius_lower=0;char reason[192]={};};
static Q actual_reference_radius_lower(const sd::Compiled&m,const sf::RealBounds&bound,const Domain&domain,const std::map<int,Edge>&edges){Q lo[2]={__builtin_huge_valq(),__builtin_huge_valq()},hi[2]={-__builtin_huge_valq(),-__builtin_huge_valq()};for(const auto&boundary:domain.boundary)for(int part=0;part<16;part++){auto curve=boundary_curve(domain,edges,boundary,sf::range(Q(part)/16,Q(part+1)/16));for(int k=0;k<2;k++){lo[k]=std::min(lo[k],sf::lower(curve.position[k]));hi[k]=std::max(hi[k],sf::upper(curve.position[k]));}}Q r2=0;for(int k=0;k<2;k++){lo[k]=mf::FP<Q>::down(lo[k]-bound.displacement[k]);hi[k]=sf::add(hi[k],bound.displacement[k]);Q d=lo[k]>0?lo[k]:hi[k]<0?-hi[k]:Q(0);r2=mf::FP<Q>::down(r2+mf::FP<Q>::down(d*d));}need(r2>0,"actual clipped reference region density-origin separation unproved");return sf::lower(mf::root(S(r2)));}
static ClipAtlas clipped_atlas(const curved::Packet&p,const sd::Compiled&m,const sf::RealBounds&bound,const Domain&domain,const std::map<int,Edge>&edges,const std::vector<Plane>&planes,const std::vector<ClippedChart>&charts,uint64_t query_budget){
    ClipAtlas out;try{
        out.reference_radius_lower=actual_reference_radius_lower(m,bound,domain,edges);
        Q scale=std::min(out.reference_radius_lower,Q(p.parameters[3]))/Q(8);
        std::array<sd::Chart,6>reference;for(int ci=0;ci<6;ci++)reference[ci]=charts[ci].reference;auto shape=shape_bounds(m);
        for(int ci=0;ci<6;ci++){
            const auto&chart=charts[ci];need(chart.complete,chart.reason);if(chart.coordinate.empty)continue;
            auto f=frame(m,reference[ci],sf::range(0,1));
            if(chart.coordinate.whole){
                int nr=std::max(1,int(ceilq(std::max(sf::mag(f.D[0]),sf::mag(f.D[1]))/scale))),nt=std::max(4,int(ceilq(std::max(sf::mag(f.E[0]),sf::mag(f.E[1]))/scale)));
                for(int i=0;i<nr;i++)for(int j=0;j<nt;j++){
                    need(out.panel.size()<query_budget/256,"fixed combined ordinary/clipped panel budget");
                    ClipPanel panel;panel.reference=reference[ci];panel.whole=true;panel.strip.chart=ci;panel.ordinary={ci,i,j,nr,nt,S(i)/S(nr),S(i+1)/S(nr),S(j)/S(nt),S(j+1)/S(nt)};out.panel.push_back(panel);
                }continue;
            }
            // Only actual selected boundary branches determine the graph
            // metric of THIS event slab. Other faces remain in ownership
            // checks and cannot inflate a metric for an absent branch.
            for(const auto&strip:chart.strip){
                auto local=frame(m,reference[ci],strip.dependent?sf::range(0,1):sf::range(sf::lower(strip.lo),sf::upper(strip.hi)));
                Q D=std::max(sf::mag(local.D[0]),sf::mag(local.D[1])),E=std::max(sf::mag(local.E[0]),sf::mag(local.E[1]));
                Q base_free=strip.dependent?D:E,free_bound=base_free,dependent_bound=strip.dependent?E:D;
                for(int index:{strip.lower,strip.upper})if(index>=0){
                    auto cone=gradient_cone(pullback(m,planes[index]),shape);S numerator,denominator;
                    for(int k=0;k<2;k++){numerator=numerator+cone.gradient[k]*(strip.dependent?local.D[k]:local.E[k]);denominator=denominator+cone.gradient[k]*(strip.dependent?local.E[k]:local.D[k]);}
                    Q minimum=sf::lower(denominator)>0?sf::lower(denominator):-sf::upper(denominator);
                    need(minimum>0,"actual selected boundary derivative separates on this event slab");
                    Q slope=sf::quotient(sf::mag(numerator),minimum);
                    free_bound=std::max(free_bound,sf::add(base_free,sf::mul(dependent_bound,slope)));
                }
                need(free_bound>0&&dependent_bound>0,"positive actual selected-branch chart metric");
                Q max_step=scale/free_bound;int neta=std::max(1,int(ceilq(dependent_bound/scale)));S left=strip.lo;
                bool pole=strip.dependent==1&&(strip.lower>=0||strip.upper>=0);
                while(left.n<strip.hi.n){Q step=max_step;if(pole){need(sf::lower(left)>0,"implicit angular graph starts at an unresolved apex event");step=std::min(step,sf::lower(left)/Q(3));}
                    Q next=mf::FP<Q>::down(left.n+step);bool last=next>=sf::lower(strip.hi);S right=last?strip.hi:S(next);
                    need(sf::lower(right)>sf::upper(left),"geometric event-distance atlas has separated panel boundaries");
                    for(int j=0;j<neta;j++){need(out.panel.size()<query_budget/256,"deterministic clipped atlas exceeds fixed positive-query budget before integration");Strip part=strip;part.lo=left;part.hi=right;ClipPanel panel;panel.strip=part;panel.eta_lo=S(j)/S(neta);panel.eta_hi=S(j+1)/S(neta);panel.reference=reference[ci];out.panel.push_back(panel);}
                    if(last)break;left=right;
                }
            }
        }
        out.complete=true;snprintf(out.reason,sizeof(out.reason),"finite_actual_selected_branch_atlas_from_metric_density_and_implicit_apex_distances");
    }catch(const std::exception&e){snprintf(out.reason,sizeof(out.reason),"%s",e.what());}return out;
}
static Certificate clip_panel_certificate(const curved::Packet&p,const sd::Compiled&m,const std::vector<Face>&faces,const Axis&a,const ClipPanel&panel,const Thermal&context,bool regular_q=true){if(panel.whole)return panel_certificate_on_chart(p,m,panel.reference,panel.ordinary,context);Certificate out;try{auto chart=panel.reference;S du=panel.strip.hi-panel.strip.lo,de=panel.eta_hi-panel.eta_lo,um=(panel.strip.lo+panel.strip.hi)/S(2),em=(panel.eta_lo+panel.eta_hi)/S(2);auto center=strip_point(m,chart,faces,a,panel.strip,um,em,false);auto zc=reference_moments(p,center.value,reference_radius(m,chart,center.s,center.t,center.value),context);Q constants[12];for(int k=0;k<12;k++)constants[k]=(zc[k]*center.jacobian*de).n;Q maxima[12]={};for(int axis=0;axis<2;axis++){S u=axis?sf::range(0,1):S(Q(.5L),Q(1.25L)),e=axis?S(Q(.5L),Q(1.25L)):sf::range(0,1),free=panel.strip.lo+du*u,eta=panel.eta_lo+de*e;auto q=strip_point(m,chart,faces,a,panel.strip,free,eta,true,regular_q);out.implicit_observed=panel.strip.lower>=0||panel.strip.upper>=0;out.contraction=std::max(out.contraction,q.contraction);out.invariance=std::max(out.invariance,q.invariance);out.forcing=std::max(out.forcing,q.forcing);out.correction_radius=std::max(out.correction_radius,q.correction_radius);out.prior_radius=std::max(out.prior_radius,q.prior_radius);out.real_root_radius=std::max(out.real_root_radius,q.real_root_radius);auto z=reference_moments(p,q.value,reference_radius(m,chart,q.s,q.t,q.value),context);for(int k=0;k<12;k++){S value=z[k]*q.jacobian*de;maxima[k]=std::max(maxima[k],sf::add(fabsq(value.n-constants[k]),value.e));}}Q C=sf::mag(cauchy_coefficient());for(int k=0;k<12;k++)out.remainder[k]=sf::mul(C,maxima[k]);out.pass=true;snprintf(out.reason,sizeof(out.reason),"actual_clipped_panel_positive_mapping_and_analytic_all12_remainder");}catch(const GraphFailure&e){out.implicit_observed=true;out.contraction=e.contraction;out.invariance=e.invariance;out.forcing=e.forcing;out.correction_radius=e.radius;out.prior_radius=e.prior;out.real_root_radius=e.real_root_radius;snprintf(out.reason,sizeof(out.reason),"%s",e.what());}catch(const sd::DomainFailure&e){out.argument_observed=true;out.argument=e.argument;out.argument_limit=e.limit;out.argument_phase=e.phase;out.argument_edge=e.edge;out.coordinate_observed=e.coordinate_observed;out.argument_coordinate[0]=e.coordinate[0];out.argument_coordinate[1]=e.coordinate[1];out.argument_anchor[0]=e.anchor[0];out.argument_anchor[1]=e.anchor[1];out.argument_parameter=e.edge_parameter;snprintf(out.reason,sizeof(out.reason),"%s",e.what());}catch(const std::exception&e){snprintf(out.reason,sizeof(out.reason),"%s",e.what());}return out;}
}
#endif
