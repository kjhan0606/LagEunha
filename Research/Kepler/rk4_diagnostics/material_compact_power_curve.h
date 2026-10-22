#ifndef MATERIAL_COMPACT_POWER_CURVE_H
#define MATERIAL_COMPACT_POWER_CURVE_H
// Canonical power-plane / SAME-source-curve intersection primitive. Include
// after source geometry/target registry. No polygon proposal can certify Ω*.
namespace compact_power {
namespace sf=shared_function;
using Q=sf::Q;using S=sf::S;using D=sf::D;
using R=long double;
inline void require(bool b,const char*s){compact_p3::require(b,s);}
struct Plane {
    uint64_t cell[2]={};
    S origin[2],normal[2],offset;
    // For canonical low-cell inside: normal dot (x-origin) - offset <=0.
    // normal=2(q_high-q_low); offset=|q_high-q_low|²+W_low-W_high.
    // q enclosures retain actual-source transport error, not an IC reset.
};
inline Plane plane(const compact_source::Site&a,R wa,const compact_source::Site&b,R wb){
    require(a.cell!=b.cell&&std::isfinite(wa)&&std::isfinite(wb),"two distinct finite weighted source sites");
    const auto&lo=a.cell<b.cell?a:b;const auto&hi=a.cell<b.cell?b:a;
    Plane p;p.cell[0]=lo.cell;p.cell[1]=hi.cell;p.offset=S(Q(a.cell<b.cell?wa:wb))-S(Q(a.cell<b.cell?wb:wa));
    S separation;
    for(int k=0;k<2;k++){
        require(std::isfinite(lo.position[k])&&std::isfinite(hi.position[k])&&std::isfinite(lo.position_error[k])&&std::isfinite(hi.position_error[k])&&lo.position_error[k]>=0&&hi.position_error[k]>=0,"finite canonical generator enclosures");
        p.origin[k]=compact_source::canonical_position(lo,k);
        S delta=compact_source::canonical_position(hi,k)-p.origin[k];p.normal[k]=S(2)*delta;separation=separation+delta*delta;
    }
    require(sf::lower(separation)>0,"distinct canonical transported generators, no tolerance merge");p.offset=p.offset+separation;return p;
}
inline D value(const Plane&p,const std::array<D,2>&x){
    D out(-p.offset);for(int k=0;k<2;k++)out=out+D(p.normal[k])*(x[k]-D(p.origin[k]));return out;
}
struct Curve {
    uint64_t edge=0,source_hash=0,source_epoch=0;
    uint64_t vertex[2]={UINT64_MAX,UINT64_MAX};
    sf::Edge scalar={};STCoefficient coefficient[6]={};
};
inline Curve curve(const compact_source::Exterior&e){
    Curve c;c.edge=e.edge;c.source_hash=e.frozen.source_hash;c.source_epoch=e.frozen.source_epoch;c.scalar.native=e.frozen.function;
    c.vertex[0]=e.vertex[0];c.vertex[1]=e.vertex[1];
    require(e.edge==e.frozen.edge&&c.scalar.native.rank>=0&&c.scalar.native.rank<=2,"same full canonical exterior scalar rank");
    for(int k=0;k<2;k++){c.scalar.A[k]=sf::exact(e.reference[0][k]);c.scalar.Dref[k]=sf::exact(e.reference[1][k])-c.scalar.A[k];}
    if(c.scalar.native.rank)for(int j=0;j<4;j++){
        S x=c.scalar.A[0]+S(j)*c.scalar.Dref[0]/S(3),y=c.scalar.A[1]+S(j)*c.scalar.Dref[1]/S(3);
        auto v=sf::velocity(D(x),D(y),c.scalar.native.eps);for(int k=0;k<2;k++)c.scalar.node[j][k]=v[k].v;
    }
    for(int j=0;j<6;j++)c.coefficient[j]=e.frozen.coefficient[j];return c;
}
inline Curve source_edge_curve(const st::Factory&f,size_t source_triangle,int local_edge,
        uint64_t original_edge,uint64_t source_hash,uint64_t source_epoch){
    require(source_triangle<f.triangle.size()&&local_edge>=0&&local_edge<3&&
        original_edge!=UINT64_MAX&&source_hash&&source_epoch,
        "known original source edge and accepted coefficient epoch");
    const auto&t=f.triangle[source_triangle];const auto&e=t.edge[local_edge];
    require(e.rank>=1&&e.rank<=2,"same full original shared source edge rank");
    Curve c;c.edge=original_edge;c.source_hash=source_hash;c.source_epoch=source_epoch;
    c.scalar=sf::edge(e,t,local_edge);
    for(int j=0;j<4;j++){
        require(e.poly[j]>=0&&size_t(e.poly[j])<f.coefficient.size(),
            "complete current shared polynomial source edge coefficient");
        c.coefficient[j]=f.coefficient[e.poly[j]];
    }
    for(int j=0;j<e.rank;j++){
        require(e.id[j]>=0&&size_t(e.id[j])<f.coefficient.size(),
            "complete current shared enriched source edge coefficient");
        c.coefficient[4+j]=f.coefficient[e.id[j]];
    }
    return c;
}
inline std::array<D,2> position(const Curve&c,S u){
    D parameter(u,S(1),S(0));std::array<D,2>out,raw;
    if(c.scalar.native.rank)raw=sf::edge_raw(c.scalar,parameter);
    for(int k=0;k<2;k++){
        S origin=sf::exact(c.coefficient[0].x[k]),delta=sf::exact(c.coefficient[3].x[k])-origin;
        S bend1=(sf::exact(c.coefficient[1].x[k])-origin)-delta/S(3);
        S bend2=(sf::exact(c.coefficient[2].x[k])-origin)-S(2)*delta/S(3);
        D x=D(origin)+parameter*D(delta)+sf::lagrange(1,parameter)*D(bend1)+sf::lagrange(2,parameter)*D(bend2);
        for(int j=0;j<c.scalar.native.rank;j++){
            D N;for(int d=0;d<2;d++)N=N+D(sf::exact(c.scalar.native.encode[j][d]))*raw[d];
            x=x+D(sf::exact(c.coefficient[4+j].x[k]))*N;
        }
        // Exact canonical vertex chart, not a numerical curve endpoint clamp.
        // Lagrange values and scalar residuals are algebraically Kronecker/0.
        // Keep the analytically evaluated derivative and all nonvertex balls.
        if(u.e==0&&(u.n==0||u.n==1))x.v=sf::exact(c.coefficient[u.n==0?0:3].x[k]);
        out[k]=x;
    }return out;
}
// Second-order interval jet of the same canonical curve. Its centered
// enclosure retains parameter correlation on nearly coincident source and
// target edges; it is intersected with the original direct enclosure.
struct J2 {S v,d,dd;J2(S x=S(),S dx=S(),S ddx=S()):v(x),d(dx),dd(ddx){}};
inline J2 operator+(J2 a,J2 b){return {a.v+b.v,a.d+b.d,a.dd+b.dd};}
inline J2 operator-(J2 a){return {-a.v,-a.d,-a.dd};}
inline J2 operator-(J2 a,J2 b){return a+(-b);}
inline J2 operator*(J2 a,J2 b){return {
    a.v*b.v,a.d*b.v+a.v*b.d,
    a.dd*b.v+S(2)*a.d*b.d+a.v*b.dd};}
inline J2 operator/(J2 a,S b){return {a.v/b,a.d/b,a.dd/b};}
inline J2 second_power(J2 a,Q exponent){
    S value=sf::mf::power(a.v,exponent);
    S first=S(exponent)*value/a.v;
    S second=S(exponent*(exponent-Q(1)))*value/(a.v*a.v);
    return {value,first*a.d,first*a.dd+second*a.d*a.d};
}
inline J2 second_lagrange(int j,J2 u){
    J2 out(S(1));for(int k=0;k<4;k++)if(j!=k)
        out=out*((J2(S(3))*u-J2(S(k)))/S(j-k));
    return out;
}
inline std::array<J2,2> position_second(const Curve&c,S parameter){
    J2 u(parameter,S(1)),raw[2];
    if(c.scalar.native.rank){
        J2 X=J2(c.scalar.A[0])+u*J2(c.scalar.Dref[0]);
        J2 Y=J2(c.scalar.A[1])+u*J2(c.scalar.Dref[1]);
        S eps=sf::exact(c.scalar.native.eps);
        J2 w=second_power(X*X+Y*Y+J2(eps*eps),Q(-.75L));
        raw[0]=-Y*w;raw[1]=X*w;
        for(int j=0;j<4;j++){
            J2 L=second_lagrange(j,u);
            for(int d=0;d<2;d++)raw[d]=raw[d]-L*J2(c.scalar.node[j][d]);
        }
    }
    std::array<J2,2>out;
    for(int k=0;k<2;k++){
        S origin=sf::exact(c.coefficient[0].x[k]);
        S delta=sf::exact(c.coefficient[3].x[k])-origin;
        S bend1=(sf::exact(c.coefficient[1].x[k])-origin)-delta/S(3);
        S bend2=(sf::exact(c.coefficient[2].x[k])-origin)-S(2)*delta/S(3);
        J2 x=J2(origin)+u*J2(delta)+second_lagrange(1,u)*J2(bend1)+
            second_lagrange(2,u)*J2(bend2);
        for(int j=0;j<c.scalar.native.rank;j++){
            J2 N;
            for(int d=0;d<2;d++)
                N=N+J2(sf::exact(c.scalar.native.encode[j][d]))*raw[d];
            x=x+J2(sf::exact(c.coefficient[4+j].x[k]))*N;
        }
        if(parameter.e==0&&(parameter.n==0||parameter.n==1))
            x.v=sf::exact(c.coefficient[parameter.n==0?0:3].x[k]);
        out[k]=x;
    }
    return out;
}
inline std::array<D,2> position_centered(const Curve&c,S parameter){
    auto direct=position(c,parameter);
    if(parameter.e==0)return direct;
    S middle(parameter.n),h=parameter-middle;
    auto center=position_second(c,middle),whole=position_second(c,parameter);
    std::array<D,2>out;
    for(int k=0;k<2;k++){
        S centered=center[k].v+center[k].d*h+whole[k].dd*h*h/S(2);
        S derivative=center[k].d+whole[k].dd*h;
        Q lo=std::max(sf::lower(centered),sf::lower(direct[k].v));
        Q hi=std::min(sf::upper(centered),sf::upper(direct[k].v));
        require(lo<=hi,"centered and direct canonical source-edge positions intersect");
        out[k].v=sf::range(lo,hi);
        lo=std::max(sf::lower(derivative),sf::lower(direct[k].x));
        hi=std::min(sf::upper(derivative),sf::upper(direct[k].x));
        require(lo<=hi,"centered and direct canonical source-edge derivatives intersect");
        out[k].x=sf::range(lo,hi);
    }
    return out;
}
inline D value_centered(const Plane&p,const Curve&c,S parameter){
    D direct=value(p,position(c,parameter));
    if(parameter.e==0)return direct;
    S middle(parameter.n),h=parameter-middle;
    D at_middle=value(p,position(c,middle));
    auto whole=position_second(c,parameter);
    S second;
    for(int k=0;k<2;k++)second=second+p.normal[k]*whole[k].dd;
    S centered=at_middle.v+at_middle.x*h+second*h*h/S(2);
    S derivative=at_middle.x+second*h;
    Q lo=std::max(sf::lower(centered),sf::lower(direct.v));
    Q hi=std::min(sf::upper(centered),sf::upper(direct.v));
    require(lo<=hi,"centered and direct canonical source-edge signed distances intersect");
    D out;out.v=sf::range(lo,hi);
    lo=std::max(sf::lower(derivative),sf::lower(direct.x));
    hi=std::min(sf::upper(derivative),sf::upper(direct.x));
    require(lo<=hi,"centered and direct canonical source-edge signed slopes intersect");
    out.x=sf::range(lo,hi);return out;
}
inline int strict_sign(S x){return sf::lower(x)>0?1:(sf::upper(x)<0?-1:0);}
struct RootBudget {
    size_t maximum_boxes=0;unsigned maximum_depth=0,maximum_newton=0;
    Q parameter_width=0; // MUST be assigned from caller's existing error ledger.
};
struct Event {
    uint64_t cell[2]={},edge=0;size_t ordinal=0;
    Q lower=0,upper=0;int derivative_sign=0;
};
enum class Unresolved { none,box_budget,depth_budget,endpoint_or_coincidence,root_width,overlapping_events };
struct Roots {
    std::vector<Event>event;size_t boxes=0,newton=0;
    bool complete=false;Unresolved unresolved=Unresolved::none;
    Q witness_lower=0,witness_upper=1;
};
inline Roots intersections(const Plane&p,const Curve&c,const RootBudget&budget,
        bool exact_source_endpoint_0=false,bool exact_source_endpoint_1=false){
    require(budget.maximum_boxes>0&&budget.maximum_boxes<=65536&&budget.maximum_depth>0&&budget.maximum_depth<=64&&budget.maximum_newton>0&&budget.maximum_newton<=64&&finiteq(budget.parameter_width)&&budget.parameter_width>0&&budget.parameter_width<1,"explicit bounded canonical event budget, no hidden geometry tolerance");
    struct Box {Q a,b;unsigned depth;};std::vector<Box>stack={{0,1,0}};Roots out;
    auto fail=[&](Unresolved reason,Q a,Q b){out.unresolved=reason;out.witness_lower=a;out.witness_upper=b;return out;};
    auto at=[&](S u){auto f=value_centered(p,c,u);require(finiteq(f.v.n)&&finiteq(f.v.e)&&finiteq(f.x.n)&&finiteq(f.x.e),"finite canonical power/curve enclosure");return f;};
    while(!stack.empty()){
        Box box=stack.back();stack.pop_back();if(out.boxes==budget.maximum_boxes)return fail(Unresolved::box_budget,box.a,box.b);out.boxes++;
        auto f=at(sf::range(box.a,box.b));
        if(strict_sign(f.v)){
            require(!(box.a==0&&exact_source_endpoint_0)&&
                !(box.b==1&&exact_source_endpoint_1),
                "certified source endpoint must remain enclosed by its target line");
            continue;
        }
        int direction=strict_sign(f.x);
        if(direction){
            if((box.a==0&&exact_source_endpoint_0)||(box.b==1&&exact_source_endpoint_1)){
                require(!(box.a==0&&exact_source_endpoint_0&&
                    box.b==1&&exact_source_endpoint_1),
                    "one strictly monotone source edge/target line cannot contain both endpoint roots");
                Event e;e.cell[0]=p.cell[0];e.cell[1]=p.cell[1];e.edge=c.edge;
                e.lower=e.upper=(box.a==0&&exact_source_endpoint_0)?Q(0):Q(1);
                e.derivative_sign=direction;out.event.push_back(e);continue;
            }
            auto left=at(S(box.a)).v,right=at(S(box.b)).v;int ls=strict_sign(left),rs=strict_sign(right);
            if(ls&&rs){
                if(ls==rs)continue;
                require(ls==-direction&&rs==direction,"monotone canonical bracket orientation");
                Q a=box.a,b=box.b;bool narrow=false;
                // Strict endpoint bracket + strict derivative certify ONE
                // root. Interval Newton contracts its enclosure; it never
                // invents existence from a sampled near-zero value.
                for(unsigned it=0;it<budget.maximum_newton;it++){
                    if(sf::upper(S(b)-S(a))<=budget.parameter_width){narrow=true;break;}
                    Q mid=(a+b)/2;auto center=at(S(mid));auto whole=at(sf::range(a,b));
                    require(strict_sign(whole.x)==direction,"retained monotonicity on root enclosure");
                    S next=S(mid)-center.v/whole.x;
                    Q na=std::max(a,sf::lower(next)),nb=std::min(b,sf::upper(next));out.newton++;
                    require(na<=nb,"nonempty interval Newton enclosure of bracketed canonical root");
                    if(nb-na<(b-a)*Q(.75L)){a=na;b=nb;continue;}
                    int sign=strict_sign(center.v);if(!sign)break;
                    if(sign==ls)a=mid;else b=mid;
                }
                narrow|=sf::upper(S(b)-S(a))<=budget.parameter_width;
                if(!narrow)return fail(Unresolved::root_width,a,b);
                Event e;e.cell[0]=p.cell[0];e.cell[1]=p.cell[1];e.edge=c.edge;e.lower=a;e.upper=b;e.derivative_sign=direction;out.event.push_back(e);continue;
            }
            // A true endpoint event, a near endpoint or uncertain q may be
            // indistinguishable here. No endpoint snap/tolerance tie. A later
            // algebraic shared-event adapter must resolve it or reject.
            return fail(Unresolved::endpoint_or_coincidence,box.a,box.b);
        }
        if(box.depth==budget.maximum_depth)return fail(Unresolved::depth_budget,box.a,box.b);
        Q mid=(box.a+box.b)/2;stack.push_back({mid,box.b,box.depth+1});stack.push_back({box.a,mid,box.depth+1});
    }
    std::sort(out.event.begin(),out.event.end(),[](const Event&a,const Event&b){return a.lower<b.lower;});
    for(size_t j=0;j<out.event.size();j++){
        out.event[j].ordinal=j;
        if(j&&out.event[j-1].upper>=out.event[j].lower)return fail(Unresolved::overlapping_events,out.event[j-1].lower,out.event[j].upper);
    }
    out.complete=true;return out;
}
} // namespace compact_power
#endif
