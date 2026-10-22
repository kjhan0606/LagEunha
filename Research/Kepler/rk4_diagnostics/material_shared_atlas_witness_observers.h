#ifndef MATERIAL_SHARED_ATLAS_WITNESS_OBSERVERS_H
#define MATERIAL_SHARED_ATLAS_WITNESS_OBSERVERS_H
// Read-only hooks in PRIVATE copied headers. Every original operation and
// rejection condition remains in place. No hook changes a field or radius.
namespace atlas_witness {
namespace sd=shared_duffy;namespace sf=shared_function;using Q=sf::Q;using S=sf::S;
struct AxisTrial {uint64_t slot,source,target,chart,trial,dependent,valid,active_mask;sd::Chart reference;Q margin;S derivative[23];};
struct Event {uint64_t slot,source,target,chart,face,boundary,complete,windows,index,count;S parameter;};
struct InverseIteration {uint64_t slot,source,target,chart,iteration;S input[2],lambda[3];Q center[2],residual_nominal[2],correction[2],step;};
struct InverseFinal {uint64_t slot,source,target,chart,inside;S input[2],output[2],position[2],F[4],J,lambda[3];Q anchor[2],center[2],L,radius,residual;};
struct RadiusPiece {uint64_t slot,source,target,index,part;int32_t first,last,canonical_edge,reserved;S lo,hi,position[2];};
struct RadiusFinal {uint64_t slot,source,target;Q reference_lo[2],reference_hi[2],r2,displacement[2];};
struct Capture {
    uint64_t slot=0,source=0,target=0,chart=0;
    std::vector<AxisTrial>axis;std::vector<Event>event;
    std::vector<InverseIteration>iteration;std::vector<InverseFinal>inverse;
    std::vector<RadiusPiece>radius_piece;std::vector<RadiusFinal>radius;
};
static thread_local Capture* current=nullptr;
struct Scope {Capture* previous;explicit Scope(Capture&capture):previous(current){current=&capture;}~Scope(){current=previous;}};
static void axis_trial(const sd::Compiled&m,const sd::Chart&chart,int dependent,const std::vector<int>&active,const std::vector<S>&derivatives,bool valid,Q margin){
    if(!current)return;AxisTrial row={};row.slot=current->slot;row.source=m.base.source;row.target=current->target;row.chart=current->chart;row.trial=current->axis.size();row.dependent=dependent;row.valid=valid;row.reference=chart;row.margin=margin;
    for(int face:active){need(face>=0&&face<23&&size_t(face)<derivatives.size(),"witness fixed plane identity");row.active_mask|=uint64_t(1)<<face;row.derivative[face]=derivatives[face];}current->axis.push_back(row);
}
template<class Roots>static void edge_events(int face,int boundary,const Roots&r){
    if(!current)return;uint64_t count=r.root.size();
    for(uint64_t i=0;i<std::max(uint64_t(1),count);i++){Event row={};row.slot=current->slot;row.source=current->source;row.target=current->target;row.chart=current->chart;row.face=face;row.boundary=boundary;row.complete=r.complete;row.windows=r.windows;row.index=count?i:UINT64_MAX;row.count=count;if(count)row.parameter=r.root[i];current->event.push_back(row);}
}
template<class Lambda>static void inverse_iteration(const sd::Compiled&m,S x,S y,int iteration,Q X,Q Y,Q fx,Q fy,Q dx,Q dy,Q step,const Lambda&lambda){
    if(!current)return;InverseIteration row={};row.slot=current->slot;row.source=m.base.source;row.target=current->target;row.chart=current->chart;row.iteration=iteration;row.input[0]=x;row.input[1]=y;row.center[0]=X;row.center[1]=Y;row.residual_nominal[0]=fx;row.residual_nominal[1]=fy;row.correction[0]=dx;row.correction[1]=dy;row.step=step;for(int i=0;i<3;i++)row.lambda[i]=lambda[i].v;current->iteration.push_back(row);
}
template<class Point,class Lambda>static void inverse_final(const sd::Compiled&m,S x,S y,Q X,Q Y,const Point&q,const sf::Inverse&inverse,const Lambda&lambda,Q residual){
    if(!current)return;InverseFinal row={};row.slot=current->slot;row.source=m.base.source;row.target=current->target;row.chart=current->chart;row.inside=inverse.inside;row.input[0]=x;row.input[1]=y;row.output[0]=inverse.X;row.output[1]=inverse.Y;row.position[0]=q.position[0];row.position[1]=q.position[1];row.J=q.J;for(int i=0;i<4;i++)row.F[i]=q.F[i];for(int i=0;i<3;i++)row.lambda[i]=lambda[i].v;row.anchor[0]=m.base.anchor[0];row.anchor[1]=m.base.anchor[1];row.center[0]=X;row.center[1]=Y;row.L=inverse.L;row.radius=inverse.radius;row.residual=residual;current->inverse.push_back(row);
}
template<class Boundary,class Point>static void radius_piece(const sd::Compiled&m,const Boundary&b,int part,const Point&q){
    if(!current)return;RadiusPiece row={};row.slot=current->slot;row.source=m.base.source;row.target=current->target;row.index=current->radius_piece.size();row.part=part;row.first=b.first;row.last=b.last;row.canonical_edge=b.edge;row.lo=b.lo;row.hi=b.hi;for(int k=0;k<2;k++)row.position[k]=q.position[k];current->radius_piece.push_back(row);
}
template<class Bounds>static void radius_final(const sd::Compiled&m,const Bounds&bounds,const Q*lo,const Q*hi,Q r2){
    if(!current)return;RadiusFinal row={};row.slot=current->slot;row.source=m.base.source;row.target=current->target;row.r2=r2;for(int k=0;k<2;k++){row.reference_lo[k]=lo[k];row.reference_hi[k]=hi[k];row.displacement[k]=bounds.displacement[k];}current->radius.push_back(row);
}
}
#endif
