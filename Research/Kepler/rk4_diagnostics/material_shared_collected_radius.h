#ifndef MATERIAL_SHARED_COLLECTED_RADIUS_H
#define MATERIAL_SHARED_COLLECTED_RADIUS_H
namespace shared_target {
// SAME reference radius. Collect repeated angular dependence before ball
// magnitudes, instead of summing independent barycentric-direction balls.
// At an origin apex D(t)=(B-A)+t(C-B), r=s*sqrt(D.D). No sampled fit.
static S collected_reference_radius(const sd::Compiled&m,const sd::Chart&chart,S s,S t,const sd::Point&point,S*ray=nullptr,S*argument=nullptr){
 bool origin=true;for(int k=0;k<2;k++)origin=origin&&m.base.vertex[chart.apex][k].n==0&&m.base.vertex[chart.apex][k].e==0;
 if(!origin)return reference_radius(m,chart,s,t,point);
 auto exact_direction=frame(m,chart,t);if(ray)for(int k=0;k<2;k++)ray[k]=exact_direction.D[k];S square=exact_direction.D[0]*exact_direction.D[0]+exact_direction.D[1]*exact_direction.D[1];if(argument)*argument=square;return s*complex_sqrt(square);
}
}
#endif
