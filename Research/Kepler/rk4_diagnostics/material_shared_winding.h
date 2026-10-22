#ifndef MATERIAL_SHARED_WINDING_H
#define MATERIAL_SHARED_WINDING_H
#include "material_shared_target_geometry.h"
namespace shared_target {
enum class Location {outside=0,inside=1,boundary=2,unproved=3};
struct Membership {Location location=Location::unproved;int winding=0,ray=-1;uint64_t crossings=0,endpoint_events=0,windows=0;int failed_edge=-1;Roots root_failure;char reason[192]={};};
static Plane ray_plane(S x,S y,int nx,int ny){Plane p;p.en[0]=Exact(nx);p.en[1]=Exact(ny);p.n[0]=S(nx);p.n[1]=S(ny);p.eb=Exact(nx)*Exact(x.n)+Exact(ny)*Exact(y.n);p.b=ball(p.eb);p.input_offset_radius=sf::add(sf::mul(Q(abs(nx)),x.e),sf::mul(Q(abs(ny)),y.e));return p;}
// Four geometrically fixed ray normals. A ray is selected before root work
// only if the WHOLE real source Lipschitz enclosure makes every edge
// transverse and vertex predicates are exact or separated. This is not a
// precision/radius search. No Cartesian inverse is used for membership.
static Membership winding(const sd::Compiled&m,const sf::RealBounds&rb,const std::map<int,Edge>&edges,S x,S y,int forced_ray=-1){
    Membership out;try{need(rb.L<1,"winding requires whole-source injectivity");const auto&t=m.base.state->triangle[m.base.source];for(int a=0;a<3;a++){const auto&v=m.base.state->coefficient[t.id[a]];if(x.e==0&&y.e==0&&x.n==Q(v.x[0])&&y.n==Q(v.x[1])){out.location=Location::boundary;snprintf(out.reason,sizeof(out.reason),"exact_original_global_vertex");return out;}}
        int normals[4][2]={{0,1},{1,0},{1,1},{1,-1}};S derivative[3];Plane transverse,along;for(int ray=0;ray<4;ray++){if(forced_ray>=0&&ray!=forced_ray)continue;int nx=normals[ray][0],ny=normals[ray][1];auto candidate=ray_plane(x,y,nx,ny);bool valid=true;for(int e=0;e<3;e++){const auto&reference=m.base.edge[e];S slope=S(nx)*reference.Dref[0]+S(ny)*reference.Dref[1];Q radius=sf::mul(Q(abs(nx)+abs(ny)),sf::mul(rb.L,std::max(sf::mag(reference.Dref[0]),sf::mag(reference.Dref[1]))));derivative[e]=S(slope.n,sf::add(slope.e,radius));valid=valid&&(sf::lower(derivative[e])>0||sf::upper(derivative[e])<0);auto projection=project(edges.at(t.edge[e].poly[1]),candidate);for(int j:{0,3}){S endpoint=projection.node[j];endpoint.e=sf::add(endpoint.e,candidate.input_offset_radius);valid=valid&&(sf::lower(endpoint)>0||sf::upper(endpoint)<0||(endpoint.n==0&&endpoint.e==0));}}if(valid){out.ray=ray;transverse=candidate;along=ray_plane(x,y,ny,-nx);break;}}
        need(out.ray>=0,"no separated finite geometric winding ray; boundary event remains explicit");for(int e=0;e<3;e++){const auto&edge=edges.at(t.edge[e].poly[1]);auto roots=isolate(edge,transverse,&derivative[e]);out.windows+=roots.evaluations;if(!roots.complete){out.failed_edge=edge.key;out.root_failure=roots;}need(roots.complete,roots.reason);auto projection=project(edge,along);bool forward=(t.endpoint[e][1]==(t.endpoint[e][0]+1)%3);int orientation=forward?1:-1,up=sf::lower(derivative[e])>0?orientation:-orientation;for(const auto&root:roots.root){S coordinate=projected(projection,root.t).v;need(sf::lower(coordinate)>0||sf::upper(coordinate)<0,"point box meets canonical boundary; boundary identity required");if(sf::upper(coordinate)<0)continue;bool include=true;if(root.endpoint>=0){out.endpoint_events++;bool start=forward?root.endpoint==0:root.endpoint==1;include=(start&&up>0)||(!start&&up<0);}if(include){out.winding+=up;out.crossings++;}}}
        need(out.winding==0||out.winding==1,"CCW injective source has winding zero or one");out.location=out.winding?Location::inside:Location::outside;snprintf(out.reason,sizeof(out.reason),"canonical_transverse_ray_and_half_open_endpoint_winding");
    }catch(const std::exception&e){snprintf(out.reason,sizeof(out.reason),"%s",e.what());}return out;
}
static int inside_or_outside(const Membership&m){need(m.location==Location::inside||m.location==Location::outside,m.reason);return m.location==Location::inside?1:0;}
}
#endif
