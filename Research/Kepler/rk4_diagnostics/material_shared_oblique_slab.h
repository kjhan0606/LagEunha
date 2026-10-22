#ifndef MATERIAL_SHARED_OBLIQUE_SLAB_H
#define MATERIAL_SHARED_OBLIQUE_SLAB_H
namespace oblique_slab {
struct Internal {uint64_t slot,source,target,chart,face,side,index,count,complete,windows;S parameter,s,t,reference[2],projection[2];char reason[192];};
static_assert(sizeof(Internal)==496,"saved internal event ABI");
struct Event {uint64_t slot,source,target,chart,index,kind;int64_t face,side,node;S u,v,s,t;char identity[192];};
struct Slab {uint64_t slot,source,target,chart,index,empty;int64_t geometric_lower,geometric_upper,lower,upper;S u0,u1,sample,geometric_vlo,geometric_vhi,selected_vlo,selected_vhi;};
struct Proof {uint64_t slot,slab,face,case_id,iterations;S derivative,left,right,root;};
struct Work {Frame29 frame;std::vector<Event>events;std::vector<Slab>slabs;std::vector<Proof>proofs;Slab pending_slab={};Proof pending_proof={};bool ordered=false,complete=false;uint64_t skipped=0;std::string reason;};
static S edge_v(const Frame29&f,int side,S u){int b=(side+1)%4;S du=f.projected[b][0]-f.projected[side][0];need(sf::lower(du)>0||sf::upper(du)<0,"oblique section boundary transverse to projection");return f.projected[side][1]+(u-f.projected[side][0])*(f.projected[b][1]-f.projected[side][1])/du;}
static S value(const sd::Compiled&m,const target::Face&p,const Frame29&f,S u,S v){S X[2];for(int k=0;k<2;k++)X[k]=f.apex[k]+u*f.transverse[k]+v*f.direction[k];auto c=target::coordinates(m,f.original,X[0],X[1]);return target::face(p,f.original,c.s,c.t).v;}
static S graph(const sd::Compiled&m,const target::Face&p,const Frame29&f,S u,S lower,S upper,Proof&proof){
 Q lo=sf::lower(lower),hi=sf::upper(upper);S d=proof.derivative;need(sf::lower(d)>0||sf::upper(d)<0,"global oblique graph derivative separated");
 for(int it=0;it<256;it++){proof.iterations=it+1;S box=sf::range(lo,hi);need(box.n>sf::upper(lower)&&box.n<sf::lower(upper),"oblique graph center inside original trapezoid section");S v=value(m,p,f,u,S(box.n)),N=S(box.n)-v/d;Q a=std::max(lo,sf::lower(N)),b=std::min(hi,sf::upper(N));need(a<=b,"oblique root inclusion");if(b-a<1e-27Q)return sf::range(a,b);if(b-a<Q(.75L)*(hi-lo)){lo=a;hi=b;continue;}bool increasing=sf::lower(d)>0;if(sf::lower(v)>0){if(increasing)hi=box.n;else lo=box.n;continue;}if(sf::upper(v)<0){if(increasing)lo=box.n;else hi=box.n;continue;}return sf::range(a,b);}
 throw std::runtime_error("fixed256 oblique graph work budget");
}
static void build(const sd::Compiled&m,const std::vector<Site>&sites,const std::vector<Event29>&nodes,const std::vector<Internal>&internal,Work&w){const auto&f=w.frame;
 auto add=[&](uint64_t kind,int64_t face,int64_t side,int64_t node,S u,S v,S s,S t,const std::string&id){Event e={};e.slot=f.slot;e.source=f.source;e.target=f.target;e.chart=f.chart;e.kind=kind;e.face=face;e.side=side;e.node=node;e.u=u;e.v=v;e.s=s;e.t=t;need(id.size()<sizeof(e.identity),"full event identity fits");snprintf(e.identity,sizeof(e.identity),"%s",id.c_str());w.events.push_back(e);};
 for(int a=0;a<4;a++)add(0,-1,a,-1,f.projected[a][0],f.projected[a][1],S(a==0||a==3?f.cut:Q(1)),S(a>=2?1:0),"corner;"+std::to_string(f.slot)+";"+std::to_string(a));
 for(const auto&e:internal)if(e.slot==f.slot){need(e.complete,"retained complete internal root family");if(e.index!=UINT64_MAX)add(1,e.face,e.side,-1,e.projection[0],e.projection[1],e.s,e.t,"internal;"+std::to_string(f.slot)+";"+std::to_string(e.side)+";"+std::to_string(e.face)+";"+std::to_string(e.index));}
 const auto&tri=m.base.state->triangle[m.base.source];int other=f.original.part?(f.original.apex+2)%3:(f.original.apex+1)%3,key=-1;
 for(int e=0;e<3;e++){int a=tri.endpoint[e][0],b=tri.endpoint[e][1];if((a==f.original.apex&&b==other)||(b==f.original.apex&&a==other))key=m.base.edge[e].native.poly[1];}need(key>=0,"canonical source-edge chart identity");
 for(const auto&e:nodes)if(e.slot==f.slot){need(e.passed,"retained used reference event passed");bool apex=true;for(int k=0;k<2;k++)apex=apex&&e.reference[k].n==f.apex[k].n&&e.reference[k].e==0;if(apex){w.skipped++;continue;}auto c=target::coordinates(m,f.original,e.reference[0],e.reference[1]);bool canonical=e.kind==1&&std::string(e.event).rfind("e;"+std::to_string(key)+";",0)==0;if(canonical)c.t=S(f.original.part?1:0);
  if(sf::upper(c.s)<f.cut||sf::lower(c.s)>1||sf::upper(c.t)<0||sf::lower(c.t)>1){w.skipped++;continue;}
  need(sf::lower(c.s)>f.cut&&sf::upper(c.s)<1,"used event radial chart-boundary ownership unresolved");need(canonical||(sf::lower(c.t)>0&&sf::upper(c.t)<1),"used power/internal chart-boundary ownership unresolved");S uv[2];project29(f,e.reference,uv);add(e.kind==2?3:2,-1,canonical?(f.original.part?2:0):-1,e.node,uv[0],uv[1],c.s,c.t,e.event);
 }
 std::sort(w.events.begin(),w.events.end(),[](const Event&a,const Event&b){return a.u.n<b.u.n;});for(size_t j=0;j<w.events.size();j++)w.events[j].index=j;for(size_t j=1;j<w.events.size();j++)need(sf::upper(w.events[j-1].u)<sf::lower(w.events[j].u),"oblique event projections ordered without coalescing");w.ordered=true;
 auto planes=target::planes(sites,f.target);std::vector<target::Face>faces;for(const auto&p:planes)faces.push_back(target::pullback(m,p));
 for(size_t j=0;j+1<w.events.size();j++){Slab z={};z.slot=f.slot;z.source=f.source;z.target=f.target;z.chart=f.chart;z.index=j;z.u0=w.events[j].u;z.u1=w.events[j+1].u;z.sample=S(((z.u0+z.u1)/S(2)).n);z.geometric_lower=z.geometric_upper=z.lower=z.upper=-1;
  std::vector<std::pair<S,int>>hit;for(int a=0;a<4;a++){int b=(a+1)%4;Q lo=std::min(f.projected[a][0].n,f.projected[b][0].n),hi=std::max(f.projected[a][0].n,f.projected[b][0].n);if(z.sample.n>lo&&z.sample.n<hi)hit.push_back({edge_v(f,a,z.sample),a});}need(hit.size()==2,"open projection slab meets exactly two trapezoid sides");if(hit[1].first.n<hit[0].first.n)std::swap(hit[0],hit[1]);z.geometric_lower=hit[0].second;z.geometric_upper=hit[1].second;z.geometric_vlo=z.selected_vlo=hit[0].first;z.geometric_vhi=z.selected_vhi=hit[1].first;need(sf::upper(z.geometric_vlo)<sf::lower(z.geometric_vhi),"positive geometric section");
  w.pending_slab=z;
  for(size_t face=0;face<planes.size();face++)if(f.active_mask&(uint64_t(1)<<face)){w.pending_proof=Proof{};auto&p=w.pending_proof;p.slot=f.slot;p.slab=j;p.face=face;p.derivative=f.derivative[face];p.left=value(m,faces[face],f,z.sample,z.geometric_vlo);p.right=value(m,faces[face],f,z.sample,z.geometric_vhi);
   if(sf::upper(p.left)<0&&sf::upper(p.right)<0)p.case_id=0;
   else if(sf::lower(p.left)>0&&sf::lower(p.right)>0){p.case_id=1;z.empty=1;}
   else{p.case_id=2;need((sf::upper(p.left)<0&&sf::lower(p.right)>0)||(sf::lower(p.left)>0&&sf::upper(p.right)<0),"oblique graph endpoint signs separated");p.root=graph(m,faces[face],f,z.sample,z.geometric_vlo,z.geometric_vhi,p);
    if(sf::lower(p.derivative)>0){if(sf::upper(p.root)<sf::lower(z.selected_vhi)){z.selected_vhi=p.root;z.upper=face;}else need(sf::lower(p.root)>sf::upper(z.selected_vhi),"distinct upper oblique selected branches");}
    else{if(sf::lower(p.root)>sf::upper(z.selected_vlo)){z.selected_vlo=p.root;z.lower=face;}else need(sf::upper(p.root)<sf::lower(z.selected_vlo),"distinct lower oblique selected branches");}
   }w.proofs.push_back(p);
  }
  if(sf::lower(z.selected_vlo)>sf::upper(z.selected_vhi))z.empty=1;else if(!z.empty)need(sf::upper(z.selected_vlo)<sf::lower(z.selected_vhi),"positive selected section without discarded sliver");w.slabs.push_back(z);
 }w.complete=true;w.reason="complete_actual_oblique_projected_event_and_selected_branch_arrangement";
}
}
#endif
