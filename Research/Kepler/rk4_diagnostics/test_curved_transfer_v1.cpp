#include <algorithm>
#include <array>
#include <vector>
#include <map>
#include <set>
#include <string>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cfloat>
#include <chrono>
#ifdef GFS_CURVED_PROOF_EXCEPTIONS
#include <stdexcept>
#endif
#include <mpi.h>
#include <omp.h>
static bool active_mpi=false;
static void need(bool ok,const char *message){if(!ok){
#ifdef GFS_CURVED_PROOF_EXCEPTIONS
throw std::runtime_error(message);
#else
fprintf(stderr,"FAIL: %s\n",message);if(active_mpi && !omp_in_parallel())MPI_Abort(MPI_COMM_WORLD,2);exit(2);
#endif
}}
#include "material_curved_packet.h"
#include "material_curved_trace_bounds.h"
#include "material_curved_partition.h"
using R=long double;using I=curved_bounds::I;using Pt=std::array<double,2>;using Key=std::pair<Pt,Pt>;
static curved_bounds::Enclosure positive_source_bounds(const curved::Record &r,const curved_bounds::Tri &tri,R eps){
    auto out=curved_bounds::enclose(r,tri,eps);if(out.J.lo>0)return out;curved_bounds::Certificate c;curved_bounds::certify(r,tri,eps,0,6,.24L,c);need(c.complete && c.min_J>0,"bounded positive source-triangle certificate");out.J=I(c.min_J,c.max_J);return out;
}
enum {MASS,PX,PY,ANGULAR,KINETIC,THERMAL,VOLUME,QX,QY,POLAR,THERMAL_RATE,POTENTIAL,NMOMENT};
struct Moment {R a[NMOMENT]={};};
static void sum(Moment &a,const Moment &b){for(int k=0;k<NMOMENT;k++)a.a[k]+=b.a[k];}
struct QuadPoint {uint64_t triangle,target,flags;double X,Y,reference_weight,rho0;};
static_assert(sizeof(QuadPoint)==56,"frozen positive quadrature record size");
struct Row {uint64_t triangle,source_index,radial_panels=0,radial_queries=0,root_panels=0,pair_polynomials=0,angular_fragments=0,points=0,unknown_points=0,budget_failures=0;Moment target[24],unknown_target[24],angular_uncertainty;R unknown_reference_measure=0;};
static void gauss(int n,std::vector<R> &x,std::vector<R> &w){x.resize(n);w.resize(n);for(int i=0;i<(n+1)/2;i++){R z=cosl(acosl(-1.L)*(i+.75L)/(n+.5L)),d=0;for(int it=0;it<100;it++){R a=1,b=z;for(int k=2;k<=n;k++){R c=((2*k-1)*z*b-(k-1)*a)/k;a=b;b=c;}d=n*(z*b-a)/(z*z-1);R v=z-b/d;if(fabsl(v-z)<4*LDBL_EPSILON){z=v;break;}z=v;}R a=1,b=z;for(int k=2;k<=n;k++){R c=((2*k-1)*z*b-(k-1)*a)/k;a=b;b=c;}d=n*(z*b-a)/(z*z-1);x[i]=(1-z)/2;x[n-1-i]=(1+z)/2;w[i]=w[n-1-i]=1/((1-z*z)*d*d);}}
static R density(const curved::Packet &p,R r){return p.parameters[6]+p.parameters[9]/(r+p.parameters[7])*.5L*(1+tanhl((r-p.parameters[4])/p.parameters[8]))*.5L*(1+tanhl((p.parameters[5]-r)/p.parameters[8]));}
static bool inside(const curved::Record &r,R x,R y){const auto &g=r.triangle.geometry;R dx=x-g.x[0][0],dy=y-g.x[0][1],ux=(R)g.x[1][0]-g.x[0][0],uy=(R)g.x[1][1]-g.x[0][1],vx=(R)g.x[2][0]-g.x[0][0],vy=(R)g.x[2][1]-g.x[0][1],det=ux*vy-uy*vx,a=dx*vy-dy*vx,b=ux*dy-uy*dx;return a>=0 && b>=0 && a+b<=det;}
static std::vector<R> radial_cuts(const curved::Record &r){
    const auto &g=r.triangle.geometry;std::vector<R> cuts;R minimum=1e99L,maximum=0;for(int a=0;a<3;a++){R v=hypotl(g.x[a][0],g.x[a][1]);cuts.push_back(v);maximum=std::max(maximum,v);int b=(a+1)%3;R x=g.x[a][0],y=g.x[a][1],dx=(R)g.x[b][0]-x,dy=(R)g.x[b][1]-y,t=std::max(0.L,std::min(1.L,-(x*dx+y*dy)/(dx*dx+dy*dy))),d=hypotl(x+t*dx,y+t*dy);minimum=std::min(minimum,d);if(t>0 && t<1)cuts.push_back(d);}
    if(inside(r,0,0))minimum=0;cuts.push_back(minimum);cuts.push_back(maximum);for(R a=1e-6L;a<maximum;a*=2)if(a>minimum)cuts.push_back(a);std::sort(cuts.begin(),cuts.end());std::vector<R> v;for(R c:cuts)if(c>=minimum && c<=maximum && (v.empty() || c>v.back()+32*LDBL_EPSILON*std::max(1.L,fabsl(c))))v.push_back(c);need(v.size()>=2,"finite source radial support");return v;
}
static int owner(const std::vector<curved_partition::Site> &site,const std::vector<R> &weight,const GfsMaterialPoint12 &p){int j=0;R best=1e99L;for(size_t k=0;k<site.size();k++){R f=site[k].x*site[k].x+site[k].y*site[k].y-2*site[k].x*p.position[0]-2*site[k].y*p.position[1]-weight[k];if(f<best){best=f;j=k;}}return j;}
static Moment integrand(const curved::Packet &packet,const GfsMaterialPoint12 &p,R rho,R w){Moment z;R dm=rho*w,x=p.position[0],y=p.position[1],vx=p.velocity[0],vy=p.velocity[1],P=packet.parameters[2]*powl(p.J,-packet.parameters[1]);z.a[MASS]=dm;z.a[PX]=dm*vx;z.a[PY]=dm*vy;z.a[ANGULAR]=dm*(x*vy-y*vx);z.a[KINETIC]=.5L*dm*(vx*vx+vy*vy);z.a[THERMAL]=P*p.J*w/(packet.parameters[1]-1);z.a[VOLUME]=p.J*w;z.a[QX]=dm*x;z.a[QY]=dm*y;z.a[POLAR]=dm*(x*x+y*y);z.a[THERMAL_RATE]=-P*p.J*w*(p.grad_v[0]+p.grad_v[3]);z.a[POTENTIAL]=-dm/sqrtl(x*x+y*y+packet.parameters[3]*packet.parameters[3]);return z;}
static std::vector<int> source_candidates(const curved::Record &r,const curved_bounds::Enclosure &bound,const std::vector<curved_partition::Site> &site,const std::vector<R> &weight,int own){
    I position[2];for(int k=0;k<2;k++){R lo=r.triangle.geometry.x[0][k],hi=lo;for(int a=1;a<3;a++){lo=std::min(lo,(R)r.triangle.geometry.x[a][k]);hi=std::max(hi,(R)r.triangle.geometry.x[a][k]);}position[k]=I(lo,hi)+bound.displacement[k];}std::vector<int> candidate={own};for(int j=0;j<(int)site.size();j++)if(j!=own){I nx=I(site[j].x)-I(site[own].x),ny=I(site[j].y)-I(site[own].y),c=(curved_bounds::square(I(site[j].x))+curved_bounds::square(I(site[j].y))-curved_bounds::square(I(site[own].x))-curved_bounds::square(I(site[own].y))+I(weight[own])-I(weight[j]))/I(2);if((nx*position[0]+ny*position[1]-c).hi>=0)candidate.push_back(j);}return candidate;
}
static Row integrate(size_t it,const curved::Packet &p,const std::vector<curved_partition::Site> &site,const std::vector<R> &weight,int nr,int na,std::vector<QuadPoint> *saved,const std::vector<R> *override_cuts=nullptr,const curved_bounds::Enclosure *fixed_bounds=nullptr,bool outgoing_only=false
#ifdef MATERIAL_CURVED_COMPILED_H
    ,const curved_compiled::Map *cached_map=nullptr,const std::vector<int> *whole_source_candidates=nullptr
#endif
    ){
    const auto &r=p.record[it];Row row={};row.triangle=it;for(size_t i=0;i<site.size();i++)if(site[i].id==r.cell)row.source_index=i;curved_bounds::Tri tri;for(int a=0;a<3;a++)for(int k=0;k<2;k++)tri[a][k]=r.triangle.geometry.x[a][k];auto bound=fixed_bounds?*fixed_bounds:positive_source_bounds(r,tri,p.parameters[3]);need(bound.J.lo>0,"interval positive source map");R vmax[2]={curved_bounds::absmax(bound.velocity[0]),curved_bounds::absmax(bound.velocity[1])},v2=vmax[0]*vmax[0]+vmax[1]*vmax[1],xmax[2]={};for(int k=0;k<2;k++){for(int a=0;a<3;a++)xmax[k]=std::max(xmax[k],fabsl(tri[a][k]));xmax[k]+=curved_bounds::absmax(bound.displacement[k]);}
    std::vector<R> rn,rw,an,aw;gauss(nr,rn,rw);gauss(na,an,aw);auto cuts=override_cuts?*override_cuts:radial_cuts(r);row.radial_panels=cuts.size()-1;R pi=acosl(-1.L);
#ifdef MATERIAL_CURVED_COMPILED_H
    curved_compiled::Map local_map;if(!cached_map){local_map=curved_compiled::compile(r,p.parameters[3]);cached_map=&local_map;}const auto &compiled=*cached_map;
#endif
    for(size_t panel=0;panel+1<cuts.size();panel++)for(int iq=0;iq<nr;iq++){
        R a=cuts[panel],b=cuts[panel+1],t=rn[iq],s=sinl(pi*t/2),radius=a+(b-a)*s*s,dr=(b-a)*pi/2*sinl(pi*t)*rw[iq],rho=density(p,radius);row.radial_queries++;
#ifdef MATERIAL_CURVED_PARTITION_RAY_CACHE
        auto part=curved_partition::partition(r,radius,p.parameters[3],site,weight,bound.native_position_envelope,1e-12L,&compiled,whole_source_candidates);
#elif defined(MATERIAL_CURVED_COMPILED_H)
        auto part=curved_partition::partition(r,radius,p.parameters[3],site,weight,bound.native_position_envelope,1e-12L,&compiled);
#else
        auto part=curved_partition::partition(r,radius,p.parameters[3],site,weight,bound.native_position_envelope,1e-12L);
#endif
        row.root_panels+=part.root_panels;row.pair_polynomials+=part.pair_polynomials;row.budget_failures+=!part.complete;
        auto angular=[&](int chart,R left,R right,int target,int flags){if(outgoing_only && target==(int)row.source_index)return;R theta0=2*atanl(left),theta1=2*atanl(right);if(theta1<=theta0)return;row.angular_fragments++;for(int j=0;j<na;j++){R theta=theta0+(theta1-theta0)*an[j],X=chart*radius*cosl(theta),Y=chart*radius*sinl(theta);if(flags==1 && !inside(r,X,Y))continue;R w=radius*dr*(theta1-theta0)*aw[j];need(w>0,"positive common-refinement quadrature weight");auto point=curved::native_point(r,(double)X,(double)Y,p.parameters[3]);int dest=target>=0?target:owner(site,weight,point);if(outgoing_only && dest==(int)row.source_index)continue;auto z=integrand(p,point,rho,w);sum(row.target[dest],z);if(flags)sum(row.unknown_target[dest],z);row.points++;row.unknown_points+=flags!=0;if(saved)saved->push_back({it,(uint64_t)dest,(uint64_t)flags,(double)X,(double)Y,(double)w,(double)rho});}};
        for(const auto &f:part.fragment)angular(f.chart,f.left,f.right,f.target,0);
        auto unknown=[&](const std::vector<std::pair<int,curved_roots::Window>> &v,int flags){for(auto entry:v){if(outgoing_only){
#ifdef MATERIAL_CURVED_PARTITION_RAY_CACHE
                const auto &ray=part.chart_ray[entry.first==1?0:1];
#elif defined(MATERIAL_CURVED_COMPILED_H)
                auto ray=curved_roots::ray(r,I(radius),entry.first,p.parameters[3],&compiled);
#else
                auto ray=curved_roots::ray(r,I(radius),entry.first,p.parameters[3]);
#endif
                bool retained_by_own=true;for(int j=0;j<(int)site.size();j++)if(j!=(int)row.source_index
#ifdef MATERIAL_CURVED_COMPILED_H
                    && (!whole_source_candidates || std::find(whole_source_candidates->begin(),whole_source_candidates->end(),j)!=whole_source_candidates->end())
#endif
                    ){auto f=curved_partition::plane(ray,site[row.source_index],site[j],weight[row.source_index],weight[j],bound.native_position_envelope);if(curved_roots::sign(curved_roots::bernstein(f,entry.second.lo,entry.second.hi),f.degree)!=-1){retained_by_own=false;break;}}if(retained_by_own)continue;}
            R width=2*(atanl(entry.second.hi)-atanl(entry.second.lo)),measure=radius*dr*width;row.unknown_reference_measure+=measure;auto &u=row.angular_uncertainty.a;u[MASS]+=measure*rho;u[PX]+=measure*rho*vmax[0];u[PY]+=measure*rho*vmax[1];u[ANGULAR]+=measure*rho*(xmax[0]*vmax[1]+xmax[1]*vmax[0]);u[KINETIC]+=measure*rho*v2/2;u[THERMAL]+=measure*p.parameters[2]*powl(bound.J.lo,1-p.parameters[1])/(p.parameters[1]-1);u[VOLUME]+=measure*bound.J.hi;u[QX]+=measure*rho*xmax[0];u[QY]+=measure*rho*xmax[1];u[POLAR]+=measure*rho*(xmax[0]*xmax[0]+xmax[1]*xmax[1]);u[POTENTIAL]+=measure*rho/p.parameters[3];
#ifdef MATERIAL_CURVED_VELOCITY_GRADIENT_BOUNDS
            u[THERMAL_RATE]+=measure*p.parameters[2]*powl(bound.J.lo,1-p.parameters[1])*curved_bounds::absmax(curved_bounds::divergence(bound));
#endif
            angular(entry.first,entry.second.lo,entry.second.hi,-1,flags);}};
        unknown(part.source_unknown,1);unknown(part.target_unknown,2);
    }return row;
}
struct Edge {size_t triangle;int component;Pt a,b;};
static void trace_audit(const curved::Packet &p,const char *dest){
    std::map<Key,std::vector<Edge>> edge;for(size_t it=0;it<p.record.size();it++){const auto &r=p.record[it];for(int j=0;j<3;j++){Pt a={r.triangle.geometry.x[j][0],r.triangle.geometry.x[j][1]},b={r.triangle.geometry.x[(j+1)%3][0],r.triangle.geometry.x[(j+1)%3][1]};edge[a<b?Key(a,b):Key(b,a)].push_back({it,r.component,a,b});}}
    std::string filename=std::string(dest)+"/whole_trace_bounds.csv";FILE *f=fopen(filename.c_str(),"w");need(f,"whole-edge trace bound ledger");fprintf(f,"triangle_left,triangle_right,length,analytic_position_jump_bound,analytic_velocity_jump_bound\n");R maxx=0,maxv=0;uint64_t count=0;
    for(const auto &entry:edge)if(entry.second.size()==2){const auto &a=entry.second[0],&b=entry.second[1];R len=hypotl((R)a.a[0]-a.b[0],(R)a.a[1]-a.b[1]);auto x1=curved_trace::edge_form(p.record[a.triangle],a.a,a.b,p.parameters[3]),x2=curved_trace::edge_form(p.record[b.triangle],a.a,a.b,p.parameters[3]);auto v1=curved_trace::edge_form(p.record[a.triangle],a.a,a.b,p.parameters[3],true),v2=curved_trace::edge_form(p.record[b.triangle],a.a,a.b,p.parameters[3],true);R dx=curved_trace::difference_bound(x1,x2,p.parameters[3]),dv=curved_trace::difference_bound(v1,v2,p.parameters[3]);fprintf(f,"%zu,%zu,%.21Lg,%.21Lg,%.21Lg\n",a.triangle,b.triangle,len,dx,dv);maxx=std::max(maxx,dx);maxv=std::max(maxv,dv);count++;}fclose(f);
    filename=std::string(dest)+"/whole_trace_bounds.json";f=fopen(filename.c_str(),"w");need(f,"trace bound summary");fprintf(f,"{\"scope\":\"entire shared-edge analytic functional mismatch enclosures, all stored coefficients; not exact C0 equality or native libm proof\",\"edges\":%llu,\"max_position_jump_bound\":%.17Lg,\"max_velocity_jump_bound\":%.17Lg,\"exact_C0_proved\":false}\n",(unsigned long long)count,maxx,maxv);fclose(f);
}
#ifndef GFS_CURVED_TRANSFER_LIBRARY
int main(int argc,char **argv){
    int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);active_mpi=true;need(provided>=MPI_THREAD_FUNNELED,"MPI FUNNELED");int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);need(argc==8 || argc==9,"packet, sites, destination, radial order, angular order, weighted, save quadrature, optional debug triangle");auto p=curved::read(argv[1]);auto site=curved_partition::read_sites(argv[2]);int nr=atoi(argv[4]),na=atoi(argv[5]),weighted=atoi(argv[6]),save=atoi(argv[7]),filter=argc==9?atoi(argv[8]):-1;need(nr>=4 && na>=4,"positive bounded integration orders");std::vector<R> weights(site.size());for(size_t j=0;j<site.size();j++)weights[j]=weighted?1e-7L*sinl(.7L*j):0;auto started=std::chrono::steady_clock::now();if(rank==0)trace_audit(p,argv[3]);
    std::vector<Row> row(p.record.size());std::vector<std::vector<QuadPoint>> point(save?p.record.size():0);
    #pragma omp parallel for schedule(static)
    for(int i=0;i<(int)p.record.size();i++)if(i%ranks==rank && (filter<0 || filter==i))row[i]=integrate(i,p,site,weights,nr,na,save?&point[i]:nullptr);
    if(save){uint64_t count=0;for(const auto &v:point)count+=v.size();char name[1024];snprintf(name,sizeof(name),"%s/positive_quadrature_rank%d.bin",argv[3],rank);FILE *f=fopen(name,"wb");need(f,"new positive transported quadrature packet");uint64_t h[6]={0x474650334f564c31ULL,1,(uint64_t)rank,(uint64_t)ranks,count,56};need(fwrite(h,sizeof(h),1,f)==1,"positive quadrature header");for(const auto &v:point)if(v.size())need(fwrite(v.data(),sizeof(QuadPoint),v.size(),f)==v.size(),"positive source-field quadrature");fclose(f);}
    std::vector<Row> gathered(rank==0?row.size()*ranks:0);MPI_Gather(row.data(),row.size()*sizeof(Row),MPI_BYTE,gathered.data(),row.size()*sizeof(Row),MPI_BYTE,0,MPI_COMM_WORLD);
    if(rank==0){
        for(size_t i=0;i<row.size();i++)row[i]=gathered[(i%ranks)*row.size()+i];Moment transfer[24][24]={},uncertainty={};uint64_t npoint=0,nunknown=0,nroot=0,npair=0,nradial=0,nfragment=0,fail=0;for(const auto &r:row){for(int j=0;j<24;j++)sum(transfer[r.source_index][j],r.target[j]);sum(uncertainty,r.angular_uncertainty);npoint+=r.points;nunknown+=r.unknown_points;nroot+=r.root_panels;npair+=r.pair_polynomials;nradial+=r.radial_queries;nfragment+=r.angular_fragments;fail+=r.budget_failures;}
        std::string path=std::string(argv[3])+"/transfer_moments.csv";FILE *f=fopen(path.c_str(),"w");need(f,"common-reference transported moments");fprintf(f,"source,target,mass,px,py,angular,K,U,volume,qx,qy,polar,thermal_rate,potential\n");for(int i=0;i<24;i++)for(int j=0;j<24;j++){fprintf(f,"%llu,%llu",(unsigned long long)site[i].id,(unsigned long long)site[j].id);for(R x:transfer[i][j].a)fprintf(f,",%.21Lg",x);fprintf(f,"\n");}fclose(f);
        path=std::string(argv[3])+"/angular_uncertainty.csv";f=fopen(path.c_str(),"w");need(f,"angular uncertainty ledger");fprintf(f,"triangle,source,radial_queries,root_panels,points,unknown_points,mass_bound,K_bound,U_bound,reference_measure\n");for(const auto &r:row)fprintf(f,"%llu,%llu,%llu,%llu,%llu,%llu,%.21Lg,%.21Lg,%.21Lg,%.21Lg\n",(unsigned long long)r.triangle,(unsigned long long)site[r.source_index].id,(unsigned long long)r.radial_queries,(unsigned long long)r.root_panels,(unsigned long long)r.points,(unsigned long long)r.unknown_points,r.angular_uncertainty.a[MASS],r.angular_uncertainty.a[KINETIC],r.angular_uncertainty.a[THERMAL],r.unknown_reference_measure);fclose(f);
        double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();printf("{\"scope\":\"actual curved source-to-Laguerre positive common-reference quadrature estimates of m,p,L,K,U,V and thermal rate; original source fields retained; angular uncertainty bounded at radial quadrature nodes only; radial topology not yet certified; no weight Newton, compression or time evolution\",\"ranks\":%d,\"OMP_threads\":%d,\"radial_order\":%d,\"angular_order\":%d,\"weighted_pattern\":%s,\"triangles\":%zu,\"radial_queries\":%llu,\"angular_fragments\":%llu,\"root_panels\":%llu,\"pair_polynomials\":%llu,\"positive_points\":%llu,\"points_in_retained_uncertain_bands\":%llu,\"budget_failures\":%llu,\"angular_mass_uncertainty_bound\":%.17Lg,\"angular_K_uncertainty_bound\":%.17Lg,\"angular_U_uncertainty_bound\":%.17Lg,\"radial_topology_proved\":false,\"exact_shared_C0_proved\":false,\"cold_remap_promoted\":false,\"seconds\":%.17g}\n",ranks,omp_get_max_threads(),nr,na,weighted?"true":"false",row.size(),(unsigned long long)nradial,(unsigned long long)nfragment,(unsigned long long)nroot,(unsigned long long)npair,(unsigned long long)npoint,(unsigned long long)nunknown,(unsigned long long)fail,uncertainty.a[MASS],uncertainty.a[KINETIC],uncertainty.a[THERMAL],seconds);need(!fail,"fixed-budget angular isolation for all radial quadrature nodes");
    }active_mpi=false;MPI_Finalize();return 0;
}
#endif
