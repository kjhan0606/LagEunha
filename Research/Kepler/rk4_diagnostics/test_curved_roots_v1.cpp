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
#include <mpi.h>
#include <omp.h>
static bool active_mpi=false;
static void need(bool ok,const char *message){if(!ok){fprintf(stderr,"FAIL: %s\n",message);if(active_mpi && !omp_in_parallel())MPI_Abort(MPI_COMM_WORLD,2);exit(2);}}
#include "material_curved_packet.h"
#include "material_curved_roots.h"
using R=long double;using I=curved_bounds::I;using W=curved_roots::Window;
struct Site {uint64_t id;R x,y;};
struct Row {uint64_t triangle,source_intervals=0,pair_polynomials=0,root_bands=0,classified_fragments=0,transferred_fragments=0,root_panels=0,budget_failures=0,native_owner_mismatches=0;R source_angle=0,target_uncertain_angle=0,source_uncertain_angle=0,max_native_plane_excess=0;};
static std::vector<Site> sites(const char *path){FILE *f=fopen(path,"r");need(f,"actual mapped site CSV");char header[512];need(fgets(header,sizeof(header),f),"site CSV header");std::vector<Site> result;unsigned long long id;double X,Y,x,y,vx,vy;while(fscanf(f,"%llu,%lf,%lf,%lf,%lf,%lf,%lf\n",&id,&X,&Y,&x,&y,&vx,&vy)==7)result.push_back({id,x,y});need(feof(f) && result.size()==24,"all mapped sites, no steering");fclose(f);return result;}
static std::vector<W> merge(std::vector<W> input){std::sort(input.begin(),input.end(),[](W a,W b){return a.lo<b.lo;});std::vector<W> result;for(auto w:input){if(!result.empty() && w.lo<=result.back().hi)result.back().hi=std::max(result.back().hi,w.hi);else result.push_back(w);}return result;}
static curved_roots::P plane(const curved_roots::Ray &ray,const Site &a,const Site &b,R wa,R wb,R uncertainty){I ax=I(b.x)-I(a.x),ay=I(b.y)-I(a.y),c=(I(b.x)*I(b.x)+I(b.y)*I(b.y)-I(a.x)*I(a.x)-I(a.y)*I(a.y)+I(wa)-I(wb))/I(2);return curved_roots::plane(ray,ax,ay,c,uncertainty);}
static uint64_t fixtures(){
    std::vector<std::vector<R>> roots={{-.85L,-.4L,-.00003L,.00003L,.3L,.8L},{-.2L,.1L,.1000001L,.5L},{.2L,.2L},{-1,0,1},{.123455L,.123457L}};
    uint64_t panels=0;for(const auto &known:roots){auto p=curved_roots::constant(I(1));for(R x:known){curved_roots::P factor(1);factor.c[0]=I(-x);factor.c[1]=I(1);p=curved_roots::mul(p,factor);}auto out=curved_roots::isolate(p,-1,1,1e-10L,42,32768);need(out.budget_complete,"root fixture bounded isolation");for(R x:known){bool contained=false;for(auto w:out.window)contained=contained || (w.lo<=x && x<=w.hi);need(contained,"simple, endpoint, close or repeated root retained");}need(out.uncertain_length<1e-6L,"fixture uncertainty bounded, no root deletion");panels+=out.panels;}
    return panels;
}
static Row process(size_t it,const curved::Packet &p,const std::vector<Site> &site){
    const auto &r=p.record[it];Row row={};row.triangle=it;R X=0,Y=0;for(int a=0;a<3;a++){X+=r.triangle.geometry.x[a][0]/3.L;Y+=r.triangle.geometry.x[a][1]/3.L;}R radius=hypotl(X,Y);need(radius>0,"positive source centroid radius");int own=-1;for(int j=0;j<(int)site.size();j++)if(site[j].id==r.cell)own=j;need(own>=0,"original cell site");
    auto bound=curved_bounds::enclose(r,curved_bounds::Tri{std::array<R,2>{X,Y},std::array<R,2>{X,Y},std::array<R,2>{X,Y}},p.parameters[3]);R uncertainty=bound.native_position_envelope;
    for(int chart:{1,-1}){
        auto ray=curved_roots::ray(r,I(radius),chart,p.parameters[3]);R source_uncertain=0;uint64_t source_panels=0;auto source=curved_roots::source_intervals(ray,source_uncertain,source_panels);row.source_uncertain_angle+=2*source_uncertain;row.root_panels+=source_panels;
        for(int weights=0;weights<2;weights++)for(auto segment:source){
            row.source_intervals++;row.source_angle+=2*(atanl(segment.hi)-atanl(segment.lo));std::vector<R> weight(site.size());for(size_t j=0;j<site.size();j++)weight[j]=weights?1e-7L*sinl(.7L*j):0;
            std::vector<int> candidate={own};for(int j=0;j<(int)site.size();j++)if(j!=own){auto f=plane(ray,site[own],site[j],weight[own],weight[j],uncertainty);if(curved_roots::sign(curved_roots::bernstein(f,segment.lo,segment.hi),f.degree)!=-1)candidate.push_back(j);}
            std::vector<W> uncertain;for(size_t i=0;i<candidate.size();i++)for(size_t j=i+1;j<candidate.size();j++){
                int a=candidate[i],b=candidate[j];auto f=plane(ray,site[a],site[b],weight[a],weight[b],uncertainty);auto cut=curved_roots::isolate(f,segment.lo,segment.hi,1e-10L,42,16384);row.pair_polynomials++;row.root_panels+=cut.panels;row.budget_failures+=!cut.budget_complete;uncertain.insert(uncertain.end(),cut.window.begin(),cut.window.end());
            }
            auto bands=merge(uncertain);row.root_bands+=bands.size();R left=segment.lo;
            auto classify=[&](R a,R b){if(b<=a)return;R mid=(a+b)/2,D=1+mid*mid;I xx=curved_roots::value(ray.position[0],I(mid))/I(D*D*D),yy=curved_roots::value(ray.position[1],I(mid))/I(D*D*D);int winner=0;R best=1e99L;
                for(int j=0;j<(int)site.size();j++){R score=site[j].x*site[j].x+site[j].y*site[j].y-site[j].x*(xx.lo+xx.hi)-site[j].y*(yy.lo+yy.hi)-weight[j];if(score<best){best=score;winner=j;}}
                for(int j=0;j<(int)site.size();j++)if(j!=winner){auto f=plane(ray,site[winner],site[j],weight[winner],weight[j],uncertainty);I v=curved_roots::value(f,I(mid));need(v.hi<0,"target label certified outside retained root bands");}
                R x=chart*radius*(1-mid*mid)/D,y=chart*radius*2*mid/D;auto native=curved::native_point(r,(double)x,(double)y,p.parameters[3]);int native_winner=0;R native_best=1e99L;
                for(int j=0;j<(int)site.size();j++){R score=site[j].x*site[j].x+site[j].y*site[j].y-2*site[j].x*native.position[0]-2*site[j].y*native.position[1]-weight[j];if(score<native_best){native_best=score;native_winner=j;}}
                row.native_owner_mismatches+=native_winner!=winner;row.classified_fragments++;row.transferred_fragments+=winner!=own;
                for(int j=0;j<(int)site.size();j++)if(j!=winner){auto f=plane(ray,site[winner],site[j],weight[winner],weight[j],uncertainty);I v=curved_roots::value(f,I(mid))/I(D*D*D);R ax=site[j].x-site[winner].x,ay=site[j].y-site[winner].y,c=(site[j].x*site[j].x+site[j].y*site[j].y-site[winner].x*site[winner].x-site[winner].y*site[winner].y+weight[winner]-weight[j])/2,truth=ax*native.position[0]+ay*native.position[1]-c;row.max_native_plane_excess=std::max({row.max_native_plane_excess,v.lo-truth,truth-v.hi});}
            };
            for(auto w:bands){classify(left,w.lo);row.target_uncertain_angle+=2*(atanl(w.hi)-atanl(w.lo));left=w.hi;}classify(left,segment.hi);
        }
    }return row;
}
int main(int argc,char **argv){
    int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);active_mpi=true;need(provided>=MPI_THREAD_FUNNELED,"MPI FUNNELED");int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);need(argc==4,"packet, advected site CSV, new output directory");auto p=curved::read(argv[1]);auto site=sites(argv[2]);auto start=std::chrono::steady_clock::now();uint64_t fixture_panels=0;if(rank==0)fixture_panels=fixtures();
    std::vector<Row> row(p.record.size());
    #pragma omp parallel for schedule(static)
    for(int it=0;it<(int)p.record.size();it++)if(it%ranks==rank)row[it]=process(it,p,site);
    std::vector<Row> gathered(rank==0?row.size()*ranks:0);MPI_Gather(row.data(),row.size()*sizeof(Row),MPI_BYTE,gathered.data(),row.size()*sizeof(Row),MPI_BYTE,0,MPI_COMM_WORLD);
    if(rank==0){
        for(size_t i=0;i<row.size();i++)row[i]=gathered[(i%ranks)*row.size()+i];std::string filename=std::string(argv[3])+"/angular_partition.csv";FILE *f=fopen(filename.c_str(),"w");need(f,"angular partition output");fprintf(f,"triangle,source_intervals,pair_polynomials,root_bands,classified_fragments,transferred_fragments,root_panels,budget_failures,native_owner_mismatches,source_angle,target_uncertain_angle,source_uncertain_angle,max_native_plane_excess\n");Row total={};R max_fraction=0;
        for(const auto &r:row){fprintf(f,"%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%llu,%.21Lg,%.21Lg,%.21Lg,%.21Lg\n",(unsigned long long)r.triangle,(unsigned long long)r.source_intervals,(unsigned long long)r.pair_polynomials,(unsigned long long)r.root_bands,(unsigned long long)r.classified_fragments,(unsigned long long)r.transferred_fragments,(unsigned long long)r.root_panels,(unsigned long long)r.budget_failures,(unsigned long long)r.native_owner_mismatches,r.source_angle,r.target_uncertain_angle,r.source_uncertain_angle,r.max_native_plane_excess);total.source_intervals+=r.source_intervals;total.pair_polynomials+=r.pair_polynomials;total.root_bands+=r.root_bands;total.classified_fragments+=r.classified_fragments;total.transferred_fragments+=r.transferred_fragments;total.root_panels+=r.root_panels;total.budget_failures+=r.budget_failures;total.native_owner_mismatches+=r.native_owner_mismatches;total.source_angle+=r.source_angle;total.target_uncertain_angle+=r.target_uncertain_angle;total.source_uncertain_angle+=r.source_uncertain_angle;total.max_native_plane_excess=std::max(total.max_native_plane_excess,r.max_native_plane_excess);max_fraction=std::max(max_fraction,(r.target_uncertain_angle+2*r.source_uncertain_angle)/r.source_angle);}
        fclose(f);double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();printf("{\"scope\":\"interval coefficient angular ray root coverage and actual curved source/target partitions at one reference radius per triangle; native rounding envelope retained; no radial integration, overlap moments, injectivity proof, weight solve or re-expression\",\"ranks\":%d,\"OMP_threads\":%d,\"triangles\":%zu,\"fixtures\":5,\"fixture_panels\":%llu,\"source_intervals\":%llu,\"pair_polynomials\":%llu,\"root_bands\":%llu,\"classified_fragments\":%llu,\"transferred_fragments\":%llu,\"root_panels\":%llu,\"budget_failures\":%llu,\"native_owner_mismatches\":%llu,\"native_plane_enclosure_excess\":%.17Lg,\"source_angle_sum_two_weight_cases\":%.17Lg,\"target_uncertain_angle_sum\":%.17Lg,\"source_uncertain_angle_sum_one_geometry\":%.17Lg,\"max_uncertain_angular_fraction_bound\":%.17Lg,\"seconds\":%.17g}\n",ranks,omp_get_max_threads(),row.size(),(unsigned long long)fixture_panels,(unsigned long long)total.source_intervals,(unsigned long long)total.pair_polynomials,(unsigned long long)total.root_bands,(unsigned long long)total.classified_fragments,(unsigned long long)total.transferred_fragments,(unsigned long long)total.root_panels,(unsigned long long)total.budget_failures,(unsigned long long)total.native_owner_mismatches,total.max_native_plane_excess,total.source_angle,total.target_uncertain_angle,total.source_uncertain_angle,max_fraction,seconds);
        need(!total.budget_failures && !total.native_owner_mismatches && total.max_native_plane_excess==0,"bounded native angular geometry coverage");need(total.transferred_fragments>0,"actual finite transport, not identity partition");
    }active_mpi=false;MPI_Finalize();return 0;
}
