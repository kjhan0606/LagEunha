#define GFS_CURVED_TRANSFER_LIBRARY
#include "test_curved_transfer_v1.cpp"
#include "material_curved_topology.h"
struct IntegrationPanel {
    curved_topology::Panel geometry;
    Row fine;
    std::vector<QuadPoint> points;
    R estimated_error=0;
    int quadrature_depth=0;
};
struct AdaptiveResult {
    Row total;
    R quadrature_error_sum=0,event_mass_bound=0,event_K_bound=0;
    uint64_t evaluated_panels=0,accepted_panels=0;
    int complete=0;
};
static std::vector<std::vector<curved_topology::Panel>> read_guard(const char *directory,size_t n){
    std::vector<std::vector<curved_topology::Panel>> out(n);
    for(int rank=0;rank<4;rank++){
        char path[1024];snprintf(path,sizeof(path),"%s/radial_panels_rank%d.csv",directory,rank);FILE *f=fopen(path,"r");need(f,"preserved radial guard partition");char header[512];need(fgets(header,sizeof(header),f),"guard header");size_t triangle;curved_topology::Panel p;int stable;
        while(fscanf(f,"%zu,%Lg,%Lg,%d,%d,%Lg,%Lg\n",&triangle,&p.low,&p.high,&p.depth,&stable,&p.mass_bound,&p.K_bound)==7){need(triangle<n && p.low<p.high,"guard panel identity");p.stable=stable;out[triangle].push_back(p);}need(feof(f),"guard complete rows");fclose(f);
    }
    for(auto &v:out){need(!v.empty(),"every source triangle has guard leaves");std::sort(v.begin(),v.end(),[](const auto &a,const auto &b){return a.low<b.low;});for(size_t j=1;j<v.size();j++)need(v[j].low==v[j-1].high,"exact serialized radial coverage");}
    return out;
}
static void add_row(Row &to,const Row &from){
    to.radial_panels+=from.radial_panels;to.radial_queries+=from.radial_queries;to.root_panels+=from.root_panels;to.pair_polynomials+=from.pair_polynomials;to.angular_fragments+=from.angular_fragments;to.points+=from.points;to.unknown_points+=from.unknown_points;to.budget_failures+=from.budget_failures;to.unknown_reference_measure+=from.unknown_reference_measure;sum(to.angular_uncertainty,from.angular_uncertainty);for(int j=0;j<24;j++)sum(to.target[j],from.target[j]);
}
static AdaptiveResult adaptive(size_t it,const curved::Packet &p,const std::vector<curved_partition::Site> &sites,const std::vector<curved_topology::Panel> &guard,std::vector<QuadPoint> &retained){
    const auto &r=p.record[it];curved_bounds::Tri tri;for(int a=0;a<3;a++)for(int k=0;k<2;k++)tri[a][k]=r.triangle.geometry.x[a][k];auto bound=positive_source_bounds(r,tri,p.parameters[3]);R Ulower=p.parameters[2]*r.triangle.geometry.area*powl(bound.J.hi,1-p.parameters[1])/(p.parameters[1]-1),Mlower=p.parameters[6]*r.triangle.geometry.area,rmax=0;for(const auto &x:tri)rmax=std::max(rmax,hypotl(x[0],x[1]));R vmax=hypotl(curved_bounds::absmax(bound.velocity[0]),curved_bounds::absmax(bound.velocity[1])),Pallow=1e-12L*Mlower*std::max(vmax,1e-30L),Lallow=1e-12L*Mlower*std::max(rmax*vmax,1e-30L);std::vector<R> weights(sites.size(),0);std::vector<IntegrationPanel> panel;AdaptiveResult result;
    auto evaluate=[&](curved_topology::Panel geometry,int depth){
        IntegrationPanel out;out.geometry=geometry;out.quadrature_depth=depth;std::vector<R> interval={geometry.low,geometry.high};auto coarse=integrate(it,p,sites,weights,8,12,nullptr,&interval,&bound);out.fine=integrate(it,p,sites,weights,16,16,&out.points,&interval,&bound);result.evaluated_panels++;
        if(geometry.stable)for(int j=0;j<24;j++){
            const auto &a=coarse.target[j].a,&b=out.fine.target[j].a;R e=0;e=std::max(e,fabsl(a[MASS]-b[MASS])/(1e-12L*Mlower));e=std::max(e,fabsl(a[KINETIC]-b[KINETIC])/(1e-10L*Ulower));e=std::max(e,fabsl(a[POTENTIAL]-b[POTENTIAL])/(1e-10L*Ulower));e=std::max(e,fabsl(a[THERMAL]-b[THERMAL])/(1e-10L*Ulower));e=std::max(e,fabsl(a[THERMAL_RATE]-b[THERMAL_RATE])/(1e-9L*Ulower));e=std::max(e,fabsl(a[VOLUME]-b[VOLUME])/(1e-12L*r.triangle.geometry.area));e=std::max(e,fabsl(a[PX]-b[PX])/Pallow);e=std::max(e,fabsl(a[PY]-b[PY])/Pallow);e=std::max(e,fabsl(a[ANGULAR]-b[ANGULAR])/Lallow);out.estimated_error=std::max(out.estimated_error,8*e);
        }
        need(!coarse.budget_failures && !out.fine.budget_failures,"bounded adaptive angular root coverage");return out;
    };
    for(const auto &g:guard)panel.push_back(evaluate(g,0));
    while(true){
        R total=0,worst=0;int selected=-1;for(size_t j=0;j<panel.size();j++){total+=panel[j].estimated_error;if(panel[j].quadrature_depth<24 && panel[j].estimated_error>worst){worst=panel[j].estimated_error;selected=j;}}
        result.quadrature_error_sum=total;if(total<=1){result.complete=1;break;}if(selected<0 || result.evaluated_panels>=guard.size()+256)break;
        auto geometry=panel[selected].geometry;int depth=panel[selected].quadrature_depth+1;panel.erase(panel.begin()+selected);R mid=(geometry.low+geometry.high)/2;if(!(mid>geometry.low && mid<geometry.high))break;auto right=geometry;right.low=mid;geometry.high=mid;panel.push_back(evaluate(geometry,depth));panel.push_back(evaluate(right,depth));
    }
    std::sort(panel.begin(),panel.end(),[](const auto &a,const auto &b){return a.geometry.low<b.geometry.low;});result.total.triangle=it;for(size_t j=0;j<sites.size();j++)if(sites[j].id==r.cell)result.total.source_index=j;result.accepted_panels=panel.size();
    for(auto &v:panel){add_row(result.total,v.fine);result.event_mass_bound+=v.geometry.mass_bound;result.event_K_bound+=v.geometry.K_bound;retained.insert(retained.end(),v.points.begin(),v.points.end());}return result;
}
int main(int argc,char **argv){
    int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);active_mpi=true;need(provided>=MPI_THREAD_FUNNELED,"MPI FUNNELED");int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);need(argc==5,"packet, mapped sites, successful guard directory, new output");auto p=curved::read(argv[1]);auto sites=curved_partition::read_sites(argv[2]);auto guard=read_guard(argv[3],p.record.size());std::vector<AdaptiveResult> row(p.record.size());std::vector<std::vector<QuadPoint>> point(p.record.size());auto start=std::chrono::steady_clock::now();
    #pragma omp parallel for schedule(dynamic,1)
    for(int it=0;it<(int)p.record.size();it++)if(it%ranks==rank)row[it]=adaptive(it,p,sites,guard[it],point[it]);
    uint64_t npoint=0;for(const auto &v:point)npoint+=v.size();char path[1024];snprintf(path,sizeof(path),"%s/positive_quadrature_rank%d.bin",argv[4],rank);FILE *f=fopen(path,"wb");need(f,"new adaptive positive packet");uint64_t header[6]={0x474650334f564c31ULL,1,(uint64_t)rank,(uint64_t)ranks,npoint,56};need(fwrite(header,sizeof(header),1,f)==1,"adaptive packet header");for(const auto &v:point)if(v.size())need(fwrite(v.data(),sizeof(QuadPoint),v.size(),f)==v.size(),"retained adaptive quadrature");fclose(f);
    std::vector<AdaptiveResult> gathered(rank==0?row.size()*ranks:0);MPI_Gather(row.data(),row.size()*sizeof(AdaptiveResult),MPI_BYTE,gathered.data(),row.size()*sizeof(AdaptiveResult),MPI_BYTE,0,MPI_COMM_WORLD);
    if(rank==0){
        for(size_t it=0;it<row.size();it++)row[it]=gathered[(it%ranks)*row.size()+it];Moment transfer[24][24]={};uint64_t passed=0,points=0,evaluated=0,accepted=0;R max_error=0,eventK=0,eventM=0,angularK=0;std::string name=std::string(argv[4])+"/adaptive_rows.csv";f=fopen(name.c_str(),"w");need(f,"adaptive quadrature ledger");fprintf(f,"triangle,complete,evaluated_panels,accepted_panels,estimated_error_ratio,event_mass_bound,event_K_bound,angular_K_bound\n");for(const auto &v:row){for(int j=0;j<24;j++)sum(transfer[v.total.source_index][j],v.total.target[j]);passed+=v.complete;points+=v.total.points;evaluated+=v.evaluated_panels;accepted+=v.accepted_panels;max_error=std::max(max_error,v.quadrature_error_sum);eventK+=v.event_K_bound;eventM+=v.event_mass_bound;angularK+=v.total.angular_uncertainty.a[KINETIC];fprintf(f,"%llu,%d,%llu,%llu,%.21Lg,%.21Lg,%.21Lg,%.21Lg\n",(unsigned long long)v.total.triangle,v.complete,(unsigned long long)v.evaluated_panels,(unsigned long long)v.accepted_panels,v.quadrature_error_sum,v.event_mass_bound,v.event_K_bound,v.total.angular_uncertainty.a[KINETIC]);}fclose(f);
        name=std::string(argv[4])+"/transfer_moments.csv";f=fopen(name.c_str(),"w");need(f,"adaptive transported moment matrix");fprintf(f,"source,target,mass,px,py,angular,K,U,volume,qx,qy,polar,thermal_rate,potential\n");for(int i=0;i<24;i++)for(int j=0;j<24;j++){fprintf(f,"%llu,%llu",(unsigned long long)sites[i].id,(unsigned long long)sites[j].id);for(R v:transfer[i][j].a)fprintf(f,",%.21Lg",v);fprintf(f,"\n");}fclose(f);
        double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();printf("{\"scope\":\"W0 adaptive positive native curved quadrature on event-guarded radial panels; Gauss error is an estimate, radial event bounds separate; no weight solve, compression, exact C0 or long evolution\",\"ranks\":%d,\"OMP_threads\":%d,\"triangles\":%zu,\"complete_triangles\":%llu,\"positive_points\":%llu,\"evaluated_panels\":%llu,\"accepted_panels\":%llu,\"max_estimated_error_ratio\":%.17Lg,\"event_mass_bound\":%.17Lg,\"event_K_bound\":%.17Lg,\"angular_K_bound\":%.17Lg,\"all_quadrature_estimates_pass\":%s,\"seconds\":%.17g}\n",ranks,omp_get_max_threads(),row.size(),(unsigned long long)passed,(unsigned long long)points,(unsigned long long)evaluated,(unsigned long long)accepted,max_error,eventM,eventK,angularK,passed==row.size()?"true":"false",seconds);
    }active_mpi=false;MPI_Finalize();return 0;
}
