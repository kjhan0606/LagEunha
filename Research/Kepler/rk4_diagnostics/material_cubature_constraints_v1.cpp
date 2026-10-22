// One fixed representative set from the latest accepted full state.
// Native preprocessing for positive cubature; not an evolution/IC/mass rebuild.
#include "material_full_initial_driver.cpp"
struct ConstraintPoint {STPoint q;std::array<STJet,ST_MAX>N;GfsMaterialPoint12 p,p0;};
static ConstraintPoint native_point(const full_evolution::Material&m,const GfsEnrichConfig&cfg,size_t s,int index){
    const auto&t=m.owned.f.triangle[s];double x[3],y[3];for(int j=0;j<3;j++){x[j]=t.vertex[j][0];y[j]=t.vertex[j][1];}GfsMomentPolygon poly={x,y,3,0};ConstraintPoint p={};p.q.boundary=-1;p.q.triangle=s;
    need(gfs_enrich_quad(poly,&cfg,&m.owned.rule[s],index,&p.q.X,&p.q.Y,&p.q.weight),"accepted-state native reference point");if(p.q.weight==0)return p;
    double ic[4];need(p.q.weight>0&&gfs_enrich_reference(&cfg,p.q.X+12,p.q.Y+12,ic,0,0),"unchanged positive reference density");p.q.rho=ic[0];
    R ux=R(x[1])-x[0],uy=R(y[1])-y[0],vx=R(x[2])-x[0],vy=R(y[2])-y[0],dx=R(p.q.X)-x[0],dy=R(p.q.Y)-y[0],det=ux*vy-uy*vx;
    p.q.l[1]=(dx*vy-dy*vx)/det;p.q.l[2]=(ux*dy-uy*dx)/det;p.q.l[0]=1-p.q.l[1]-p.q.l[2];
    need(st_shape(t,p.q,p.N.data())&&st_point(t,m.owned.f.coefficient.data(),p.N.data(),&p.p)&&st_point(t,m.generating.data(),p.N.data(),&p.p0),"actual accepted/current generating fields");return p;
}
int main(int argc,char**argv){
    int provided;MPI_Init_thread(&argc,&argv,MPI_THREAD_FUNNELED,&provided);int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);
    try{
        {
        need(argc==5&&ranks==4,"saved_initial immutable_manifest accepted_prefix NEW_output; original four owners");double begin=MPI_Wtime();
        full_evolution::Material m;m.saved=argv[1];m.immutable=full_evolution::textfile(argv[2]);m.gamma=double(1.666667f);m.P0=1e-6;m.eps=.047;
        m.owned=material_full::load_basis(m.saved+"/full_conditioned_basis_r"+std::to_string(rank)+".bin");md::Ownership layout(m.owned.global,m.owned.owner,MPI_COMM_WORLD);full_evolution::thermal_import(m,rank);full_evolution::checkpoint(m,layout,std::string(argv[3])+std::to_string(rank)+".bin",false);full_evolution::sync(m,layout);
        need(m.ledger.steps>=4&&m.ledger.stage==0,"latest actually accepted state only");auto geometry=gp::geometry(m.owned.f);GfsEnrichConfig cfg;config_rule(cfg,m.gamma,m.P0,m.eps);
        std::set<size_t>selected;for(double target:{0.,.047,.5,2.,8.,12.}){size_t best=0;double error=1e99;for(size_t s=0;s<m.owned.f.triangle.size();s++){const auto&t=m.owned.f.triangle[s];double x=0,y=0;for(int j=0;j<3;j++){x+=t.vertex[j][0]/3;y+=t.vertex[j][1]/3;}double e=fabs(hypot(x,y)-target);if(e<error){error=e;best=s;}}selected.insert(best);}
        std::map<size_t,R>delta;for(size_t s:selected)delta.emplace(m.owned.f.triangle[s].cell,0);
        std::map<size_t,std::array<R,2>>cell;
        for(size_t s=0;s<m.owned.f.triangle.size();s++)if(delta.count(m.owned.f.triangle[s].cell)){
            const auto&t=m.owned.f.triangle[s];double x[3],y[3];for(int j=0;j<3;j++){x[j]=t.vertex[j][0];y[j]=t.vertex[j][1];}GfsMomentPolygon poly={x,y,3,0};auto&sum=cell[t.cell];
            for(int index=0;index<gfs_enrich_nquad(poly,&cfg,&m.owned.rule[s]);index++){auto p=native_point(m,cfg,s,index);if(p.q.weight){sum[0]+=p.q.weight;sum[1]+=R(p.q.weight)*p.p.J;}}
        }
        for(auto&entry:delta){auto&sum=cell.at(entry.first);need(sum[0]>0,"positive whole-source accepted anchor measure");entry.second=double((geometry.volume[entry.first]-sum[1])/sum[0]);}
        std::string dest=argv[4];FILE*manifest=fopen((dest+"/constraints_rank"+std::to_string(rank)+".jsonl").c_str(),"wx");need(manifest,"new fixed representative manifest");
        uint64_t executed=0;
        for(size_t s:selected){
            double start=MPI_Wtime();const auto&t=m.owned.f.triangle[s];std::vector<ConstraintPoint>point;std::vector<uint32_t>slot;double x[3],y[3];for(int j=0;j<3;j++){x[j]=t.vertex[j][0];y[j]=t.vertex[j][1];}GfsMomentPolygon poly={x,y,3,0};
            for(int index=0;index<gfs_enrich_nquad(poly,&cfg,&m.owned.rule[s]);index++){auto p=native_point(m,cfg,s,index);if(p.q.weight){point.push_back(std::move(p));slot.push_back(index);}}
            size_t rows=1+t.n*(t.n+1)/2+t.n*(t.n-1)/2+2+4*t.n+1,n=point.size();std::vector<R>A(rows*n),target(rows),absolute(rows);std::vector<double>w(n);
            for(size_t q=0;q<n;q++){
                const auto&p=point[q];double a=p.p.J+double(delta.at(t.cell)),a0=p.p0.J+double(m.generating_delta[t.cell][0]);need(a>0&&a0>0,"immutable anchor thermal context");
                double lr=log1p((a-a0)/a0),P=m.P0*exp(-m.gamma*lr),U=m.P0*a0/(m.gamma-1)*exp((1-m.gamma)*lr);
                double z=p.p.position[0]*p.p.position[0]+p.p.position[1]*p.p.position[1]+m.eps*m.eps,g=pow(z,-1.5),acc[2]={-p.p.position[0]*g,-p.p.position[1]*g};
                size_t r=0;auto put=[&](R value){A[r*n+q]=value;r++;};put(1);
                for(int j=0;j<t.n;j++)for(int k=j;k<t.n;k++)put(R(p.q.rho)*p.N[j].v*p.N[k].v);
                for(int j=0;j<t.n;j++)for(int k=j+1;k<t.n;k++)put(R(p.N[j].x)*p.N[k].y-R(p.N[j].y)*p.N[k].x);
                put(U);put(P);
                for(int d=0;d<2;d++)for(int j=0;j<t.n;j++){R derivative=d==0?R(p.N[j].x)*p.p.F[3]-R(p.N[j].y)*p.p.F[2]:R(p.N[j].y)*p.p.F[0]-R(p.N[j].x)*p.p.F[1];put(-R(P)*derivative);}
                for(int d=0;d<2;d++)for(int j=0;j<t.n;j++)put(R(p.q.rho)*p.N[j].v*acc[d]);
                put(-R(p.q.rho)/sqrt(z));need(r==rows,"explicit mass/determinant/heat/gradient/potential/gravity constraints");w[q]=p.q.weight;
            }
            for(size_t r=0;r<rows;r++)for(size_t q=0;q<n;q++){R z=A[r*n+q]*w[q];target[r]+=z;absolute[r]+=fabsl(z);}
            std::string name="constraints_r"+std::to_string(rank)+"_triangle"+std::to_string(m.owned.source_triangle[s])+".bin";FILE*out=fopen((dest+"/"+name).c_str(),"wbx");need(out,"new native representative constraints");uint64_t h[8]={0x4d46435542434e31ULL,1,m.owned.source_triangle[s],rows,n,sizeof(R),uint64_t(t.n),uint64_t(m.owned.f.cell[t.cell])};material_full::write(out,h,8);material_full::write(out,slot.data(),n);material_full::write(out,w.data(),n);material_full::write(out,A.data(),A.size());material_full::write(out,target.data(),rows);material_full::write(out,absolute.data(),rows);material_full::durable(out);
            fprintf(manifest,"{\"triangle\":%llu,\"cell\":%llu,\"rows\":%zu,\"positive_points\":%zu,\"accepted_time\":%.17g,\"source_delta\":%.21Lg,\"file\":\"%s\",\"native_constraint_seconds\":%.9g,\"weights_changed\":false,\"candidate_rule_committed\":false}\n",(unsigned long long)m.owned.source_triangle[s],(unsigned long long)m.owned.f.cell[t.cell],rows,n,m.ledger.time,delta.at(t.cell),name.c_str(),MPI_Wtime()-start);fflush(manifest);executed+=n;
        }
        material_full::durable(manifest);uint64_t global=executed;MPI_Allreduce(MPI_IN_PLACE,&global,1,MPI_UINT64_T,MPI_SUM,MPI_COMM_WORLD);if(!rank)printf("{\"scope\":\"fixed representative native accepted-state cubature feasibility constraints\",\"positive_triangle_points_exported\":%llu,\"seconds\":%.9g,\"new_RK_step\":false,\"mass_operator_changed\":false,\"rule_committed\":false}\n",(unsigned long long)global,MPI_Wtime()-begin);
        } // Destroy Ownership's duplicated communicator while MPI is active.
        MPI_Finalize();return 0;
    }catch(const std::exception&e){fprintf(stderr,"rank%d CONSTRAINT_PREPROCESS_FAIL %s\n",rank,e.what());MPI_Abort(MPI_COMM_WORLD,2);}return 2;
}
