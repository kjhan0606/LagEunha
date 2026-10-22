#ifndef MATERIAL_SHARED_TRACE_RHS_H
#define MATERIAL_SHARED_TRACE_RHS_H
namespace shared_trace_rhs {
using R=long double;using Vec=std::vector<R>;
struct Factor {
    int n=0;R L[ST_MAX][ST_MAX]={};
    Factor(const R*M,int count):n(count){for(int i=0;i<n;i++)for(int j=0;j<=i;j++){R value=M[i*ST_MAX+j];for(int k=0;k<j;k++)value-=L[i][k]*L[j][k];need(i!=j||value>0,"local shared-trace mass Cholesky SPD");L[i][j]=i==j?sqrtl(value):value/L[j][j];}}
    void solve(R*x)const{for(int i=0;i<n;i++){for(int j=0;j<i;j++)x[i]-=L[i][j]*x[j];x[i]/=L[i][i];}for(int i=n-1;i>=0;i--){for(int j=i+1;j<n;j++)x[i]-=L[j][i]*x[j];x[i]/=L[i][i];}}
};
static R dot(const Vec&a,const Vec&b){R sum=0;for(size_t i=0;i<a.size();i++)sum+=a[i]*b[i];return sum;}
struct Mass {
    const st::Factory&f;Vec local,value,partition;std::vector<int> row,column;std::vector<Factor> factor;
    Mass(const st::Factory&factory,const Vec&raw):f(factory),local(raw.size()),partition(f.dof.size()){
        std::vector<std::map<int,R>> rows(f.dof.size());
        for(size_t s=0;s<f.triangle.size();s++){const auto&t=f.triangle[s];R T[ST_MAX][ST_MAX]={},temp[ST_MAX][ST_MAX]={};for(int a=0;a<t.n;a++)T[a][a]=1;
            for(int a=0;a<2;a++){int r=t.bubble_offset+a;for(int b=0;b<t.n;b++)T[r][b]=0;for(int b=0;b<2;b++){T[r][t.bubble_offset+b]=t.bubble_encode[a][b];T[r][9]-=R(t.bubble_encode[a][b])*t.bubble_shift[b];}}
            const R*M=raw.data()+s*ST_MAX*ST_MAX;R*dest=local.data()+s*ST_MAX*ST_MAX;for(int a=0;a<t.n;a++)for(int b=0;b<t.n;b++)for(int k=0;k<t.n;k++)temp[a][b]+=T[a][k]*M[k*ST_MAX+b];
            for(int a=0;a<t.n;a++)for(int b=0;b<=a;b++){R z=0;for(int k=0;k<t.n;k++)z+=temp[a][k]*T[b][k];dest[a*ST_MAX+b]=dest[b*ST_MAX+a]=z;}
            factor.emplace_back(dest,t.n);for(int a=0;a<t.n;a++){partition[t.id[a]]++;for(int b=0;b<t.n;b++)rows[t.id[a]][t.id[b]]+=dest[a*ST_MAX+b];}
        }
        row.push_back(0);for(auto&r:rows){for(auto&v:r){column.push_back(v.first);value.push_back(v.second);}row.push_back(value.size());}for(auto&v:partition){need(v>0,"every canonical DOF has support");v=1/v;}
    }
    Vec apply(const Vec&x)const{Vec y(x.size());for(size_t i=0;i<x.size();i++)for(int j=row[i];j<row[i+1];j++)y[i]+=value[j]*x[column[j]];return y;}
    Vec precondition(const Vec&r)const{Vec z(r.size());for(size_t s=0;s<f.triangle.size();s++){const auto&t=f.triangle[s];R local[ST_MAX];for(int a=0;a<t.n;a++)local[a]=r[t.id[a]]*partition[t.id[a]];factor[s].solve(local);for(int a=0;a<t.n;a++)z[t.id[a]]+=partition[t.id[a]]*local[a];}return z;}
};
struct Solution {Vec x;int iterations=0;R true_relative_residual=0;};
static Solution solve(const Mass&M,const Vec&b){
    Solution out;out.x.resize(b.size());Vec r=b,z=M.precondition(r),p=z;R rz=dot(r,z),norm=dot(b,b);need(norm>0&&rz>0,"positive sparse mass solve residual metric");
    for(size_t it=0;it<4*b.size();it++){auto Ap=M.apply(p);R den=dot(p,Ap);need(den>0,"sparse mass operator positive");R alpha=rz/den;for(size_t i=0;i<b.size();i++){out.x[i]+=alpha*p[i];r[i]-=alpha*Ap[i];}out.iterations=it+1;
        if(dot(r,r)<norm*1e-32L){auto exact=M.apply(out.x);for(size_t i=0;i<b.size();i++)exact[i]=b[i]-exact[i];out.true_relative_residual=sqrtl(dot(exact,exact)/norm);if(out.true_relative_residual<1e-16L)return out;r=std::move(exact);z=M.precondition(r);p=z;rz=dot(r,z);continue;}
        z=M.precondition(r);R next=dot(r,z);need(next>0,"positive additive element preconditioner");R beta=next/rz;for(size_t i=0;i<b.size();i++)p[i]=z[i]+beta*p[i];rz=next;
    }need(false,"fixed-work sparse mass solve did not meet original equation accuracy");return out;
}
struct Boundary {R work=0,torque=0,force[2]={0,0};};
static Boundary external(const st::Factory&f,const STPhysics&physics,Vec&element,int rank,int ranks){
    std::map<int,std::vector<std::pair<int,int>>> edges;for(size_t s=0;s<f.triangle.size();s++)for(int a=0;a<3;a++)edges[f.triangle[s].edge[a].poly[1]].push_back({s,a});
    std::vector<double>x,w;st::gauss(64,x,w);Boundary out;
    for(auto&entry:edges)if(entry.second.size()==1){int s=entry.second[0].first,a=entry.second[0].second;if(s%ranks!=rank)continue;const auto&t=f.triangle[s];const auto&e=t.edge[a];int b=(a+1)%3;double dx=t.vertex[b][0]-t.vertex[a][0],dy=t.vertex[b][1]-t.vertex[a][1],length=hypot(dx,dy),nx=dy/length,ny=-dx/length;
        for(int j=0;j<64;j++){STPoint q={};q.boundary=a;q.triangle=s;q.X=e.A[0]+x[j]*e.D[0];q.Y=e.A[1]+x[j]*e.D[1];q.l[t.endpoint[a][0]]=1-x[j];q.l[t.endpoint[a][1]]=x[j];STJet N[ST_MAX];GfsMaterialPoint12 p;need(st_shape(t,q,N)&&st_point(t,f.coefficient.data(),N,&p)&&st_trace(e,f.coefficient.data(),x[j],p.position,p.velocity),"external boundary original shared state");double P=physics.P0*pow(p.J,-physics.gamma),sigma[3]={-P,0,-P},traction[2],ratio;need(gfs_material12_piola_traction(&p,sigma,nx,ny,traction,&ratio),"native external Piola traction");R weight=R(length)*w[j];for(int d=0;d<2;d++){out.force[d]+=weight*traction[d];for(int k=0;k<t.n;k++)element[(s*2+d)*ST_MAX+k]+=weight*N[k].v*traction[d];}out.work+=weight*(R(p.velocity[0])*traction[0]+R(p.velocity[1])*traction[1]);out.torque+=weight*(R(p.position[0])*traction[1]-R(p.position[1])*traction[0]);}
    }
    R fields[4]={out.work,out.torque,out.force[0],out.force[1]};MPI_Allreduce(MPI_IN_PLACE,fields,4,MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);out.work=fields[0];out.torque=fields[1];out.force[0]=fields[2];out.force[1]=fields[3];return out;
}
static std::array<Vec,2> forces(const st::Factory&f,const std::vector<STPoint>&point,const std::vector<STOutput>&output,const STPhysics&physics,int rank,int ranks,Boundary&boundary){
    Vec element(f.triangle.size()*2*ST_MAX);for(size_t i=0;i<point.size();i++){const auto&q=point[i];need(output[i].valid,"valid native stress point before sparse assembly");for(int d=0;d<2;d++)for(int a=0;a<f.triangle[q.triangle].n;a++)element[(q.triangle*2+d)*ST_MAX+a]+=output[i].force[d][a];}
    boundary=external(f,physics,element,rank,ranks);MPI_Allreduce(MPI_IN_PLACE,element.data(),element.size(),MPI_LONG_DOUBLE,MPI_SUM,MPI_COMM_WORLD);std::array<Vec,2> out={Vec(f.dof.size()),Vec(f.dof.size())};
    for(size_t s=0;s<f.triangle.size();s++)for(int d=0;d<2;d++)for(int a=0;a<f.triangle[s].n;a++)out[d][f.triangle[s].id[a]]+=element[(s*2+d)*ST_MAX+a];return out;
}
}
#endif
