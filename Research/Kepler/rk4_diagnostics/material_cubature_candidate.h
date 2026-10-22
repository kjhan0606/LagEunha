#ifndef MATERIAL_CUBATURE_CANDIDATE_H
#define MATERIAL_CUBATURE_CANDIDATE_H
// Same native binary80 physical contractions as the accepted 24-source audit.
// This validates a rule candidate; full-cell/stage oracle checks still follow.
struct CubatureCandidate {
    bool pass=false;std::array<R,16>diagnostic{};
    std::vector<R>residual,rounding_bound;
};
static CubatureCandidate check_cubature_candidate(const full_evolution::Material&m,size_t s,
        const std::vector<STCoefficient>&k1,const R*fixed,size_t n,size_t rows,const std::vector<R>&A,
        const std::vector<R>&target,const std::vector<R>&absolute,const std::vector<uint32_t>&index,const std::vector<double>&weight){
    CubatureCandidate out;const auto&t=m.owned.f.triangle[s];size_t nn=index.size();if(!nn||nn!=weight.size()||4*nn>n)return out;
    for(size_t j=0;j<nn;j++)if(index[j]>=n||(j&&index[j-1]>=index[j])||!(weight[j]>0)||!std::isfinite(weight[j]))return out;
    std::vector<R>error(rows),bound(rows);R maximum=0;
    for(size_t r=0;r<rows;r++){R sum=0,abs_sum=0;for(size_t q=0;q<nn;q++){R z=A[r*n+index[q]]*weight[q];sum+=z;abs_sum+=fabsl(z);}error[r]=sum-target[r];bound[r]=(n+nn+2)*LDBL_EPSILON*(absolute[r]+abs_sum);if(absolute[r]>0)maximum=std::max(maximum,fabsl(error[r])/absolute[r]);else if(error[r]!=0)return out;}
    R E[ST_MAX][ST_MAX]={},Q[ST_MAX][ST_MAX]={},C[ST_MAX][ST_MAX]={},CB[ST_MAX][ST_MAX]={};size_t row=1;
    for(int a=0;a<t.n;a++)for(int b=a;b<t.n;b++,row++){E[a][b]=E[b][a]=error[row];Q[a][b]=Q[b][a]=target[row]+error[row];}
    for(int a=0;a<t.n;a++)for(int b=a+1;b<t.n;b++,row++){C[a][b]=error[row];C[b][a]=-error[row];CB[a][b]=CB[b][a]=fabsl(error[row])+bound[row];}
    size_t gradient=row+2,gravity=gradient+2*t.n,potential=gravity+2*t.n;need(potential+1==rows,"complete original moment system");R U=target[row],oldP=target[row+1]/target[0],newP=(target[row+1]+error[row+1])/(target[0]+error[0]);if(!(U>0&&target[0]>0))return out;
    R x[ST_MAX][2]={},v[ST_MAX][2]={},G[2][ST_MAX]={},dG[2][ST_MAX]={};const auto&origin=m.owned.f.coefficient[t.id[0]];
    for(int a=0;a<t.n;a++)for(int d=0;d<2;d++){const auto&c=m.owned.f.coefficient[t.id[a]];x[a][d]=double(c.x[d]-(a<10?origin.x[d]:0));v[a][d]=double(c.v[d]-(a<10?origin.v[d]:0));}
    R oldK=0,newK=0,dK=0,dL=0,dWork=0,dJ=0,dJdot=0,Jbound=0;
    for(int a=0;a<t.n;a++)for(int b=0;b<t.n;b++){
        const auto&ca=m.owned.f.coefficient[t.id[a]];const auto&cb=m.owned.f.coefficient[t.id[b]];
        for(int d=0;d<2;d++){R vv=R(ca.v[d])*cb.v[d];oldK+=.5L*vv*fixed[a*ST_MAX+b];newK+=.5L*vv*Q[a][b];dK+=.5L*vv*E[a][b];dWork+=R(ca.v[d])*k1[t.id[b]].v[d]*E[a][b];}
        dL+=(R(ca.x[0])*cb.v[1]-R(ca.x[1])*cb.v[0])*E[a][b];
        dJ+=x[a][0]*x[b][1]*C[a][b];dJdot+=(v[a][0]*x[b][1]+x[a][0]*v[b][1])*C[a][b];Jbound+=fabsl(x[a][0]*x[b][1])*CB[a][b];dG[0][a]+=x[b][1]*C[a][b];dG[1][b]+=x[a][0]*C[a][b];
    }
    size_t cstart=1+t.n*(t.n+1)/2,r=cstart;
    for(int a=0;a<t.n;a++)for(int b=a+1;b<t.n;b++,r++){R c=target[r];G[0][a]+=x[b][1]*c;G[0][b]-=x[a][1]*c;G[1][b]+=x[a][0]*c;G[1][a]-=x[b][0]*c;}
    R gerror=0,gscale=0,force_error=0,force_scale=0;
    for(int d=0;d<2;d++)for(int a=0;a<t.n;a++){gerror=std::max(gerror,fabsl(dG[d][a]));gscale=std::max(gscale,fabsl(G[d][a]));R f=-target[gradient+d*t.n+a]-oldP*G[d][a]+target[gravity+d*t.n+a],df=-error[gradient+d*t.n+a]-(newP-oldP)*G[d][a]-newP*dG[d][a]+error[gravity+d*t.n+a];force_scale=std::max(force_scale,fabsl(f));force_error=std::max(force_error,fabsl(df));}
    if(!(oldK>0&&gscale>0&&force_scale>0))return out;
    out.diagnostic={maximum,dK/U,(newK-oldK)/U,dL/(2*oldK),dWork/U,error[row]/U,dJ/target[0],dJdot/target[0],Jbound/target[0],gerror/gscale,(newP-oldP)/m.P0,force_error/force_scale,error[potential],target[potential],U,error[potential]/U};
    for(R z:out.diagnostic)if(!std::isfinite(z))return out;
    out.pass=fabsl(dK)/U<1e-8L&&fabsl(newK-oldK)/U<1e-8L&&fabsl(dL)/(2*oldK)<1e-11L&&fabsl(dWork)/U<1e-8L&&fabsl(error[row])/U<1e-8L&&fabsl(error[potential])/U<1e-8L&&Jbound/target[0]<1e-10L&&fabsl(dJdot)/target[0]<1e-8L&&gerror/gscale<1e-11L&&fabsl(newP-oldP)/m.P0<1e-11L&&force_error/force_scale<1e-11L;
    out.residual=std::move(error);out.rounding_bound=std::move(bound);return out;
}
#endif
