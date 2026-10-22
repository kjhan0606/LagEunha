#ifndef GFS_MATERIAL_P3_REFERENCE_TRACE_H
#define GFS_MATERIAL_P3_REFERENCE_TRACE_H
/* Reference trace rank/constraints. Geometry and IC only; no force inputs. */
struct P3Row {std::vector<std::pair<int,Real>> entries;int face=-1,trace=-1;};
struct P3Face {
    int triangle[2],cell[2],endpoint[2],node[4],selected[6]={-1,-1,-1,-1,-1,-1},residual_rank=0;
    Real combination[2][2]={{1,0},{0,1}};double length,nx,ny;Vec gram;
    Real rank_uncertainty=0,omitted_function_max=0,condition=0,quadrature_difference=0;
};
static void p3_face_shape(const P3Mesh &m,const GfsEnrichConfig &cfg,const std::vector<Pair> &scale,const P3Face &f,int side,double q,double N[12],double Nx[12],double Ny[12]){
    const auto &a=m.linear.node[f.endpoint[0]],&b=m.linear.node[f.endpoint[1]];double x=(1-q)*a[0]+q*b[0],y=(1-q)*a[1]+q*b[1],v[2],B[4];velocity(cfg,x+12,y+12,v,B);const auto &t=m.linear.tri[f.triangle[side]];
    p3_shape_stable(m.triangle[f.triangle[side]],cfg,x,y,scale[t.component].data(),N,Nx,Ny);
}
static void p3_face_functions(double q,const double N[12],Real phi[6]){
    double a=1-q,b=q;phi[0]=.5*a*(3*a-1)*(3*a-2);phi[1]=4.5*a*b*(3*a-1);phi[2]=4.5*a*b*(3*b-1);phi[3]=.5*b*(3*b-1)*(3*b-2);phi[4]=N[10];phi[5]=N[11];
}
static Vec p3_face_gram(const P3Mesh &m,const GfsEnrichConfig &cfg,const std::vector<Pair> &scale,const P3Face &f,int nq,Real *amplitude=nullptr){
    std::vector<double> q(nq),w(nq);gauss(nq,q.data(),w.data());Vec G(36);
    for(int j=0;j<nq;j++){
        double N[12],Nx[12],Ny[12];p3_face_shape(m,cfg,scale,f,0,q[j],N,Nx,Ny);Real phi[6];p3_face_functions(q[j],N,phi);
        for(int a=0;a<6;a++)for(int b=0;b<=a;b++){Real z=f.length*w[j]*phi[a]*phi[b];G[a*6+b]+=z;if(a!=b)G[b*6+a]+=z;}
    }
    if(amplitude){*amplitude=0;for(int side=0;side<2;side++)for(int a=0;a<10;a++)for(int k=0;k<2;k++)*amplitude=std::max(*amplitude,(Real)fabs(m.triangle[f.triangle[side]].v0[a][k]));}
    return G;
}
static std::array<Real,3> p3_residual_gram(const Vec &G){
    Vec P(16);for(int i=0;i<4;i++)for(int j=0;j<4;j++)P[i*4+j]=G[i*6+j];Factor factor(P,4);Vec coupling[2]={Vec(4),Vec(4)},projected[2];
    for(int k=0;k<2;k++){for(int i=0;i<4;i++)coupling[k][i]=G[i*6+4+k];projected[k]=factor.solve(coupling[k]);}
    std::array<Real,3> H={G[4*6+4],G[4*6+5],G[5*6+5]};for(int i=0;i<4;i++){H[0]-=coupling[0][i]*projected[0][i];H[1]-=coupling[0][i]*projected[1][i];H[2]-=coupling[1][i]*projected[1][i];}return H;
}
static void p3_trace_rank(const P3Mesh &m,const GfsEnrichConfig &cfg,const std::vector<Pair> &scale,P3Face &f){
    Real amplitude=0;auto coarse=p3_face_gram(m,cfg,scale,f,48),fine=p3_face_gram(m,cfg,scale,f,96,&amplitude);auto H=p3_residual_gram(fine),Hc=p3_residual_gram(coarse);
    Real uncertainty=0,norm=0;for(int i=0;i<3;i++)uncertainty=std::max(uncertainty,fabsl(H[i]-Hc[i]));for(auto a:fine)norm=std::max(norm,fabsl(a));f.quadrature_difference=uncertainty;
    int component=m.linear.tri[f.triangle[0]].component;Real noise=256*DBL_EPSILON*amplitude/std::min(scale[component][0],scale[component][1]);Real signal=sqrtl(std::max(fabsl(H[0]),fabsl(H[2]))/f.length);
    uncertainty=32*uncertainty+f.length*(2*signal*noise+noise*noise)+2048*LDBL_EPSILON*norm;f.rank_uncertainty=uncertainty;
    Real hi=.5L*(H[0]+H[2]+hypotl(H[0]-H[2],2*H[1])),lo=hi>0?(H[0]*H[2]-H[1]*H[1])/hi:0;
    need(hi>100*uncertainty,"P3 trace leading residual rank uncertain");
    if(lo>100*uncertainty)f.residual_rank=2;
    else{
        Real u[2]={H[1],hi-H[0]};if(hypotl(u[0],u[1])<32*LDBL_EPSILON*hi){u[0]=H[0]>=H[2]?1:0;u[1]=H[0]>=H[2]?0:1;}Real len=hypotl(u[0],u[1]);u[0]/=len;u[1]/=len;
        double q[96],w[96];gauss(96,q,w);Real omission=0;for(int j=0;j<96;j++){double N[12],Nx[12],Ny[12];p3_face_shape(m,cfg,scale,f,0,q[j],N,Nx,Ny);omission=std::max(omission,fabsl(-u[1]*N[10]+u[0]*N[11]));}
        f.omitted_function_max=omission;need(omission<=2*noise && fabsl(lo)<=uncertainty,"P3 trace small residual rank uncertain: no silent truncation");f.residual_rank=1;f.combination[0][0]=u[0];f.combination[0][1]=u[1];
    }
    int n=4+f.residual_rank;f.gram.resize(n*n);
    auto coeff=[&](int i,int a)->Real{if(i<4)return a==i?1:0;return a>=4?f.combination[i-4][a-4]:0;};
    for(int i=0;i<n;i++)for(int j=0;j<n;j++)for(int a=0;a<6;a++)for(int b=0;b<6;b++)f.gram[i*n+j]+=coeff(i,a)*fine[a*6+b]*coeff(j,b);
    Factor factor(f.gram,n);auto e=extrema(f.gram,n);need(e[0]>0,"P3 retained full trace Gram SPD");f.condition=e[1]/e[0];
}
static std::vector<P3Face> p3_physical_faces(const P3Mesh &m,const GfsEnrichConfig &cfg,const std::vector<Pair> &scale){
    std::vector<P3Face> out;
    for(const auto &entry:m.linear.edge){
        if(entry.second.size()!=2)continue;const auto &left=m.linear.tri[entry.second[0]],&right=m.linear.tri[entry.second[1]];if(left.cell==right.cell)continue;need(left.component==right.component,"P3 face connected-component ownership");P3Face f;
        f.triangle[0]=entry.second[0];f.triangle[1]=entry.second[1];f.cell[0]=left.cell;f.cell[1]=right.cell;f.endpoint[0]=entry.first.first;f.endpoint[1]=entry.first.second;
        auto a=m.linear.node[f.endpoint[0]],b=m.linear.node[f.endpoint[1]];f.length=hypot(b[0]-a[0],b[1]-a[1]);bool normal=false;
        for(int e=0;e<3;e++){int va=left.vertex[e],vb=left.vertex[(e+1)%3];if(std::min(va,vb)==f.endpoint[0] && std::max(va,vb)==f.endpoint[1]){auto x=m.linear.node[va],y=m.linear.node[vb];f.nx=(y[1]-x[1])/f.length;f.ny=-(y[0]-x[0])/f.length;normal=true;}}
        need(normal,"P3 physical outward normal");
        for(int side=0;side<2;side++){
            const auto &t=m.linear.tri[f.triangle[side]];int ia=-1,ib=-1;for(int j=0;j<3;j++){if(t.vertex[j]==f.endpoint[0])ia=j;if(t.vertex[j]==f.endpoint[1])ib=j;}need(ia>=0 && ib>=0,"P3 edge topology");
            for(int j=0;j<4;j++){
                int found=-1;for(int a=0;a<10;a++){int index[3];gfs_material_p3_index(a,index);if(index[ia]==3-j && index[ib]==j)found=m.connectivity[f.triangle[side]][a];}need(found>=0,"P3 edge interpolation node");if(side==0)f.node[j]=found;else need(f.node[j]==found,"P3 exact shared polynomial trace topology");
            }
        }
        p3_trace_rank(m,cfg,scale,f);out.push_back(f);
    }return out;
}
#endif
