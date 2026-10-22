#ifndef MATERIAL_SHARED_REAL_BOUND_H
#define MATERIAL_SHARED_REAL_BOUND_H
#include "material_shared_function.h"
namespace shared_function {
struct CubicBound {Q value=0,gradient[2]={};};
static CubicBound cubic_bound(const Map&m,const S*nodal){
    // Exact P3 nodal->Bernstein conversion. Convex-hull bounds of the
    // derivative control net preserve affine cancellation algebraically.
    int at[4][4]={};for(int a=0;a<10;a++){int index[3];gfs_material_p3_index(a,index);at[index[1]][index[2]]=a;}S b[4][4]={};int vertex[3]={at[0][0],at[3][0],at[0][3]};for(int k=0;k<3;k++)b[k==1?3:0][k==2?3:0]=nodal[vertex[k]];
    for(int k=0;k<3;k++){int j=(k+1)%3,index1[3]={},index2[3]={};index1[k]=2;index1[j]=1;index2[k]=1;index2[j]=2;int p=at[index1[1]][index1[2]],q=at[index2[1]][index2[2]];b[index1[1]][index1[2]]=S(3)*nodal[p]-S(Q(1.5L))*nodal[q]-S(5)*nodal[vertex[k]]/S(6)+nodal[vertex[j]]/S(3);b[index2[1]][index2[2]]=S(3)*nodal[q]-S(Q(1.5L))*nodal[p]+nodal[vertex[k]]/S(3)-S(5)*nodal[vertex[j]]/S(6);}
    S boundary;for(int i=0;i<=3;i++)for(int j=0;i+j<=3;j++)if((i||j)&&i!=3&&j!=3&&!(i==1&&j==1))boundary=boundary+b[i][j];b[1][1]=(S(27)*nodal[at[1][1]]-nodal[vertex[0]]-nodal[vertex[1]]-nodal[vertex[2]]-S(3)*boundary)/S(6);
    CubicBound out;for(int i=0;i<=3;i++)for(int j=0;i+j<=3;j++)out.value=std::max(out.value,mag(b[i][j]));for(int i=0;i<=2;i++)for(int j=0;i+j<=2;j++)for(int k=0;k<2;k++){S derivative=S(3)*((b[i+1][j]-b[i][j])*m.grad[1][k]+(b[i][j+1]-b[i][j])*m.grad[2][k]);out.gradient[k]=std::max(out.gradient[k],mag(derivative));}return out;
}
struct RemainderBound {Q value[2]={},gradient[2][2]={},series_domain=0;};
static Q quotient(Q numerator,Q denominator){need(denominator>0,"positive analytic Taylor bound denominator");return mf::FP<Q>::up(numerator/denominator);}
static RemainderBound taylor_bound(const S*anchor,const Q*delta,double eps){
    RemainderBound out;S z=anchor[0]*anchor[0]+anchor[1]*anchor[1]+exact(eps)*exact(eps);Q zl=lower(z),L=quotient(mul(Q(2),add(mul(mag(anchor[0]),delta[0]),mul(mag(anchor[1]),delta[1]))),zl),W=quotient(add(mul(delta[0],delta[0]),mul(delta[1],delta[1])),zl),s=add(L,W);need(s<Q(.75L),"analytic Taylor remainder whole real source domain");out.series_domain=s;Q d=mf::FP<Q>::down(Q(1)-s),s2=mul(s,s),s3=mul(s2,s),s4=mul(s3,s),tail=quotient(mul(Q(.56396484375L),s4),d),dtail=mul(mul(Q(.56396484375L),s3),add(quotient(Q(4),d),quotient(s,mf::FP<Q>::down(d*d))));Q c2=Q(.65625L),c3=Q(.6015625L),L2=mul(L,L),W2=mul(W,W);
    Q A=add(add(mul(c2,W2),mul(c3,add(add(mul(mul(Q(3),L2),W),mul(mul(Q(3),L),W2)),mul(W2,W)))),tail),B=add(add(mul(c2,add(mul(mul(Q(2),L),W),W2)),mul(c3,s3)),tail),AL=add(mul(c3,add(mul(mul(Q(6),L),W),mul(Q(3),W2))),dtail),AQ=add(add(mul(mul(Q(2),c2),W),mul(c3,add(add(mul(Q(3),L2),mul(mul(Q(6),L),W)),mul(Q(3),W2)))),dtail),BL=add(add(mul(mul(Q(2),c2),W),mul(mul(Q(3),c3),s2)),dtail),BQ=add(add(mul(mul(Q(2),c2),add(L,W)),mul(mul(Q(3),c3),s2)),dtail),omega=mag(mf::power(S(zl),Q(-.75L)));
    for(int c=0;c<2;c++){Q ja=mag(anchor[1-c]),jd=delta[1-c];out.value[c]=mul(omega,add(mul(ja,A),mul(jd,B)));for(int j=0;j<2;j++){Q ld=quotient(mul(Q(2),mag(anchor[j])),zl),qd=quotient(mul(Q(2),delta[j]),zl),dA=add(mul(AL,ld),mul(AQ,qd)),dB=add(mul(BL,ld),mul(BQ,qd));out.gradient[c][j]=mul(omega,add(add(mul(ja,dA),mul(jd,dB)),c!=j?B:Q(0)));}}return out;
}
static RealBounds tight_real_bounds(const Map&m){
    // SAME canonical function and stored coefficients as the frozen full
    // field proof. Bounds retain v0-I3(v0) as a small Taylor residual rather
    // than separating two large terms. No fit or guard/order/precision scan.
    RealBounds out;const auto&t=m.state->triangle[m.source];Q delta[2]={};for(int a=0;a<3;a++)for(int k=0;k<2;k++)delta[k]=std::max(delta[k],mag(m.vertex[a][k]-m.vertex[0][k]));auto q=taylor_bound(m.vertex[0],delta,m.eps);S qnode[2][10];for(int a=0;a<10;a++){auto X=reference_node(m,a);auto values=remainder(t.vertex[0],X[0],X[1],m.eps);for(int c=0;c<2;c++)qnode[c][a]=values[c];}Q Rvalue[2],Rgradient[2][2];for(int c=0;c<2;c++){auto interpolation=cubic_bound(m,qnode[c]);Rvalue[c]=add(q.value[c],interpolation.value);for(int j=0;j<2;j++)Rgradient[c][j]=add(q.gradient[c][j],interpolation.gradient[j]);}
    for(int k=0;k<2;k++){S p[10];for(int a=0;a<10;a++){p[a]=exact(m.state->coefficient[t.id[a]].x[k])-reference_node(m,a)[k];if(a==9)for(int c=0;c<2;c++)p[a]=p[a]-m.bulk[0][k][c]*exact(t.bubble_shift[c]);}auto polynomial=cubic_bound(m,p);out.displacement[k]=polynomial.value;for(int j=0;j<2;j++)out.FminusI[2*k+j]=polynomial.gradient[j];for(int c=0;c<2;c++){Q coefficient=mag(m.bulk[0][k][c]);out.displacement[k]=add(out.displacement[k],mul(coefficient,Rvalue[c]));for(int j=0;j<2;j++)out.FminusI[2*k+j]=add(out.FminusI[2*k+j],mul(coefficient,Rgradient[c][j]));}}
    for(int e=0;e<3;e++){const auto&edge=m.edge[e];Q change[2]={mag(edge.Dref[0]),mag(edge.Dref[1])};auto qedge=taylor_bound(edge.A,change,m.eps);D parameter(S(Q(.5L),Q(.5L)),S(1),S(0));S nodal[4][2];for(int a=0;a<4;a++){auto value=remainder(edge.native.A,edge.A[0]+S(a)*edge.Dref[0]/S(3),edge.A[1]+S(a)*edge.Dref[1]/S(3),m.eps);for(int c=0;c<2;c++)nodal[a][c]=value[c];}Q residual[2]={qedge.value[0],qedge.value[1]},derivative[2]={};for(int c=0;c<2;c++){for(int j=0;j<2;j++)derivative[c]=add(derivative[c],mul(qedge.gradient[c][j],change[j]));for(int a=0;a<4;a++){auto L=lagrange(a,parameter);residual[c]=add(residual[c],mul(mag(L.v),mag(nodal[a][c])));derivative[c]=add(derivative[c],mul(mag(L.x),mag(nodal[a][c])));}}
        for(int k=0;k<2;k++)for(int c=0;c<2;c++){Q coefficient=mag(m.extension[0][k][e][c]);out.displacement[k]=add(out.displacement[k],mul(coefficient,residual[c]));for(int j=0;j<2;j++){Q ds=mag(m.grad[t.endpoint[e][0]][j]+m.grad[t.endpoint[e][1]][j]),b=mag(m.grad[t.endpoint[e][1]][j]),gradient=add(mul(mul(Q(4),ds),residual[c]),mul(add(b,ds),derivative[c]));out.FminusI[2*k+j]=add(out.FminusI[2*k+j],mul(coefficient,gradient));}}
    }for(int k=0;k<2;k++)out.L=std::max(out.L,add(out.FminusI[2*k],out.FminusI[2*k+1]));return out;
}
}
#endif
