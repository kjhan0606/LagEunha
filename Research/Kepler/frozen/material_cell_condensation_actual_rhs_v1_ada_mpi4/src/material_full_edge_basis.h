// Exact same canonical basis arithmetic; immutable64-point Gauss rule cached once.
static void cached_edge_basis(STEdge&e,shared_trace::Pt A,shared_trace::Pt B,double eps,const std::vector<double>&x,const std::vector<double>&w){
    e.eps=eps;for(int k=0;k<2;k++){e.A[k]=A[k];e.D[k]=B[k]-A[k];}
    for(int j=0;j<4;j++){double X=e.A[0]+(j/3.)*e.D[0],Y=e.A[1]+(j/3.)*e.D[1],Dq[4];need(gfs_material_p3_radial_remainder(e.A,X,Y,eps,e.qnode[j],Dq),"canonical edge original IC remainder");}
    // Exact algebraic rank, no sampled rank truncation. For a nonradial
    // straight edge the functions (a+b t)(r(t)^2+eps^2)^(-3/4) modulo P3
    // are independent. A radial edge has one scalar residual direction.
    __float128 det=__float128(A[0])*__float128(B[1])-__float128(A[1])*__float128(B[0]);e.rank=det==0?1:2;
    R G[2][2]={};for(int j=0;j<64;j++){STJet r[2];need(st_edge_raw(e,STJet(x[j]),r),"finite canonical residual");for(int a=0;a<2;a++)for(int b=0;b<2;b++)G[a][b]+=R(w[j])*r[a].v*r[b].v;}
    int p=G[0][0]>=G[1][1]?0:1,o=1-p;need(G[p][p]>0,"nonzero canonical residual norm");double s0=sqrtl(G[p][p]);e.encode[0][p]=1/s0;e.decode[p][0]=s0;
    if(e.rank==1){shared_trace::Pt V=fabs(A[0])+fabs(A[1])>fabs(B[0])+fabs(B[1])?A:B;double component[2]={-V[1],V[0]};need(component[p]!=0,"radial canonical pivot");e.decode[o][0]=s0*(component[o]/component[p]);}
    else {double beta=G[p][o]/G[p][p];R variance=0;for(int j=0;j<64;j++){STJet r[2];need(st_edge_raw(e,STJet(x[j]),r),"canonical second residual");R z=R(r[o].v)-R(beta)*r[p].v;variance+=R(w[j])*z*z;}need(variance>0,"independent canonical second residual norm");double s1=sqrtl(variance);e.encode[1][p]=-beta/s1;e.encode[1][o]=1/s1;e.decode[o][0]=beta*s0;e.decode[o][1]=s1;}
}
