/* Native fixed-reference hybrid INITIAL-RHS reference. Not a runtime solver. */
#define GFS_MATERIAL_SHARED_LIBRARY
#include "test_material_shared_v1.cpp"
#include "gfs_material_trace.h"
using Real=long double;
using Vec=std::vector<Real>;
struct Factor {
    int n=0;Vec L;Real minpivot=0,maxpivot=0;
    explicit Factor(const Vec &A,int size):n(size),L(n*n),minpivot(LDBL_MAX){
        for(int i=0;i<n;i++)for(int j=0;j<=i;j++){
            Real a=A[i*n+j];for(int k=0;k<j;k++)a-=L[i*n+k]*L[j*n+k];
            if(i==j){need(a>0 && std::isfinite(a),"hybrid SPD factor");minpivot=std::min(minpivot,a);maxpivot=std::max(maxpivot,a);L[i*n+j]=sqrtl(a);}
            else L[i*n+j]=a/L[j*n+j];
        }
    }
    Vec lower(Vec b)const{for(int i=0;i<n;i++){for(int j=0;j<i;j++)b[i]-=L[i*n+j]*b[j];b[i]/=L[i*n+i];}return b;}
    Vec solve(Vec b)const{b=lower(b);for(int i=n-1;i>=0;i--){for(int j=i+1;j<n;j++)b[i]-=L[j*n+i]*b[j];b[i]/=L[i*n+i];}return b;}
};
struct Groups {
    std::vector<int> p;
    explicit Groups(int n):p(n){std::iota(p.begin(),p.end(),0);}
    int root(int a){return p[a]==a?a:p[a]=root(p[a]);}
    bool join(int a,int b){a=root(a);b=root(b);if(a==b)return false;p[b]=a;return true;}
};
struct BrokenCell {
    std::vector<int> ids;int offset=0;Vec M,x,y,vx,vy,Fx,Fy,oldx,oldy;
    int local(int id)const{auto it=std::find(ids.begin(),ids.end(),id);need(it!=ids.end(),"cell trace ownership");return it-ids.begin();}
};
struct Constraint {int plus,minus,face,trace;};
struct PhysicalFace {
    int triangle[2],cell[2],endpoint[2],dof[2][4],row[4]={-1,-1,-1,-1};
    int active[4]={0,1,2,3},nt=4;double length,nx,ny;
};
static std::array<Real,2> extrema(Vec A,int n){
    // Native symmetric Jacobi reference eigensolve, not a production algorithm.
    Real norm=0;for(auto x:A)norm=std::max(norm,fabsl(x));
    for(int sweep=0;sweep<100;sweep++){
        Real off=0;
        for(int p=0;p<n;p++)for(int q=p+1;q<n;q++){
            Real apq=A[p*n+q];off=std::max(off,fabsl(apq));if(fabsl(apq)<=32*LDBL_EPSILON*norm)continue;
            Real tau=(A[q*n+q]-A[p*n+p])/(2*apq),t=copysignl(1.L,tau)/(fabsl(tau)+sqrtl(1+tau*tau));
            Real c=1/sqrtl(1+t*t),s=t*c,app=A[p*n+p],aqq=A[q*n+q];
            for(int k=0;k<n;k++)if(k!=p && k!=q){Real kp=A[k*n+p],kq=A[k*n+q];A[k*n+p]=A[p*n+k]=c*kp-s*kq;A[k*n+q]=A[q*n+k]=s*kp+c*kq;}
            A[p*n+p]=app-t*apq;A[q*n+q]=aqq+t*apq;A[p*n+q]=A[q*n+p]=0;
        }
        if(off<=64*LDBL_EPSILON*norm)break;need(sweep!=99,"Schur eigensolve convergence");
    }
    Real lo=LDBL_MAX,hi=0;for(int i=0;i<n;i++){lo=std::min(lo,A[i*n+i]);hi=std::max(hi,A[i*n+i]);}return {lo,hi};
}
static void face_shape(const SharedCaseView &s,const PhysicalFace &f,int side,double q,double N[5],double Nx[5],double Ny[5]){
    const auto &t=s.mesh.tri[f.triangle[side]];const auto &a=s.mesh.node[f.endpoint[0]],&b=s.mesh.node[f.endpoint[1]];
    double X=(1-q)*a[0]+q*b[0],Y=(1-q)*a[1]+q*b[1],v[2],B[4];velocity(s.cfg,X+12,Y+12,v,B);
    gfs_material_triangle_shape(&t.t,X,Y,v,B,s.scale[t.component].data(),true,N,Nx,Ny);
}
static Pair evaluate(const SharedCaseView &s,const PhysicalFace &f,int side,double q,const Vec &x,const Vec &y,const std::vector<BrokenCell> &cell){
    double N[5],Nx[5],Ny[5];face_shape(s,f,side,q,N,Nx,Ny);const auto &t=s.mesh.tri[f.triangle[side]];auto ids=indices(s.mesh,t,true);Pair v={0,0};const auto &c=cell[t.cell];
    for(int a=0;a<5;a++){int j=c.offset+c.local(ids[a]);v[0]+=(double)(N[a]*x[j]);v[1]+=(double)(N[a]*y[j]);}return v;
}
static void hybrid_case(const SharedCaseView &s){
    if(s.rank!=0)return;
    const auto &mesh=s.mesh;int oldn=s.x.size(),nc=mesh.cell.size(),nb=0;
    std::vector<BrokenCell> cell(nc);
    for(int c=0;c<nc;c++){
        auto &a=cell[c];std::set<int> ids;
        for(const auto &t:mesh.tri)if(t.cell==c){auto v=indices(mesh,t,true);ids.insert(v.begin(),v.end());}
        a.ids.assign(ids.begin(),ids.end());int n=a.ids.size();a.offset=nb;nb+=n;a.M.resize(n*n);
        for(int i:a.ids){a.x.push_back(s.x[i]);a.y.push_back(s.y[i]);a.vx.push_back(s.vx[i]);a.vy.push_back(s.vy[i]);a.oldx.push_back(s.old_ax[i]);a.oldy.push_back(s.old_ay[i]);a.Fx.push_back(s.cell[c].explicit_fx[i]);a.Fy.push_back(s.cell[c].explicit_fy[i]);}
        for(int i=0;i<n;i++)for(int j=0;j<n;j++)a.M[i*n+j]=s.cell[c].M[a.ids[i]*oldn+a.ids[j]];
    }
    Groups groups(nb);std::vector<PhysicalFace> face;std::vector<Constraint> rows,raw;
    for(const auto &edge:mesh.edge){
        if(edge.second.size()!=2)continue;const auto &t=mesh.tri[edge.second[0]],&u=mesh.tri[edge.second[1]];if(t.cell==u.cell)continue;
        PhysicalFace f={};f.triangle[0]=edge.second[0];f.triangle[1]=edge.second[1];f.cell[0]=t.cell;f.cell[1]=u.cell;f.endpoint[0]=edge.first.first;f.endpoint[1]=edge.first.second;
        const auto &a=mesh.node[f.endpoint[0]],&b=mesh.node[f.endpoint[1]];
        f.length=hypot(b[0]-a[0],b[1]-a[1]);double ex=t.t.x[2][0]-t.t.x[1][0],ey=t.t.x[2][1]-t.t.x[1][1];f.nx=ey/f.length;f.ny=-ex/f.length;
        if(a[0]==0 && b[0]==0){f.nt=3;f.active[2]=2;}
        else if(a[1]==0 && b[1]==0){f.nt=3;f.active[2]=3;}
        else need(a[0]*b[1]-a[1]*b[0]!=0,"non-axis radial face requires general trace rank implementation");
        for(int side=0;side<2;side++){
            const auto &tri=mesh.tri[f.triangle[side]];const auto &c=cell[tri.cell];
            int ids[4]={f.endpoint[0],f.endpoint[1],(int)mesh.node.size()+2*tri.component,(int)mesh.node.size()+2*tri.component+1};
            for(int k=0;k<4;k++)f.dof[side][k]=c.offset+c.local(ids[k]);
        }
        int fi=face.size();for(int j=0;j<f.nt;j++){
            int k=f.active[j];Constraint b={f.dof[0][k],f.dof[1][k],fi,k};raw.push_back(b);
            if(groups.join(b.plus,b.minus)){f.row[k]=rows.size();rows.push_back(b);}
        }
        face.push_back(f);
    }
    int nr=rows.size();std::set<int> roots;for(int i=0;i<nb;i++)roots.insert(groups.root(i));
    int nullity=roots.size();need(nullity==nb-nr,"independent physical trace rank");
    int old_inclusion=0,disconnected=0;for(const auto &b:raw){
        int owner[2]={-1,-1},gid[2]={-1,-1};for(int c=0;c<nc;c++)for(int k=0;k<(int)cell[c].ids.size();k++)for(int side=0;side<2;side++)if(cell[c].offset+k==(side?b.minus:b.plus)){owner[side]=c;gid[side]=cell[c].ids[k];}
        if(gid[0]!=gid[1])old_inclusion++;if(mesh.tri[face[b.face].triangle[0]].component!=mesh.tri[face[b.face].triangle[1]].component)disconnected++;
    }
    std::vector<Factor> factors;Vec Fx(nb),Fy(nb),x(nb),y(nb),vx(nb),vy(nb),oldx(nb),oldy(nb),bx(nb),by(nb);
    for(const auto &c:cell){
        int n=c.ids.size(),o=c.offset;factors.emplace_back(c.M,n);auto ax=factors.back().solve(c.Fx),ay=factors.back().solve(c.Fy);
        for(int i=0;i<n;i++){Fx[o+i]=c.Fx[i];Fy[o+i]=c.Fy[i];bx[o+i]=ax[i];by[o+i]=ay[i];x[o+i]=c.x[i];y[o+i]=c.y[i];vx[o+i]=c.vx[i];vy[o+i]=c.vy[i];oldx[o+i]=c.oldx[i];oldy[o+i]=c.oldy[i];}
    }
    // Build the Schur matrix as an explicit sum of local positive Gram forms.
    std::vector<Vec> lower(nb,Vec(nr)),minvbt(nb,Vec(nr));
    for(int c=0;c<nc;c++)for(int j=0;j<nr;j++){
        int n=cell[c].ids.size(),o=cell[c].offset;Vec b(n);if(rows[j].plus>=o && rows[j].plus<o+n)b[rows[j].plus-o]=1;if(rows[j].minus>=o && rows[j].minus<o+n)b[rows[j].minus-o]=-1;
        auto l=factors[c].lower(b),v=factors[c].solve(b);for(int i=0;i<n;i++){lower[o+i][j]=l[i];minvbt[o+i][j]=v[i];}
    }
    Vec S(nr*nr),rx(nr),ry(nr);
    for(int i=0;i<nr;i++){
        rx[i]=-bx[rows[i].plus]+bx[rows[i].minus];ry[i]=-by[rows[i].plus]+by[rows[i].minus];
        for(int j=0;j<=i;j++){Real a=0;for(int k=0;k<nb;k++)a+=lower[k][i]*lower[k][j];S[i*nr+j]=S[j*nr+i]=a;}
    }
    Factor schur(S,nr);auto eigen=extrema(S,nr);auto lx=schur.solve(rx),ly=schur.solve(ry);Vec ax=bx,ay=by,loadx(nb),loady(nb);
    for(int j=0;j<nr;j++){
        loadx[rows[j].plus]+=lx[j];loadx[rows[j].minus]-=lx[j];loady[rows[j].plus]+=ly[j];loady[rows[j].minus]-=ly[j];
        for(int k=0;k<nb;k++){ax[k]+=minvbt[k][j]*lx[j];ay[k]+=minvbt[k][j]*ly[j];}
    }
    Real constraint=0,initial_constraint=0,old_delta=0,force_scale=0,equation=0,old_null=0,bubble=0;
    Vec oldr_x(nb),oldr_y(nb),groupx(nb),groupy(nb);std::vector<int> groupn(nb);
    for(int c=0;c<nc;c++){
        int n=cell[c].ids.size(),o=cell[c].offset;
        for(int i=0;i<n;i++){
            Real ex=-Fx[o+i]-loadx[o+i],ey=-Fy[o+i]-loady[o+i],ox=-Fx[o+i],oy=-Fy[o+i];
            for(int j=0;j<n;j++){ex+=cell[c].M[i*n+j]*ax[o+j];ey+=cell[c].M[i*n+j]*ay[o+j];ox+=cell[c].M[i*n+j]*oldx[o+j];oy+=cell[c].M[i*n+j]*oldy[o+j];}
            equation=std::max({equation,fabsl(ex),fabsl(ey)});oldr_x[o+i]=ox;oldr_y[o+i]=oy;int g=groups.root(o+i);groupx[g]+=ox;groupy[g]+=oy;groupn[g]++;
            force_scale=std::max({force_scale,fabsl(Fx[o+i]),fabsl(Fy[o+i]),fabsl(loadx[o+i]),fabsl(loady[o+i])});old_delta=std::max({old_delta,fabsl(ax[o+i]-oldx[o+i]),fabsl(ay[o+i]-oldy[o+i])});
        }
    }
    for(int i=0;i<nb;i++)if(groupn[i])old_null=std::max({old_null,fabsl(groupx[i])/groupn[i],fabsl(groupy[i])/groupn[i]});
    for(const auto &b:raw){constraint=std::max({constraint,fabsl(ax[b.plus]-ax[b.minus]),fabsl(ay[b.plus]-ay[b.minus])});initial_constraint=std::max({initial_constraint,fabsl(x[b.plus]-x[b.minus]),fabsl(y[b.plus]-y[b.minus]),fabsl(vx[b.plus]-vx[b.minus]),fabsl(vy[b.plus]-vy[b.minus])});}
    // The cell-center P1 basis is a physical-boundary bubble.
    for(int c=0;c<nc;c++){
        int center=-1;for(const auto &t:mesh.tri)if(t.cell==c){center=t.vertex[0];break;}int k=cell[c].offset+cell[c].local(center);bubble=std::max({bubble,fabsl(loadx[k]),fabsl(loady[k]),fabsl(oldr_x[k]),fabsl(oldr_y[k])});
    }
    /* Direct physical-face quadrature and local ledgers are added below. */
    Vec directx(nb),directy(nb),work(nc),torque(nc),physical_work(nc),physical_torque(nc);
    Real gram_residual=0,gram_condition=0,traction_max=0,trace_v=0,trace_a=0,trace_x=0,pair_work=0,pair_torque=0,pair_momentum=0,geometry_measure=0,finite_trace=0,traction_difference_square=0,face_measure=0;
    Vec movedx=x,movedy=y;const Real h=1.e-5L;for(int k=0;k<nb;k++){movedx[k]+=h*vx[k]+.5L*h*h*ax[k];movedy[k]+=h*vy[k]+.5L*h*h*ay[k];}
    std::string sample_path=std::string(s.dest)+"/"+s.name+"_traction.csv";FILE *samples=fopen(sample_path.c_str(),"w");need(samples,"new traction samples");fprintf(samples,"face,q,x,y,tx_reference,ty_reference,tx_pressure,ty_pressure,physical_length_ratio\n");
    for(size_t fi=0;fi<face.size();fi++){
        const auto &f=face[fi];double G[4][4]={{0}},L[4][4]={{0}},mux[4]={0},muy[4]={0},eta_x[4]={0},eta_y[4]={0};
        for(int j=0;j<f.nt;j++){int k=f.active[j];if(f.row[k]>=0){mux[j]=(double)lx[f.row[k]];muy[j]=(double)ly[f.row[k]];}eta_x[j]=mux[j];eta_y[j]=muy[j];}
        for(int q=0;q<s.cfg.na;q++){
            double z=s.cfg.an[q],N[5],Nx[5],Ny[5];face_shape(s,f,0,z,N,Nx,Ny);double phi[4]={1-z,z,N[3],N[4]},w=f.length*s.cfg.aw[q];
            for(int i=0;i<f.nt;i++)for(int j=0;j<f.nt;j++)G[i][j]+=w*phi[f.active[i]]*phi[f.active[j]];
        }
        need(gfs_material_trace_factor(f.nt,G,L),"actual independent physical trace Gram positive");gfs_material_trace_solve(f.nt,L,eta_x);gfs_material_trace_solve(f.nt,L,eta_y);
        Vec gv(f.nt*f.nt);for(int i=0;i<f.nt;i++)for(int j=0;j<f.nt;j++)gv[i*f.nt+j]=G[i][j];auto ge=extrema(gv,f.nt);need(ge[0]>0,"trace Gram spectrum positive");gram_condition=std::max(gram_condition,ge[1]/ge[0]);
        for(int i=0;i<f.nt;i++){
            Real ex=-mux[i],ey=-muy[i];for(int j=0;j<f.nt;j++){ex+=G[i][j]*eta_x[j];ey+=G[i][j]*eta_y[j];}gram_residual=std::max({gram_residual,fabsl(ex),fabsl(ey)});
        }
        for(int q=0;q<s.cfg.na;q++){
            double z=s.cfg.an[q],N[5],Nx[5],Ny[5];face_shape(s,f,0,z,N,Nx,Ny);double phi[4]={1-z,z,N[3],N[4]},w=f.length*s.cfg.aw[q],tx=0,ty=0;
            for(int j=0;j<f.nt;j++){tx+=phi[f.active[j]]*eta_x[j];ty+=phi[f.active[j]]*eta_y[j];}
            traction_max=std::max(traction_max,(Real)hypot(tx,ty));auto vl=evaluate(s,f,0,z,vx,vy,cell),vr=evaluate(s,f,1,z,vx,vy,cell),al=evaluate(s,f,0,z,ax,ay,cell),ar=evaluate(s,f,1,z,ax,ay,cell),xl=evaluate(s,f,0,z,x,y,cell),xr=evaluate(s,f,1,z,x,y,cell),ml=evaluate(s,f,0,z,movedx,movedy,cell),mr=evaluate(s,f,1,z,movedx,movedy,cell);
            for(int k=0;k<2;k++){trace_v=std::max(trace_v,(Real)fabs(vl[k]-vr[k]));trace_a=std::max(trace_a,(Real)fabs(al[k]-ar[k]));trace_x=std::max(trace_x,(Real)fabs(xl[k]-xr[k]));finite_trace=std::max(finite_trace,(Real)fabs(ml[k]-mr[k]));}
            pair_work+=w*((vl[0]-vr[0])*tx+(vl[1]-vr[1])*ty);pair_torque+=w*((xl[0]-xr[0])*ty-(xl[1]-xr[1])*tx);
            Pair physical[2];double face_impulse[2]={0,0};
            traction_difference_square+=w*(pow(tx+s.cfg.field.ic.pressure*f.nx,2)+pow(ty+s.cfg.field.ic.pressure*f.ny,2));face_measure+=w;
            for(int side=0;side<2;side++){
                const auto &t=mesh.tri[f.triangle[side]];auto ids=indices(mesh,t,true);const auto &c=cell[t.cell];double ns[5],dx[5],dy[5];face_shape(s,f,side,z,ns,dx,dy);int sign=side?-1:1;
                for(int a=0;a<5;a++){int k=c.offset+c.local(ids[a]);directx[k]+=sign*w*ns[a]*tx;directy[k]+=sign*w*ns[a]*ty;}
                const auto &v=side?vr:vl,&p=side?xr:xl;work[t.cell]+=sign*w*(v[0]*tx+v[1]*ty);torque[t.cell]+=sign*w*(p[0]*ty-p[1]*tx);
                double ptx=-s.cfg.field.ic.pressure*f.nx,pty=-s.cfg.field.ic.pressure*f.ny;physical_work[t.cell]+=sign*w*(v[0]*ptx+v[1]*pty);physical_torque[t.cell]+=sign*w*(p[0]*pty-p[1]*ptx);
                // Independently evaluate the deformed tangent and measure ratio.
                double ex=mesh.node[f.endpoint[1]][0]-mesh.node[f.endpoint[0]][0],ey=mesh.node[f.endpoint[1]][1]-mesh.node[f.endpoint[0]][1],tangent[2]={0,0};
                for(int a=0;a<5;a++){int k=c.offset+c.local(ids[a]);double d=dx[a]*ex+dy[a]*ey;tangent[0]+=(double)(movedx[k]*d);tangent[1]+=(double)(movedy[k]*d);}
                double ratio=hypot(tangent[0],tangent[1])/f.length;need(ratio>0 && std::isfinite(ratio),"finite material face length positive");physical[side]={(double)(sign*tx/ratio),(double)(sign*ty/ratio)};
                double partition=ns[0]+ns[1]+ns[2];for(int k=0;k<2;k++)face_impulse[k]+=w*ratio*physical[side][k]*partition;
                geometry_measure=std::max({geometry_measure,(Real)fabs(physical[side][0]*ratio-sign*tx),(Real)fabs(physical[side][1]*ratio-sign*ty)});
                if(side==0)fprintf(samples,"%zu,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",fi,q,xl[0],xl[1],tx,ty,ptx,pty,ratio);
            }
            pair_momentum=std::max(pair_momentum,(Real)hypot(face_impulse[0],face_impulse[1]));
        }
    }
    fclose(samples);
    Real direct_error=0;for(int i=0;i<nb;i++)direct_error=std::max({direct_error,fabsl(directx[i]-loadx[i]),fabsl(directy[i]-loady[i])});
    // Independent point reconstruction against the actual initial Euler source.
    // P is spatially constant here, so the exact initial acceleration is g.
    Vec mass(nc),source_work(nc),source_torque(nc),error_square(nc),gravity_square(nc),error_max(nc),volume(nc),curvature_error_L1(nc),curvature_error_max(nc),curvature_exact_max(nc),curvature_numerical_max(nc);
    Vec curvature_fd_error[3]={Vec(nc),Vec(nc),Vec(nc)};const Real curvature_fd_step[3]={2.e-5L,1.e-5L,5.e-6L};
    Real finite_min_J=LDBL_MAX;
    for(const auto &t:mesh.tri){
        int c=t.cell;const auto &bc=cell[c];double px[3],py[3];for(int a=0;a<3;a++){px[a]=t.t.x[a][0];py[a]=t.t.x[a][1];}
        GfsMomentPolygon poly={px,py,3,0};GfsEnrichBasis qb={};qb.generator[0]=qb.generator[1]=12;GfsEnrichQuadrature qr;
        need(gfs_enrich_quadrature_init(poly,&s.cfg,&qb,&qr),"independent hybrid source quadrature");auto ids=indices(mesh,t,true);
        for(int q=0;q<gfs_enrich_nquad(poly,&s.cfg,&qr);q++){
            double X,Y,w,N[5],Nx[5],Ny[5],iv[2],B[4],u[4];need(gfs_enrich_quad(poly,&s.cfg,&qr,q,&X,&Y,&w),"source assessment point");if(w==0)continue;
            velocity(s.cfg,X+12,Y+12,iv,B);gfs_material_triangle_shape(&t.t,X,Y,iv,B,s.scale[t.component].data(),true,N,Nx,Ny);need(gfs_enrich_reference(&s.cfg,X+12,Y+12,u,0,0),"source assessment rho");
            Real aa[2]={0,0},vv[2]={0,0},diva=0,Fzero[4]={0},Vgrad[4]={0},Agrad[4]={0};GfsMaterialMap moved={};moved.n=5;
            for(int a=0;a<5;a++){
                int k=bc.offset+bc.local(ids[a]);aa[0]+=N[a]*ax[k];aa[1]+=N[a]*ay[k];diva+=Nx[a]*ax[k]+Ny[a]*ay[k];vv[0]+=N[a]*vx[k];vv[1]+=N[a]*vy[k];moved.x[0][a]=(double)movedx[k];moved.x[1][a]=(double)movedy[k];
                Real coords[2]={x[k]-(a<3?t.t.x[0][0]:0),y[k]-(a<3?t.t.x[0][1]:0)},vel[2]={vx[k],vy[k]},acc[2]={ax[k],ay[k]};
                for(int d=0;d<2;d++){Fzero[2*d]+=coords[d]*Nx[a];Fzero[2*d+1]+=coords[d]*Ny[a];Vgrad[2*d]+=vel[d]*Nx[a];Vgrad[2*d+1]+=vel[d]*Ny[a];Agrad[2*d]+=acc[d]*Nx[a];Agrad[2*d+1]+=acc[d]*Ny[a];}
            }
            GfsMaterialPoint mp;need(gfs_material_point(&moved,N,Nx,Ny,&mp),"finite-map positive volume Jacobian");finite_min_J=std::min(finite_min_J,(Real)mp.J);
            Real gx=0,gy=0,divg=0;if(s.source_mode==1){Real r2=(Real)X*X+(Real)Y*Y,z=r2+(Real)s.cfg.field.ic.eps*s.cfg.field.ic.eps,f=powl(z,-1.5L);gx=-X*f;gy=-Y*f;divg=-2*f+3*r2*f/z;}if(s.source_mode==2){gx=-X;gy=-Y;divg=-2;}
            Real dx=aa[0]-gx,dy=aa[1]-gy,dm=(Real)w*u[0];mass[c]+=dm;error_square[c]+=dm*(dx*dx+dy*dy);gravity_square[c]+=dm*(gx*gx+gy*gy);error_max[c]=std::max(error_max[c],sqrtl(dx*dx+dy*dy));source_work[c]+=dm*(vv[0]*dx+vv[1]*dy);source_torque[c]+=dm*(X*dy-Y*dx);
            // At initial div(v)=0: (D^2 P/Dt^2)/P=-gamma*(div(a)-tr(B^2)).
            // This is an operator Taylor coefficient, not an evolved pressure.
            Real trB2=(Real)B[0]*B[0]+2.L*B[1]*B[2]+(Real)B[3]*B[3],exact=-s.cfg.field.ic.gamma*(divg-trB2),numerical=-s.cfg.field.ic.gamma*(diva-trB2),ce=numerical-exact;
            volume[c]+=w;curvature_error_L1[c]+=w*fabsl(ce);curvature_error_max[c]=std::max(curvature_error_max[c],fabsl(ce));curvature_exact_max[c]=std::max(curvature_exact_max[c],fabsl(exact));curvature_numerical_max[c]=std::max(curvature_numerical_max[c],fabsl(numerical));
            // Independent finite-map central difference. Stable determinant
            // increments avoid subtracting nearly equal unit Jacobians.
            Real Jzero=Fzero[0]*Fzero[3]-Fzero[1]*Fzero[2];need(Jzero>0,"initial independent map determinant");
            for(int ih=0;ih<3;ih++){
                Real hq=curvature_fd_step[ih],fd=0;
                for(int sign:{-1,1}){
                    Real D[4];for(int j=0;j<4;j++)D[j]=sign*hq*Vgrad[j]+.5L*hq*hq*Agrad[j];
                    Real dJ=(Fzero[3]*D[0]+Fzero[0]*D[3]-Fzero[2]*D[1]-Fzero[1]*D[2]+D[0]*D[3]-D[1]*D[2])/Jzero;need(dJ>-1,"two-sided finite material map determinant positive");
                    fd+=expm1l(-s.cfg.field.ic.gamma*log1pl(dJ));
                }
                fd/=hq*hq;curvature_fd_error[ih][c]=std::max(curvature_fd_error[ih][c],fabsl(fd-numerical));
            }
        }
    }
    Real maxE=0,maxL=0,maxPhysE=0,maxPhysL=0,totalE=0,totalL=0,totalU=0,totalK=0,maxPdot=0,maxEidentity=0,maxLidentity=0;
    printf("{\"case\":\"%s\",\"ranks\":%d,\"cells_count\":%d,\"broken_scalar_dofs\":%d,\"old_scalar_dofs\":%d,\"physical_faces\":%zu,\"raw_trace_rows\":%zu,\"trace_rank\":%d,\"conforming_nullity\":%d,\"extra_conforming_modes\":%d,\"old_space_not_in_kernel_rows\":%d,\"disconnected_couplings\":%d,\"Schur_min_eigenvalue\":%.17Lg,\"Schur_condition\":%.17Lg,\"Schur_min_pivot\":%.17Lg,\"trace_Gram_condition_max\":%.17Lg,\"trace_Gram_solve_absolute\":%.17Lg,\"initial_coefficient_jump\":%.17Lg,\"acceleration_constraint\":%.17Lg,\"acceleration_vs_old_max\":%.17Lg,\"old_residual_kernel_projection_relative\":%.17Lg,\"bubble_load_relative\":%.17Lg,\"momentum_equation_relative\":%.17Lg,\"direct_traction_load_relative\":%.17Lg,\"traction_magnitude_over_P_max\":%.17Lg,\"position_trace_jump\":%.17Lg,\"velocity_trace_jump\":%.17Lg,\"acceleration_trace_jump\":%.17Lg,\"finite_map_position_jump\":%.17Lg,\"finite_map_measure_error\":%.17Lg,\"pair_momentum_absolute\":%.17Lg,\"pair_work_absolute\":%.17Lg,\"pair_torque_absolute\":%.17Lg,\"cells\":[",
        s.name,s.ranks,nc,nb,oldn,face.size(),raw.size(),nr,nullity,nullity-oldn,old_inclusion,disconnected,eigen[0],eigen[1]/eigen[0],schur.minpivot,gram_condition,gram_residual,initial_constraint,constraint,old_delta,old_null/force_scale,bubble/force_scale,equation/force_scale,direct_error/force_scale,traction_max/s.cfg.field.ic.pressure,trace_x,trace_v,trace_a,finite_trace,geometry_measure,pair_momentum,pair_work,pair_torque);
    for(int c=0;c<nc;c++){
        const auto &a=cell[c];const auto &ca=s.cell[c];int n=a.ids.size(),o=a.offset;Real K=0,kdot=0,Ldot=0;
        for(int i=0;i<n;i++)for(int j=0;j<n;j++){Real m=a.M[i*n+j];K+=.5L*m*(vx[o+i]*vx[o+j]+vy[o+i]*vy[o+j]);kdot+=m*(vx[o+i]*ax[o+j]+vy[o+i]*ay[o+j]);Ldot+=m*(x[o+i]*ay[o+j]-y[o+i]*ax[o+j]+vx[o+i]*vy[o+j]-vy[o+i]*vx[o+j]);}
        Real E=kdot+ca.Udot-work[c]-ca.external_work-ca.gpower,Lr=Ldot-torque[c]-ca.external_torque-ca.gtorque;
        Real PE=kdot+ca.Udot-ca.work-ca.gpower,PL=Ldot-ca.torque-ca.gtorque;
        Real dwork=work[c]+ca.external_work-ca.work,dtorque=torque[c]+ca.external_torque-ca.torque;
        maxE=std::max(maxE,fabsl(E)/ca.U);maxL=std::max(maxL,fabsl(Lr)/(2*K));maxPhysE=std::max(maxPhysE,fabsl(PE)/ca.U);maxPhysL=std::max(maxPhysL,fabsl(PL)/(2*K));totalE+=E;totalL+=Lr;totalK+=K;totalU+=ca.U;maxPdot=std::max(maxPdot,(Real)(s.cfg.field.ic.gamma*ca.divmax));maxEidentity=std::max(maxEidentity,fabsl(PE-E-dwork)/ca.U);maxLidentity=std::max(maxLidentity,fabsl(PL-Lr-dtorque)/(2*K));
        if(c)printf(",");printf("{\"id\":%llu,\"U\":%.17g,\"K\":%.17Lg,\"numerical_energy_residual_over_U\":%.17Lg,\"numerical_angular_residual_over_2K\":%.17Lg,\"physical_energy_residual_over_U\":%.17Lg,\"physical_angular_residual_over_2K\":%.17Lg,\"traction_difference_work_over_U\":%.17Lg,\"traction_difference_torque_over_2K\":%.17Lg,\"acceleration_error_RMS\":%.17Lg,\"gravity_RMS\":%.17Lg,\"acceleration_error_max\":%.17Lg,\"independent_source_projection_work_over_U\":%.17Lg,\"independent_source_projection_torque_over_2K\":%.17Lg,\"pressure_material_curvature_error_L1_volume\":%.17Lg,\"pressure_material_curvature_error_max\":%.17Lg,\"pressure_material_curvature_exact_max\":%.17Lg,\"pressure_material_curvature_numerical_max\":%.17Lg,\"pressure_curvature_finite_map_error_max\":[%.17Lg,%.17Lg,%.17Lg]}",(unsigned long long)mesh.cell[c].id,ca.U,K,E/ca.U,Lr/(2*K),PE/ca.U,PL/(2*K),dwork/ca.U,dtorque/(2*K),sqrtl(error_square[c]/mass[c]),sqrtl(gravity_square[c]/mass[c]),error_max[c],source_work[c]/ca.U,source_torque[c]/(2*K),curvature_error_L1[c]/volume[c],curvature_error_max[c],curvature_exact_max[c],curvature_numerical_max[c],curvature_fd_error[0][c],curvature_fd_error[1][c],curvature_fd_error[2][c]);
    }
    printf("],\"finite_map_min_volume_J\":%.17Lg,\"max_numerical_energy_residual_over_U\":%.17Lg,\"max_numerical_angular_residual_over_2K\":%.17Lg,\"max_physical_energy_residual_over_U\":%.17Lg,\"max_physical_angular_residual_over_2K\":%.17Lg,\"global_numerical_energy_residual_over_U\":%.17Lg,\"global_numerical_angular_residual_over_2K\":%.17Lg,\"pressure_derivative_over_P_max\":%.17Lg,\"traction_difference_L2_over_P\":%.17Lg,\"energy_traction_difference_identity_over_U\":%.17Lg,\"angular_traction_difference_identity_over_2K\":%.17Lg}",finite_min_J,maxE,maxL,maxPhysE,maxPhysL,totalE/totalU,totalL/(2*totalK),maxPdot,sqrtl(traction_difference_square/face_measure)/s.cfg.field.ic.pressure,maxEidentity,maxLidentity);
    std::string state_path=std::string(s.dest)+"/"+s.name+"_state.bin";FILE *state=fopen(state_path.c_str(),"wb");need(state,"new broken hybrid state");uint64_t header[3]={0x4746485942523031ULL,(uint64_t)nb,(uint64_t)nr};fwrite(header,sizeof(header),1,state);
    for(const auto *v:{&x,&y,&vx,&vy,&ax,&ay,&lx,&ly})for(auto a:*v){double d=(double)a;need(fwrite(&d,8,1,state)==1,"hybrid state export");}fclose(state);fflush(stdout);
}
#ifndef GFS_MATERIAL_HYBRID_LIBRARY
int main(int argc,char **argv){
    MPI_Init(&argc,&argv);int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);need(argc==3 || argc==4,"polygon packet, new output, optional fine_quadrature");
    bool fine=argc==4;if(fine)need(std::string(argv[3])=="fine_quadrature","known quadrature option");shared_hook=hybrid_case;
    std::set<uint64_t> ids={17894,20502,39878,44945,67927,67934,67935,67942,67943,67992,67993,68000,68001,68002,68525,68526,68527,68534,68535,68584,68585,68592,68593,68600};auto mesh=mesh_from(load(argv[1],ids));
    if(rank==0)printf("{\"scope\":\"independently solved native material hybrid initial RHS; actual physical trace constraints; dual numerical traction is reported separately from physical pressure traction; no evolution or shock claim\",\"quadrature\":[%d,%d],\"cases\":[",fine?64:32,fine?28:20);
    for(int test=0;test<4;test++){
        auto cfg=config(fine?64:32,fine?28:20);const char *name=test==0?"gravity_off":test==1?"linear_gravity":test==2?"Kepler":"Kepler_perturb";int source=test==0?0:test==1?2:1;
        if(rank==0 && test)printf(",\n");run_case(mesh,cfg,name,test==3?.013:0,test==3?-.007:0,rank,ranks,argv[2],source);
    }
    if(rank==0)printf("]}\n");MPI_Finalize();return 0;
}
#endif
