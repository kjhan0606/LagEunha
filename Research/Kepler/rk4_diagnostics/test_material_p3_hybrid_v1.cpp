/* Fixed P3 plus one h/2 level. Native INITIAL-RHS accuracy reference only. */
#define GFS_MATERIAL_HYBRID_LIBRARY
#include "test_material_hybrid_v1.cpp"
#include "material_p3_mesh.h"
#include "material_p3_trace.h"
#include <chrono>
#include <sys/resource.h>
struct P3Point {
    int triangle;double X,Y,w,rho,N[12],Nx[12],Ny[12],v[2],B[4];
};
struct P3CellData {
    int n=0,offset=0,component=-1;std::vector<int> node;std::vector<double> raw;
    std::vector<double> force[3][2];Vec M;
    double mass=0,Kdirect=0,U=0,Udot=0,external_work=0,external_torque=0,external_force[2]={0},physical_work=0,physical_torque=0;
    int local(int id)const{auto i=std::lower_bound(node.begin(),node.end(),id);need(i!=node.end() && *i==id,"P3 physical-cell node ownership");return i-node.begin();}
};
struct P3Assembly {
    std::vector<P3CellData> cell;std::vector<Pair> scale;std::vector<P3Point> point,owned_point;
    int dofs=0;Real min_mass_pivot=LDBL_MAX,max_mass_condition=0;double seconds=0;
};
static void p3_local_ids(const P3Mesh &m,const std::vector<P3CellData> &cell,int triangle,int out[12]){
    const auto &c=cell[m.linear.tri[triangle].cell];for(int a=0;a<10;a++)out[a]=c.local(m.connectivity[triangle][a]);out[10]=c.n-2;out[11]=c.n-1;
}
static void p3_source(int mode,Real X,Real Y,Real eps,Real g[2],Real D[4]){
    for(int k=0;k<2;k++)g[k]=0;for(int k=0;k<4;k++)D[k]=0;
    if(mode==1){g[0]=-X;g[1]=-Y;D[0]=D[3]=-1;}
    if(mode==2){Real z=X*X+Y*Y+eps*eps,f=powl(z,-1.5L);g[0]=-X*f;g[1]=-Y*f;D[0]=-f+3*X*X*f/z;D[1]=D[2]=3*X*Y*f/z;D[3]=-f+3*Y*Y*f/z;}
}
static P3Assembly p3_assemble(const P3Mesh &m,const GfsEnrichConfig &cfg,int rank,int ranks){
    auto start=std::chrono::steady_clock::now();P3Assembly result;int nc=m.linear.cell.size();result.cell.resize(nc);
    for(int c=0;c<nc;c++){
        auto &a=result.cell[c];a.node=m.cell_nodes[c];a.n=a.node.size()+2;a.offset=result.dofs;result.dofs+=a.n;a.raw.resize(a.n*a.n);for(auto &mode:a.force)for(auto &v:mode)v.resize(a.n);
    }
    std::vector<P3Point> local_points;
    for(size_t it=0;it<m.linear.tri.size();it++){
        const auto &t=m.linear.tri[it];auto &c=result.cell[t.cell];c.component=t.component;if(t.cell%ranks!=rank)continue;int ids[12];p3_local_ids(m,result.cell,it,ids);
        double px[3],py[3];for(int a=0;a<3;a++){px[a]=t.t.x[a][0];py[a]=t.t.x[a][1];}GfsMomentPolygon poly={px,py,3,0};GfsEnrichBasis qb={};qb.generator[0]=qb.generator[1]=12;GfsEnrichQuadrature qr;need(gfs_enrich_quadrature_init(poly,&cfg,&qb,&qr),"P3 operator positive volume quadrature");
        for(int q=0;q<gfs_enrich_nquad(poly,&cfg,&qr);q++){
            P3Point p={};p.triangle=it;need(gfs_enrich_quad(poly,&cfg,&qr,q,&p.X,&p.Y,&p.w),"P3 operator point");if(p.w==0)continue;double u[4],one[2]={1,1};need(gfs_enrich_reference(&cfg,p.X+12,p.Y+12,u,0,0),"P3 original input field");p.rho=u[0];velocity(cfg,p.X+12,p.Y+12,p.v,p.B);p3_shape_stable(m.triangle[it],cfg,p.X,p.Y,one,p.N,p.Nx,p.Ny);
            double dm=p.w*p.rho;for(int a=0;a<12;a++)for(int b=a;b<12;b++){double z=dm*p.N[a]*p.N[b];c.raw[ids[a]*c.n+ids[b]]+=z;if(a!=b)c.raw[ids[b]*c.n+ids[a]]+=z;}
            for(int mode=0;mode<3;mode++){
                Real g[2],D[4];p3_source(mode,p.X,p.Y,cfg.field.ic.eps,g,D);
                for(int a=0;a<12;a++){c.force[mode][0][ids[a]]+=p.w*cfg.field.ic.pressure*p.Nx[a]+(double)(dm*p.N[a]*g[0]);c.force[mode][1][ids[a]]+=p.w*cfg.field.ic.pressure*p.Ny[a]+(double)(dm*p.N[a]*g[1]);}
            }
            c.mass+=dm;c.Kdirect+=.5*dm*(p.v[0]*p.v[0]+p.v[1]*p.v[1]);c.U+=p.w*cfg.field.ic.pressure/(cfg.field.ic.gamma-1);c.Udot-=p.w*cfg.field.ic.pressure*(p.B[0]+p.B[3]);local_points.push_back(p);
        }
    }
    // Actual physical boundaries; internal material triangulation edges are omitted.
    for(const auto &edge:m.linear.edge){
        if(edge.second.size()==2 && m.linear.tri[edge.second[0]].cell==m.linear.tri[edge.second[1]].cell)continue;
        for(int triangle:edge.second){
            const auto &t=m.linear.tri[triangle];if(t.cell%ranks!=rank)continue;auto &c=result.cell[t.cell];int ids[12];p3_local_ids(m,result.cell,triangle,ids);int va=-1,vb=-1;
            for(int e=0;e<3;e++){int a=t.vertex[e],b=t.vertex[(e+1)%3];if(std::min(a,b)==edge.first.first && std::max(a,b)==edge.first.second){va=a;vb=b;}}
            need(va>=0 && vb>=0,"P3 boundary orientation");auto a=m.linear.node[va],b=m.linear.node[vb];double length=hypot(b[0]-a[0],b[1]-a[1]),tx=-cfg.field.ic.pressure*(b[1]-a[1])/length,ty=cfg.field.ic.pressure*(b[0]-a[0])/length;
            for(int q=0;q<cfg.na;q++){
                double z=cfg.an[q],w=length*cfg.aw[q],X=(1-z)*a[0]+z*b[0],Y=(1-z)*a[1]+z*b[1],v[2],B[4],N[12],Nx[12],Ny[12],one[2]={1,1};velocity(cfg,X+12,Y+12,v,B);p3_shape_stable(m.triangle[triangle],cfg,X,Y,one,N,Nx,Ny);
                c.physical_work+=w*(v[0]*tx+v[1]*ty);c.physical_torque+=w*(X*ty-Y*tx);
                if(edge.second.size()==1){
                    for(int mode=0;mode<3;mode++)for(int j=0;j<12;j++){c.force[mode][0][ids[j]]+=w*N[j]*tx;c.force[mode][1][ids[j]]+=w*N[j]*ty;}
                    c.external_work+=w*(v[0]*tx+v[1]*ty);c.external_torque+=w*(X*ty-Y*tx);c.external_force[0]+=w*tx;c.external_force[1]+=w*ty;
                }
            }
        }
    }
    for(auto &c:result.cell){
        all_sum(c.raw);for(auto &mode:c.force)for(auto &v:mode)all_sum(v);std::vector<double> a={c.mass,c.Kdirect,c.U,c.Udot,c.external_work,c.external_force[0],c.external_force[1],c.physical_work,c.physical_torque,c.external_torque};all_sum(a);c.mass=a[0];c.Kdirect=a[1];c.U=a[2];c.Udot=a[3];c.external_work=a[4];c.external_force[0]=a[5];c.external_force[1]=a[6];c.physical_work=a[7];c.physical_torque=a[8];c.external_torque=a[9];
    }
    result.scale.resize(m.linear.components,{0,0});std::vector<double> mass(m.linear.components);
    for(const auto &c:result.cell){mass[c.component]+=c.mass;for(int k=0;k<2;k++)result.scale[c.component][k]+=c.raw[(c.n-2+k)*c.n+c.n-2+k];}
    for(int k=0;k<m.linear.components;k++)for(int d=0;d<2;d++){need(result.scale[k][d]>0,"P3 input residual mass norm");result.scale[k][d]=sqrt(result.scale[k][d]/mass[k]);}
    for(auto &c:result.cell){
        std::vector<double> norm(c.n,1);for(int d=0;d<2;d++)norm[c.n-2+d]=result.scale[c.component][d];c.M.resize(c.n*c.n);
        for(int i=0;i<c.n;i++){for(int j=0;j<c.n;j++)c.M[i*c.n+j]=(Real)c.raw[i*c.n+j]/norm[i]/norm[j];for(auto &mode:c.force)for(auto &v:mode)v[i]/=norm[i];}
    }
    need(local_points.size()*sizeof(P3Point)<INT32_MAX,"bounded P3 point gather size");int bytes=local_points.size()*sizeof(P3Point);std::vector<int> counts(ranks),offsets(ranks);MPI_Gather(&bytes,1,MPI_INT,counts.data(),1,MPI_INT,0,MPI_COMM_WORLD);
    if(rank==0){size_t sum=0;for(int i=0;i<ranks;i++){offsets[i]=sum;sum+=counts[i];}need(sum<INT32_MAX && sum%sizeof(P3Point)==0,"bounded P3 global point packet");result.point.resize(sum/sizeof(P3Point));}
    MPI_Gatherv(local_points.data(),bytes,MPI_BYTE,result.point.data(),counts.data(),offsets.data(),MPI_BYTE,0,MPI_COMM_WORLD);
    if(rank==0)for(size_t it=0;it<m.triangle.size();it++){
        auto add=[&](double X,double Y){P3Point p={};p.triangle=it;p.X=X;p.Y=Y;double u[4],one[2]={1,1};need(gfs_enrich_reference(&cfg,X+12,Y+12,u,0,0),"P3 boundary diagnostic density");p.rho=u[0];velocity(cfg,X+12,Y+12,p.v,p.B);p3_shape_stable(m.triangle[it],cfg,X,Y,one,p.N,p.Nx,p.Ny);result.point.push_back(p);};
        for(int j=0;j<10;j++){auto x=m.node[m.connectivity[it][j]];add(x[0],x[1]);}
        for(int e=0;e<3;e++)for(double q:{1./6,.5,5./6}){const auto &x=m.triangle[it].geometry.x[e],&y=m.triangle[it].geometry.x[(e+1)%3];add((1-q)*x[0]+q*y[0],(1-q)*x[1]+q*y[1]);}
    }
    if(rank==0)std::sort(result.point.begin(),result.point.end(),[](const P3Point &a,const P3Point &b){return std::tie(a.triangle,a.X,a.Y,a.w)<std::tie(b.triangle,b.X,b.Y,b.w);});
    result.owned_point.swap(local_points);
    result.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();return result;
}
struct P3System {
    std::vector<P3Face> faces;std::vector<P3Row> rows,raw;
    std::vector<Factor> factor;std::vector<Vec> minvbt;Vec schur;
    int scalar_nullity=0;Real residual_rank_margin=LDBL_MAX,mass_condition=0,schur_condition=0,schur_min_eigen=0;double seconds=0;
};
static P3System p3_system(const P3Mesh &m,const GfsEnrichConfig &cfg,P3Assembly &a){
    auto start=std::chrono::steady_clock::now();P3System out;out.faces=p3_physical_faces(m,cfg,a.scale);Groups groups(a.dofs);std::vector<Vec> residual_q;int nc=a.cell.size();
    for(size_t fi=0;fi<out.faces.size();fi++){
        auto &f=out.faces[fi];const auto &left=a.cell[f.cell[0]],&right=a.cell[f.cell[1]];
        for(int j=0;j<4+f.residual_rank;j++){
            P3Row row;row.face=fi;row.trace=j;bool independent=false;
            if(j<4){int p=left.offset+left.local(f.node[j]),n=right.offset+right.local(f.node[j]);row.entries={{p,1},{n,-1}};independent=groups.join(p,n);}
            else{
                Vec b(2*nc);for(int k=0;k<2;k++){Real v=f.combination[j-4][k];if(v==0)continue;row.entries.push_back({left.offset+left.n-2+k,v});row.entries.push_back({right.offset+right.n-2+k,-v});b[2*f.cell[0]+k]=v;b[2*f.cell[1]+k]=-v;}
                Real initial=0;for(auto v:b)initial+=v*v;
                for(int pass=0;pass<2;pass++)for(const auto &q:residual_q){Real d=0;for(size_t k=0;k<b.size();k++)d+=b[k]*q[k];for(size_t k=0;k<b.size();k++)b[k]-=d*q[k];}
                Real norm=0;for(auto v:b)norm+=v*v;norm=sqrtl(norm/initial);Real uncertainty=4096*DBL_EPSILON*(1+out.faces.size());
                if(norm>100*uncertainty){independent=true;Real length=sqrtl(initial)*norm;for(auto &v:b)v/=length;residual_q.push_back(b);out.residual_rank_margin=std::min(out.residual_rank_margin,norm/uncertainty);}
                else need(norm<=uncertainty,"P3 global trace rank uncertain");
            }
            out.raw.push_back(row);if(independent){f.selected[j]=out.rows.size();out.rows.push_back(row);}
        }
    }
    int nr=out.rows.size(),nb=a.dofs;out.scalar_nullity=nb-nr;std::vector<Vec> lower(nb,Vec(nr));out.minvbt.resize(nb,Vec(nr));
    for(auto &c:a.cell){
        out.factor.emplace_back(c.M,c.n);auto eigen=extrema(c.M,c.n);need(eigen[0]>0,"P3 mass positive spectrum");out.mass_condition=std::max(out.mass_condition,eigen[1]/eigen[0]);
        for(int j=0;j<nr;j++){
            Vec b(c.n);bool touches=false;for(const auto &e:out.rows[j].entries)if(e.first>=c.offset && e.first<c.offset+c.n){b[e.first-c.offset]=e.second;touches=true;}if(!touches)continue;
            auto l=out.factor.back().lower(b),v=out.factor.back().solve(b);for(int i=0;i<c.n;i++){lower[c.offset+i][j]=l[i];out.minvbt[c.offset+i][j]=v[i];}
        }
    }
    out.schur.resize(nr*nr);for(int i=0;i<nr;i++)for(int j=0;j<=i;j++){Real v=0;for(int k=0;k<nb;k++)v+=lower[k][i]*lower[k][j];out.schur[i*nr+j]=out.schur[j*nr+i]=v;}
    auto eig=extrema(out.schur,nr);need(eig[0]>0,"P3 Schur independent rank SPD");out.schur_min_eigen=eig[0];out.schur_condition=eig[1]/eig[0];out.seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();return out;
}
static Pair p3_field(const P3Mesh &m,const P3Assembly &a,int triangle,const double N[12],const Vec &x,const Vec &y){
    int ids[12];p3_local_ids(m,a.cell,triangle,ids);const auto &c=a.cell[m.linear.tri[triangle].cell];Pair out={0,0};for(int j=0;j<12;j++){out[0]+=(double)(N[j]*x[c.offset+ids[j]]);out[1]+=(double)(N[j]*y[c.offset+ids[j]]);}return out;
}
struct P3Solution {Vec x,y,vx,vy,ax,ay,lx,ly,loadx,loady,error_ax,error_ay;Real equation=0,constraint=0,initial_constraint=0,force_scale=0,linear_correction=0;};
static P3Solution p3_solve(const P3Mesh &m,const GfsEnrichConfig &cfg,const P3Assembly &a,const P3System &s,int mode,double boostx,double boosty,const P3Solution *input_state=nullptr){
    int nb=a.dofs,nr=s.rows.size();P3Solution out;for(auto *v:{&out.x,&out.y,&out.vx,&out.vy,&out.ax,&out.ay,&out.loadx,&out.loady})v->resize(nb);
    if(input_state){out.x=input_state->x;out.y=input_state->y;out.vx=input_state->vx;out.vy=input_state->vy;}
    for(size_t ic=0;ic<a.cell.size();ic++){
        const auto &c=a.cell[ic];Vec fx(c.force[mode][0].begin(),c.force[mode][0].end()),fy(c.force[mode][1].begin(),c.force[mode][1].end());auto ax=s.factor[ic].solve(fx),ay=s.factor[ic].solve(fy);
        for(int i=0;i<c.n;i++){
            int k=c.offset+i;out.ax[k]=ax[i];out.ay[k]=ay[i];if(!input_state && i<c.n-2){auto x=m.node[c.node[i]];double v[2],B[4];velocity(cfg,x[0]+12,x[1]+12,v,B);out.x[k]=x[0];out.y[k]=x[1];out.vx[k]=v[0]+boostx;out.vy[k]=v[1]+boosty;}
        }
        if(!input_state){out.vx[c.offset+c.n-2]=a.scale[c.component][0];out.vy[c.offset+c.n-1]=a.scale[c.component][1];}
    }
    Vec rx(nr),ry(nr);for(int j=0;j<nr;j++)for(auto e:s.rows[j].entries){rx[j]-=e.second*out.ax[e.first];ry[j]-=e.second*out.ay[e.first];}Factor factor(s.schur,nr);out.lx=factor.solve(rx);out.ly=factor.solve(ry);
    for(int j=0;j<nr;j++){
        for(auto e:s.rows[j].entries){out.loadx[e.first]+=e.second*out.lx[j];out.loady[e.first]+=e.second*out.ly[j];}
        for(int i=0;i<nb;i++){out.ax[i]+=s.minvbt[i][j]*out.lx[j];out.ay[i]+=s.minvbt[i][j]*out.ly[j];}
    }
    for(const auto &row:s.raw){Real ax=0,ay=0,x=0,y=0,vx=0,vy=0;for(auto e:row.entries){ax+=e.second*out.ax[e.first];ay+=e.second*out.ay[e.first];x+=e.second*out.x[e.first];y+=e.second*out.y[e.first];vx+=e.second*out.vx[e.first];vy+=e.second*out.vy[e.first];}out.constraint=std::max({out.constraint,fabsl(ax),fabsl(ay)});out.initial_constraint=std::max({out.initial_constraint,fabsl(x),fabsl(y),fabsl(vx),fabsl(vy)});}
    for(const auto &c:a.cell)for(int i=0;i<c.n;i++){
        Real rx=-c.force[mode][0][i]-out.loadx[c.offset+i],ry=-c.force[mode][1][i]-out.loady[c.offset+i];for(int j=0;j<c.n;j++){rx+=c.M[i*c.n+j]*out.ax[c.offset+j];ry+=c.M[i*c.n+j]*out.ay[c.offset+j];}
        out.equation=std::max({out.equation,fabsl(rx),fabsl(ry)});out.force_scale=std::max({out.force_scale,(Real)fabs(c.force[mode][0][i]),(Real)fabs(c.force[mode][1][i]),fabsl(out.loadx[c.offset+i]),fabsl(out.loady[c.offset+i])});
    }
    // Read-only residual correction probe of the SAME saddle equations.
    // It is never applied to a, lambda, heat, or the exported state.
    out.error_ax.resize(nb);out.error_ay.resize(nb);
    for(size_t ic=0;ic<a.cell.size();ic++){
        const auto &c=a.cell[ic];Vec rx(c.n),ry(c.n);
        for(int i=0;i<c.n;i++){rx[i]=c.force[mode][0][i]+out.loadx[c.offset+i];ry[i]=c.force[mode][1][i]+out.loady[c.offset+i];for(int j=0;j<c.n;j++){rx[i]-=c.M[i*c.n+j]*out.ax[c.offset+j];ry[i]-=c.M[i*c.n+j]*out.ay[c.offset+j];}}
        auto x=s.factor[ic].solve(rx),y=s.factor[ic].solve(ry);for(int i=0;i<c.n;i++){out.error_ax[c.offset+i]=x[i];out.error_ay[c.offset+i]=y[i];}
    }
    Vec crx(nr),cry(nr);for(int j=0;j<nr;j++)for(auto e:s.rows[j].entries){crx[j]-=e.second*(out.ax[e.first]+out.error_ax[e.first]);cry[j]-=e.second*(out.ay[e.first]+out.error_ay[e.first]);}auto clx=factor.solve(crx),cly=factor.solve(cry);
    for(int i=0;i<nb;i++){for(int j=0;j<nr;j++){out.error_ax[i]+=s.minvbt[i][j]*clx[j];out.error_ay[i]+=s.minvbt[i][j]*cly[j];}out.linear_correction=std::max({out.linear_correction,fabsl(out.error_ax[i]),fabsl(out.error_ay[i])});}
    return out;
}
struct P3CellMeasure {
    Real volume=0,mass=0,Kdirect=0,initial_v_error=0,initial_B_error=0,error_a2=0,gravity2=0,error_a_max=0,error_grad_L1=0,error_grad_max=0,curvature_L1=0,curvature_max=0,chi_L1=0,chi_max=0,chi_x=0,chi_y=0,gpower=0,gtorque=0,Udot=0,Pdot=0,source_work=0,source_torque=0,fd_error=0,fd_chi_error=0,min_fd_J=LDBL_MAX,linear_chi_error=0;
};
static std::vector<P3CellMeasure> p3_measure_volume(const P3Mesh &m,const GfsEnrichConfig &cfg,const P3Assembly &a,const P3Solution &s,int mode,bool finite_map,double bx,double by){
    std::vector<P3CellMeasure> measure(a.cell.size());
    for(const auto &p:a.point){
        const auto &t=m.linear.tri[p.triangle];const auto &c=a.cell[t.cell];auto &z=measure[t.cell];int ids[12];p3_local_ids(m,a.cell,p.triangle,ids);Real aa[2]={0,0},vv[2]={0,0},Da[4]={0},Dv[4]={0},F0[4]={0},delta_div=0;
        for(int j=0;j<12;j++){
            int k=c.offset+ids[j];Real norm=j<10?1:a.scale[c.component][j-10],N=p.N[j]/norm,Nx=p.Nx[j]/norm,Ny=p.Ny[j]/norm;aa[0]+=N*s.ax[k];aa[1]+=N*s.ay[k];vv[0]+=N*s.vx[k];vv[1]+=N*s.vy[k];Da[0]+=Nx*s.ax[k];Da[1]+=Ny*s.ax[k];Da[2]+=Nx*s.ay[k];Da[3]+=Ny*s.ay[k];Dv[0]+=Nx*s.vx[k];Dv[1]+=Ny*s.vx[k];Dv[2]+=Nx*s.vy[k];Dv[3]+=Ny*s.vy[k];
            Real x=s.x[k]-(j<10?t.t.x[0][0]:0),y=s.y[k]-(j<10?t.t.x[0][1]:0);F0[0]+=Nx*x;F0[1]+=Ny*x;F0[2]+=Nx*y;F0[3]+=Ny*y;
            delta_div+=Nx*s.error_ax[k]+Ny*s.error_ay[k];
        }
        Real g[2],Dg[4];p3_source(mode,p.X,p.Y,cfg.field.ic.eps,g,Dg);Real ex=aa[0]-g[0],ey=aa[1]-g[1],dm=(Real)p.w*p.rho,graderror=0;for(int k=0;k<4;k++)graderror+=(Da[k]-Dg[k])*(Da[k]-Dg[k]);graderror=sqrtl(graderror);
        z.volume+=p.w;z.mass+=dm;z.error_a2+=dm*(ex*ex+ey*ey);z.gravity2+=dm*(g[0]*g[0]+g[1]*g[1]);z.error_a_max=std::max(z.error_a_max,hypotl(ex,ey));z.error_grad_L1+=p.w*graderror;z.error_grad_max=std::max(z.error_grad_max,graderror);
        z.Kdirect+=.5L*dm*(powl((Real)p.v[0]+bx,2)+powl((Real)p.v[1]+by,2));z.initial_v_error=std::max({z.initial_v_error,fabsl(vv[0]-p.v[0]-bx),fabsl(vv[1]-p.v[1]-by)});for(int k=0;k<4;k++)z.initial_B_error=std::max(z.initial_B_error,fabsl(Dv[k]-p.B[k]));
        Real exact=-cfg.field.ic.gamma*(Dg[0]+Dg[3]-((Real)p.B[0]*p.B[0]+2.L*p.B[1]*p.B[2]+(Real)p.B[3]*p.B[3]));Real numerical=-cfg.field.ic.gamma*(Da[0]+Da[3]-(Dv[0]*Dv[0]+2*Dv[1]*Dv[2]+Dv[3]*Dv[3]));Real C=numerical-exact,omega2=powl((Real)p.X*p.X+(Real)p.Y*p.Y+(Real)cfg.field.ic.eps*cfg.field.ic.eps,-1.5L),chi=fabsl(C)/omega2;
        z.curvature_L1+=p.w*fabsl(C);z.curvature_max=std::max(z.curvature_max,fabsl(C));z.chi_L1+=p.w*chi;if(chi>z.chi_max){z.chi_max=chi;z.chi_x=p.X;z.chi_y=p.Y;}
        z.linear_chi_error=std::max(z.linear_chi_error,cfg.field.ic.gamma*fabsl(delta_div)/omega2);
        z.Udot-=p.w*cfg.field.ic.pressure*(Dv[0]+Dv[3]);z.Pdot=std::max(z.Pdot,(Real)cfg.field.ic.gamma*fabsl(Dv[0]+Dv[3]));z.gpower+=dm*(vv[0]*g[0]+vv[1]*g[1]);z.gtorque+=dm*(p.X*g[1]-p.Y*g[0]);z.source_work+=dm*(vv[0]*ex+vv[1]*ey);z.source_torque+=dm*(p.X*ey-p.Y*ex);
        if(finite_map){
            Real h=5.e-6L,J0=F0[0]*F0[3]-F0[1]*F0[2],fd=0;need(J0>0,"P3 finite-map initial determinant positive");
            for(int sign:{-1,1}){Real D[4];for(int j=0;j<4;j++)D[j]=sign*h*Dv[j]+.5L*h*h*Da[j];Real dJ=(F0[3]*D[0]+F0[0]*D[3]-F0[2]*D[1]-F0[1]*D[2]+D[0]*D[3]-D[1]*D[2])/J0;need(dJ>-1,"P3 two-sided finite-map determinant positive");fd+=expm1l(-cfg.field.ic.gamma*log1pl(dJ));z.min_fd_J=std::min(z.min_fd_J,J0*(1+dJ));}
            fd/=h*h;z.fd_error=std::max(z.fd_error,fabsl(fd-numerical));z.fd_chi_error=std::max(z.fd_chi_error,fabsl(fd-numerical)/omega2);
        }
    }return measure;
}
struct P3FaceMeasure {
    Vec directx,directy,work,torque;Real load_error=0,position_jump=0,velocity_jump=0,acceleration_jump=0,pair_work=0,pair_torque=0,traction_difference_square=0,face_measure=0;
};
static P3FaceMeasure p3_measure_faces(const P3Mesh &m,const GfsEnrichConfig &cfg,const P3Assembly &a,const P3System &sys,const P3Solution &s,const char *dest,const char *name,int level){
    P3FaceMeasure out;out.directx.resize(a.dofs);out.directy.resize(a.dofs);out.work.resize(a.cell.size());out.torque.resize(a.cell.size());double q[96],wq[96];gauss(96,q,wq);
    std::string path=std::string(dest)+"/level"+std::to_string(level)+"_"+name+"_traction.csv";FILE *file=fopen(path.c_str(),"w");need(file,"new P3 traction samples");fprintf(file,"face,q,x,y,tx,ty,ax,ay\n");
    for(size_t fi=0;fi<sys.faces.size();fi++){
        const auto &f=sys.faces[fi];int n=4+f.residual_rank;Vec mux(n),muy(n);for(int j=0;j<n;j++)if(f.selected[j]>=0){mux[j]=s.lx[f.selected[j]];muy[j]=s.ly[f.selected[j]];}Factor factor(f.gram,n);auto ex=factor.solve(mux),ey=factor.solve(muy);
        for(int iq=0;iq<96;iq++){
            double N[2][12],Nx[12],Ny[12];for(int side=0;side<2;side++)p3_face_shape(m,cfg,a.scale,f,side,q[iq],N[side],Nx,Ny);Real original[6],phi[6];p3_face_functions(q[iq],N[0],original);for(int j=0;j<4;j++)phi[j]=original[j];for(int j=0;j<f.residual_rank;j++)phi[4+j]=f.combination[j][0]*original[4]+f.combination[j][1]*original[5];
            Real tx=0,ty=0,w=f.length*wq[iq];for(int j=0;j<n;j++){tx+=phi[j]*ex[j];ty+=phi[j]*ey[j];}Pair v[2],x[2],acc[2];
            for(int side=0;side<2;side++){
                int it=f.triangle[side];const auto &c=a.cell[f.cell[side]];v[side]=p3_field(m,a,it,N[side],s.vx,s.vy);x[side]=p3_field(m,a,it,N[side],s.x,s.y);acc[side]=p3_field(m,a,it,N[side],s.ax,s.ay);int ids[12];p3_local_ids(m,a.cell,it,ids);Real sign=side?-1:1;
                for(int j=0;j<12;j++){out.directx[c.offset+ids[j]]+=sign*w*N[side][j]*tx;out.directy[c.offset+ids[j]]+=sign*w*N[side][j]*ty;}out.work[f.cell[side]]+=sign*w*(v[side][0]*tx+v[side][1]*ty);out.torque[f.cell[side]]+=sign*w*(x[side][0]*ty-x[side][1]*tx);
            }
            for(int k=0;k<2;k++){out.position_jump=std::max(out.position_jump,(Real)fabs(x[0][k]-x[1][k]));out.velocity_jump=std::max(out.velocity_jump,(Real)fabs(v[0][k]-v[1][k]));out.acceleration_jump=std::max(out.acceleration_jump,(Real)fabs(acc[0][k]-acc[1][k]));}
            out.pair_work+=w*((v[0][0]-v[1][0])*tx+(v[0][1]-v[1][1])*ty);out.pair_torque+=w*((x[0][0]-x[1][0])*ty-(x[0][1]-x[1][1])*tx);out.traction_difference_square+=w*(powl(tx+cfg.field.ic.pressure*f.nx,2)+powl(ty+cfg.field.ic.pressure*f.ny,2));out.face_measure+=w;
            fprintf(file,"%zu,%d,%.17g,%.17g,%.17Lg,%.17Lg,%.17g,%.17g\n",fi,iq,x[0][0],x[0][1],tx,ty,acc[0][0],acc[0][1]);
        }
    }
    fclose(file);for(int i=0;i<a.dofs;i++)out.load_error=std::max({out.load_error,fabsl(out.directx[i]-s.loadx[i]),fabsl(out.directy[i]-s.loady[i])});return out;
}
static Real p3_internal_jump(const P3Mesh &m,const GfsEnrichConfig &cfg,const P3Assembly &a,const P3Solution &s){
    Real error=0;
    for(const auto &entry:m.linear.edge)if(entry.second.size()==2 && m.linear.tri[entry.second[0]].cell==m.linear.tri[entry.second[1]].cell){
        P3Face f;for(int k=0;k<2;k++){f.triangle[k]=entry.second[k];f.cell[k]=m.linear.tri[entry.second[k]].cell;}f.endpoint[0]=entry.first.first;f.endpoint[1]=entry.first.second;
        for(double q:{.13,.37,.89}){Pair acc[2];for(int side=0;side<2;side++){double N[12],Nx[12],Ny[12];p3_face_shape(m,cfg,a.scale,f,side,q,N,Nx,Ny);acc[side]=p3_field(m,a,f.triangle[side],N,s.ax,s.ay);}for(int k=0;k<2;k++)error=std::max(error,(Real)fabs(acc[0][k]-acc[1][k]));}
    }return error;
}
static void p3_report_case(const P3Mesh &m,const GfsEnrichConfig &cfg,const P3Assembly &a,const P3System &sys,int mode,double bx,double by,int level,const char *name,const char *dest){
    auto start=std::chrono::steady_clock::now();auto s=p3_solve(m,cfg,a,sys,mode,bx,by);bool finite=level==1 && mode==2 && bx==0 && by==0;auto measure=p3_measure_volume(m,cfg,a,s,mode,finite,bx,by);auto face=p3_measure_faces(m,cfg,a,sys,s,dest,name,level);Real internal_jump=p3_internal_jump(m,cfg,a,s);
    std::vector<std::set<int>> boundary(a.cell.size());for(const auto &edge:m.linear.edge){
        if(edge.second.size()==2 && m.linear.tri[edge.second[0]].cell==m.linear.tri[edge.second[1]].cell)continue;
        for(int it:edge.second){const auto &t=m.linear.tri[it];int opposite=-1;for(int j=0;j<3;j++)if(t.vertex[j]!=edge.first.first && t.vertex[j]!=edge.first.second)opposite=j;need(opposite>=0,"P3 physical-cell bubble boundary");for(int j=0;j<10;j++){int index[3];gfs_material_p3_index(j,index);if(index[opposite]==0)boundary[t.cell].insert(m.connectivity[it][j]);}}
    }
    int bubbles=0;Real bubble_load=0,bubble_virtual=0,linear_chi=0;
    for(size_t ic=0;ic<a.cell.size();ic++){const auto &c=a.cell[ic];linear_chi=std::max(linear_chi,measure[ic].linear_chi_error);for(int i=0;i<c.n-2;i++)if(!boundary[ic].count(c.node[i])){
        bubbles++;int k=c.offset+i;Real rx=-c.force[mode][0][i],ry=-c.force[mode][1][i];for(int j=0;j<c.n;j++){rx+=c.M[i*c.n+j]*s.ax[c.offset+j];ry+=c.M[i*c.n+j]*s.ay[c.offset+j];}bubble_virtual=std::max({bubble_virtual,fabsl(rx),fabsl(ry)});bubble_load=std::max({bubble_load,fabsl(s.loadx[k]),fabsl(s.loady[k]),fabsl(face.directx[k]),fabsl(face.directy[k])});
    }}
    printf("{\"physical_cell_polynomial_bubbles\":%d,\"bubble_traction_load_relative\":%.17Lg,\"bubble_virtual_force_relative\":%.17Lg,\"linear_residual_correction_acceleration_max\":%.17Lg,\"linear_residual_correction_chi_max\":%.17Lg,",bubbles,bubble_load/s.force_scale,bubble_virtual/s.force_scale,s.linear_correction,linear_chi);
    Real maxE=0,maxL=0,maxPE=0,maxPL=0,maxC=0,maxCL1=0,maxchi=0,maxchiL1=0,totalE=0,totalL=0,totalU=0,totalK=0,totalmass=0,maxPdot=0,maxv=0,maxB=0,maxFD=0,maxFDchi=0,minFDJ=1;
    printf("\"case\":\"%s\",\"finite_map_applied\":%s,\"momentum_equation_relative\":%.17Lg,\"raw_constraint_max\":%.17Lg,\"initial_constraint_max\":%.17Lg,\"direct_traction_load_relative\":%.17Lg,\"position_trace_jump\":%.17Lg,\"velocity_trace_jump\":%.17Lg,\"acceleration_trace_jump\":%.17Lg,\"internal_triangle_acceleration_jump\":%.17Lg,\"pair_work_absolute\":%.17Lg,\"pair_torque_absolute\":%.17Lg,\"cells\":[",name,finite?"true":"false",s.equation/s.force_scale,s.constraint,s.initial_constraint,face.load_error/s.force_scale,face.position_jump,face.velocity_jump,face.acceleration_jump,internal_jump,face.pair_work,face.pair_torque);
    for(size_t ic=0;ic<a.cell.size();ic++){
        const auto &c=a.cell[ic];const auto &z=measure[ic];Real K=0,Kdot=0,Ldot=0;
        for(int i=0;i<c.n;i++)for(int j=0;j<c.n;j++){int I=c.offset+i,J=c.offset+j;Real m=c.M[i*c.n+j];K+=.5L*m*(s.vx[I]*s.vx[J]+s.vy[I]*s.vy[J]);Kdot+=m*(s.vx[I]*s.ax[J]+s.vy[I]*s.ay[J]);Ldot+=m*(s.x[I]*s.ay[J]-s.y[I]*s.ax[J]);}
        Real external=c.external_work+bx*c.external_force[0]+by*c.external_force[1],E=Kdot+z.Udot-face.work[ic]-external-z.gpower,L=Ldot-face.torque[ic]-c.external_torque-z.gtorque,PE=Kdot+z.Udot-c.physical_work-z.gpower,PL=Ldot-c.physical_torque-z.gtorque;
        totalE+=E;totalL+=L;totalU+=c.U;totalK+=K;totalmass+=c.mass;maxE=std::max(maxE,fabsl(E)/c.U);maxL=std::max(maxL,fabsl(L)/(2*K));maxPE=std::max(maxPE,fabsl(PE)/c.U);maxPL=std::max(maxPL,fabsl(PL)/(2*K));maxC=std::max(maxC,z.curvature_max);maxCL1=std::max(maxCL1,z.curvature_L1/z.volume);maxchi=std::max(maxchi,z.chi_max);maxchiL1=std::max(maxchiL1,z.chi_L1/z.volume);maxPdot=std::max(maxPdot,z.Pdot);maxv=std::max(maxv,z.initial_v_error);maxB=std::max(maxB,z.initial_B_error);maxFD=std::max(maxFD,z.fd_error);maxFDchi=std::max(maxFDchi,z.fd_chi_error);if(finite)minFDJ=std::min(minFDJ,z.min_fd_J);
        if(ic)printf(",");printf("{\"id\":%llu,\"mass\":%.17g,\"U\":%.17g,\"K\":%.17Lg,\"initial_kinetic_integral_error_over_U\":%.17Lg,\"numerical_energy_residual_over_U\":%.17Lg,\"numerical_angular_residual_over_2K\":%.17Lg,\"physical_energy_residual_over_U\":%.17Lg,\"physical_angular_residual_over_2K\":%.17Lg,\"independent_source_work_over_U\":%.17Lg,\"independent_source_torque_over_2K\":%.17Lg,\"acceleration_error_RMS\":%.17Lg,\"gravity_RMS\":%.17Lg,\"acceleration_error_max\":%.17Lg,\"acceleration_gradient_error_L1_volume\":%.17Lg,\"acceleration_gradient_error_max\":%.17Lg,\"pressure_curvature_error_max\":%.17Lg,\"pressure_curvature_error_L1_volume\":%.17Lg,\"chi_max\":%.17Lg,\"chi_L1_volume\":%.17Lg,\"chi_max_x\":%.17Lg,\"chi_max_y\":%.17Lg}",(unsigned long long)m.linear.cell[ic].id,c.mass,c.U,K,(K-z.Kdirect)/c.U,E/c.U,L/(2*K),PE/c.U,PL/(2*K),z.source_work/c.U,z.source_torque/(2*K),sqrtl(z.error_a2/z.mass),sqrtl(z.gravity2/z.mass),z.error_a_max,z.error_grad_L1/z.volume,z.error_grad_max,z.curvature_max,z.curvature_L1/z.volume,z.chi_max,z.chi_L1/z.volume,z.chi_x,z.chi_y);
    }
    double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    printf("],\"max_numerical_energy_residual_over_U\":%.17Lg,\"max_numerical_angular_residual_over_2K\":%.17Lg,\"max_physical_energy_residual_over_U\":%.17Lg,\"max_physical_angular_residual_over_2K\":%.17Lg,\"global_numerical_energy_residual_over_U\":%.17Lg,\"global_numerical_angular_residual_over_2K\":%.17Lg,\"pressure_derivative_over_P_max\":%.17Lg,\"initial_velocity_error\":%.17Lg,\"initial_gradient_error\":%.17Lg,\"curvature_error_max\":%.17Lg,\"max_cell_curvature_error_L1\":%.17Lg,\"chi_max\":%.17Lg,\"max_cell_chi_L1\":%.17Lg,\"traction_difference_L2_over_P\":%.17Lg,\"finite_map_curvature_difference_max\":%.17Lg,\"finite_map_curvature_chi_difference_max\":%.17Lg,\"finite_map_min_J\":%.17Lg,\"total_mass\":%.17Lg,\"total_U\":%.17Lg,\"total_K\":%.17Lg,\"solve_and_diagnostics_seconds\":%.17g}",maxE,maxL,maxPE,maxPL,totalE/totalU,totalL/(2*totalK),maxPdot,maxv,maxB,maxC,maxCL1,maxchi,maxchiL1,sqrtl(face.traction_difference_square/face.face_measure)/cfg.field.ic.pressure,maxFD,maxFDchi,minFDJ,totalmass,totalU,totalK,seconds);
    std::string path=std::string(dest)+"/level"+std::to_string(level)+"_"+name+"_state.bin";FILE *file=fopen(path.c_str(),"wb");need(file,"new P3 hybrid state");uint64_t header[3]={0x4746503348594231ULL,(uint64_t)a.dofs,(uint64_t)sys.rows.size()};fwrite(header,sizeof(header),1,file);for(const auto *v:{&s.x,&s.y,&s.vx,&s.vy,&s.ax,&s.ay,&s.lx,&s.ly})for(auto value:*v){double d=(double)value;need(fwrite(&d,8,1,file)==1,"P3 hybrid state export");}fclose(file);fflush(stdout);
}
static void p3_level(const SharedMesh &linear,int level,const GfsEnrichConfig &cfg,int rank,int ranks,const char *dest){
    auto mesh=p3_build(linear,cfg);auto assembly=p3_assemble(mesh,cfg,rank,ranks);if(rank!=0)return;auto system=p3_system(mesh,cfg,assembly);Real gram_condition=0,gram_uncertainty=0,omission=0;int rank1=0,rank2=0;
    std::string remainder_path=std::string(dest)+"/level"+std::to_string(level)+"_remainder.csv";FILE *remainder=fopen(remainder_path.c_str(),"w");need(remainder,"new P3 stable-remainder oracle queries");fprintf(remainder,"triangle,ax,ay,bx,by,cx,cy,X,Y,eps,Rx,Ry,Rxx,Rxy,Ryx,Ryy\n");
    for(size_t it=0;it<mesh.triangle.size();it++)if(it%17==0 || it+1==mesh.triangle.size()){
        const auto &t=mesh.triangle[it];double X=.2*t.geometry.x[0][0]+.3*t.geometry.x[1][0]+.5*t.geometry.x[2][0],Y=.2*t.geometry.x[0][1]+.3*t.geometry.x[1][1]+.5*t.geometry.x[2][1],N[12],Nx[12],Ny[12],one[2]={1,1};p3_shape_stable(t,cfg,X,Y,one,N,Nx,Ny);
        fprintf(remainder,"%zu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",it,t.geometry.x[0][0],t.geometry.x[0][1],t.geometry.x[1][0],t.geometry.x[1][1],t.geometry.x[2][0],t.geometry.x[2][1],X,Y,cfg.field.ic.eps,N[10],N[11],Nx[10],Ny[10],Nx[11],Ny[11]);
    }fclose(remainder);
    for(const auto &f:system.faces){gram_condition=std::max(gram_condition,f.condition);gram_uncertainty=std::max(gram_uncertainty,f.rank_uncertainty);omission=std::max(omission,f.omitted_function_max);if(f.residual_rank==1)rank1++;else rank2++;}
    size_t matrix_bytes=0;for(const auto &c:assembly.cell)matrix_bytes+=c.M.size()*sizeof(Real)+c.raw.size()*sizeof(double);size_t schur_bytes=system.schur.size()*sizeof(Real);struct rusage use;getrusage(RUSAGE_SELF,&use);
    printf("{\"level\":%d,\"triangles\":%zu,\"physical_cells\":%zu,\"global_polynomial_nodes\":%zu,\"broken_scalar_dofs\":%d,\"physical_faces\":%zu,\"trace_rows\":%zu,\"raw_trace_rows\":%zu,\"conforming_scalar_dofs\":%d,\"residual_rank1_faces\":%d,\"residual_rank2_faces\":%d,\"mass_condition_max\":%.17Lg,\"Schur_condition\":%.17Lg,\"Schur_min_eigenvalue\":%.17Lg,\"trace_Gram_condition_max\":%.17Lg,\"trace_rank_uncertainty_max\":%.17Lg,\"omitted_trace_function_max\":%.17Lg,\"assembly_seconds\":%.17g,\"trace_and_Schur_seconds\":%.17g,\"mass_storage_bytes\":%zu,\"Schur_storage_bytes\":%zu,\"cached_point_bytes\":%zu,\"rank0_maxrss_kib\":%ld,\"cases\":[",level,linear.tri.size(),linear.cell.size(),mesh.node.size(),assembly.dofs,system.faces.size(),system.rows.size(),system.raw.size(),system.scalar_nullity,rank1,rank2,system.mass_condition,system.schur_condition,system.schur_min_eigen,gram_condition,gram_uncertainty,omission,assembly.seconds,system.seconds,matrix_bytes,schur_bytes,assembly.point.size()*sizeof(P3Point),use.ru_maxrss);
    for(int test=0;test<4;test++){if(test)printf(",\n");p3_report_case(mesh,cfg,assembly,system,test==0?0:test==1?1:2,test==3?.013:0,test==3?-.007:0,level,test==0?"gravity_off":test==1?"linear_gravity":test==2?"Kepler":"Kepler_perturb",dest);}printf("]}");fflush(stdout);
}
#ifndef GFS_MATERIAL_P3_HYBRID_LIBRARY
int main(int argc,char **argv){
    MPI_Init(&argc,&argv);int rank,ranks;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&ranks);need(argc==3 || argc==4,"native polygon packet, new output, optional fine_quadrature");bool fine=argc==4;if(fine)need(std::string(argv[3])=="fine_quadrature","known P3 quadrature option");
    std::set<uint64_t> ids={17894,20502,39878,44945,67927,67934,67935,67942,67943,67992,67993,68000,68001,68002,68525,68526,68527,68534,68535,68584,68585,68592,68593,68600};auto mesh=mesh_from(load(argv[1],ids));auto cfg=config(fine?64:32,fine?28:20);
    if(rank==0)printf("{\"scope\":\"fixed P3 material-map INITIAL-RHS with one h/2 internal refinement; generic physical trace rank; no changed generator IC, remap, shock or time-evolution claim\",\"ranks\":%d,\"quadrature\":[%d,%d],\"levels\":[",ranks,cfg.na,cfg.nr);
    p3_level(mesh,0,cfg,rank,ranks,argv[2]);auto refined=p3_refine_once(mesh);if(rank==0)printf(",\n");p3_level(refined,1,cfg,rank,ranks,argv[2]);if(rank==0)printf("]}\n");MPI_Finalize();return 0;
}
#endif
