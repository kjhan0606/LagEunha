#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <string>
#include "gfs_enriched_euler.h"
#include "gfs_ale_flux.h"
static void need(bool a,const char *s){if(!a){fprintf(stderr,"FAIL %s\n",s);exit(1);}}
static void gauss(int n,double *x,double *w){
    for(int i=0;i<(n+1)/2;i++){
        long double z=cosl(acosl(-1.L)*(i+.75L)/(n+.5L)),d=0;
        for(int it=0;it<100;it++){
            long double a=1,b=z;for(int k=2;k<=n;k++){long double c=((2*k-1)*z*b-(k-1)*a)/k;a=b;b=c;}
            d=n*(z*b-a)/(z*z-1);long double next=z-b/d;
            if(fabsl(next-z)<4*LDBL_EPSILON){z=next;break;}z=next;
        }
        long double a=1,b=z;for(int k=2;k<=n;k++){long double c=((2*k-1)*z*b-(k-1)*a)/k;a=b;b=c;}
        d=n*(z*b-a)/(z*z-1);double ww=1/((1-z*z)*d*d);
        x[i]=(1-z)*.5L;x[n-1-i]=(1+z)*.5L;w[i]=w[n-1-i]=ww;
    }
}
struct Cell{uint64_t id;double gx,gy;std::vector<double>x,y;};
static std::map<uint64_t,Cell> load(const char *path,const std::set<uint64_t> &selected){
    std::map<uint64_t,Cell> out;std::string root(path);
    FILE *f=fopen((root+"/dg_project_r0.bin").c_str(),"rb");need(f,"project packet");
    uint64_t h[8];need(fread(h,sizeof(h),1,f)==1 && h[0]==0x4746534447503031ULL && h[1]==1,"project header");
    std::vector<char> record(h[6]);
    for(uint64_t n=0;n<h[2];n++){
        need(fread(record.data(),record.size(),1,f)==1,"project row");
        uint64_t id;memcpy(&id,record.data(),8);if(!selected.count(id))continue;
        Cell c={};c.id=id;memcpy(&c.gx,record.data()+8,8);memcpy(&c.gy,record.data()+16,8);out[id]=c;
    }fclose(f);
    f=fopen((root+"/dg_polygons_r0.bin").c_str(),"rb");need(f,"polygon packet");
    uint64_t hp[4];need(fread(hp,sizeof(hp),1,f)==1 && hp[0]==0x4746534447504f31ULL,"polygon header");
    for(uint64_t n=0;n<hp[2];n++){
        uint64_t id;int tag[2];need(fread(&id,8,1,f)==1 && fread(tag,8,1,f)==1,"polygon row");
        std::vector<double> xy(2*tag[0]);need(fread(xy.data(),8,xy.size(),f)==xy.size(),"polygon data");
        if(!out.count(id))continue;need(tag[1]==1,"native polygon valid");
        for(int k=0;k<tag[0];k++){out[id].x.push_back(xy[2*k]);out[id].y.push_back(xy[2*k+1]);}
    }fclose(f);need(out.size()==selected.size(),"all selected cells");return out;
}
static GfsEnrichConfig config(int na,int nr){
    GfsEnrichConfig c={};c.na=na;c.nr=nr;gauss(na,c.an,c.aw);gauss(nr,c.rn,c.rw);
    c.field.ic.kind=2;c.field.ic.gamma=(double)1.666667f;c.field.ic.pressure=1.e-6;
    c.field.ic.cx=c.field.ic.cy=12;c.field.ic.rin=2;c.field.ic.rout=8;c.field.ic.eps=.047;
    return c;
}
static GfsEnrichState boost(const GfsEnrichState &s,double bx,double by){
    GfsEnrichState t=s;for(int k=0;k<GFS_ENRICH_N;k++){
        t.c[1][k]=s.c[1][k]+bx*s.c[0][k];t.c[2][k]=s.c[2][k]+by*s.c[0][k];
        t.c[3][k]=s.c[3][k]+bx*s.c[1][k]+by*s.c[2][k]+.5*(bx*bx+by*by)*s.c[0][k];
    }return t;
}
static void transform(double u[4],double bx,double by){
    u[3]+=bx*u[1]+by*u[2]+.5*(bx*bx+by*by)*u[0];u[1]+=bx*u[0];u[2]+=by*u[0];
}
struct Assessment{
    double minp=DBL_MAX,maxp=-DBL_MAX,perror=0,field_error=0,conservation=0,integral[4]={};int nonpositive=0;
};
static void assess_point(const GfsEnrichBasis &b,const GfsEnrichConfig &cfg,
    const GfsEnrichState &s,double x,double y,double w,Assessment &a,double bx=0,double by=0){
    double u[4],v[4];need(gfs_enrich_evaluate(&b,&cfg,&s,x,y,u),"evaluate");
    need(gfs_enrich_reference(&cfg,b.generator[0]+x,b.generator[1]+y,v,0,0),"reference");transform(v,bx,by);
    double p=gfs_dg_pressure(u,cfg.field.ic.gamma);a.minp=fmin(a.minp,p);a.maxp=fmax(a.maxp,p);
    a.perror=fmax(a.perror,fabs(p/cfg.field.ic.pressure-1));if(!(p>0))a.nonpositive++;
    for(int c=0;c<4;c++){
        a.field_error=fmax(a.field_error,fabs(u[c]-v[c])/b.scale[c]);
        if(w>0)a.integral[c]+=w*u[c];
    }
}
/* Positive line rule segmented about the point closest to the physical
 * origin. Same reference origin and r+1e-6 scale as volume integration. */
template<class Fun> static void faces(const Cell &cell,const GfsEnrichConfig &cfg,Fun fun){
    for(size_t e=0;e<cell.x.size();e++){
        size_t j=(e+1)%cell.x.size();double ex=cell.x[j]-cell.x[e],ey=cell.y[j]-cell.y[e],len=hypot(ex,ey);
        if(len==0)continue;
        double ox=cfg.field.ic.cx-cell.gx-cell.x[e],oy=cfg.field.ic.cy-cell.gy-cell.y[e];
        double t0=fmax(0.,fmin(1.,(ox*ex+oy*ey)/(len*len)));
        double d=hypot(cell.x[e]+t0*ex+cell.gx-cfg.field.ic.cx,cell.y[e]+t0*ey+cell.gy-cfg.field.ic.cy);
        for(int side=0;side<2;side++){
            double delta=(side?1:0)-t0;if(delta==0)continue;
            double hi=fmin(1.,(d+1.e-6)/(fabs(delta)*len)),lo=0;
            for(int seg=0;seg<GFS_ENRICH_SEGMENTS;seg++){
                for(int k=0;k<cfg.na;k++){
                    double s=lo+(hi-lo)*cfg.an[k],t=t0+delta*s;
                    fun(cell.x[e]+t*ex,cell.y[e]+t*ey,ey/len,-ex/len,len*fabs(delta)*(hi-lo)*cfg.aw[k]);
                }
                if(hi==1)break;lo=hi;hi=fmin(1.,4*hi);
            }
        }
    }
}
static Assessment assess(const Cell &cell,const GfsEnrichBasis &b,const GfsEnrichConfig &cfg,
    const GfsEnrichState &s,const double integral[4],double bx=0,double by=0){
    Assessment a;GfsMomentPolygon p={cell.x.data(),cell.y.data(),(int)cell.x.size(),0};
    GfsEnrichQuadrature r;need(gfs_enrich_quadrature_init(p,&cfg,&b,&r),"assessment rule");
    for(int n=0;n<gfs_enrich_nquad(p,&cfg,&r);n++){
        double x,y,w;need(gfs_enrich_quad(p,&cfg,&r,n,&x,&y,&w),"assessment quadrature");
        if(w>0)assess_point(b,cfg,s,x,y,w,a,bx,by);
    }
    for(size_t k=0;k<cell.x.size();k++)assess_point(b,cfg,s,cell.x[k],cell.y[k],0,a,bx,by);
    faces(cell,cfg,[&](double x,double y,double nx,double ny,double w){assess_point(b,cfg,s,x,y,0,a,bx,by);});
    double target[4];for(int c=0;c<4;c++)target[c]=integral[c];transform(target,bx,by);
    for(int c=0;c<4;c++)a.conservation=fmax(a.conservation,fabs(a.integral[c]-target[c])/(b.geometry.area*b.scale[c]));
    return a;
}
static void analytic_rhs(const GfsEnrichConfig &cfg,double x,double y,double bx,double by,
    double gx,double gy,double out[4]){
    double u[4],dx[4],dy[4];need(gfs_enrich_reference(&cfg,x,y,u,dx,dy),"analytic gradients");
    transform(u,bx,by);transform(dx,bx,by);transform(dy,bx,by);
    double vx=u[1]/u[0],vy=u[2]/u[0],p=gfs_dg_pressure(u,cfg.field.ic.gamma);
    double vxx=(dx[1]-vx*dx[0])/u[0],vxy=(dy[1]-vx*dy[0])/u[0];
    double vyx=(dx[2]-vy*dx[0])/u[0],vyy=(dy[2]-vy*dy[0])/u[0];
    double px=(cfg.field.ic.gamma-1)*(dx[3]-vx*dx[1]-vy*dx[2]+.5*(vx*vx+vy*vy)*dx[0]);
    double py=(cfg.field.ic.gamma-1)*(dy[3]-vx*dy[1]-vy*dy[2]+.5*(vx*vx+vy*vy)*dy[0]);
    out[0]=-dx[1]-dy[2];
    out[1]=-(dx[1]*vx+u[1]*vxx+px+dy[1]*vy+u[1]*vyy)+u[0]*gx;
    out[2]=-(dx[2]*vx+u[2]*vxx+dy[2]*vy+u[2]*vyy+py)+u[0]*gy;
    out[3]=-((dx[3]+px)*vx+(u[3]+p)*vxx+(dy[3]+py)*vy+(u[3]+p)*vyy)+u[1]*gx+u[2]*gy;
}
struct RhsResult{double error=0,weaknorm=0,strongnorm=0,poly_dot=0,enrich_dot=0,mean[4]={};};
static RhsResult rhs(const Cell &cell,const GfsEnrichBasis &b,const GfsEnrichConfig &cfg,
    const GfsEnrichState &s,double bx,double by){
    RhsResult result;double weak[4][GFS_ENRICH_N]={},strong[4][GFS_ENRICH_N]={};
    double gram[GFS_ENRICH_N][GFS_ENRICH_N]={};
    GfsMomentPolygon p={cell.x.data(),cell.y.data(),(int)cell.x.size(),0};GfsEnrichQuadrature r;
    need(gfs_enrich_quadrature_init(p,&cfg,&b,&r),"RHS rule");
    for(int n=0;n<gfs_enrich_nquad(p,&cfg,&r);n++){
        double x,y,w;need(gfs_enrich_quad(p,&cfg,&r,n,&x,&y,&w),"RHS quadrature");if(w==0)continue;
        double X=cell.gx+x-cfg.field.ic.cx,Y=cell.gy+y-cfg.field.ic.cy;
        double inv=pow(X*X+Y*Y+cfg.field.ic.eps*cfg.field.ic.eps,-1.5),gx=-X*inv,gy=-Y*inv;
        need(gfs_enrich_volume_point(&b,&cfg,&s,x,y,gx,gy,w,weak),"full Euler volume/source");
        double q[GFS_ENRICH_N],udot[4];need(gfs_enrich_q(&b,&cfg,GFS_ENRICH_N,x,y,q,0,0),"RHS tests");
        analytic_rhs(cfg,cell.gx+x,cell.gy+y,bx,by,gx,gy,udot);
        for(int k=0;k<GFS_ENRICH_N;k++){
            double test=k?q[k]-b.mean[k]:1;
            for(int c=0;c<4;c++)strong[c][k]+=w*test*udot[c];
            for(int j=0;j<=k;j++)gram[k][j]+=w*test*(j?q[j]-b.mean[j]:1);
        }
    }
    faces(cell,cfg,[&](double x,double y,double nx,double ny,double w){
        double u[4];need(gfs_enrich_evaluate(&b,&cfg,&s,x,y,u),"face evaluate");
        GfsAlePrimitive v={u[0],gfs_dg_pressure(u,cfg.field.ic.gamma),u[1]/u[0],u[2]/u[0]};
        GfsAleFlux f=gfs_ale_flux(v,v,nx,ny,0,0,cfg.field.ic.gamma);need(f.valid,"actual full HLLC face flux");
        need(gfs_enrich_face_point(&b,&cfg,x,y,f.flux,w,weak),"full Euler face weak form");
    });
    double speed=sqrt(2*b.scale[3]/b.scale[0]);double rate=speed/fmin(b.geometry.l00,b.geometry.l11);
    for(int c=0;c<4;c++)for(int k=0;k<GFS_ENRICH_N;k++){
        double scale=b.geometry.area*b.scale[c]*rate;
        result.error=fmax(result.error,fabs(weak[c][k]-strong[c][k])/scale);
        result.weaknorm=fmax(result.weaknorm,fabs(weak[c][k])/scale);
        result.strongnorm=fmax(result.strongnorm,fabs(strong[c][k])/scale);
    }
    int active[GFS_ENRICH_N],nactive=0;
    for(int k=0;k<GFS_ENRICH_N;k++)if(b.L[k][k]>0)active[nactive++]=k;
    double L[GFS_ENRICH_N][GFS_ENRICH_N]={};
    for(int i=0;i<nactive;i++)for(int j=0;j<=i;j++){
        int a=active[i],bb=active[j];double v=gram[a][bb];
        for(int k=0;k<j;k++)v=fma(-L[i][k],L[j][k],v);
        need(j<i || v>0,"actual mean-free Gram has positive active pivots");
        L[i][j]=j==i?sqrt(v):v/L[j][j];
    }
    for(int c=0;c<4;c++){
        double z[GFS_ENRICH_N];for(int i=0;i<nactive;i++){
            z[i]=weak[c][active[i]];for(int j=0;j<i;j++)z[i]-=L[i][j]*z[j];z[i]/=L[i][i];
        }
        for(int i=nactive-1;i>=0;i--){for(int j=i+1;j<nactive;j++)z[i]-=L[j][i]*z[j];z[i]/=L[i][i];}
        for(int i=0;i<nactive;i++){
            int k=active[i];double v=fabs(z[i])/(b.scale[c]*rate);
            if(k>0 && k<6)result.poly_dot=fmax(result.poly_dot,v);
            if(k>=6)result.enrich_dot=fmax(result.enrich_dot,v);
        }result.mean[c]=weak[c][0];
    }return result;
}
static double gradient_contract(){
    GfsEnrichConfig c=config(16,12);double maxerror=0;
    const double xy[5][2]={{12.001,12.002},{12.02,12.01},{13.7,12.3},{12.1,16.8},{20.3,12.1}};
    for(auto &p:xy){
        double u[4],dx[4],dy[4];need(gfs_enrich_reference(&c,p[0],p[1],u,dx,dy),"derivative reference");
        long double h=1.e-4L*hypotl((long double)p[0]-12,(long double)p[1]-12);
        for(int d=0;d<2;d++){
            long double v[4][4];int offsets[4]={-2,-1,1,2};
            for(int n=0;n<4;n++){
                long double x=(long double)p[0]-12+(d?0:offsets[n]*h),y=(long double)p[1]-12+(d?offsets[n]*h:0);
                long double r=hypotl(x,y),rho=.001L+2/(r+1.e-6L)*.5L*(1+tanhl((r-2)/.25L))*.5L*(1+tanhl((8-r)/.25L));
                long double omega=powl(r*r+(long double)c.field.ic.eps*c.field.ic.eps,-.75L),vx=-y*omega,vy=x*omega;
                v[n][0]=rho;v[n][1]=rho*vx;v[n][2]=rho*vy;
                v[n][3]=(long double)c.field.ic.pressure/(c.field.ic.gamma-1)+.5L*rho*(vx*vx+vy*vy);
            }
            for(int k=0;k<4;k++){
                long double deriv=(v[0][k]-8*v[1][k]+8*v[2][k]-v[3][k])/(12*h);
                maxerror=fmax(maxerror,fabsl(deriv-(d?dy[k]:dx[k]))/(1+fabs(d?dy[k]:dx[k])));
            }
        }
    }return maxerror;
}
#ifndef GFS_ENRICH_TEST_LIBRARY
int main(int argc,char **argv){
    need(argc>=2,"native polygon directory argument");
    std::set<uint64_t> selected={67943,68000,68527,68584,67935,68001,68526,68592};
    if(argc>2){selected.clear();for(int k=2;k<argc;k++)selected.insert(strtoull(argv[k],0,10));}
    auto cells=load(argv[1],selected);
    double grad_error=gradient_contract();need(grad_error<1.e-8,"independent analytic reference derivative gate");
    printf("{\"scope\":\"selected actual native polygons; enriched representation and static-cell full Euler RHS, no evolution or moving-weight solver\",\n\"gradient_error\":%.17g,\n\"cells\":[\n",grad_error);
    bool first=true;
    for(auto &entry:cells){
        Cell &c=entry.second;GfsMomentPolygon polygon={c.x.data(),c.y.data(),(int)c.x.size(),0};
        for(auto count:std::vector<std::pair<int,int>>{{24,12},{48,20},{64,28}}){
            auto cfg=config(count.first,count.second);GfsEnrichBasis b={};GfsEnrichState s={};double integral[4];
            need(gfs_enrich_build(polygon,&cfg,c.gx,c.gy,&b,&s,integral),"enriched basis build");
            auto check=config(64,28);Assessment a=assess(c,b,check,s,integral);
            auto boosted=boost(s,.013,-.007);Assessment ab=assess(c,b,check,boosted,integral,.013,-.007);
            RhsResult eq=rhs(c,b,check,s,0,0),pert=rhs(c,b,check,boosted,.013,-.007);
            if(!first)printf(",\n");first=false;
            printf("{\"id\":%llu,\"angular_n\":%d,\"radial_n\":%d,\"rank\":%d,\"minimum_norm\":%.17g,\"omitted_norm\":%.17g,\"orthogonality_error\":%.17g,\"area_error\":%.17g,\"min_P\":%.17g,\"max_relative_P_error\":%.17g,\"relative_field_error\":%.17g,\"independent_integral_error\":%.17g,\"nonpositive_points\":%d,\"boost_max_relative_P_error\":%.17g,\"eq_rhs_error\":%.17g,\"eq_rhs_norm\":%.17g,\"perturbed_rhs_error\":%.17g,\"perturbed_rhs_norm\":%.17g,\"perturbed_polynomial_coefficient_dot\":%.17g,\"perturbed_enrichment_coefficient_dot\":%.17g,\"integral\":[%.17g,%.17g,%.17g,%.17g]}",
                (unsigned long long)c.id,cfg.na,cfg.nr,b.rank,b.minimum_norm,b.omitted_norm,b.orthogonality_error,b.area_error,
                a.minp,a.perror,a.field_error,a.conservation,a.nonpositive,ab.perror,eq.error,eq.weaknorm,pert.error,pert.weaknorm,
                pert.poly_dot,pert.enrich_dot,integral[0],integral[1],integral[2],integral[3]);fflush(stdout);
        }
    }printf("\n]}\n");return 0;
}
#endif
