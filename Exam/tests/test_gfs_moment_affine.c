#include <math.h>
#include <stdio.h>
#include "../gfs_moment.h"

static int failures;

static void check_close(const char *name,double got,double want,double tol){
    double scale=fmax(1.0,fabs(want));
    if(!isfinite(got) || fabs(got-want)>tol*scale){
        fprintf(stderr,"FAIL %-34s got=%.17g want=%.17g tol=%.3g\n",
                name,got,want,tol);
        failures++;
    }
}

static void make_uniform_affine_cell(GfsMomentState *s,
    double a,double b,double rho,double cx,double cy,double ox,double oy,
    double vx,double vy,const double B[4],double u){
    const double g=1.0/sqrt(3.0);
    gfs_moment_zero(s);
    for(int ix=-1;ix<=1;ix+=2)for(int iy=-1;iy<=1;iy+=2){
        double dx=ix*a*g,dy=iy*b*g;
        double dm=rho*a*b;
        double x=cx+dx-ox,y=cy+dy-oy;
        double qx=vx+B[0]*dx+B[1]*dy;
        double qy=vy+B[2]*dx+B[3]*dy;
        gfs_moment_accumulate(s,dm,x,y,qx,qy,dm*u);
    }
}

static void test_affine_primitive(void){
    const double a=.6,b=.35,rho=1.2,cx=.3,cy=-.4,ox=-.2,oy=.15;
    const double vx=.7,vy=-1.1,u=2.3;
    const double B[4]={.13,-2.0,.8,-.07};
    const double mass=4*a*b*rho;
    const double qxx=mass*a*a/3.0,qyy=mass*b*b/3.0;
    GfsMomentState s;
    GfsMomentPrimitive p;
    make_uniform_affine_cell(&s,a,b,rho,cx,cy,ox,oy,vx,vy,B,u);
    if(!gfs_moment_primitive(&s,&p)){
        fprintf(stderr,"FAIL affine state rejected\n");failures++;return;
    }
    check_close("mass",s.a[GFM_M],mass,2e-14);
    check_close("centroid x",p.c[0],cx-ox,2e-14);
    check_close("centroid y",p.c[1],cy-oy,2e-14);
    check_close("mean vx",p.v[0],vx,2e-14);
    check_close("mean vy",p.v[1],vy,2e-14);
    for(int k=0;k<4;k++){
        char label[32];snprintf(label,sizeof(label),"affine B[%d]",k);
        check_close(label,p.B[k],B[k],3e-13);
    }
    check_close("affine kinetic energy",p.Kaffine,
        .5*(qxx*(B[0]*B[0]+B[2]*B[2])+
            qyy*(B[1]*B[1]+B[3]*B[3])),3e-13);
    check_close("mean kinetic energy",p.Kmean,.5*mass*(vx*vx+vy*vy),3e-13);
    check_close("thermal energy",p.U,mass*u,3e-13);

    double intrinsic=B[2]*qxx-B[1]*qyy;
    double Lwant=(cx*vy-cy*vx)*mass+intrinsic;
    check_close("orbital plus intrinsic Lz",gfs_moment_Lz(&s,ox,oy),Lwant,3e-13);

    GfsMomentState shifted=gfs_moment_shift(s,.41,-.27);
    GfsMomentPrimitive ps;
    if(!gfs_moment_primitive(&shifted,&ps)){
        fprintf(stderr,"FAIL shifted affine state rejected\n");failures++;return;
    }
    check_close("shifted centroid x",ps.c[0],p.c[0]-.41,3e-13);
    check_close("shifted centroid y",ps.c[1],p.c[1]+.27,3e-13);
    for(int k=0;k<4;k++){
        char label[32];snprintf(label,sizeof(label),"shifted B[%d]",k);
        check_close(label,ps.B[k],p.B[k],3e-13);
    }
    check_close("origin-shift Lz",gfs_moment_Lz(&shifted,ox+.41,oy-.27),Lwant,3e-13);

    GfsMomentState boosted=gfs_moment_boost(s,.9,-.35);
    GfsMomentPrimitive pb;
    if(!gfs_moment_primitive(&boosted,&pb)){
        fprintf(stderr,"FAIL boosted affine state rejected\n");failures++;return;
    }
    check_close("boost affine energy",pb.Kaffine,p.Kaffine,3e-13);
    check_close("boost thermal energy",pb.U,p.U,3e-13);
    for(int k=0;k<4;k++){
        char label[32];snprintf(label,sizeof(label),"boosted B[%d]",k);
        check_close(label,pb.B[k],p.B[k],3e-13);
    }
    check_close("boosted Lz",gfs_moment_Lz(&boosted,ox,oy),
        Lwant+mass*(cx*(-.35)-cy*.9),3e-13);
}

static void test_keplerian_shear_split(void){
    const double a=.5,b=.5,rho=1.0,cx=1.3,cy=.0,ox=.0,oy=.0;
    const double omega=1.7,u=1.0;
    const double B[4]={0,-omega,-.5*omega,0};
    GfsMomentState s;GfsMomentPrimitive p;
    make_uniform_affine_cell(&s,a,b,rho,cx,cy,ox,oy,0,0,B,u);
    if(!gfs_moment_primitive(&s,&p)){
        fprintf(stderr,"FAIL Kepler shear state rejected\n");failures++;return;
    }
    double I=p.S[0]+p.S[2];
    double spin_rate=.5*(p.B[2]-p.B[1]);
    double shear_rate=.5*(p.B[2]+p.B[1]);
    check_close("Kepler local spin rate",spin_rate,.25*omega,3e-13);
    check_close("Kepler symmetric shear",shear_rate,-.75*omega,3e-13);
    check_close("Kepler intrinsic Lz",gfs_moment_Lz(&s,ox,oy),
        spin_rate*I,3e-13);
    check_close("Kepler affine energy",p.Kaffine,
        .5*(p.S[0]*omega*omega+
            p.S[2]*.25*omega*omega),3e-13);
}

static void test_pair_face_conservation(void){
    const double ox[2]={-.7,1.2},oy[2]={.4,-.9};
    const double xf=.35,yf=.62,area=.8;
    const double flux[4]={.17,-.23,.41,1.3};
    GfsMomentState rhs[2];
    gfs_moment_zero(&rhs[0]);gfs_moment_zero(&rhs[1]);
    gfs_moment_face_rhs(&rhs[0],xf-ox[0],yf-oy[0],flux,area);
    gfs_moment_face_rhs(&rhs[1],xf-ox[1],yf-oy[1],flux,-area);
    GfsMomentState global0=gfs_moment_shift(rhs[0],-ox[0],-oy[0]);
    GfsMomentState global1=gfs_moment_shift(rhs[1],-ox[1],-oy[1]);
    for(int k=0;k<GFM_NSTATE;k++){
        char label[32];snprintf(label,sizeof(label),"paired global moment[%d]",k);
        check_close(label,global0.a[k]+global1.a[k],0,2e-13);
    }
}

int main(void){
    test_affine_primitive();
    test_keplerian_shear_split();
    test_pair_face_conservation();
    if(failures){fprintf(stderr,"%d affine moment checks failed\n",failures);return 1;}
    puts("PASS affine moment recovery, Kepler spin/shear split, and pair-face conservation");
    return 0;
}
