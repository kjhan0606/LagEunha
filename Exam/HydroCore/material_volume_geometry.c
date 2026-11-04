/* Actual native geometry exporter for a finite material-volume algebra proof.
 * Not a hydrodynamic solver. Optional displacement is a directional geometry
 * probe, with velocity and weight rate supplied by the independent verifier. */
#define NC 1024
#define Voro3D_FindVC Voro3D_FindVCStrict
#define main archived_geometry_audit_main
#include "geometry_audit.c"
#undef main
#include "../gfs_velocity_fit.h"
static void vector(const double *a,int n){
    printf("[");for(int i=0;i<n;i++)printf("%s%.17g",i?",":"",a[i]);printf("]");
}
int main(int argc,char **argv){
    if(argc!=4&&argc!=6)return 2;
    int dim=atoi(argv[1]),side=atoi(argv[2]);double fraction=atof(argv[3]),epsilon=0;
    if(dim<1||dim>3||side<3||!isfinite(fraction))return 2;
    int n=1;for(int a=0;a<dim;a++)n*=side;if(n>NC)return 2;
    int primitive=dim==1?16:3;
    for(int i=0;i<n;i++){
        int q=i,tile=0,stride=1,cell[3]={0};
        for(int a=0;a<dim;a++){cell[a]=q%side;q/=side;tile+=(cell[a]%primitive)*stride;stride*=primitive;}
        for(int a=0;a<3;a++)pos[i][a]=a<dim?(cell[a]+.5+.07*sin(1.3*(tile+1)*(a+2)))/side:0;
        site_weight[i]=fraction/(side*side)*sin(.79*(tile+1));
    }
    double mean=0;for(int i=0;i<n;i++)mean+=site_weight[i]/n;
    for(int i=0;i<n;i++)site_weight[i]-=mean;
    if(argc==6){
        FILE *in=fopen(argv[4],"rb");int header[2];epsilon=atof(argv[5]);
        if(!in||!isfinite(epsilon)||fread(header,sizeof(header),1,in)!=1||header[0]!=dim||header[1]!=n)return 2;
        for(int i=0;i<n;i++)if(fread(vel[i],sizeof(double),dim,in)!=(size_t)dim)return 2;
        if(fread(site_weight_rate,sizeof(double),n,in)!=(size_t)n||fgetc(in)!=EOF||fclose(in))return 2;
    }
    double V[NC];tess(dim,n,epsilon,V,1);
    GfsVelocityFit fit[NC];double inverse[NC][3][3];
    for(int i=0;i<n;i++)gfs_velocity_fit_init(fit+i,dim);
    for(int k=0;k<nf;k++){
        double zero[3]={0};if(!gfs_velocity_fit_add_face(fit+faces[k].i,faces[k].r,zero,faces[k].area))return 2;
    }
    for(int i=0;i<n;i++){
        for(int a=0;a<dim;a++)for(int b=0;b<dim;b++)fit[i].R[a][b]=a==b;
        if(!gfs_velocity_fit_solve(fit+i,inverse[i],NULL))return 2;
    }
    printf("{\"dimension\":%d,\"side\":%d,\"cells\":%d,\"weight_fraction\":%.17g,\"epsilon\":%.17g,\"volume\":",dim,side,n,fraction,epsilon);vector(V,n);
    printf(",\"positions\":[");
    for(int i=0;i<n;i++){double x[3];for(int a=0;a<dim;a++)x[a]=pos[i][a]+epsilon*vel[i][a];if(i)printf(",");vector(x,dim);}
    printf("],\"weights\":[");for(int i=0;i<n;i++)printf("%s%.17g",i?",":"",site_weight[i]+epsilon*site_weight_rate[i]);
    printf("],\"native_inverse_fit\":[");
    for(int i=0;i<n;i++){if(i)printf(",");printf("[");for(int a=0;a<dim;a++){if(a)printf(",");vector(inverse[i][a],dim);}printf("]");}
    printf("],\"faces\":[");
    for(int k=0;k<nf;k++){
        Face *f=faces+k;GfsVariationalFace b;if(!gfs_variational_face(dim,f->r,f->c,f->area,&b))return 2;
        printf("%s{\"i\":%d,\"j\":%d,\"area\":%.17g,\"r\":",k?",":"",f->i,f->j,f->area);vector(f->r,dim);
        printf(",\"c\":");vector(f->c,dim);printf(",\"bi\":");vector(b.bi,dim);printf(",\"bj\":");vector(b.bj,dim);printf("}");
    }
    printf("]}\n");return 0;
}
