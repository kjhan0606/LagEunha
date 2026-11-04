/* Original native geometry with degree-three-exact positive face quadrature.
 * Diagnostic moments/perturbations only, no PDE or evolved moment state.
 * The Duffy map is the existing gfs_dg_quad construction generalized to a
 * three-dimensional face triangle. No production DG closure is enabled. */
#define NC 256
#define Voro3D_FindVC Voro3D_FindVCStrict
static void capture_face(int,int,int,const double (*)[3]);
#define GFS_CAPTURE_GEOMETRY_FACE capture_face
#define main archived_geometry_audit_main
#include "geometry_audit.c"
#undef main

#define MAXQUAD 131072
typedef struct {int face;double x[3],weight;} Quad;
static Quad quad[MAXQUAD];static int nq,moments_only;
static double measure[NC],first[NC][3],second[NC][3][3];
static void point(int k,int dim,const double x[3],double w){
    if(w==0)return;
    if(!(w>0)||!isfinite(w))exit(2);
    Face *f=faces+k;int i=f->i;double distance=sqrt(dot(f->r,f->r,dim));
    double radial=dot(x,f->r,dim)/distance;
    measure[i]+=w*radial/dim;
    for(int a=0;a<dim;a++){
        first[i][a]+=w*radial*x[a]/(dim+1);
        for(int b=0;b<dim;b++)second[i][a][b]+=w*radial*x[a]*x[b]/(dim+2);
    }
    if(!moments_only){
        if(nq>=MAXQUAD)exit(2);
        quad[nq].face=k;quad[nq].weight=w;
        for(int a=0;a<3;a++)quad[nq].x[a]=x[a];nq++;
    }
}
static void capture_face(int k,int dim,int count,const double (*v)[3]){
    if(dim==1){point(k,dim,v[0],1);return;}
    if(dim==2){
        if(count!=2)exit(2);
        double a=1/sqrt(3.),length=faces[k].area;
        for(int sign=-1;sign<=1;sign+=2){
            double x[3];for(int d=0;d<3;d++)x[d]=.5*((1-sign*a)*v[0][d]+(1+sign*a)*v[1][d]);
            point(k,dim,x,.5*length);
        }
        return;
    }
    const double node[3]={.5*(1-sqrt(.6)),.5,.5*(1+sqrt(.6))};
    const double weight[3]={5./18,4./9,5./18};
    for(int t=1;t<count-1;t++){
        double a[3],b[3],cross[3];for(int d=0;d<3;d++){a[d]=v[t][d]-v[0][d];b[d]=v[t+1][d]-v[0][d];}
        cross[0]=a[1]*b[2]-a[2]*b[1];cross[1]=a[2]*b[0]-a[0]*b[2];cross[2]=a[0]*b[1]-a[1]*b[0];
        double jacobian=sqrt(dot(cross,cross,3));
        for(int u=0;u<3;u++)for(int z=0;z<3;z++){
            double x[3],s=node[u],r=node[z];
            for(int d=0;d<3;d++)x[d]=(1-s)*(1-r)*v[0][d]+s*v[t][d]+(1-s)*r*v[t+1][d];
            point(k,dim,x,jacobian*(1-s)*weight[u]*weight[z]);
        }
    }
}
static void vector(const double *a,int d){
    printf("[");for(int k=0;k<d;k++)printf("%s%.17g",k?",":"",a[k]);printf("]");
}
int main(int argc,char **argv){
    if(argc!=4&&argc!=6)return 2;
    int dim=atoi(argv[1]),side=atoi(argv[2]);double epsilon=atof(argv[3]);
    if(dim<1||dim>3||side<3||!isfinite(epsilon))return 2;
    int n=1;for(int a=0;a<dim;a++)n*=side;if(n>NC)return 2;
    int primitive=dim==1?16:3;
    for(int i=0;i<n;i++){
        int q=i,tile=0,stride=1,cell[3]={0};
        for(int a=0;a<dim;a++){cell[a]=q%side;q/=side;tile+=(cell[a]%primitive)*stride;stride*=primitive;}
        for(int a=0;a<3;a++)pos[i][a]=a<dim?(cell[a]+.5+.07*sin(1.3*(tile+1)*(a+2)))/side:0;
    }
    if(argc==6){
        FILE *in=fopen(argv[4],"rb");int header[2];moments_only=atoi(argv[5]);
        if(!in||fread(header,sizeof(header),1,in)!=1||header[0]!=dim||header[1]!=n)return 2;
        for(int i=0;i<n;i++)if(fread(vel[i],sizeof(double),dim,in)!=(size_t)dim)return 2;
        if(fread(site_weight_rate,sizeof(double),n,in)!=(size_t)n||fgetc(in)!=EOF||fclose(in))return 2;
    }else if(epsilon!=0)return 2;
    double V[NC];tess(dim,n,epsilon,V,1);
    for(int i=0;i<n;i++)if(fabs(measure[i]-V[i])>1e-11*V[i])return 2;
    printf("{\"scope\":\"original native geometry and degree-three-exact positive face quadrature; no PDE\",\"dimension\":%d,\"side\":%d,\"cells\":%d,\"epsilon\":%.17g,\"positions\":[",dim,side,n,epsilon);
    for(int i=0;i<n;i++){double x[3];for(int a=0;a<dim;a++)x[a]=pos[i][a]+epsilon*vel[i][a];if(i)printf(",");vector(x,dim);}
    printf("],\"volume\":");vector(V,n);
    printf(",\"centroid\":[");
    for(int i=0;i<n;i++){double c[3];for(int a=0;a<dim;a++)c[a]=pos[i][a]+epsilon*vel[i][a]+first[i][a]/V[i];if(i)printf(",");vector(c,dim);}
    printf("],\"central_second_moment\":[");
    for(int i=0;i<n;i++){if(i)printf(",");printf("[");for(int a=0;a<dim;a++){if(a)printf(",");double row[3];for(int b=0;b<dim;b++)row[b]=second[i][a][b]-first[i][a]*first[i][b]/V[i];vector(row,dim);}printf("]");}
    printf("]");
    if(!moments_only){
        printf(",\"faces\":[");int cursor=0;
        for(int k=0;k<nf;k++){
            Face *f=faces+k;
            printf("%s{\"i\":%d,\"j\":%d,\"area\":%.17g,\"r\":",k?",":"",f->i,f->j,f->area);vector(f->r,dim);
            printf(",\"c\":");vector(f->c,dim);printf(",\"quad\":[");int first_point=1;
            while(cursor<nq&&quad[cursor].face==k){if(!first_point)printf(",");first_point=0;printf("[%.17g",quad[cursor].weight);for(int a=0;a<dim;a++)printf(",%.17g",quad[cursor].x[a]);printf("]");cursor++;}printf("]}");
        }
        if(cursor!=nq)return 2;printf("]");
    }
    printf(",\"accepted_steps\":0}\n");return 0;
}
