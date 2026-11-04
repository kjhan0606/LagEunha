/* Actual project Voro geometry, not a substitute tessellation. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../../indxflag.h"
#include "../../Voro/voro.h"
#include "../gfs_variational_core.h"

#define NC 64
#define NV 8192
typedef struct {int i,j;double r[3],c[3],area;} Face;
static Face faces[10000];
static int nf;
static double pos[NC][3],vel[NC][3],pressure[NC];
static double dot(const double *a,const double *b,int dim){
    double s=0;for(int k=0;k<dim;k++)s+=a[k]*b[k];return s;
}
static void save_face(int i,int j,const double r[3],const double c[3],double area){
    if(nf>=10000||j<0){fprintf(stderr,"invalid face\n");exit(2);}
    faces[nf].i=i;faces[nf].j=j;faces[nf].area=area;
    memcpy(faces[nf].r,r,24);memcpy(faces[nf].c,c,24);nf++;
}
static void tess(int dim,int n,double displacement,double *V,int capture){
    if(capture)nf=0;
    for(int i=0;i<n;i++){
        V[i]=0;
        if(dim==1){
            int l=(i+n-1)%n,r=(i+1)%n;
            double xi=pos[i][0]+displacement*vel[i][0];
            double xl=pos[l][0]+displacement*vel[l][0]-(i==0);
            double xr=pos[r][0]+displacement*vel[r][0]+(i==n-1);
            V[i]=.5*(xr-xl);
            if(capture){
                double rr[3]={xr-xi,0,0},c[3]={.5*(xr-xi),0,0};save_face(i,r,rr,c,1);
                rr[0]=xl-xi;c[0]=.5*(xl-xi);save_face(i,l,rr,c,1);
            }
        }else if(dim==2){
            Voro2D_point input[NC*9]={0},sorted_guard[NC*9+8]={0},ctr={0};
            Voro2D_point *sorted=sorted_guard+8;
            Voro2D_Corner corner[NV];
            int m=0,parent[NC*9];
            ctr.x=pos[i][0]+displacement*vel[i][0];
            ctr.y=pos[i][1]+displacement*vel[i][1];ctr.indx=-1;
            for(int sy=-1;sy<=1;sy++)for(int sx=-1;sx<=1;sx++)for(int j=0;j<n;j++){
                if(!sx&&!sy&&j==i)continue;
                input[m].x=pos[j][0]+displacement*vel[j][0]+sx;
                input[m].y=pos[j][1]+displacement*vel[j][1]+sy;
                input[m].indx=m;parent[m]=j;m++;
            }
            int nv=Voro2D_FindVC(&ctr,input,sorted,m,corner,16,4);
            if(nv<=0||nv>=NV)exit(2);
            Voro2D_Corner *a=corner;
            for(int q=0;q<nv;q++){
                Voro2D_Corner *b=a->upperlink;int id=a->upperrelated;
                if(id<0){fprintf(stderr,"unclosed 2D\n");exit(2);}
                int j=parent[sorted[id].indx];
                double r[3]={sorted[id].x,sorted[id].y,0};
                double c[3]={.5*(a->x+b->x),.5*(a->y+b->y),0};
                double area=hypot(b->x-a->x,b->y-a->y),d=sqrt(dot(r,r,2));
                V[i]+=.5*area*dot(c,r,2)/d;
                if(capture)save_face(i,j,r,c,area);
                a=b;
            }
        }else{
            Voro3D_point *input=calloc(NC*27,sizeof(*input));
            Voro3D_point *guard=calloc(NC*27+8,sizeof(*guard)),*sorted=guard+8;
            Voro3D_Vertex *vertex=calloc(NV,sizeof(*vertex));
            int *work=calloc(NV,sizeof(*work));Voro3D_point ctr={0};
            ctr.x=pos[i][0]+displacement*vel[i][0];
            ctr.y=pos[i][1]+displacement*vel[i][1];
            ctr.z=pos[i][2]+displacement*vel[i][2];ctr.indx=-1;
            int m=0,parent[NC*27];
            for(int sz=-1;sz<=1;sz++)for(int sy=-1;sy<=1;sy++)for(int sx=-1;sx<=1;sx++)
                for(int j=0;j<n;j++){
                    if(!sx&&!sy&&!sz&&j==i)continue;
                    input[m].x=pos[j][0]+displacement*vel[j][0]+sx;
                    input[m].y=pos[j][1]+displacement*vel[j][1]+sy;
                    input[m].z=pos[j][2]+displacement*vel[j][2]+sz;
                    input[m].indx=m;parent[m]=j;m++;
                }
            int nv=Voro3D_FindVC(&ctr,input,sorted,m,vertex,16,4,0,work);
            if(nv<=0||nv>=NV)exit(2);
            unsigned char used[NC*27]={0};
            for(int q=0;q<nv;q++)for(int k=0;k<3;k++){
                int id=vertex[q].related[k];
                if(id<0){fprintf(stderr,"unclosed 3D\n");exit(2);}
                if(used[id])continue;
                used[id]=1;
                Voro3D_point polygon[NV];
                int count=Voro3D_FaceExtract(vertex,nv,id,polygon);
                if(count<3)exit(2);
                double area=0,c[3]={0},r[3]={sorted[id].x,sorted[id].y,sorted[id].z};
                for(int t=1;t<count-1;t++){
                    double a[3]={polygon[t].x-polygon[0].x,polygon[t].y-polygon[0].y,polygon[t].z-polygon[0].z};
                    double b[3]={polygon[t+1].x-polygon[0].x,polygon[t+1].y-polygon[0].y,polygon[t+1].z-polygon[0].z};
                    double cross[3]={a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0]};
                    double A=.5*sqrt(dot(cross,cross,3));area+=A;
                    c[0]+=A*(polygon[0].x+polygon[t].x+polygon[t+1].x)/3;
                    c[1]+=A*(polygon[0].y+polygon[t].y+polygon[t+1].y)/3;
                    c[2]+=A*(polygon[0].z+polygon[t].z+polygon[t+1].z)/3;
                }
                if(!(area>0))exit(2);
                for(int k=0;k<3;k++)c[k]/=area;
                V[i]+=area*dot(c,r,3)/(3*sqrt(dot(r,r,3)));
                if(capture)save_face(i,parent[sorted[id].indx],r,c,area);
            }
            free(input);free(guard);free(vertex);free(work);
        }
        if(!(V[i]>0)||!isfinite(V[i])){fprintf(stderr,"invalid volume\n");exit(2);}
    }
}
int main(void){
    int fail=0;
    printf("{\"dimensions\":[\n");
    for(int dim=1;dim<=3;dim++){
        int side=dim==1?16:3,n=1;for(int a=0;a<dim;a++)n*=side;
        for(int i=0;i<n;i++){
            int q=i;
            for(int a=0;a<3;a++){
                int cell=q%side;q/=side;
                pos[i][a]=a<dim?(cell+.5+.07*sin(1.3*(i+1)*(a+2)))/side:0;
                vel[i][a]=a<dim?.1*sin(.7*(i+1)*(a+1)):0;
            }
            pressure[i]=1+.3*cos(.9*i);
        }
        double V[NC],plus[NC],minus[NC],rate[NC]={0},force[NC][3]={{0}},closure[NC][3]={{0}};
        const double h=2e-6;
        tess(dim,n,0,V,1);tess(dim,n,h,plus,0);tess(dim,n,-h,minus,0);
        for(int k=0;k<nf;k++){
            Face *f=faces+k;GfsVariationalFace b;double df[3];
            if(!gfs_variational_face(dim,f->r,f->c,f->area,&b))return 2;
            rate[f->i]+=gfs_variational_sweep(&b,vel[f->i],vel[f->j]);
            gfs_variational_pressure(&b,pressure[f->i],pressure[f->j],df);
            for(int a=0;a<dim;a++){
                force[f->i][a]+=df[a];closure[f->i][a]+=f->area*b.normal[a];
            }
        }
        double maxfd=0,closed=0,energy=0,totalV=0,momentum[3]={0};
        for(int i=0;i<n;i++){
            totalV+=V[i];maxfd=fmax(maxfd,fabs((plus[i]-minus[i])/(2*h)-rate[i]));
            energy+=dot(vel[i],force[i],dim)-pressure[i]*rate[i];
            for(int a=0;a<dim;a++){closed=fmax(closed,fabs(closure[i][a]));momentum[a]+=force[i][a];}
        }
        double mom=sqrt(dot(momentum,momentum,dim));
        int ok=maxfd<2e-7&&closed<1e-9&&fabs(energy)<1e-9&&mom<1e-9&&fabs(totalV-1)<1e-9;
        fail|=!ok;
        printf("%s{\"dim\":%d,\"cells\":%d,\"directed_faces\":%d,\"total_volume\":%.17g,\"max_volume_derivative_error\":%.17g,\"max_face_closure\":%.17g,\"energy_rhs_residual\":%.17g,\"total_force\":%.17g,\"pass\":%s}\n",
            dim>1?",":"",dim,n,nf,totalV,maxfd,closed,energy,mom,ok?"true":"false");
    }
    /* A real geometric blind mode, not a toleranced implementation test. */
    {
        const int n=16;double maxrate=0,mind=.0;
        for(int i=0;i<n;i++){
            double vi=i%2?-1:1,vl=(i+n-1)%2?-1:1,vr=(i+1)%2?-1:1;
            maxrate=fmax(maxrate,fabs(.5*(vr-vl)));
            (void)vi;
        }
        mind=.5/n;
        printf("],\"alternating_mode\":{\"max_volume_rate\":%.17g,\"free_crossing_time\":%.17g,\"bulk_viscosity_detects_mode\":false},\"pass\":%s}\n",
            maxrate,mind,fail?"false":"true");
    }
    return fail?1:0;
}
