/* Actual native geometry and point walk versus exhaustive power minimum.
 * Exhaustive O(N^2) work exists ONLY in this small read-only audit. */
#include <stdio.h>
#include <stdlib.h>
#include "../gfs_fv_native.h"
static long double potential(const GfsFvMesh *m,int j,const double x[3]){
    long double value=(long double)m->cell[0].W-m->cell[j].W;
    for(int a=0;a<m->dim;a++){
        long double dx=(long double)x[a]-m->cell[j].q[a];
        if(gfs_fv_boundary_type(m->bc[a],0)==GFS_FV_PERIODIC)dx-=m->L[a]*floorl(dx/m->L[a]+.5L);
        value+=dx*dx;
    }return value;
}
static int check(const GfsFvMesh *m,const double x[3],int *maximum,long long *count){
    int owner=-1,other=-1,hops=0;double rel[3];
    int code=gfs_fv_point_owner(m,0,x,&owner,rel,&hops);if(code){fprintf(stderr,"walk failure dim=%d code=%d point=(%.17g %.17g %.17g)\n",m->dim,code,x[0],x[1],x[2]);return 1;}
    if(hops>*maximum)*maximum=hops;
    if(gfs_fv_point_owner(m,m->n-1,x,&other,rel,&hops)||owner!=other){fprintf(stderr,"seed-dependent owner dim=%d owners=%d,%d\n",m->dim,owner,other);return 1;}
    long double best=potential(m,0,x);for(int j=1;j<m->n;j++)best=fminl(best,potential(m,j,x));
    if(potential(m,owner,x)-best>512*DBL_EPSILON){fprintf(stderr,"incorrect exhaustive minimum\n");return 1;}
    (*count)++;return 0;
}
int main(void){
    int bad=0,maximum=0;long long count=0;
    for(int dim=1;dim<=3;dim++)for(int fixture=0;fixture<2;fixture++){
        GfsFvMesh m={0};m.dim=dim;m.n=1;m.gamma=1.4;for(int a=0;a<3;a++)m.L[a]=1;for(int a=0;a<dim;a++)m.n*=4;m.cell=calloc(m.n,sizeof(*m.cell));
        for(int i=0;i<m.n;i++){
            int cell=i;for(int a=0;a<dim;a++){m.cell[i].q[a]=(cell%4+.5)/4;cell/=4;if(fixture==0)m.cell[i].q[a]+=.01*sin(.7*(i+1)*(a+1));}
            int x=i%4;m.cell[i].W=fixture==1?(x==1||x==2?1.2/16:0):.01/16*cos(.9*i);m.cell[i].id=i;
        }
        if(gfs_fv_geometry(&m))return 2;
        for(int i=0;i<m.n;i++){
            bad|=check(&m,m.cell[i].q,&maximum,&count);
            for(int k=0;k<m.geometry[i].faces;k++){
                double x[3];for(int a=0;a<3;a++)x[a]=m.cell[i].q[a]+m.geometry[i].face[k].centroid[a];
                bad|=check(&m,x,&maximum,&count);
            }
        }
        for(int k=0;k<50;k++){
            double x[3]={0};for(int a=0;a<dim;a++)x[a]=fmod((k+.13)*(a==0?.6180339887:a==1?.4142135623:.7320508075),1.);
            bad|=check(&m,x,&maximum,&count);
            double shifted[3];for(int a=0;a<3;a++)shifted[a]=x[a]+(a<dim?(a==0?3:a==1?-2:1):0);
            int i,j,h;double rel[3];if(gfs_fv_point_owner(&m,0,x,&i,rel,&h)||gfs_fv_point_owner(&m,0,shifted,&j,rel,&h)||i!=j)bad=1;
        }
        double vertex[3]={.5,dim>=2?.5:0,dim==3?.5:0};bad|=check(&m,vertex,&maximum,&count);
        vertex[0]=0;bad|=check(&m,vertex,&maximum,&count);
        int owners[64];double rel[3];int hops;
        for(int i=0;i<m.n;i++)if(gfs_fv_point_owner(&m,i,m.cell[i].q,owners+i,rel,&hops))bad=1;
        for(int i=0;i<m.n;i++)m.cell[i].W+=.03125;
        if(gfs_fv_geometry(&m))return 2;
        for(int i=0;i<m.n;i++){int owner;if(gfs_fv_point_owner(&m,i,m.cell[i].q,&owner,rel,&hops)||owner!=owners[i])bad=1;}
        gfs_fv_geometry_free(&m);free(m.cell);
    }
    /* Fixed-wall points outside the physical box must fail, never wrap. */
    GfsFvMesh m={0};m.dim=1;m.n=4;m.gamma=1.4;m.L[0]=1;m.bc[0]=GFS_FV_REFLECT;m.cell=calloc(4,sizeof(*m.cell));for(int i=0;i<4;i++){m.cell[i].q[0]=(i+.5)/4;m.cell[i].id=i;}
    if(gfs_fv_geometry(&m))return 2;double outside[3]={-.1,0,0},relative[3];int owner,hops;
    if(gfs_fv_point_owner(&m,0,outside,&owner,relative,&hops)==0)bad=1;gfs_fv_geometry_free(&m);free(m.cell);
    printf("{\"scope\":\"native D1/D2/D3 point walk, true weighted cells, face ties, periodic shifts, gauge and boundary rejection; no PDE\",\"queries\":%lld,\"maximum_hops\":%d,\"pass\":%s}\n",count,maximum,bad?"false":"true");return bad?1:0;
}
