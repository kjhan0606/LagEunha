/* Native power-cell point location by real face adjacency, not a second
 * tessellation or an all-particle search. Ties use the least cell ID in
 * the connected roundoff-bounded containing-cell fan (fixed cap64).
 * Strict descent has at most N distinct owners; failure never falls back
 * to evaluating an unrelated cell field. */
#include <float.h>
#include <math.h>
#include "gfs_fv_native.h"
static void point_relative(const GfsFvMesh *m,int k,const double point[3],long double x[3]){
    for(int a=0;a<3;a++){
        x[a]=a<m->dim?(long double)point[a]-m->cell[k].q[a]:0;
        if(a<m->dim&&gfs_fv_boundary_type(m->bc[a],0)==GFS_FV_PERIODIC){
            long double L=m->L[a];x[a]-=L*floorl(x[a]/L+.5L);
        }
    }
}
static long double power_difference(const GfsFvMesh *m,int k,
        const GfsFvNativeFace *f,const long double x[3],long double *bound){
    long double r2=0,dot=0,size=0,dw=(long double)m->cell[k].W-m->cell[f->j].W;
    for(int a=0;a<m->dim;a++){
        long double r=f->r[a];r2+=r*r;dot+=r*x[a];size+=fabsl(r*x[a]);
    }
    /* Same double-input, long-double, roundoff-bounded power half-space
     * arithmetic as the strict native geometry. Algebraic W differences
     * do not modify or recenter stored weights. */
    *bound=128*DBL_EPSILON*(r2+2*size+fabsl(dw));
    return r2-2*dot+dw;
}
static int containing(const GfsFvMesh *m,int k,const double point[3],
        int *next,int *touching){
    long double x[3],worst=0;point_relative(m,k,point,x);*next=-1;*touching=0;
    const GfsFvGeometry *g=m->geometry+k;
    for(int f=0;f<g->faces;f++){
        const GfsFvNativeFace *face=g->face+f;long double bound,difference=power_difference(m,k,face,x,&bound);
        if(!isfinite(difference)||!isfinite(bound))return -1;
        if(difference < -bound){
            if(face->mirror)return -2;
            if(*next<0||difference<worst||(difference==worst&&face->j<*next)){*next=face->j;worst=difference;}
        }else if(fabsl(difference)<=bound)*touching=1;
    }
    return *next<0;
}
static long double reduced_potential(const GfsFvMesh *m,int seed,int k,const double point[3]){
    long double x[3],value=(long double)m->cell[seed].W-m->cell[k].W;point_relative(m,k,point,x);
    for(int a=0;a<m->dim;a++)value+=x[a]*x[a];return value;
}
int gfs_fv_point_owner(const GfsFvMesh *m,int seed,const double point[3],
        int *owner,double relative[3],int *hops){
    if(seed<0||seed>=m->n||!m->geometry)return 1;
    for(int a=0;a<m->dim;a++)if(!isfinite(point[a]))return 1;
    int current=seed,touch=0,next=-1;*hops=0;
    for(;;){
        int state=containing(m,current,point,&next,&touch);
        if(state<0)return state==-2?3:2;
        if(state)break;
        if(next==current||++*hops>m->n)return 4;
        long double before=reduced_potential(m,seed,current,point),after=reduced_potential(m,seed,next,point);
        if(!(after<before))return 5; /* numerical inconsistency/cycle */
        current=next;
    }
    if(touch){
        int queue[64],count=1;queue[0]=current;int best=current;
        for(int at=0;at<count;at++){
            int k=queue[at];long double x[3];point_relative(m,k,point,x);
            const GfsFvGeometry *g=m->geometry+k;
            for(int f=0;f<g->faces;f++){
                const GfsFvNativeFace *face=g->face+f;if(face->mirror)continue;
                long double bound,difference=power_difference(m,k,face,x,&bound);
                if(fabsl(difference)>bound)continue;
                int j=face->j,visited=0;for(int v=0;v<count;v++)if(queue[v]==j){visited=1;break;}
                if(visited)continue;
                int candidate,ties,state=containing(m,j,point,&candidate,&ties);
                if(state<0)return 2;if(!state)continue;
                if(count==64)return 6;queue[count++]=j;(*hops)++;
                if(m->cell[j].id<m->cell[best].id)best=j;
            }
        }
        current=best;
    }
    int state=containing(m,current,point,&next,&touch);if(state!=1)return 7;
    long double x[3];point_relative(m,current,point,x);
    for(int a=0;a<3;a++)relative[a]=(double)x[a];*owner=current;return 0;
}
