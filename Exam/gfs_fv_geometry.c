/* Actual project interval/Voro2D/strict Voro3D geometry for the new FV
 * branch. The linked-cell search follows the existing native g3 adapter;
 * its security radius includes the actual weight range. No jitter, mesh
 * steering, finite face-area deletion, or substitute tessellation. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <omp.h>
#include "../indxflag.h"
#include "../Voro/voro.h"
#include "gfs_fv_native.h"
#define FV_VERTEX_CAP 8192
typedef struct {
    Voro2D_point *p2,*s2;Voro2D_Corner *corner;
    Voro3D_point *p3,*s3,*polygon;Voro3D_Vertex *vertex;
    int *parent,*mirror,*work;unsigned char *used;double (*rel)[3];int capacity;
} NativeWork;
static void *fv_memory(void *p,size_t bytes){
    void *out=realloc(p,bytes?bytes:1);if(!out){fprintf(stderr,"FV geometry allocation failed\n");exit(12);}return out;
}
static double fv_dot(const double a[3],const double b[3],int d){double s=0;for(int k=0;k<d;k++)s+=a[k]*b[k];return s;}
static double fv_wrap(double x,double lo,double L){return lo+(x-lo-L*floor((x-lo)/L));}
static void fv_grow(NativeWork *w,int n,int d){
    if(n<=w->capacity)return;int cap=2*n+64;
    if(d==2){w->p2=fv_memory(w->p2,cap*sizeof(*w->p2));w->s2=fv_memory(w->s2,(cap+8)*sizeof(*w->s2));memset(w->s2,0,8*sizeof(*w->s2));}
    if(d==3){w->p3=fv_memory(w->p3,cap*sizeof(*w->p3));w->s3=fv_memory(w->s3,(cap+8)*sizeof(*w->s3));memset(w->s3,0,8*sizeof(*w->s3));}
    w->parent=fv_memory(w->parent,cap*sizeof(int));w->mirror=fv_memory(w->mirror,cap*sizeof(int));
    w->rel=fv_memory(w->rel,cap*sizeof(*w->rel));w->used=fv_memory(w->used,cap);w->capacity=cap;
}
static GfsFvNativeFace *fv_face(GfsFvGeometry *g){
    if(g->faces==g->face_capacity){g->face_capacity=2*g->face_capacity+16;g->face=fv_memory(g->face,g->face_capacity*sizeof(*g->face));}
    GfsFvNativeFace *f=g->face+g->faces++;memset(f,0,sizeof(*f));f->point=g->points;return f;
}
static void fv_point(GfsFvGeometry *g,GfsFvNativeFace *f,const double x[3],double weight){
    if(!(weight>0))return;
    if(g->points==g->point_capacity){g->point_capacity=2*g->point_capacity+64;g->point=fv_memory(g->point,g->point_capacity*sizeof(*g->point));}
    GfsFvPoint *p=g->point+g->points++;memcpy(p->x,x,3*sizeof(double));p->weight=weight;f->count++;
}
static void fv_finish_face(GfsFvGeometry *g,GfsFvNativeFace *f,int dim){
    f->distance=sqrt(fv_dot(f->r,f->r,dim));
    for(int a=0;a<dim;a++)f->normal[a]=f->r[a]/f->distance;
    for(int k=0;k<f->count;k++){
        GfsFvPoint *p=g->point+f->point+k;f->area+=p->weight;
        for(int a=0;a<dim;a++)f->centroid[a]+=p->weight*p->x[a];
    }
    for(int a=0;a<dim;a++)f->centroid[a]/=f->area;
}
static int fv_cell_moments(GfsFvGeometry *g,int dim){
    g->V=0;double first[3]={0},closure[3]={0},area=0;
    for(int k=0;k<g->faces;k++){
        GfsFvNativeFace *f=g->face+k;
        if(!(f->area>0)||!(f->distance>0))return 1;
        area+=f->area;for(int a=0;a<dim;a++)closure[a]+=f->area*f->normal[a];
        for(int q=0;q<f->count;q++){
            GfsFvPoint *p=g->point+f->point+q;double radial=fv_dot(p->x,f->normal,dim);
            g->V+=p->weight*radial/dim;
            for(int a=0;a<dim;a++)first[a]+=p->weight*radial*p->x[a]/(dim+1);
        }
    }
    if(!(g->V>0)||!isfinite(g->V))return 1;
    for(int a=0;a<dim;a++){
        if(fabs(closure[a])>4096*DBL_EPSILON*area)return 1;
        g->centroid[a]=first[a]/g->V;
    }
    return 0;
}
static int fv_intervals(GfsFvMesh *m){
    int n=m->n;
    for(int i=0;i<n;i++){
        GfsFvGeometry *g=m->geometry+i;g->faces=g->points=0;
        for(int sign=-1;sign<=1;sign+=2){
            int j=i+sign;double xj,point[3]={0};int wall=0;
            if(j<0||j>=n){
                if(m->bc[0]==GFS_FV_PERIODIC){j=(j+n)%n;xj=m->cell[j].q[0]+sign*m->L[0];}
                else{j=i;wall=1;xj=2*(m->lo[0]+(sign>0?m->L[0]:0))-m->cell[i].q[0];}
            }else xj=m->cell[j].q[0];
            double r=xj-m->cell[i].q[0];if(!(sign*r>0))return 1;
            point[0]=.5*r+(m->cell[i].W-m->cell[j].W)/(2*r);
            GfsFvNativeFace *f=fv_face(g);f->j=j;f->mirror=wall;f->r[0]=r;
            fv_point(g,f,point,1);fv_finish_face(g,f,1);
        }
        if(fv_cell_moments(g,1))return 1;
    }
    return 0;
}
static int fv_native_cell(GfsFvMesh *m,int i,NativeWork *w,const int nc[3],
        const double cw[3],const int *head,const int *next,const double (*x)[3],double Wmax,int *max_candidates){
    const int dim=m->dim;int ci[3]={0},Rlimit=2;
    double cwmin=HUGE_VAL;
    for(int a=0;a<dim;a++){ci[a]=(int)floor((x[i][a]-m->lo[a])/cw[a]);if(ci[a]>=nc[a])ci[a]=nc[a]-1;if(ci[a]<0)ci[a]=0;Rlimit=Rlimit>2*nc[a]+4?Rlimit:2*nc[a]+4;cwmin=fmin(cwmin,cw[a]);}
    GfsFvGeometry *g=m->geometry+i;
    for(int R=1;R<=Rlimit;R++){
        int count=0;
        for(int dz=dim==3?-R:0;dz<=(dim==3?R:0);dz++)for(int dy=-R;dy<=R;dy++)for(int dx=-R;dx<=R;dx++){
            int delta[3]={dx,dy,dz},cell[3]={0},mirror=0,skip=0;double shift[3]={0};
            for(int a=0;a<dim;a++){
                int c=ci[a]+delta[a];
                if(m->bc[a]==GFS_FV_PERIODIC){int wrap=(int)floor((double)c/nc[a]);cell[a]=c-wrap*nc[a];shift[a]=wrap*m->L[a];}
                else{if(c<0){cell[a]=-1-c;mirror|=1<<a;}else if(c>=nc[a]){cell[a]=2*nc[a]-1-c;mirror|=1<<(a+3);}else cell[a]=c;if(cell[a]<0||cell[a]>=nc[a])skip=1;}
            }
            if(skip)continue;
            int b=(cell[2]*nc[1]+cell[1])*nc[0]+cell[0];
            for(int j=head[b];j>=0;j=next[j]){
                double rel[3]={0};int identical=mirror==0;
                for(int a=0;a<dim;a++){
                    double q=x[j][a]+shift[a];if(shift[a]!=0)identical=0;
                    if(mirror&(1<<a))q=2*m->lo[a]-q;
                    if(mirror&(1<<(a+3)))q=2*(m->lo[a]+m->L[a])-q;
                    rel[a]=q-x[i][a];
                }
                if(identical&&j==i)continue;
                fv_grow(w,count+1,dim);w->parent[count]=j;w->mirror[count]=mirror;
                memcpy(w->rel[count],rel,sizeof(rel));
                if(dim==2){Voro2D_point *p=w->p2+count;memset(p,0,sizeof(*p));p->x=rel[0];p->y=rel[1];p->w2=m->cell[j].W;p->indx=count;}
                else{Voro3D_point *p=w->p3+count;memset(p,0,sizeof(*p));p->x=rel[0];p->y=rel[1];p->z=rel[2];p->w2=m->cell[j].W;p->indx=count;}
                count++;
            }
        }
        if(count>*max_candidates)*max_candidates=count;
        if(count<dim+1)continue;
        double radius2=0;int bad=0,vertices;
        if(dim==2){
            Voro2D_point center={0};center.w2=m->cell[i].W;center.indx=-1;
            vertices=Voro2D_FindVCPowerStrict(&center,w->p2,w->s2+8,count,w->corner,FV_VERTEX_CAP,2*R*cwmin);
            if(vertices<=0||vertices>=FV_VERTEX_CAP)return 1;
            Voro2D_Corner *p=w->corner;
            for(int k=0;k<vertices;k++){if(p->upperrelated<0)bad=1;radius2=fmax(radius2,p->x*p->x+p->y*p->y);p=p->upperlink;}
        }else{
            Voro3D_point center={0};center.w2=m->cell[i].W;center.indx=-1;
            vertices=Voro3D_FindVCPowerStrict(&center,w->p3,w->s3+8,count,w->vertex,FV_VERTEX_CAP,2*R*cwmin,0,w->work);
            if(vertices<=0||vertices>=FV_VERTEX_CAP)return 1;
            for(int k=0;k<vertices;k++)if(w->vertex[k].status==Active){
                Voro3D_Vertex *p=w->vertex+k;radius2=fmax(radius2,p->x*p->x+p->y*p->y+p->z*p->z);
                for(int a=0;a<3;a++)if(p->related[a]<0)bad=1;
            }
        }
        /* Omitted neighbors are at least R*cwmin away. For every cell point
         * |x|<=radius their radical half-space cannot cut this cell if
         * d>radius+sqrt(radius^2+max(Wmax-Wi,0)). Gauge invariant. */
        double radius=sqrt(radius2),needed=radius+sqrt(radius2+fmax(0,Wmax-m->cell[i].W));
        if(bad||needed>R*cwmin)continue;
        g->faces=g->points=0;
        if(dim==2){
            const double node[3]={-sqrt(.6),0,sqrt(.6)},weight[3]={5./9,8./9,5./9};
            Voro2D_Corner *a=w->corner;
            for(int k=0;k<vertices;k++){
                Voro2D_Corner *b=a->upperlink;int id=w->s2[8+a->upperrelated].indx;
                if(id<0||id>=count)return 1;
                double A=hypot(b->x-a->x,b->y-a->y);
                if(A>0){GfsFvNativeFace *f=fv_face(g);f->j=w->parent[id];f->mirror=w->mirror[id];memcpy(f->r,w->rel[id],sizeof(f->r));
                    for(int q=0;q<3;q++){double v[3]={.5*((1-node[q])*a->x+(1+node[q])*b->x),.5*((1-node[q])*a->y+(1+node[q])*b->y),0};fv_point(g,f,v,.5*A*weight[q]);}
                    fv_finish_face(g,f,dim);}
                a=b;
            }
        }else{
            const double node[3]={.5*(1-sqrt(.6)),.5,.5*(1+sqrt(.6))},weight[3]={5./18,4./9,5./18};
            memset(w->used,0,count);
            for(int k=0;k<vertices;k++)if(w->vertex[k].status==Active)for(int a=0;a<3;a++){
                int related=w->vertex[k].related[a];if(related<0||related>=count)return 1;if(w->used[related])continue;w->used[related]=1;
                int id=w->s3[8+related].indx;if(id<0||id>=count)return 1;
                int size=Voro3D_FaceExtract(w->vertex,vertices,related,w->polygon);if(size<3)return 1;
                GfsFvNativeFace *f=fv_face(g);f->j=w->parent[id];f->mirror=w->mirror[id];memcpy(f->r,w->rel[id],sizeof(f->r));
                for(int t=1;t<size-1;t++){
                    double v0[3]={w->polygon[0].x,w->polygon[0].y,w->polygon[0].z};
                    double v1[3]={w->polygon[t].x,w->polygon[t].y,w->polygon[t].z};
                    double v2[3]={w->polygon[t+1].x,w->polygon[t+1].y,w->polygon[t+1].z};
                    double u[3],v[3];for(int b=0;b<3;b++){u[b]=v1[b]-v0[b];v[b]=v2[b]-v0[b];}
                    double cross[3]={u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0]},jac=sqrt(fv_dot(cross,cross,3));
                    for(int s=0;s<3;s++)for(int r=0;r<3;r++){
                        double p[3];for(int b=0;b<3;b++)p[b]=(1-node[s])*(1-node[r])*v0[b]+node[s]*v1[b]+(1-node[s])*node[r]*v2[b];
                        fv_point(g,f,p,jac*(1-node[s])*weight[s]*weight[r]);
                    }
                }
                if(!f->count){g->faces--;continue;}fv_finish_face(g,f,dim);
            }
        }
        return fv_cell_moments(g,dim);
    }
    return 1;
}
int gfs_fv_geometry(GfsFvMesh *m){
    double start=omp_get_wtime();int bad=0,d=m->dim,n=m->n;
    if(d<1||d>3||n<2)return 1;
    for(int i=0;i<n;i++)for(int a=0;a<d;a++)if(!isfinite(m->cell[i].q[a])||!isfinite(m->cell[i].W)||
        (m->bc[a]!=GFS_FV_PERIODIC&&!(m->cell[i].q[a]>m->lo[a]&&m->cell[i].q[a]<m->lo[a]+m->L[a])))return 1;
    if(!m->geometry)m->geometry=calloc(n,sizeof(*m->geometry));
    if(d==1)bad=fv_intervals(m);
    else{
        double box=1;for(int a=0;a<d;a++)box*=m->L[a];double h=pow(box/n,1./d),cw[3]={1,1,1};int nc[3]={1,1,1},total=1;
        for(int a=0;a<d;a++){nc[a]=(int)floor(m->L[a]/(2*h));if(nc[a]<1)nc[a]=1;cw[a]=m->L[a]/nc[a];total*=nc[a];}
        int *head=fv_memory(NULL,total*sizeof(int)),*next=fv_memory(NULL,n*sizeof(int));double (*x)[3]=fv_memory(NULL,n*sizeof(*x));double Wmax=-HUGE_VAL;
        for(int i=0;i<total;i++)head[i]=-1;
        for(int i=0;i<n;i++){
            int c[3]={0};for(int a=0;a<3;a++){x[i][a]=a<d?fv_wrap(m->cell[i].q[a],m->lo[a],m->L[a]):0;if(a<d){c[a]=(int)floor((x[i][a]-m->lo[a])/cw[a]);if(c[a]>=nc[a])c[a]=nc[a]-1;}}
            int b=(c[2]*nc[1]+c[1])*nc[0]+c[0];next[i]=head[b];head[b]=i;Wmax=fmax(Wmax,m->cell[i].W);
        }
        int maximum=0;
#pragma omp parallel if(n>1024) reduction(max:bad,maximum)
        {
            NativeWork w={0};
            if(d==2)w.corner=fv_memory(NULL,FV_VERTEX_CAP*sizeof(*w.corner));
            else{w.vertex=fv_memory(NULL,FV_VERTEX_CAP*sizeof(*w.vertex));w.polygon=fv_memory(NULL,FV_VERTEX_CAP*sizeof(*w.polygon));w.work=fv_memory(NULL,FV_VERTEX_CAP*sizeof(int));}
#pragma omp for schedule(dynamic,32)
            for(int i=0;i<n;i++)if(fv_native_cell(m,i,&w,nc,cw,head,next,x,Wmax,&maximum))bad=1;
            free(w.p2);free(w.s2);free(w.corner);free(w.p3);free(w.s3);free(w.polygon);free(w.vertex);free(w.parent);free(w.mirror);free(w.work);free(w.rel);free(w.used);
        }
        if(maximum>m->maximum_candidates)m->maximum_candidates=maximum;
        free(head);free(next);free(x);
    }
    if(!bad){
        long double volume=0;double box=1;for(int a=0;a<d;a++)box*=m->L[a];
        for(int i=0;i<n;i++)volume+=m->geometry[i].V;
        double error=fabs((double)volume/box-1);m->maximum_volume_closure=fmax(m->maximum_volume_closure,error);if(error>5e-10)bad=1;
    }
    m->geometry_seconds+=omp_get_wtime()-start;return bad;
}
void gfs_fv_geometry_free(GfsFvMesh *m){
    if(m->geometry)for(int i=0;i<m->n;i++){free(m->geometry[i].face);free(m->geometry[i].point);}free(m->geometry);m->geometry=NULL;
}
