/* Static raw geometry exchange and local common-segment reconstruction.
 * int32 N,nface; owner[N]; double rows[][14]: raw i,j,xi,yi,dx,dy,ax,ay,bx,by,epoch, reference area,cx,cy.
 * Distributed tessellation, HG and time integration are outside this harness. */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define W 14
static int cmp(const void *aa,const void *bb){
    const double *a=aa,*b=bb;
    for(int k=0;k<2;k++) if(a[k]!=b[k]) return a[k]<b[k]?-1:1;
    return 0;
}
static double *reciprocal(double *local,int n,const double *f){
    double key[W]={0}; key[0]=f[1];key[1]=f[0];
    int lo=0,hi=n;
    while(lo<hi){int m=(lo+hi)/2;if(cmp(local+m*W,key)<0)lo=m+1;else hi=m;}
    for(int k=lo;k<n && cmp(local+k*W,key)==0;k++){
        double *g=local+k*W;
        if(fabs(g[4]+f[4])+fabs(g[5]+f[5])<=1e-9*hypot(f[4],f[5])) return g;
    }
    return NULL;
}
int main(int argc,char **argv){
    MPI_Init(&argc,&argv);
    int rank,np;MPI_Comm_rank(MPI_COMM_WORLD,&rank);MPI_Comm_size(MPI_COMM_WORLD,&np);
    if(argc<2)MPI_Abort(MPI_COMM_WORLD,2);
    char path[4096];snprintf(path,sizeof(path),"%s/part_%d_%d.bin",argv[1],np,rank);
    FILE *fp=fopen(path,"rb");if(!fp)MPI_Abort(MPI_COMM_WORLD,3);
    int n,nf;if(fread(&n,4,1,fp)!=1 || fread(&nf,4,1,fp)!=1 || n<1 || nf<0)MPI_Abort(MPI_COMM_WORLD,4);
    int *owner=malloc(n*sizeof(int));double *local=malloc((size_t)(nf+1)*W*sizeof(double));
    if(fread(owner,4,n,fp)!=(size_t)n || fread(local,sizeof(double)*W,nf,fp)!=(size_t)nf)MPI_Abort(MPI_COMM_WORLD,5);
    fclose(fp);
    int *sc=calloc(np,sizeof(int)),*rc=calloc(np,sizeof(int)),*sd=calloc(np,sizeof(int)),*rd=calloc(np,sizeof(int)),*next=calloc(np,sizeof(int));
    long long remote=0;
    for(int k=0;k<nf;k++){
        double *f=local+k*W;int i=(int)f[0],j=(int)f[1];
        if(i<0||i>=n||j<0||j>=n||owner[i]!=rank||owner[j]<0||owner[j]>=np)MPI_Abort(MPI_COMM_WORLD,6);
        sc[owner[j]]+=W;remote+=(owner[j]!=rank);
    }
    if(argc>2 && !strcmp(argv[2],"stale") && rank==0 && nf) local[10]=0;
    MPI_Alltoall(sc,1,MPI_INT,rc,1,MPI_INT,MPI_COMM_WORLD);
    for(int p=1;p<np;p++){sd[p]=sd[p-1]+sc[p-1];rd[p]=rd[p-1]+rc[p-1];}
    int nr=(rd[np-1]+rc[np-1])/W;
    memcpy(next,sd,np*sizeof(int));
    double *send=malloc((size_t)(nf+1)*W*sizeof(double)),*recv=malloc((size_t)(nr+1)*W*sizeof(double));
    for(int k=0;k<nf;k++){double *f=local+k*W;int p=owner[(int)f[1]];memcpy(send+next[p],f,W*sizeof(double));next[p]+=W;}
    if(argc>2 && !strcmp(argv[2],"corrupt") && rank==0){
        for(int k=0;k<nf;k++){double *f=send+k*W;if(f[11]>1e-6 && owner[(int)f[1]]!=rank){f[8]=f[6];f[9]=f[7];break;}}
    }
    MPI_Alltoallv(send,sc,sd,MPI_DOUBLE,recv,rc,rd,MPI_DOUBLE,MPI_COMM_WORLD);
    qsort(recv,nr,W*sizeof(double),cmp);
    long long stats[6]={nf,nr,remote,0,0,0},total[6];double maxima[3]={0},global[3];
    for(int k=0;k<nf;k++){
        double *f=local+k*W,*g=reciprocal(recv,nr,f);
        if(f[10]!=1 || (g && g[10]!=1)){stats[3]++;continue;}
        double dx=f[4],dy=f[5],d=hypot(dx,dy),tx=-dy/d,ty=dx/d,mx=.5*dx,my=.5*dy;
        double s1=(f[6]-mx)*tx+(f[7]-my)*ty,s2=(f[8]-mx)*tx+(f[9]-my)*ty;
        double lo=fmin(s1,s2),hi=fmax(s1,s2),area=0,cx=mx+tx*.5*(lo+hi),cy=my+ty*.5*(lo+hi);
        if(!g)stats[4]++;
        else{
            double s3=(g[6]+dx-mx)*tx+(g[7]+dy-my)*ty,s4=(g[8]+dx-mx)*tx+(g[9]+dy-my)*ty;
            lo=fmax(lo,fmin(s3,s4));hi=fmin(hi,fmax(s3,s4));
            if(hi>lo){area=hi-lo;cx=mx+tx*.5*(lo+hi);cy=my+ty*.5*(lo+hi);}
        }
        double da=fabs(area-f[11]),dc=area>1e-14?hypot(cx-f[12],cy-f[13]):0;
        maxima[0]=fmax(maxima[0],da);maxima[1]=fmax(maxima[1],dc);
        if(!isfinite(area)||!isfinite(dc)||da>1e-12*fmax(f[11],.01)||dc>1e-12)stats[5]++;
    }
    MPI_Allreduce(stats,total,6,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(maxima,global,3,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    if(!rank)printf("{\"ranks\":%d,\"sent\":%lld,\"received\":%lld,\"remote\":%lld,\"stale\":%lld,\"missing\":%lld,\"mismatch\":%lld,\"max_area_error\":%.17g,\"max_centroid_error\":%.17g,\"unused\":%.17g}\n",np,total[0],total[1],total[2],total[3],total[4],total[5],global[0],global[1],global[2]);
    int bad=total[0]!=total[1] || total[3] || total[5];
    free(owner);free(local);free(send);free(recv);free(sc);free(rc);free(sd);free(rd);free(next);
    MPI_Finalize();return bad?1:0;
}
