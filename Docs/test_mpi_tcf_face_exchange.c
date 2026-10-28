/* Static canonical TCF fixture exchange. Not production MPI hydro evolution.
 * Fixture: int32 N,nface; int32 owner[N]; double face[nface][22].
 * Columns 0..20 are the CPU face dump; 21 is a fixture generation stamp.
 */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define W 22
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
        if(fabs(g[6]+f[6])+fabs(g[7]+f[7])<2e-6) return g;
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
    if(argc>2 && !strcmp(argv[2],"stale") && rank==0 && nf) local[21]=0;
    MPI_Alltoall(sc,1,MPI_INT,rc,1,MPI_INT,MPI_COMM_WORLD);
    for(int p=1;p<np;p++){sd[p]=sd[p-1]+sc[p-1];rd[p]=rd[p-1]+rc[p-1];}
    int nr=(rd[np-1]+rc[np-1])/W;
    memcpy(next,sd,np*sizeof(int));
    double *send=malloc((size_t)(nf+1)*W*sizeof(double)),*recv=malloc((size_t)(nr+1)*W*sizeof(double));
    for(int k=0;k<nf;k++){double *f=local+k*W;int p=owner[(int)f[1]];memcpy(send+next[p],f,W*sizeof(double));next[p]+=W;}
    MPI_Alltoallv(send,sc,sd,MPI_DOUBLE,recv,rc,rd,MPI_DOUBLE,MPI_COMM_WORLD);
    qsort(local,nf,W*sizeof(double),cmp);
    long long stats[6]={nf,nr,remote,0,0,0},total[6];double maxima[3]={0},global[3];
    for(int k=0;k<nr;k++){
        double *f=recv+k*W,*g=reciprocal(local,nf,f);
        if(f[21]!=1 || (g && g[21]!=1)){stats[3]++;continue;}
        if(!g){stats[4]++;continue;}
        double da=fabs(f[4]-g[4]),dd=fabs(f[5]-g[5]);
        double fx=-f[4]*f[14]*f[6],fy=-f[4]*f[14]*f[7];
        double gx=-g[4]*g[14]*g[6],gy=-g[4]*g[14]*g[7];
        double df=hypot(fx+gx,fy+gy),fs=hypot(fx,fy)+hypot(gx,gy);
        double h=.5*f[4]*f[14]*f[9],gh=.5*g[4]*g[14]*g[9],dh=fabs(h-gh);
        maxima[0]=fmax(maxima[0],da);maxima[1]=fmax(maxima[1],df);maxima[2]=fmax(maxima[2],dh);
        if(da>1e-5*(fabs(f[4])+fabs(g[4]))+1e-18 || dd>1e-5*(fabs(f[5])+fabs(g[5]))+1e-18 || df>1e-5*fs+1e-18 || dh>1e-5*(fabs(h)+fabs(gh))+1e-18)stats[5]++;
    }
    MPI_Allreduce(stats,total,6,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(maxima,global,3,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    if(!rank)printf("{\"ranks\":%d,\"sent\":%lld,\"received\":%lld,\"remote\":%lld,\"stale\":%lld,\"missing\":%lld,\"mismatch\":%lld,\"max_area_difference\":%.17g,\"max_force_residual\":%.17g,\"max_heat_difference\":%.17g}\n",np,total[0],total[1],total[2],total[3],total[4],total[5],global[0],global[1],global[2]);
    int bad=total[0]!=total[1] || total[3] || total[4] || total[5];
    free(owner);free(local);free(send);free(recv);free(sc);free(rc);free(sd);free(rd);free(next);
    MPI_Finalize();return bad?1:0;
}
