/* Create a new Q2 comparison IC from a preserved full-stride checkpoint. */
#include <errno.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "eunha.h"
#include "Exam/gfs_moment.h"
#include "Exam/gfs_thermal64.h"

static void fail(const char *s){perror(s);exit(1);}
int main(int argc,char **argv){
    if(argc!=4){fprintf(stderr,"input output scale\n");return 2;}
    char *end;double scale=strtod(argv[3],&end);
    if(*end || !(scale>0) || !isfinite(scale))return 2;
    if(sizeof(treevorostressrk4particletype)!=576 || sizeof(PosType)!=8 ||
       offsetof(treevorostressrk4particletype,ie)!=116 ||
       offsetof(treevorostressrk4particletype,pressure)!=108) return 3;
    FILE *in=fopen(argv[1],"rb");if(!in)fail("input");
    int h[4];if(fread(h,sizeof(h),1,in)!=1 || h[0]!=-1 || h[1]<=0 || h[2]!=576 || h[3]!=1)return 4;
    size_t n=(size_t)h[1];treevorostressrk4particletype *p=malloc(n*sizeof(*p));if(!p)fail("allocation");
    if(fread(p,sizeof(*p),n,in)!=n)return 5;
    unsigned char footer[52];if(fread(footer,sizeof(footer),1,in)!=1 || fgetc(in)!=EOF)return 5;
    fclose(in);
    double time;memcpy(&time,footer,sizeof(time));if(time!=30.0)return 6;
    for(size_t i=0;i<n;i++){
        GfsMomentPayload moment;memcpy(&moment,&p[i].rk4,sizeof(moment));
        if(GFS_T64_TAGGED(p[i]) || moment.tag==GFS_MOMENT_TAG || p[i].w2!=0 ||
           !(p[i].ie>0) || !(p[i].pressure>0) || !(p[i].csound>0))return 7;
        if(scale!=1.){
            p[i].ie*=scale;p[i].pressure*=scale;p[i].csound*=sqrt(scale);
            p[i].stress.K*=scale;p[i].stress.dPdx*=scale;p[i].stress.dPdy*=scale;
            if(!(p[i].ie>0) || !isfinite(p[i].ie) || !(p[i].pressure>0) ||
               !isfinite(p[i].pressure) || !(p[i].csound>0) || !isfinite(p[i].csound))return 8;
        }
    }
    FILE *out=fopen(argv[2],"wbx");if(!out)fail("output");
    if(fwrite(h,sizeof(h),1,out)!=1 || fwrite(p,sizeof(*p),n,out)!=n ||
       fwrite(footer,sizeof(footer),1,out)!=1 || fclose(out))fail("write");
    free(p);printf("scaled %zu particles at t=30 by %.17g; offsets ie=%zu pressure=%zu csound=%zu K=%zu dPdx=%zu dPdy=%zu\n",
        n,scale,offsetof(treevorostressrk4particletype,ie),offsetof(treevorostressrk4particletype,pressure),
        offsetof(treevorostressrk4particletype,csound),offsetof(treevorostressrk4particletype,stress.K),
        offsetof(treevorostressrk4particletype,stress.dPdx),offsetof(treevorostressrk4particletype,stress.dPdy));return 0;
}
