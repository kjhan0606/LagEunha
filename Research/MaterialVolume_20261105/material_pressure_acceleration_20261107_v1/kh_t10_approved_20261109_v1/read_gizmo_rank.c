/* Read-only exporter for this archived, single-rank, gas-only GIZMO build.
 * Layout is taken from its own allvars.h, not independently guessed offsets.
 * Restart stream order is All, NumPart, P[], N_gas, SphP[]. */
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include "allvars.h"

static void read_exact(FILE *f, void *p, size_t size, size_t n)
{
    if(fread(p, size, n, f)!=n) {fprintf(stderr,"truncated restart\n"); exit(2);}
}
int main(int argc, char **argv)
{
    if(argc!=3 && argc!=4) return 2;
    const double target=argc==4?strtod(argv[3],NULL):1.35;
    if(!isfinite(target) || target<=0) return 2;
    FILE *f=fopen(argv[1],"rb"); if(!f) return 2;
    struct global_data_all_processes a;
    int n, ng;
    read_exact(f,&a,sizeof(a),1); read_exact(f,&n,sizeof(n),1);
    fprintf(stderr,"restart time=%.17g target=%.17g Ti=%lld TIMEBASE=%lld n=%d sizes=%zu,%zu,%zu\n",
            a.Time,target,(long long)a.Ti_Current,(long long)TIMEBASE,n,sizeof(a),sizeof(struct particle_data),sizeof(struct gas_cell_data));
    if(sizeof(struct particle_data)!=288 || sizeof(struct gas_cell_data)!=312 ||
       n<1 || n>1000000 || a.TotNumPart<n || a.TotN_gas!=a.TotNumPart ||
       fabs(a.Time-target)>1e-12 ||
       !(a.Ti_Current>=TIMEBASE || a.Time>a.TimeMax) || a.ComovingIntegrationOn) return 3;
    struct particle_data *p=malloc((size_t)n*sizeof(*p));
    struct gas_cell_data *s=malloc((size_t)n*sizeof(*s));
    if(!p || !s) return 4;
    read_exact(f,p,sizeof(*p),n); read_exact(f,&ng,sizeof(ng),1);
    if(ng!=n) return 3;
    read_exact(f,s,sizeof(*s),n); fclose(f);
    FILE *out=fopen(argv[2],"wx"); if(!out) return 5;
    for(int i=0;i<n;i++) {
        if(p[i].Type!=0 || p[i].Ti_current!=a.Ti_Current) return 6;
        double row[10]={p[i].Pos[0],p[i].Pos[1],p[i].Vel[0],p[i].Vel[1],
                       p[i].Mass,s[i].InternalEnergy,s[i].Density,(double)p[i].ID,
                       s[i].InternalEnergyPred,(double)p[i].TimeBin};
        for(int k=0;k<10;k++) if(!isfinite(row[k])) return 7;
        if(fwrite(row,sizeof(double),10,out)!=10) return 8;
    }
    fclose(out);
    printf("{\"time\":%.17g,\"cells\":%d,\"global_bytes\":%zu,\"particle_bytes\":%zu,\"gas_bytes\":%zu,\"Ti_Current\":%lld}\n",
           a.Time,n,sizeof(a),sizeof(*p),sizeof(*s),(long long)a.Ti_Current);
    free(p); free(s); return 0;
}
