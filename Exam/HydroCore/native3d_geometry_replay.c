/* Replay a captured native 3D geometry failure through the same adapter.
 * The captured coordinates are actual failed-stage positions, not jittered
 * substitutes. No hydrodynamic time step is taken. */
#include "../NativeU64/exam3d_gfs.c"
_Static_assert(sizeof(g3face)==80,"captured geometry face ABI");
static double mismatch(g3sim *s,double *ratio,double *centroid_moment_ratio){
    double max=0;*ratio=0;*centroid_moment_ratio=0;
    for(long i=0;i<s->n;i++){
        g3face *f=g3_faces(s,i);
        for(int k=0;k<s->fcnt[i];k++){
            int j=f[k].j;if(i>j)continue;
            double sh[3]={f[k].rel[0]-s->x[j]+s->x[i],f[k].rel[1]-s->y[j]+s->y[i],f[k].rel[2]-s->z[j]+s->z[i]};
            for(int a=0;a<3;a++)sh[a]=s->L[a]*nearbyint(sh[a]/s->L[a]);
            double reverse[3]={-sh[0],-sh[1],-sh[2]},Ai=0,Aj=0,Ci[3]={0},Cj[3]={0};
            for(int l=0;l<s->fcnt[i];l++)if(g3_u_match(s,i,f+l,j,sh)){
                double A=g3_u_area(f+l);Ai+=A;for(int a=0;a<3;a++)Ci[a]+=A*f[l].c[a];
            }
            g3face *other=g3_faces(s,j);
            for(int l=0;l<s->fcnt[j];l++)if(g3_u_match(s,j,other+l,i,reverse)){
                double A=g3_u_area(other+l);Aj+=A;for(int a=0;a<3;a++)Cj[a]+=A*(f[k].rel[a]+other[l].c[a]);
            }
            double h=fmax(cbrt(s->vol[i]),cbrt(s->vol[j])),error=fabs(Ai-Aj);
            max=fmax(max,error);*ratio=fmax(*ratio,error/(2048*DBL_EPSILON*h*h));
            double moment=0;for(int a=0;a<3;a++)moment+=pow(Ci[a]-Cj[a],2);
            *centroid_moment_ratio=fmax(*centroid_moment_ratio,sqrt(moment)/(4096*DBL_EPSILON*h*h*h));
        }
    }
    return max;
}
int main(int argc,char **argv){
    if(argc!=3)return 2;
    FILE *in=fopen(argv[1],"rb");uint64_t header[4];g3sim S={0},*s=&S;
    if(!in||fread(header,sizeof(header),1,in)!=1||header[0]!=UINT64_C(0x473347454f464149)||header[3]!=1||header[1]>10000000||header[2]>100000000)return 2;
    s->n=header[1];if(fread(&s->gamma,sizeof(double),1,in)!=1)return 2;
    g3_alloc_particles(s,s->n);s->unified=1;s->nthreads=1;s->tb=g3_alloc(sizeof(g3fbuf));
    s->tb[0].cap=header[2];s->tb[0].f=g3_alloc(header[2]*sizeof(g3face));
    for(int a=0;a<3;a++){s->bc[a]=BC_PERIODIC;s->L[a]=s->hi[a]=atof(argv[2]);if(!(s->L[a]>0))return 2;}
    for(long i=0;i<s->n;i++){
        double q[9];int nf;
        if(fread(q,sizeof(q),1,in)!=1||fread(s->id+i,sizeof(int64_t),1,in)!=1||fread(&nf,sizeof(int),1,in)!=1||nf<0||s->tb[0].n+nf>s->tb[0].cap)return 2;
        s->x[i]=q[0];s->y[i]=q[1];s->z[i]=q[2];s->vx[i]=q[3];s->vy[i]=q[4];s->vz[i]=q[5];s->m[i]=q[6];s->ie[i]=q[7];s->vol[i]=q[8];
        s->fcnt[i]=nf;s->foff[i]=s->tb[0].n;
        if(fread(s->tb[0].f+s->tb[0].n,sizeof(g3face),nf,in)!=(size_t)nf)return 2;s->tb[0].n+=nf;
    }
    if(fgetc(in)!=EOF)return 2;fclose(in);
    double old_ratio,new_ratio,old_centroid,new_centroid,old=mismatch(s,&old_ratio,&old_centroid);
    printf("{\"scope\":\"captured original geometry\",\"max_area_mismatch\":%.17g,\"max_over_roundoff_bound\":%.17g,\"centroid_moment_over_bound\":%.17g}\n",old,old_ratio,old_centroid);
    if(g3_tessellate(s))return 3;
    double corrected=mismatch(s,&new_ratio,&new_centroid);int audit=g3_u_geometry_audit(s);
    printf("{\"scope\":\"same actual coordinates, strict half spaces and intersections\",\"max_area_mismatch\":%.17g,\"max_over_roundoff_bound\":%.17g,\"centroid_moment_over_bound\":%.17g,\"geometry_audit_exit\":%d,\"jitter\":%ld}\n",corrected,new_ratio,new_centroid,audit,s->n_jitter);
    return new_ratio>1||new_centroid>1||audit?1:0;
}
