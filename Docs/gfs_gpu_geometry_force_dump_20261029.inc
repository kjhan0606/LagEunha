/* Read-only native validation output. At most four force calls per rank;
 * normal runs allocate/write nothing. Fixed-size arrays use native little
 * endian double/float/int32/int64 layout, validated by the Python reader. */
static void gfs_dump_geometry_force(SimParameters *simpar, GPUContext *ctx,
                                   int nbp, int n_total, int n_faces, int call)
{
    const char *setting=getenv("GFS_GEOM_FORCE_DUMP_MAXCALL");
    int maximum=setting ? atoi(setting) : 0;
    if(maximum<=0 || call>maximum) return;
    if(maximum>4 || sizeof(int)!=4 || sizeof(long long)!=8 ||
       sizeof(double)!=8 || sizeof(float)!=4) {
        fprintf(stderr,"[GEOM_FORCE_DUMP] unsupported count or native layout\n");
        MPI_Abort(MPI_COMM(simpar),94);
    }
    size_t nf=(size_t)n_faces;
    int *offset=(int*)malloc(((size_t)nbp+1)*sizeof(int));
    double *corners=(double*)malloc(4*nf*sizeof(double));
    int *indices=(int*)malloc(4*nf*sizeof(int));
    if(!offset || !corners || !indices) MPI_Abort(MPI_COMM(simpar),94);
    gpu_download_face_csr(ctx,nbp,n_faces,offset,
        corners,corners+nf,corners+2*nf,corners+3*nf,
        indices,indices+nf,indices+2*nf,indices+3*nf);
    char name[128];
    snprintf(name,sizeof(name),"geom_force_rank%03d_call%03d.bin",MYID(simpar),call);
    FILE *out=fopen(name,"wbx");
    if(!out) {
        fprintf(stderr,"[GEOM_FORCE_DUMP] cannot create %s\n",name);
        MPI_Abort(MPI_COMM(simpar),94);
    }
    const char magic[8]={'G','F','S','F','G','0','1',0};
    unsigned long long header[6]={(unsigned long long)nbp,(unsigned long long)n_total,
        (unsigned long long)n_faces,(unsigned long long)call,
        (unsigned long long)MYID(simpar),1ull};
    int okay=fwrite(magic,1,8,out)==8 && fwrite(header,sizeof(header),1,out)==1;
#define GFS_GF_WRITE(ptr,type,count) do { size_t n_=(size_t)(count); \
    if(fwrite((ptr),sizeof(type),n_,out)!=n_) okay=0; } while(0)
    GFS_GF_WRITE(ctx->h_parts.indx,long long,n_total);
    GFS_GF_WRITE(ctx->h_parts.x,double,n_total);
    GFS_GF_WRITE(ctx->h_parts.y,double,n_total);
    GFS_GF_WRITE(offset,int,nbp+1);
    GFS_GF_WRITE(corners,double,4*nf);
    GFS_GF_WRITE(indices,int,4*nf);
    GFS_GF_WRITE(ctx->h_parts.ax_out,double,nbp);
    GFS_GF_WRITE(ctx->h_parts.ay_out,double,nbp);
    GFS_GF_WRITE(ctx->h_parts.die_out,float,nbp);
    GFS_GF_WRITE(ctx->h_parts.dK_out,float,nbp);
    GFS_GF_WRITE(ctx->h_parts.dt_out,float,nbp);
#undef GFS_GF_WRITE
    if(fclose(out)!=0) okay=0;
    free(offset);free(corners);free(indices);
    if(!okay) MPI_Abort(MPI_COMM(simpar),94);
    fprintf(stderr,"[GEOM_FORCE_DUMP] rank=%d call=%d np=%d faces=%d file=%s\n",
            MYID(simpar),call,nbp,n_faces,name);
}
