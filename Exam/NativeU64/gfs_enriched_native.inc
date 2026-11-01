/* Included after native Gauss-rule constructor. Opt-in read-only audit. */
#include "gfs_enriched_coupling_native.inc"
static void gfs_enriched_projection_output(SimParameters *simpar,int np,const GfsDgConfig *field){
    GfsEnrichConfig cfg[2]={0};
    int na=(int)gfs_moment_env("GFS_ENRICH_ANGULAR",8),nr=(int)gfs_moment_env("GFS_ENRICH_RADIAL",8);
    int ca=(int)gfs_moment_env("GFS_ENRICH_CHECK_ANGULAR",16),cr=(int)gfs_moment_env("GFS_ENRICH_CHECK_RADIAL",16);
    if(na<4 || na>GFS_ENRICH_NA || ca<4 || ca>GFS_ENRICH_NA ||
       nr<4 || nr>GFS_ENRICH_NR || cr<4 || cr>GFS_ENRICH_NR)MPI_Abort(MPI_COMM(simpar),99);
    for(int k=0;k<2;k++){
        cfg[k].field=*field;cfg[k].na=k?ca:na;cfg[k].nr=k?cr:nr;
        gfs_dg_gauss_rule(cfg[k].na,cfg[k].an,cfg[k].aw);
        gfs_dg_gauss_rule(cfg[k].nr,cfg[k].rn,cfg[k].rw);
    }
    GfsEnrichAudit *audit=calloc((size_t)(np?np:1),sizeof(*audit));if(!audit)MPI_Abort(MPI_COMM(simpar),99);
    int bad=gpu_enriched_project_initial(gpu_get_context(),np,cfg,audit),total=0;
    MPI_Allreduce(&bad,&total,1,MPI_INT,MPI_SUM,MPI_COMM(simpar));
    char path[160];snprintf(path,sizeof(path),"enriched_project_r%d.bin",MYID(simpar));
    FILE *out=fopen(path,"wx");if(!out)MPI_Abort(MPI_COMM(simpar),99);
    unsigned long long header[8]={0x474653454e523031ULL,2,(unsigned long long)np,GFS_ENRICH_N,ER_NMETRICS,1528,
        (unsigned long long)(na+256*nr),(unsigned long long)(ca+256*cr)};
    fwrite(header,sizeof(header),1,out);
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(simpar);
    int nonpositive[2]={0},all_nonpositive[2];
    for(int i=0;i<np;i++){
        unsigned long long id=(size_t)PINDX((treevorork4particletype*)(p+i));double xy[2]={p[i].x,p[i].y};
        GfsEnrichBasis *b=&audit[i].basis;
        fwrite(&id,sizeof(id),1,out);fwrite(xy,sizeof(double),2,out);
        fwrite(&b->geometry.area,sizeof(double),6,out);fwrite(b->generator,sizeof(double),2,out);
        fwrite(b->scale,sizeof(double),4,out);fwrite(b->L,sizeof(double),100,out);fwrite(b->mean,sizeof(double),10,out);
        fwrite(&b->minimum_norm,sizeof(double),4,out);int br[2]={b->rank,b->status};fwrite(br,sizeof(br),1,out);
        fwrite(audit[i].state.c,sizeof(double),40,out);fwrite(audit[i].integral,sizeof(double),4,out);
        fwrite(audit[i].metrics,sizeof(double),ER_NMETRICS,out);
        int flags[4]={audit[i].status,audit[i].nonpositive_rho,audit[i].nonpositive_p,0};fwrite(flags,sizeof(flags),1,out);
        if(audit[i].nonpositive_rho)nonpositive[0]++;if(audit[i].nonpositive_p)nonpositive[1]++;
    }
    if(ferror(out) || fclose(out))MPI_Abort(MPI_COMM(simpar),99);
    MPI_Allreduce(nonpositive,all_nonpositive,2,MPI_INT,MPI_SUM,MPI_COMM(simpar));
    if(total){fprintf(stderr,"[ENRICH_PROJECT] geometry/representation failed cells=%d\n",total);MPI_Abort(MPI_COMM(simpar),99);}
    const char *coupled=getenv("GFS_ENRICH_COUPLED_RHS");
    if(coupled && atoi(coupled)>0){
        if(all_nonpositive[0] || all_nonpositive[1])MPI_Abort(MPI_COMM(simpar),99);
        gfs_enriched_coupling_output(simpar,np,cfg+1,audit);
    }
    free(audit);
    if(MYID(simpar)==0)fprintf(stderr,"[ENRICH_PROJECT_ONLY] completed; nonpositive_rho_cells=%d nonpositive_pressure_cells=%d; no evolution, limiter or physics pass implied\n",
        all_nonpositive[0],all_nonpositive[1]);
    MPI_Barrier(MPI_COMM(simpar));MPI_Finalize();exit(0);
}
