/* Native MPI allgather is limited to this initial RHS audit. Future runtime
 * must carry the evolving full state through normal migration/halo/restart.
 */
typedef struct {unsigned long long id;int old;} GfsCoupledSort;
static int gfs_ec_compare(const void *a,const void *b){
    unsigned long long x=((const GfsCoupledSort*)a)->id,y=((const GfsCoupledSort*)b)->id;
    return x<y?-1:(x>y?1:0);
}
static void gfs_ec_mpi_sum(void *user,double *v,int n){
    SimParameters *sp=user;MPI_Allreduce(MPI_IN_PLACE,v,n,MPI_DOUBLE,MPI_SUM,MPI_COMM(sp));
}
static void gfs_enriched_coupling_output(SimParameters *sp,int np,
    const GfsEnrichConfig *cfg,const GfsEnrichAudit *local){
    int nr=NID(sp),*counts=calloc(nr,sizeof(int)),*offset=calloc(nr,sizeof(int));
    int *bytes=calloc(nr,sizeof(int)),*b_offset=calloc(nr,sizeof(int));
    if(!counts || !offset || !bytes || !b_offset)MPI_Abort(MPI_COMM(sp),99);
    MPI_Allgather(&np,1,MPI_INT,counts,1,MPI_INT,MPI_COMM(sp));size_t n=0;
    for(int r=0;r<nr;r++){
        if((n+(size_t)counts[r])*sizeof(GfsEnrichAudit)>2147483647ULL || (n+(size_t)counts[r])*40>2147483647ULL)MPI_Abort(MPI_COMM(sp),99);
        offset[r]=(int)n;b_offset[r]=(int)(n*sizeof(GfsEnrichAudit));bytes[r]=(int)(counts[r]*sizeof(GfsEnrichAudit));n+=counts[r];
    }
    unsigned long long *local_ids=malloc((size_t)(np?np:1)*sizeof(*local_ids));
    unsigned long long *ids=malloc(n*sizeof(*ids)),*ordered_ids=malloc(n*sizeof(*ordered_ids));
    GfsEnrichAudit *all=malloc(n*sizeof(*all)),*ordered=malloc(n*sizeof(*ordered));
    GfsCoupledSort *order=malloc(n*sizeof(*order));
    double *z=calloc(n,sizeof(double)),*rhs=calloc(n*40,sizeof(double)),*geom=calloc(n*40,sizeof(double));
    if(!local_ids || !ids || !ordered_ids || !all || !ordered || !order || !z || !rhs || !geom)MPI_Abort(MPI_COMM(sp),99);
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(sp);
    for(int i=0;i<np;i++)local_ids[i]=(size_t)PINDX((treevorork4particletype*)(p+i));
    MPI_Allgatherv(local_ids,np,MPI_UNSIGNED_LONG_LONG,ids,counts,offset,MPI_UNSIGNED_LONG_LONG,MPI_COMM(sp));
    MPI_Allgatherv(local,np*(int)sizeof(*local),MPI_BYTE,all,bytes,b_offset,MPI_BYTE,MPI_COMM(sp));
    for(size_t i=0;i<n;i++){order[i].id=ids[i];order[i].old=(int)i;}
    qsort(order,n,sizeof(*order),gfs_ec_compare);
    for(size_t i=0;i<n;i++){
        if(i && order[i].id==order[i-1].id)MPI_Abort(MPI_COMM(sp),99);
        ordered_ids[i]=order[i].id;ordered[i]=all[order[i].old];
    }
    double perturb[2]={0,0};const char *s=getenv("GFS_ENRICH_PERTURB");
    if(s && *s){
        char *end=NULL;perturb[0]=strtod(s,&end);if(end==s || *end!=',')MPI_Abort(MPI_COMM(sp),99);
        const char *second=end+1;perturb[1]=strtod(second,&end);
        if(end==second || *end || !isfinite(perturb[0]) || !isfinite(perturb[1]))MPI_Abort(MPI_COMM(sp),99);
    }
    double report[EC_NREPORT]={0};
    if(gpu_enriched_coupled_rhs(gpu_get_context(),np,(int)n,MYID(sp),ordered_ids,ordered,cfg,perturb,
        gfs_ec_mpi_sum,sp,z,rhs,geom,report)){
        fprintf(stderr,"[ENRICH_COUPLED] native Euler RHS/weight-rate solve failed\n");MPI_Abort(MPI_COMM(sp),99);
    }
    double mass_link=0,max_mdot=0;
    for(size_t i=0;i<n;i++){
        double mu=ordered[i].integral[0]/ordered[i].basis.geometry.area;
        mass_link=fmax(mass_link,fabs(rhs[i*40]+mu*z[i])/mu);
        max_mdot=fmax(max_mdot,fabs(rhs[i*40])/ordered[i].integral[0]);
    }
    if(MYID(sp)==0){
        FILE *out=fopen("enriched_coupled.bin","wx");if(!out)MPI_Abort(MPI_COMM(sp),99);
        unsigned long long h[8]={0x4746534543523031ULL,1,n,10,4,688,(unsigned long long)(cfg->na+256*cfg->nr),0};
        fwrite(h,sizeof(h),1,out);fwrite(perturb,sizeof(double),2,out);
        for(size_t i=0;i<n;i++){
            double mu=ordered[i].integral[0]/ordered[i].basis.geometry.area;
            fwrite(ordered_ids+i,sizeof(*ordered_ids),1,out);fwrite(ordered[i].basis.generator,sizeof(double),2,out);
            fwrite(&mu,sizeof(double),1,out);fwrite(ordered[i].integral,sizeof(double),1,out);fwrite(z+i,sizeof(double),1,out);
            fwrite(rhs+i*40,sizeof(double),40,out);fwrite(geom+i*40,sizeof(double),40,out);
        }
        if(ferror(out) || fclose(out))MPI_Abort(MPI_COMM(sp),99);
        out=fopen("enriched_coupled.json","wx");if(!out)MPI_Abort(MPI_COMM(sp),99);
        const char *name[EC_NREPORT]={"quadrature_points","canonical_pairs","Newton_iterations","CG_iterations","backtracks",
            "scaled_root_residual","total_mass_rate","weighted_gauge_rate","Jacobian_symmetry_error","Rayleigh_over_mu_bound",
            "secant_over_mu_bound","HLL_points","max_normalized_rhs_minus_geometry"};
        fprintf(out,"{\"scope\":\"initial native moving-cell full Euler moment RHS and screened weight-rate solve; no evolution/dynamic-weight geometry/restart\",\"cells\":%zu,\"ranks\":%d,\"perturb\":[%.17g,%.17g],\"mass_weight_rate_error\":%.17g,\"max_abs_mdot_over_mass\":%.17g",
            n,nr,perturb[0],perturb[1],mass_link,max_mdot);
        for(int k=0;k<EC_NREPORT;k++)fprintf(out,",\"%s\":%.17g",name[k],report[k]);fprintf(out,"}\n");
        if(ferror(out) || fclose(out))MPI_Abort(MPI_COMM(sp),99);
        fprintf(stderr,"[ENRICH_COUPLED] cells=%zu Newton=%g CG=%g mass_weight_error=%.9g symmetry=%.9g Rayleigh=%.9g; no time advancement\n",
            n,report[EC_NEWTON],report[EC_CG],mass_link,report[EC_SYMMETRY],report[EC_RAYLEIGH]);
    }
    free(counts);free(offset);free(bytes);free(b_offset);free(local_ids);free(ids);free(ordered_ids);
    free(all);free(ordered);free(order);free(z);free(rhs);free(geom);
}
