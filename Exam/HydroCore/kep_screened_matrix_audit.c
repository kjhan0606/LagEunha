/* Read-only dimensional native RHS/matrix snapshot. Historical states are
 * operator inputs only and are never converted into new trajectories. */
#define main archived_native_main
#include "fv_native_evolution.c"
#undef main
#include "../gfs_fv_native.c"
#include <unistd.h>
#ifdef GFS_NATIVE_IC0
static int ic0_algebra(GfsFvWorkspace *w,const GfsFvMesh *m,const char *prefix){
    double *x=calloc(m->n,sizeof(double)),*y=calloc(m->n,sizeof(double)),*ax=calloc(m->n,sizeof(double)),*csr=calloc(m->n,sizeof(double));
    double *inverse_x=calloc(m->n,sizeof(double)),*inverse_y=calloc(m->n,sizeof(double));
    if(!x||!y||!ax||!csr||!inverse_x||!inverse_y)return 4;
    if(!w->ic_applications&&prepare_ic0(w,m))return 3;
    double matrix_error=0,matrix_backward_error=0,constant_mode_error=0;
    for(int mode=0;mode<2;mode++){
        for(int i=0;i<m->n;i++){x[i]=mode?1:sin(1.7*(i+1));csr[i]=w->diagonal[i]*x[i];y[i]=fabs(csr[i]);}
        multiply(w,m,x,ax);double scale=1;
        for(int i=0;i<m->n;i++)for(int p=w->ic_offset[i];p<w->ic_offset[i+1];p++){
            int j=w->ic_column[p];double a=w->ic_matrix[p];csr[i]+=a*x[j];csr[j]+=a*x[i];
            y[i]+=fabs(a*x[j]);y[j]+=fabs(a*x[i]);
        }
        for(int i=0;i<m->n;i++)scale=fmax(scale,fabs(ax[i]));
        for(int i=0;i<m->n;i++){
            matrix_error=fmax(matrix_error,fabs(ax[i]-csr[i])/scale);
            matrix_backward_error=fmax(matrix_backward_error,fabs(ax[i]-csr[i])/(1+y[i]));
            if(mode)constant_mode_error=fmax(constant_mode_error,fabs(ax[i]-m->cell[i].mu)/(1+m->cell[i].mu));
        }
    }
    for(int i=0;i<m->n;i++){x[i]=sin(1.7*(i+1));y[i]=cos(.83*(i+.4));w->r[i]=x[i];}
    if(ic0_precondition(w,m))return 3;memcpy(inverse_x,w->pre,m->n*sizeof(double));
    memcpy(w->r,y,m->n*sizeof(double));if(ic0_precondition(w,m))return 3;memcpy(inverse_y,w->pre,m->n*sizeof(double));
    double xMy=0,yMx=0,xMx=0;for(int i=0;i<m->n;i++){xMy+=x[i]*inverse_y[i];yMx+=y[i]*inverse_x[i];xMx+=x[i]*inverse_x[i];}
    double symmetry=fabs(xMy-yMx)/(1+fabs(xMy)+fabs(yMx));
    char path[4096];snprintf(path,sizeof(path),"%s_factor.csv",prefix);FILE *out=fopen(path,"wx");if(!out)return 4;
    fprintf(out,"i,j,Aij,Lij\n");
    for(int i=0;i<m->n;i++)for(int p=w->ic_offset[i];p<w->ic_offset[i+1];p++)fprintf(out,"%d,%d,%.17g,%.17g\n",i,w->ic_column[p],w->ic_matrix[p],w->ic_lower[p]);
    if(fclose(out))return 4;snprintf(path,sizeof(path),"%s_pivots.csv",prefix);out=fopen(path,"wx");if(!out)return 4;
    fprintf(out,"i,Aii,Di\n");for(int i=0;i<m->n;i++)fprintf(out,"%d,%.17g,%.17g\n",i,w->diagonal[i],w->ic_pivot[i]);if(fclose(out))return 4;
    /* A constant probe cancels O(a_ij) CSR terms to O(mu). The forward
     * discrepancy then scales with a_ij/mu even for the same coefficients.
     * Keep that observation and test coefficient-scale backward error plus
     * the native exact constant mode. This is not a root stopping rule. */
    int pass=matrix_backward_error<1e-12&&constant_mode_error<1e-12&&symmetry<1e-12&&xMx>0&&isfinite(xMx)&&w->ic_min_pivot>0;
    snprintf(path,sizeof(path),"%s_algebra.json",prefix);out=fopen(path,"wx");if(!out)return 4;
    fprintf(out,"{\"scope\":\"original native edge matvec vs assembled symmetric CSR, IC0 preconditioner bilinear probes\",\"matrix_action_relative_error\":%.17g,\"matrix_action_backward_error\":%.17g,\"native_constant_mode_relative_error\":%.17g,\"preconditioner_symmetry_relative_error\":%.17g,\"positive_quadratic_form\":%.17g,\"minimum_pivot\":%.17g,\"stored_lower_plus_diagonal_entries\":%d,\"shift_or_fill_or_reordering\":false,\"pass\":%s}\n",matrix_error,matrix_backward_error,constant_mode_error,symmetry,xMx,w->ic_min_pivot,w->ic_nnz+m->n,pass?"true":"false");
    if(fclose(out))return 4;free(x);free(y);free(ax);free(csr);free(inverse_x);free(inverse_y);return pass?0:1;
}
#endif
#ifndef GFS_MATRIX_AUDIT_ENTRY
#define GFS_MATRIX_AUDIT_ENTRY main
#endif
int GFS_MATRIX_AUDIT_ENTRY(int argc,char **argv){
    if(argc!=3&&argc!=4)return 2;int repetitions=argc==4?atoi(argv[3]):0;if(argc==4&&(repetitions<1||repetitions>10000))return 2;
    FILE *in=fopen(argv[1],"rb");FvCheckpoint h;
    if(!in||fread(&h,sizeof(h),1,in)!=1||h.magic!=magic||h.cell_bytes!=sizeof(GfsFvCell)||h.count<2||h.count>1000000||h.dim<1||h.dim>3||(h.version&~diagnostic_trial_flag)<1||(h.version&~diagnostic_trial_flag)>8)return 2;
    GfsFvMesh m={0};m.dim=h.dim;m.n=h.count;m.gamma=h.gamma;memcpy(m.bc,h.bc,sizeof(m.bc));memcpy(m.lo,h.lo,sizeof(m.lo));memcpy(m.L,h.L,sizeof(m.L));
    m.cell=calloc(m.n,sizeof(*m.cell));if(!m.cell||fread(m.cell,sizeof(*m.cell),m.n,in)!=(size_t)m.n||fgetc(in)!=EOF||fclose(in))return 2;
    GfsFvCell *copy=malloc(m.n*sizeof(*copy));memcpy(copy,m.cell,m.n*sizeof(*copy));
    double mark=omp_get_wtime();if(gfs_fv_geometry(&m))return 3;double geometry_seconds=omp_get_wtime()-mark;
    GfsFvWorkspace *w=gfs_fv_workspace_create();GfsFvRate *rate=calloc(m.n,sizeof(*rate));GfsFvRhsStats stats;
    mark=omp_get_wtime();int code=gfs_fv_native_rhs(w,&m,rate,&stats);double rhs_seconds=omp_get_wtime()-mark;
    if(code==1||code==2)return 3; /* No screened system exists for a failed reconstruction/fan. */
    double *zero=calloc(m.n,sizeof(double)),*R0=calloc(m.n,sizeof(double));if(residual(w,&m,zero,R0,0))return 3;
    char path[4096];snprintf(path,sizeof(path),"%s_matrix.csv",argv[2]);FILE *out=fopen(path,"wx");if(!out)return 4;
    fprintf(out,"i,j,aij\n");for(int k=0;k<w->ne;k++){FvEdge *e=w->edge+k;if(e->j>=0)fprintf(out,"%d,%d,%.17g\n",e->i,e->j,e->jacobian);}
    if(fclose(out))return 4;
    snprintf(path,sizeof(path),"%s_rows.csv",argv[2]);out=fopen(path,"wx");if(!out)return 4;
    fprintf(out,"id,mu,diagonal,b,scale,z_candidate,volume,qx,qy,qz,qdotx,qdoty,qdotz,z_low\n");
    const double *candidate=code>=10?w->direction:w->z;double maximum=0;
    for(int i=0;i<m.n;i++){
        maximum=fmax(maximum,w->diagonal[i]/m.cell[i].mu);
        const double *low=solution_low(w,candidate);
        fprintf(out,"%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",i,m.cell[i].mu,w->diagonal[i],-R0[i],w->scale[i],candidate[i],m.geometry[i].V,m.cell[i].q[0],m.cell[i].q[1],m.cell[i].q[2],w->velocity[i][0],w->velocity[i][1],w->velocity[i][2],low?low[i]:0);
    }if(fclose(out))return 4;
    if(residual(w,&m,candidate,w->trial,0))return 3;multiply(w,&m,candidate,w->Ap);
    double representation=0;for(int i=0;i<m.n;i++)representation=fmax(representation,fabs(w->trial[i]-w->Ap[i]-R0[i])/w->scale[i]);
    double initial_true_norm=norm(w,w->trial);
#ifdef GFS_NATIVE_IC0
    if(ic0_algebra(w,&m,argv[2]))return 3;
#endif
    if(repetitions){
        double solve_time=0,whole_rhs=0,assembly=0,factor=0,apply=0,max_true=0;long long iterations=0;int failures=0;
        FILE *quiet=fopen("/dev/null","w");int saved=dup(STDERR_FILENO);if(!quiet||saved<0)return 4;fflush(stderr);if(dup2(fileno(quiet),STDERR_FILENO)<0)return 4;fclose(quiet);
        for(int k=0;k<repetitions;k++){
            GfsFvRhsStats sample;double t=omp_get_wtime();int status=gfs_fv_native_rhs(w,&m,rate,&sample);whole_rhs+=omp_get_wtime()-t;
            failures+=status!=0;solve_time+=sample.solve_seconds;assembly+=sample.preconditioner_assembly_seconds;factor+=sample.preconditioner_factor_seconds;apply+=sample.preconditioner_apply_seconds;iterations+=sample.cg;max_true=fmax(max_true,sample.root_residual);
        }
        fflush(stderr);if(dup2(saved,STDERR_FILENO)<0)return 4;close(saved);
        snprintf(path,sizeof(path),"%s_benchmark.json",argv[2]);out=fopen(path,"wx");if(!out)return 4;
        fprintf(out,"{\"scope\":\"frozen repeated RHS after warm allocation; geometry and I/O excluded, symbolic and numeric IC0 rebuilt every solve; passive PCG history disabled\",\"repetitions\":%d,\"failures\":%d,\"mean_solve_seconds_including_all_setup_apply_and_true_checks\":%.17g,\"mean_whole_RHS_seconds\":%.17g,\"mean_assembly_and_symbolic_seconds\":%.17g,\"mean_numeric_factor_seconds\":%.17g,\"mean_triangular_apply_seconds\":%.17g,\"mean_iterations\":%.17g,\"maximum_true_scaled_residual\":%.17g,\"OMP_threads\":%d}\n",repetitions,failures,solve_time/repetitions,whole_rhs/repetitions,assembly/repetitions,factor/repetitions,apply/repetitions,(double)iterations/repetitions,max_true,omp_get_max_threads());
        if(fclose(out))return 4;
    }
    int unchanged=!memcmp(copy,m.cell,m.n*sizeof(*copy));
    double pair=fabs(stats.kinetic_rate-stats.kinetic_pressure_work+stats.kinetic_boundary_flux);
    printf("{\"scope\":\"readonly historical state under selected dimensional RHS; no trajectory\",\"source_checkpoint_version\":%llu,\"operator_checkpoint_version\":%llu,\"source_time\":%.17g,\"dim\":%d,\"cells\":%d,\"native_RHS_status\":%d,\"PCG_iterations\":%d,\"true_scaled_residual\":%.17g,\"native_flux_vs_matrix_representation_error\":%.17g,\"maximum_diagonal_over_mu\":%.17g,\"full_RHS_ledgers_available\":%s,\"kinetic_pair_defect\":%.17g,\"entropy_balance_defect\":%.17g,\"geometry_seconds\":%.17g,\"first_RHS_seconds_including_workspace_allocation\":%.17g,\"solve_seconds\":%.17g,\"conserved_states_and_birth_parameters_unchanged\":%s,\"accepted_evolution_steps\":0}\n",(unsigned long long)(h.version&~diagnostic_trial_flag),(unsigned long long)checkpoint_version,h.time,m.dim,m.n,code,stats.cg,initial_true_norm,representation,maximum,code?"false":"true",pair,stats.entropy_balance_defect,geometry_seconds,rhs_seconds,stats.solve_seconds,unchanged?"true":"false");
    free(zero);free(R0);free(copy);free(rate);gfs_fv_workspace_free(w);gfs_fv_geometry_free(&m);free(m.cell);return unchanged&&representation<1e-10?0:1;
}
