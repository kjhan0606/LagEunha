/* Exact failed-system precision controls. No accepted evolution step. */
#define GFS_MATRIX_AUDIT_ENTRY archived_matrix_entry
#include "kep_screened_matrix_audit.c"
#undef GFS_MATRIX_AUDIT_ENTRY
#ifndef GFS_NATIVE_COMPENSATED_SOLUTION
#error This precision audit needs the compensated solution backend
#endif
int main(int argc,char **argv){
    if(argc!=4)return 2;
    FILE *in=fopen(argv[1],"rb");FvCheckpoint h;
    if(!in||fread(&h,sizeof(h),1,in)!=1||h.magic!=magic||h.cell_bytes!=sizeof(GfsFvCell)||
       (h.version&~diagnostic_trial_flag)!=8||h.count<2||h.count>1000000)return 2;
    GfsFvMesh m={0};m.dim=h.dim;m.n=h.count;m.gamma=h.gamma;memcpy(m.bc,h.bc,sizeof(m.bc));memcpy(m.lo,h.lo,sizeof(m.lo));memcpy(m.L,h.L,sizeof(m.L));
    m.cell=calloc(m.n,sizeof(*m.cell));if(!m.cell||fread(m.cell,sizeof(*m.cell),m.n,in)!=(size_t)m.n||fgetc(in)!=EOF||fclose(in)||gfs_fv_geometry(&m))return 2;
    GfsFvCell *saved=malloc(m.n*sizeof(*saved));memcpy(saved,m.cell,m.n*sizeof(*saved));
    GfsFvWorkspace *w=gfs_fv_workspace_create();GfsFvRate *rate=calloc(m.n,sizeof(*rate));GfsFvRhsStats stats;
    int status=gfs_fv_native_rhs(w,&m,rate,&stats);if(status)return 3;
    double conservation=0;for(int c=0;c<5;c++){
        long double value=stats.boundary_flux[c],scale=1;
        for(int i=0;i<m.n;i++){value+=rate[i].dQ[c];scale+=fabsl(rate[i].dQ[c]);}
        conservation=fmax(conservation,(double)(fabsl(value)/scale));
    }
    double *high=calloc(m.n,sizeof(double)),*low=calloc(m.n,sizeof(double)),*R=calloc(m.n,sizeof(double));
    FILE *ref=fopen(argv[3],"rb");if(!ref||!high||!low||!R||fread(high,sizeof(double),m.n,ref)!=(size_t)m.n||fread(low,sizeof(double),m.n,ref)!=(size_t)m.n||fgetc(ref)!=EOF||fclose(ref))return 4;
    if(residual_parts(w,&m,high,low,R,0))return 4;double independent=norm(w,R),forward_error=0;
    for(int i=0;i<m.n;i++){
        FvWide delta=wide_add((FvWide){w->z[i],w->z_low[i]},(FvWide){-high[i],-low[i]});
        forward_error=fmax(forward_error,fabs(delta.hi+delta.lo)/(1+fabs(high[i])));
    }
    if(residual_parts(w,&m,high,NULL,R,0))return 4;double reference_low_removed=norm(w,R);
    if(residual_parts(w,&m,w->z,NULL,R,0))return 4;double solver_low_removed=norm(w,R);
    high[m.n/4]+=1e-6;if(residual_parts(w,&m,high,low,R,0))return 4;double corrupted=norm(w,R);high[m.n/4]-=1e-6;
    low[0]=NAN;int rejects_nan=residual_parts(w,&m,high,low,R,0)||!isfinite(norm(w,R));
    low[0]=INFINITY;int rejects_inf=residual_parts(w,&m,high,low,R,0)||!isfinite(norm(w,R));
    FvWide sum=wide_sum(1,ldexp(1,-54)),product=wide_product(1+DBL_EPSILON,1-DBL_EPSILON);
    int arithmetic=sum.hi==1&&sum.lo==ldexp(1,-54)&&product.hi==1&&product.lo==-DBL_EPSILON*DBL_EPSILON;
    int unchanged=!memcmp(saved,m.cell,m.n*sizeof(*saved));
    double kinetic=fabs(stats.kinetic_rate-stats.kinetic_pressure_work+stats.kinetic_boundary_flux);
    int pass=stats.root_residual<=1e-12&&independent<=1e-12&&reference_low_removed>1e-12&&solver_low_removed>1e-12&&
             corrupted>1e-12&&rejects_nan&&rejects_inf&&arithmetic&&unchanged&&conservation<1e-12&&kinetic<1e-10&&fabs(stats.entropy_balance_defect)<1e-10&&forward_error<1e-15;
    char path[4096];snprintf(path,sizeof(path),"%s_precision.json",argv[2]);FILE *out=fopen(path,"wx");if(!out)return 4;
    fprintf(out,"{\"scope\":\"frozen exact late-blast precision and native RHS controls, no evolution\",\"pass\":%s,\"actual_PCG_iterations\":%d,\"native_true_residual\":%.17g,\"independent_60digit_hi_lo_native_residual\":%.17g,\"solver_vs_independent_forward_error\":%.17g,\"reference_without_low_residual\":%.17g,\"solver_without_low_residual\":%.17g,\"corrupted_vector_residual\":%.17g,\"NaN_rejected\":%s,\"Inf_rejected\":%s,\"TwoSum_FMA_exact_controls\":%s,\"conservative_pair_assembly_relative_defect\":%.17g,\"kinetic_pair_defect\":%.17g,\"entropy_balance_defect\":%.17g,\"mass_link_defect_using_stored_rate\":%.17g,\"maximum_z_materialization_abs_error\":%.17g,\"maximum_mu_z_materialization_abs_error\":%.17g,\"maximum_edge_z_materialization_scaled_error\":%.17g,\"conserved_state_and_birth_parameters_unchanged\":%s,\"accepted_evolution_steps\":0}\n",pass?"true":"false",stats.cg,stats.root_residual,independent,forward_error,reference_low_removed,solver_low_removed,corrupted,rejects_nan?"true":"false",rejects_inf?"true":"false",arithmetic?"true":"false",conservation,kinetic,stats.entropy_balance_defect,stats.mass_link_defect,stats.maximum_z_materialization_error,stats.maximum_mu_z_materialization_error,stats.maximum_edge_z_materialization_scaled_error,unchanged?"true":"false");
    if(fclose(out))return 4;free(high);free(low);free(R);free(saved);free(rate);gfs_fv_workspace_free(w);gfs_fv_geometry_free(&m);free(m.cell);return pass?0:1;
}
