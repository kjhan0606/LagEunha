/* Saved strong-shock failure: geometry, global physical-field point
 * ownership, and the actual screened matrix. READ ONLY, no evolution. */
#define main archived_native_main
#include "fv_native_evolution.c"
#undef main
#include "../gfs_fv_native.c"
int main(int argc,char **argv){
    if(argc!=3)return 2;FILE *in=fopen(argv[1],"rb");FvCheckpoint h;
    if(!in||fread(&h,sizeof(h),1,in)!=1||h.magic!=magic||h.dim!=1||h.count>1000000||
       (h.version&~diagnostic_trial_flag)!=checkpoint_version)return 2;
    GfsFvMesh m={0};m.dim=1;m.n=h.count;m.gamma=h.gamma;memcpy(m.bc,h.bc,sizeof(m.bc));memcpy(m.lo,h.lo,sizeof(m.lo));memcpy(m.L,h.L,sizeof(m.L));m.cell=calloc(m.n,sizeof(*m.cell));
    if(!m.cell||fread(m.cell,sizeof(*m.cell),m.n,in)!=(size_t)m.n||fgetc(in)!=EOF||fclose(in)||gfs_fv_geometry(&m))return 2;
    GfsFvWorkspace *w=gfs_fv_workspace_create();allocate(w,m.n);GfsFvRhsStats s={0};
    if(reconstruct(w,&m)||build_points(w,&m,&s))return 3;
    if(residual(w,&m,w->z,w->R,1))return 3;
    double min_gap=HUGE_VAL,max_ratio=0;int min_i=-1,outside=0;double *diagonal=calloc(m.n,sizeof(double));
    for(int i=0;i<m.n;i++)diagonal[i]=m.cell[i].mu;
    char path[4096];snprintf(path,sizeof(path),"%s_matrix.csv",argv[2]);FILE *matrix=fopen(path,"wx");if(!matrix)return 4;
    fprintf(matrix,"i,j,aij\n");
    for(int e=0;e<w->ne;e++){
        FvEdge *edge=w->edge+e;if(edge->j<0)continue;double jac=0;
        for(int q=0;q<edge->count;q++){FvFanPoint *p=w->point+edge->offset+q;jac+=p->weight*p->inv2d*p->density;}
        diagonal[edge->i]+=jac;diagonal[edge->j]+=jac;fprintf(matrix,"%d,%d,%.17g\n",edge->i,edge->j,jac);
    }if(fclose(matrix))return 4;
    snprintf(path,sizeof(path),"%s_rhs.csv",argv[2]);FILE *rhs=fopen(path,"wx");if(!rhs)return 4;
    fprintf(rhs,"id,b,scale\n");
    for(int i=0;i<m.n;i++){
        double speed=sqrt(dot3(w->primitive[i].u,w->primitive[i].u,m.dim))+sqrt(m.gamma*w->primitive[i].p/w->primitive[i].rho);
        double scale=m.cell[i].mu*sqrt(m.cell[i].m0/m.cell[i].mu)*speed;
        fprintf(rhs,"%d,%.17g,%.17g\n",i,-w->R[i],scale);
    }if(fclose(rhs))return 4;
    snprintf(path,sizeof(path),"%s_cells.csv",argv[2]);FILE *out=fopen(path,"wx");if(!out)return 4;
    fprintf(out,"id,q,centroid,V,mass,rho,P,mean_v,own_extrapolated_qdot,actual_point_owner,global_reconstruction_qdot,mu,diagonal_over_mu,chi,raw_bulk_fraction,Eparallel,Eperpendicular\n");
    double *global=calloc(m.n,sizeof(double));int *owner=calloc(m.n,sizeof(int));
    for(int i=0;i<m.n;i++){
        double x=m.cell[i].q[0];int found=-1;
        /* Independent actual-interval location for this read-only audit.
         * This O(N^2) diagnostic is NOT the proposed production locator. */
        for(int j=0;j<m.n;j++){
            double center=m.cell[j].q[0]+m.geometry[j].centroid[0],left=center-.5*m.geometry[j].V,right=center+.5*m.geometry[j].V;
            if(x>=left&&x<=right){found=j;break;}
        }
        if(found<0)return 5;owner[i]=found;outside+=found!=i;
        double rel[3]={x-m.cell[found].q[0],0,0};global[i]=trace(w,&m,found,rel).u[0];
        max_ratio=fmax(max_ratio,diagonal[i]/m.cell[i].mu);
        if(i+1<m.n&&m.cell[i+1].q[0]-x<min_gap){min_gap=m.cell[i+1].q[0]-x;min_i=i;}
        double own=w->primitive[i].u[0]-dot3(w->gradient[i][1],m.geometry[i].centroid,1);
        fprintf(out,"%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",i,x,x+m.geometry[i].centroid[0],m.geometry[i].V,m.cell[i].Q[0],w->primitive[i].rho,w->primitive[i].p,w->primitive[i].u[0],own,found,global[i],m.cell[i].mu,diagonal[i]/m.cell[i].mu,w->stress_fraction[i],w->bulk_fraction[i],w->parallel_energy[i],w->perpendicular_energy[i]);
    }if(fclose(out))return 4;
    double own0=w->primitive[min_i].u[0]-dot3(w->gradient[min_i][1],m.geometry[min_i].centroid,1);
    double own1=w->primitive[min_i+1].u[0]-dot3(w->gradient[min_i+1][1],m.geometry[min_i+1].centroid,1);
    int status=solve(w,&m,&s),paired_cells=0;
#ifdef GFS_NATIVE_BLOCK_JACOBI
    const char *preconditioner="mutual_strongest_edge_2x2_block_Jacobi";
    for(int i=0;i<m.n;i++)paired_cells+=w->partner[i]>=0;
#else
    const char *preconditioner="diagonal_Jacobi";
#endif
    printf("{\"scope\":\"readonly same-state point ownership and matrix audit; global-field velocities are NOT evolved\",\"time\":%.17g,\"cells\":%d,\"diagnostic_trial_not_restartable\":%s,\"preconditioner\":\"%s\",\"paired_cells\":%d,\"generators_outside_own_cell\":%d,\"closest_pair\":[%d,%d],\"minimum_generator_gap\":%.17g,\"pair_mean_velocity_difference\":%.17g,\"pair_own_extrapolated_qdot_difference\":%.17g,\"pair_selected_qdot_difference\":%.17g,\"pair_global_physical_qdot_difference\":%.17g,\"pair_point_owners\":[%d,%d],\"maximum_diagonal_over_mu\":%.17g,\"screened_solve_status\":%d,\"actual_PCG_iterations\":%d,\"PCG_failure_reason\":%d,\"recursive_preconditioned_relative_residual\":%.17g,\"true_residual\":%.17g}\n",h.time,m.n,(h.version&diagnostic_trial_flag)?"true":"false",preconditioner,paired_cells,outside,min_i,min_i+1,min_gap,w->primitive[min_i+1].u[0]-w->primitive[min_i].u[0],own1-own0,w->velocity[min_i+1][0]-w->velocity[min_i][0],global[min_i+1]-global[min_i],owner[min_i],owner[min_i+1],max_ratio,status,w->pcg_iterations,w->pcg_failure,w->pcg_last_relative,s.root_residual);
    free(global);free(owner);free(diagonal);free(m.cell);gfs_fv_workspace_free(w);gfs_fv_geometry_free(&m);return 0;
}
