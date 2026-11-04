/* Read-only D1 geometry failure audit. The exhaustive halfspace oracle is
 * deliberately limited to small saved states and is never an evolver. */
#define main archived_native_main
#include "fv_native_evolution.c"
#undef main
int main(int argc,char **argv){
    if(argc!=3)return 2;
    FILE *in=fopen(argv[1],"rb");FvCheckpoint h;
    if(!in||fread(&h,sizeof(h),1,in)!=1||h.magic!=magic||h.dim!=1||h.count<2||h.count>4096||h.cell_bytes!=sizeof(GfsFvCell)||(h.version&~diagnostic_trial_flag)!=checkpoint_version)return 2;
    GfsFvMesh m={0};m.dim=1;m.n=h.count;m.gamma=h.gamma;
    memcpy(m.bc,h.bc,sizeof(m.bc));memcpy(m.lo,h.lo,sizeof(m.lo));memcpy(m.L,h.L,sizeof(m.L));
    m.cell=calloc(m.n,sizeof(*m.cell));
    if(!m.cell||fread(m.cell,sizeof(*m.cell),m.n,in)!=(size_t)m.n||fgetc(in)!=EOF||fclose(in))return 2;
    GfsFvCell *saved=malloc(m.n*sizeof(*saved));memcpy(saved,m.cell,m.n*sizeof(*saved));
    int geometry=gfs_fv_geometry(&m),rhs=-1;GfsFvRhsStats stats={0};
    GfsFvWorkspace *w=gfs_fv_workspace_create();GfsFvRate *rate=calloc(m.n,sizeof(*rate));
    if(!geometry)rhs=gfs_fv_native_rhs(w,&m,rate,&stats);
    int lower=gfs_fv_boundary_type(m.bc[0],0),upper=gfs_fv_boundary_type(m.bc[0],1);
    int periodic=lower==GFS_FV_PERIODIC;
    long double boundary_left=m.lo[0],boundary_right=m.lo[0]+m.L[0];
    if(lower==GFS_FV_MOVING_TRANSMIT)boundary_left=1.5L*m.cell[0].q[0]-.5L*m.cell[1].q[0];
    if(upper==GFS_FV_MOVING_TRANSMIT)boundary_right=1.5L*m.cell[m.n-1].q[0]-.5L*m.cell[m.n-2].q[0];
    FILE *out=fopen(argv[2],"wx");if(!out)return 4;
    fprintf(out,"id,q,W,mass,mu,adjacent_width,exhaustive_width,left_owner,right_owner,qdot,zdot,analytic_Vdot\n");
    long double min_adj=HUGE_VALL,min_true=HUGE_VALL;double min_gap=HUGE_VAL,min_wall=HUGE_VAL,max_error=0;
    int min_i=-1,min_true_i=-1,hidden=0;
    for(int i=0;i<m.n;i++){
        double qi=m.cell[i].q[0],Wi=m.cell[i].W;
        if(i+1<m.n)min_gap=fmin(min_gap,m.cell[i+1].q[0]-qi);
        if(lower!=GFS_FV_MOVING_TRANSMIT&&!periodic)min_wall=fmin(min_wall,qi-m.lo[0]);
        if(upper!=GFS_FV_MOVING_TRANSMIT&&!periodic)min_wall=fmin(min_wall,m.lo[0]+m.L[0]-qi);
        long double left=periodic?-HUGE_VALL:boundary_left,right=periodic?HUGE_VALL:boundary_right;
        int left_id=-1,right_id=-1;
        for(int shift=periodic?-1:0;shift<=(periodic?1:0);shift++)for(int j=0;j<m.n;j++){
            if(i==j&&shift==0)continue;
            long double qj=(long double)m.cell[j].q[0]+shift*m.L[0],r=qj-qi;
            if(r==0){if(m.cell[j].W>Wi){left=1;right=0;}continue;}
            long double plane=.5L*((long double)qi+qj)+((long double)Wi-m.cell[j].W)/(2*r);
            if(r<0&&plane>left){left=plane;left_id=j;}
            if(r>0&&plane<right){right=plane;right_id=j;}
        }
        long double own_left,own_right;
        if(i==0&&!periodic)own_left=boundary_left;
        else{int j=(i+m.n-1)%m.n;long double qj=(long double)m.cell[j].q[0]-(i==0?m.L[0]:0),r=qj-qi;own_left=.5L*((long double)qi+qj)+((long double)Wi-m.cell[j].W)/(2*r);}
        if(i==m.n-1&&!periodic)own_right=boundary_right;
        else{int j=(i+1)%m.n;long double qj=(long double)m.cell[j].q[0]+(i==m.n-1?m.L[0]:0),r=qj-qi;own_right=.5L*((long double)qi+qj)+((long double)Wi-m.cell[j].W)/(2*r);}
        long double adjacent=own_right-own_left,exact=right-left;
        if(adjacent<min_adj){min_adj=adjacent;min_i=i;}
        if(exact<min_true){min_true=exact;min_true_i=i;}
        hidden+=exact<=0;
        if(!geometry)max_error=fmax(max_error,fabs(m.geometry[i].V-(double)exact));
        fprintf(out,"%d,%.17g,%.17g,%.17g,%.17g,%.21Lg,%.21Lg,%d,%d,%.17g,%.17g,%.17g\n",i,qi,Wi,m.cell[i].Q[0],m.cell[i].mu,adjacent,exact,left_id,right_id,rhs?NAN:rate[i].qdot[0],rhs?NAN:rate[i].z,rhs?NAN:rate[i].Vdot);
    }
    if(fclose(out))return 4;
    char margin[64];if(isfinite(min_wall))snprintf(margin,sizeof(margin),"%.17g",min_wall);else strcpy(margin,"null");
    printf("{\"scope\":\"readonly saved D1 state with independent exhaustive halfspaces; no evolution\",\"time\":%.17g,\"diagnostic_trial_not_restartable\":%s,\"cells\":%d,\"native_geometry_status\":%d,\"native_rhs_status\":%d,\"minimum_ordered_generator_gap\":%.17g,\"minimum_fixed_wall_generator_margin\":%s,\"minimum_adjacent_cell\":%d,\"minimum_adjacent_width\":%.21Lg,\"minimum_exhaustive_cell\":%d,\"minimum_exhaustive_width\":%.21Lg,\"hidden_cells\":%d,\"maximum_native_exhaustive_width_difference_when_geometry_valid\":%.17g,\"conserved_states_and_birth_parameters_unchanged\":%s}\n",h.time,(h.version&diagnostic_trial_flag)?"true":"false",m.n,geometry,rhs,min_gap,margin,min_i,min_adj,min_true_i,min_true,hidden,max_error,memcmp(saved,m.cell,m.n*sizeof(*saved))?"false":"true");
    free(saved);free(rate);gfs_fv_workspace_free(w);gfs_fv_geometry_free(&m);free(m.cell);return 0;
}
