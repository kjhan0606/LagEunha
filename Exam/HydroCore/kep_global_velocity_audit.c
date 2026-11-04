/* Old accepted blast state is read ONLY to verify the new global physical
 * field query. No continuation, state conversion, or timestep is taken. */
#define main archived_native_main
#include "fv_native_evolution.c"
#undef main
#include "../gfs_fv_native.c"
int main(int argc,char **argv){
    if(argc!=2)return 2;FILE *in=fopen(argv[1],"rb");FvCheckpoint h;
    if(!in||fread(&h,sizeof(h),1,in)!=1||h.magic!=magic||h.version!=3||h.dim!=1||h.count>1000000)return 2;
    GfsFvMesh m={0};m.dim=1;m.n=h.count;m.gamma=h.gamma;memcpy(m.bc,h.bc,sizeof(m.bc));memcpy(m.lo,h.lo,sizeof(m.lo));memcpy(m.L,h.L,sizeof(m.L));m.cell=calloc(m.n,sizeof(*m.cell));
    if(!m.cell||fread(m.cell,sizeof(*m.cell),m.n,in)!=(size_t)m.n||fgetc(in)!=EOF||fclose(in)||gfs_fv_geometry(&m))return 2;
    GfsFvCell *saved=malloc(m.n*sizeof(*saved));memcpy(saved,m.cell,m.n*sizeof(*saved));
    GfsFvWorkspace *w=gfs_fv_workspace_create();allocate(w,m.n);if(reconstruct(w,&m))return 3;
    int i=0;double gap=HUGE_VAL,unchanged=0;
    for(int k=0;k+1<m.n;k++)if(m.cell[k+1].q[0]-m.cell[k].q[0]<gap){i=k;gap=m.cell[k+1].q[0]-m.cell[k].q[0];}
    for(int k=0;k<m.n;k++)if(w->velocity_owner[k]==k){
        double own=w->primitive[k].u[0]-dot3(w->gradient[k][1],m.geometry[k].centroid,1);
        unchanged=fmax(unchanged,fabs(own-w->velocity[k][0])/(1+fabs(own)));
    }
    int owner=w->velocity_owner[i];double actual=w->velocity[i+1][0]-w->velocity[i][0],expected=w->gradient[owner][1][0]*gap;
    double own0=w->primitive[i].u[0]-dot3(w->gradient[i][1],m.geometry[i].centroid,1),own1=w->primitive[i+1].u[0]-dot3(w->gradient[i+1][1],m.geometry[i+1].centroid,1);
    int pass=owner==w->velocity_owner[i+1]&&fabs(actual-expected)<1e-12&&unchanged<1e-15&&!memcmp(saved,m.cell,m.n*sizeof(*saved));
    printf("{\"scope\":\"readonly old accepted state: physical owner/affine gap identity only, not a new trajectory\",\"pair\":[%d,%d],\"owner\":%d,\"own_extrapolated_gap_rate\":%.17g,\"global_physical_gap_rate\":%.17g,\"affine_same_cell_expected\":%.17g,\"inside_own_cell_velocity_error\":%.17g,\"conserved_states_and_birth_parameters_unchanged\":%s,\"maximum_hops\":%d,\"pass\":%s}\n",i,i+1,owner,own1-own0,actual,expected,unchanged,memcmp(saved,m.cell,m.n*sizeof(*saved))?"false":"true",w->query_max_hops,pass?"true":"false");
    free(saved);free(m.cell);gfs_fv_workspace_free(w);gfs_fv_geometry_free(&m);return pass?0:1;
}
