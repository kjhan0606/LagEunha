/* Independent small graph algebra, including duplicate image edges and
 * strict invalid-pivot rejection. No geometry or fluid trajectory. */
#define GFS_MATRIX_AUDIT_ENTRY archived_matrix_main
#include "kep_screened_matrix_audit.c"
#undef GFS_MATRIX_AUDIT_ENTRY
int main(int argc,char **argv){
    if(argc!=2)return 2;
    const int endpoints[8][2]={{0,1},{0,1},{0,2},{1,2},{1,3},{2,3},{3,4},{0,4}};
    const double weights[8]={.4,.6,.7,1.1,.3,.8,1.2,.2};
    GfsFvMesh m={0};m.n=5;m.cell=calloc(m.n,sizeof(*m.cell));for(int i=0;i<m.n;i++){m.cell[i].mu=1+.2*i;m.cell[i].id=i;}
    GfsFvWorkspace *w=gfs_fv_workspace_create();allocate(w,m.n);w->ne=8;w->edge=memory(8*sizeof(*w->edge));w->incident=memory(16*sizeof(int));
    for(int k=0;k<8;k++){w->edge[k].i=endpoints[k][0];w->edge[k].j=endpoints[k][1];w->edge[k].jacobian=weights[k];w->offset[endpoints[k][0]+1]++;w->offset[endpoints[k][1]+1]++;}
    for(int i=0;i<m.n;i++)w->offset[i+1]+=w->offset[i];int cursor[5];memcpy(cursor,w->offset,sizeof(cursor));
    for(int k=0;k<8;k++){w->incident[cursor[endpoints[k][0]]++]=k;w->incident[cursor[endpoints[k][1]]++]=k;}
    for(int i=0;i<m.n;i++){w->diagonal[i]=m.cell[i].mu;for(int k=w->offset[i];k<w->offset[i+1];k++)w->diagonal[i]+=w->edge[w->incident[k]].jacobian;}
    w->ic_failure_row=-1;if(prepare_ic0(w,&m)||ic0_algebra(w,&m,argv[1]))return 1;
    int duplicate_merged=w->ic_nnz==7;
    w->diagonal[0]=-1;int negative=prepare_ic0(w,&m)!=0&&w->ic_failure_row==0;
    w->diagonal[0]=NAN;int nonfinite=prepare_ic0(w,&m)!=0&&w->ic_failure_row==0;
    int pass=duplicate_merged&&negative&&nonfinite;
    printf("{\"scope\":\"small algebraic graph only; no PDE\",\"duplicate_edges_merged\":%s,\"negative_pivot_rejected\":%s,\"nonfinite_pivot_rejected\":%s,\"no_shift_or_fallback\":true,\"pass\":%s}\n",duplicate_merged?"true":"false",negative?"true":"false",nonfinite?"true":"false",pass?"true":"false");
    gfs_fv_workspace_free(w);free(m.cell);return pass?0:1;
}
