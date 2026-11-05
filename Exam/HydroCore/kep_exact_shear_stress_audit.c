/* Analytic incompressible shear cell means on saved native geometry.
 * Manufactured operator input only, never an evolved/restart checkpoint. */
#define GFS_STRESS_AUDIT_ENTRY archived_stress_audit_entry
#include "kep_stress_switch_audit.c"
#undef GFS_STRESS_AUDIT_ENTRY
int main(int argc,char **argv){
    if(argc!=3)return 2;FILE *in=fopen(argv[1],"rb");FvCheckpoint h;
    if(!in||fread(&h,sizeof(h),1,in)!=1||h.magic!=magic||h.version!=9||h.dim!=2||h.cell_bytes!=sizeof(GfsFvCell)||h.count<2||h.count>1000000)return 2;
    GfsFvMesh m={0};m.dim=2;m.n=h.count;m.gamma=h.gamma;memcpy(m.bc,h.bc,sizeof(m.bc));memcpy(m.lo,h.lo,sizeof(m.lo));memcpy(m.L,h.L,sizeof(m.L));
    m.cell=calloc(m.n,sizeof(*m.cell));if(!m.cell||fread(m.cell,sizeof(*m.cell),m.n,in)!=(size_t)m.n||fgetc(in)!=EOF||fclose(in)||gfs_fv_geometry(&m))return 2;
    const double wave=2*acos(-1.);
    for(int i=0;i<m.n;i++){
        GfsFvGeometry *g=m.geometry+i;long double integral=0;
        /* div(0,-cos(k y)/k)=sin(k y). Every 2D face is a straight
         * segment. Its cosine integral is an exact sinc expression. */
        for(int k=0;k<g->faces;k++){
            GfsFvNativeFace *f=g->face+k;double y=m.cell[i].q[1]+f->centroid[1];
            double half_phase=.5*wave*f->area*f->normal[0];
            double sinc=half_phase==0?1:sin(half_phase)/half_phase;
            integral-=(long double)f->normal[1]*f->area*cos(wave*y)*sinc/wave;
        }
        double u=(double)(integral/g->V);GfsFvCell *c=m.cell+i;
        c->Q[0]=g->V;c->Q[1]=g->V*u;c->Q[2]=c->Q[3]=0;
        c->Q[4]=g->V*1e-6/(m.gamma-1)+.5*g->V*u*u;
    }
    int status=inspect(&m,"analytic_sine_shear_cell_means_on_saved_geometry",argv[2],-1);
    gfs_fv_geometry_free(&m);free(m.cell);return status;
}
