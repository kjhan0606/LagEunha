/* Native constrained-wall fixtures and read-only saved-state cause check.
 * No accepted trajectory is produced by this executable. */
#define main archived_native_main
#include "fv_native_evolution.c"
#undef main
#include "../gfs_fv_native.c"
static void affine(int d,int walls,const double x[3],double u[3],double G[3][3]){
    memset(G,0,9*sizeof(double));
    for(int a=0;a<d;a++){
        u[a]=a<walls?0:.1*(a+1);
        if(a<walls)G[a][a]=-.2*(a+1);
        else for(int b=walls;b<d;b++)G[a][b]=.1*(a+1)*(b+1);
        for(int b=0;b<d;b++)u[a]+=G[a][b]*x[b];
    }
}
static int inspect(GfsFvMesh *m,const char *name,int walls,int nonlinear){
    GfsFvCell *copy=malloc(m->n*sizeof(*copy));memcpy(copy,m->cell,m->n*sizeof(*copy));
    GfsFvWorkspace *w=gfs_fv_workspace_create();GfsFvRate *rate=calloc(m->n,sizeof(*rate));GfsFvRhsStats s;
    int rhs=gfs_fv_native_rhs(w,m,rate,&s);if(rhs){fprintf(stderr,"wall audit RHS=%d\n",rhs);return 2;}
    double wall_error=0,mean_error=0,affine_error=0,equations=0,data_identity=0;int wall_cells=0,affine_cells=0;
    for(int i=0;i<m->n;i++){
        GfsFvGeometry *g=m->geometry+i;int interior=1;
        for(int a=0;a<m->dim;a++)if(m->cell[i].q[a]>.8||(a>=walls&&m->cell[i].q[a]<.2))interior=0;
        double x[3],u[3],G[3][3];for(int a=0;a<3;a++)x[a]=m->cell[i].q[a]+g->centroid[a];
        affine(m->dim,walls,x,u,G);
        if(!nonlinear&&interior){
            affine_cells++;
            for(int a=0;a<m->dim;a++)for(int b=0;b<m->dim;b++){
                affine_error=fmax(affine_error,fabs(w->raw_velocity_gradient[i][a][b]-G[a][b]));
                affine_error=fmax(affine_error,fabs(w->gradient[i][a+1][b]-G[a][b]));
            }
            double at_q[3],ignored[3][3];affine(m->dim,walls,m->cell[i].q,at_q,ignored);
            for(int a=0;a<m->dim;a++)affine_error=fmax(affine_error,fabs(w->velocity[i][a]-at_q[a]));
        }
        GfsFvPrimitive center=trace(w,m,i,g->centroid);
        for(int a=0;a<m->dim;a++)mean_error=fmax(mean_error,fabs(center.u[a]-w->primitive[i].u[a]));
        if(w->wall_mask[i])wall_cells++;
        for(int k=0;k<g->faces;k++){
            GfsFvNativeFace *f=g->face+k;
            for(int a=0;a<m->dim;a++)if(mirrored(f->mirror,a)&&face_bc(m,f,a)==GFS_FV_REFLECT){
                for(int q=0;q<f->count;q++){
                    GfsFvPrimitive p=trace(w,m,i,g->point[f->point+q].x);
                    wall_error=fmax(wall_error,fabs(p.u[a]));
                }
                double r[3]={0};r[a]=w->wall_relative[i][a];
                wall_error=fmax(wall_error,fabs(physical_velocity(w,m,i,a,r,NULL)));
            }
        }
        for(int a=0;a<m->dim;a++)if(!mirrored(w->wall_mask[i],a)){
            double normal[3]={0},scale=1;
            for(int k=0;k<g->faces;k++){
                GfsFvNativeFace *f=g->face+k;GfsFvPrimitive other;double r[3];neighbor(m,w,i,f,&other,r);
                double r2=dot3(r,r,m->dim),omega=f->area/pow(sqrt(r2),m->dim+1);
                double e=other.u[a]-w->primitive[i].u[a]-dot3(w->raw_velocity_gradient[i][a],r,m->dim);
                for(int b=0;b<m->dim;b++){normal[b]+=omega*e*r[b];scale+=fabs(omega*(other.u[a]-w->primitive[i].u[a])*r[b]);}
            }
            for(int b=0;b<m->dim;b++)equations=fmax(equations,fabs(normal[b])/scale);
        }
        data_identity=fmax(data_identity,fabs(w->data_energy[i]-w->parallel_energy[i]-w->perpendicular_energy[i]-w->cross_energy[i])/(1+w->data_energy[i]));
    }
    double kinetic=fabs(s.kinetic_rate-s.kinetic_pressure_work+s.kinetic_boundary_flux);
    int unchanged=!memcmp(copy,m->cell,m->n*sizeof(*copy));
    int pass=wall_cells>0&&wall_error<1e-12&&mean_error<1e-12&&equations<1e-12&&data_identity<1e-12&&kinetic<1e-10&&fabs(s.entropy_balance_defect)<1e-10&&unchanged;
    if(!nonlinear)pass&=affine_cells>0&&affine_error<1e-12;
    printf("{\"case\":\"%s\",\"dim\":%d,\"wall_cells\":%d,\"affine_checked_cells\":%d,\"normal_wall_trace_error\":%.17g,\"cell_mean_error\":%.17g,\"raw_limited_gradient_and_generator_affine_error\":%.17g,\"free_row_normal_equation_error\":%.17g,\"data_identity_with_cross_term_error\":%.17g,\"kinetic_pair_defect\":%.17g,\"entropy_balance_defect\":%.17g,\"states_and_birth_parameters_unchanged\":%s,\"pass\":%s}\n",name,m->dim,wall_cells,affine_cells,wall_error,mean_error,affine_error,equations,data_identity,kinetic,s.entropy_balance_defect,unchanged?"true":"false",pass?"true":"false");
    free(copy);free(rate);gfs_fv_workspace_free(w);return pass?0:1;
}
static int fixtures(void){
    int failed=0;
    for(int d=1;d<=3;d++)for(int walls=1;walls<=d;walls++)for(int fixture=0;fixture<3;fixture++){
        int nonlinear=fixture==1;
        GfsFvMesh m={0};m.dim=d;m.gamma=1.4;m.n=1;
        for(int a=0;a<3;a++){m.L[a]=1;if(a<d){m.n*=6;m.bc[a]=a<walls?gfs_fv_boundary_pair(GFS_FV_REFLECT,GFS_FV_OUTFLOW):GFS_FV_PERIODIC;}}
        m.cell=calloc(m.n,sizeof(*m.cell));
        for(int i=0;i<m.n;i++){
            int index=i;for(int a=0;a<d;a++){m.cell[i].q[a]=(index%6+.5)/6;index/=6;}m.cell[i].id=i;
            if(fixture)m.cell[i].W=.005/36*sin(.7*i);
        }
        if(gfs_fv_geometry(&m))return 2;
        for(int i=0;i<m.n;i++){
            GfsFvCell *c=m.cell+i;double x[3],u[3],G[3][3];for(int a=0;a<3;a++)x[a]=c->q[a]+m.geometry[i].centroid[a];
            affine(d,walls,x,u,G);if(nonlinear)for(int a=0;a<d;a++)u[a]+=.05*sin((a+1)*(.3+i));
            double V=m.geometry[i].V;c->Q[0]=V;c->Q[4]=V/(m.gamma-1);c->m0=V;c->mu=V/pow(V,2./d);c->W0=c->W;
            for(int a=0;a<d;a++){c->Q[a+1]=V*u[a];c->Q[4]+=.5*V*u[a]*u[a];}
        }
        failed|=inspect(&m,nonlinear?"weighted_nonaffine_wall_corner":fixture?"weighted_compatible_affine_wall_corner":"compatible_affine_wall_corner",walls,nonlinear);
        gfs_fv_geometry_free(&m);free(m.cell);
    }
    /* Real strip cells touching both x walls. A nonzero normal mean makes
     * an affine no-penetration field infeasible and must be rejected. */
    for(int d=2;d<=3;d++)for(int nonzero=0;nonzero<2;nonzero++){
        GfsFvMesh m={0};m.dim=d;m.gamma=1.4;m.n=d==2?4:16;m.bc[0]=GFS_FV_REFLECT;
        for(int a=0;a<3;a++)m.L[a]=1;m.cell=calloc(m.n,sizeof(*m.cell));
        for(int i=0;i<m.n;i++){m.cell[i].q[0]=.5;m.cell[i].q[1]=(i%4+.5)/4;m.cell[i].q[2]=d==3?(i/4+.5)/4:0;m.cell[i].id=i;}
        if(gfs_fv_geometry(&m))return 2;
        for(int i=0;i<m.n;i++){double V=m.geometry[i].V;m.cell[i].Q[0]=V;m.cell[i].Q[1]=V*nonzero;m.cell[i].Q[4]=V/(m.gamma-1)+.5*V*nonzero;}
        GfsFvWorkspace *w=gfs_fv_workspace_create();allocate(w,m.n);int status=reconstruct(w,&m);
        int pass=nonzero?status!=0:status==0;
        printf("{\"case\":\"opposite_walls\",\"dim\":%d,\"nonzero_normal_mean\":%s,\"infeasibility_rejected_when_required\":%s,\"pass\":%s}\n",d,nonzero?"true":"false",pass?"true":"false",pass?"true":"false");
        failed|=!pass;gfs_fv_workspace_free(w);gfs_fv_geometry_free(&m);free(m.cell);
    }
    return failed;
}
int main(int argc,char **argv){
    if(argc==2&&!strcmp(argv[1],"fixtures"))return fixtures();
    if(argc!=2)return 2;FILE *in=fopen(argv[1],"rb");FvCheckpoint h;
    if(!in||fread(&h,sizeof(h),1,in)!=1||h.magic!=magic||h.version!=4||h.dim!=1||h.problem!=BLAST||h.count>4096)return 2;
    GfsFvMesh m={0};m.dim=1;m.gamma=h.gamma;m.n=h.count;memcpy(m.bc,h.bc,sizeof(m.bc));memcpy(m.lo,h.lo,sizeof(m.lo));memcpy(m.L,h.L,sizeof(m.L));m.cell=calloc(m.n,sizeof(*m.cell));
    if(fread(m.cell,sizeof(*m.cell),m.n,in)!=(size_t)m.n||fgetc(in)!=EOF||fclose(in)||gfs_fv_geometry(&m))return 2;
    int code=inspect(&m,"readonly_old_blast_new_wall_field",1,1);
    gfs_fv_geometry_free(&m);free(m.cell);return code;
}
