/* Opt-in periodic candidate. One O(Nfaces) sweep, no eigensolve or mesh build.
 * This bounds the positive acoustic part, not the full geometric Hessian. */
#include "../gfs_acoustic_cfl.h"
double gpu_unified_acoustic_timestep(void *opaque,double gamma,double courant)
{
    double begin=MPI_Wtime();
    SimParameters *simpar=opaque;GPUContext *ctx=gpu_get_context();
    int n=VORO_NP(simpar),nt=ctx->h_parts.n_total;
    if(!ctx->initialized||NID(simpar)!=1||!gfs_unified_periodic_mode()||n<1)return 0;
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(simpar);
    FaceCSR *f=&ctx->h_faces;ParticleSoA *a=&ctx->h_parts;
    int nf=ctx->d_faces.n_faces_total;
    gpu_download_face_csr(ctx,n,nf,f->face_offset,f->c1x,f->c1y,f->c2x,f->c2y,
                         f->neighbor_idx,f->kp_idx,f->km_idx,f->is_ghost);
    if(f->face_offset[0]!=0||f->face_offset[n]!=nf||nf<1)return 0;
    const int *alias=gpu_unified_owner_aliases(ctx,n,nt,
        SIMBOX(simpar).x.max-SIMBOX(simpar).x.min,
        SIMBOX(simpar).y.max-SIMBOX(simpar).y.min);
    double max_row=0,max_column=0;int bad=0;
#pragma omp parallel for schedule(static) reduction(max:max_row,max_column) reduction(|:bad)
    for(int i=0;i<n;i++){
        GfsAcousticBound b;
        if(!gfs_acoustic_bound_init(&b,2,p[i].mass,p[i].weight_thermal.V,p[i].weight_thermal.U,gamma)){bad=1;continue;}
        for(int k=f->face_offset[i];k<f->face_offset[i+1];k++){
            int raw=f->neighbor_idx[k];if(raw<0||f->is_ghost[k]){bad=1;continue;}
            int j=alias[raw];GfsVariationalFace face;
            double r[3]={a->x[raw]-p[i].x,a->y[raw]-p[i].y,0};
            double c[3]={.5*(f->c1x[k]+f->c2x[k]),.5*(f->c1y[k]+f->c2y[k]),0};
            double area=hypot(f->c2y[k]-f->c1y[k],f->c1x[k]-f->c2x[k]);
            if(!gfs_variational_face(2,r,c,area,&face)||
               !gfs_acoustic_bound_face(&b,&face,p[j].mass,p[j].weight_thermal.V,p[j].weight_thermal.U,gamma))bad=1;
        }
        gfs_acoustic_bound_finish(&b);max_row=fmax(max_row,b.row);max_column=fmax(max_column,b.column);
    }
    double dt=bad?0:gfs_acoustic_bound_timestep(courant,max_row,max_column);
    fprintf(stderr,"[UNIFIED_ACOUSTIC_CFL] positive_part_only=1 row=%.17g column=%.17g lambda_bound=%.17g candidate_dt=%.17g seconds=%.9g\n",max_row,max_column,max_row*max_column,dt,MPI_Wtime()-begin);
    return dt;
}
