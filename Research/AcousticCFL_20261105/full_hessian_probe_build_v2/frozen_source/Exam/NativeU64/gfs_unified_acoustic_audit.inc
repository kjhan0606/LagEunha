/* Frozen periodic diagnostic of the positive acoustic Hessian part
 * M^(-1/2) B^T diag(gamma*P/V) B M^(-1/2). The pressure-weighted second
 * derivative of geometry is not included, so this is not the full Hessian.
 * Power iteration never changes particles or selects the actual timestep. */
static int gfs_unified_acoustic_audit(SimParameters *simpar,GPUContext *ctx,
        const int *alias,double gamma)
{
    if(!gfs_unified_periodic_mode())return 1;
    int n=VORO_NP(simpar),nf=ctx->d_faces.n_faces_total;FaceCSR *f=&ctx->h_faces;ParticleSoA *a=&ctx->h_parts;
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(simpar);
    double *bi=calloc(2*nf,sizeof(double)),*bj=calloc(2*nf,sizeof(double)),*diag=calloc(2*n,sizeof(double));
    double *C=calloc(n,sizeof(double)),*rootm=calloc(n,sizeof(double)),*z=calloc(2*n,sizeof(double)),*w=calloc(n,sizeof(double)),*F=calloc(2*n,sizeof(double));
    if(!bi||!bj||!diag||!C||!rootm||!z||!w||!F)return 1;
    double min_dt=HUGE_VAL,min_pair_over_h=HUGE_VAL,row_bound=0,column_bound=0;
    for(int i=0;i<n;i++){
        double V=p[i].weight_thermal.V,P=(gamma-1)*p[i].weight_thermal.U/V;
        C[i]=gamma*P/V;rootm[i]=sqrt(p[i].mass);min_dt=fmin(min_dt,p[i].dt);
        z[2*i]=sin(.31*(double)PINDX(p+i));z[2*i+1]=cos(.47*(double)PINDX(p+i));
        for(int k=f->face_offset[i];k<f->face_offset[i+1];k++){
            int raw=f->neighbor_idx[k];if(raw<0||f->is_ghost[k])return 1;
            double r[3]={a->x[raw]-p[i].x,a->y[raw]-p[i].y,0},c[3]={.5*(f->c1x[k]+f->c2x[k]),.5*(f->c1y[k]+f->c2y[k]),0};
            GfsVariationalFace b;
            if(!gfs_variational_face(2,r,c,hypot(f->c2y[k]-f->c1y[k],f->c1x[k]-f->c2x[k]),&b))return 1;
            min_pair_over_h=fmin(min_pair_over_h,b.separation/sqrt(V));
            for(int d=0;d<2;d++){bi[2*k+d]=b.bi[d];bj[2*k+d]=b.bj[d];diag[2*i+d]-=b.bj[d];}
        }
    }
    for(int i=0;i<n;i++){
        double row=(fabs(diag[2*i])+fabs(diag[2*i+1]))/rootm[i],col[2]={sqrt(C[i])*fabs(diag[2*i]),sqrt(C[i])*fabs(diag[2*i+1])};
        for(int k=f->face_offset[i];k<f->face_offset[i+1];k++){
            int j=alias[f->neighbor_idx[k]];
            row+=(fabs(bj[2*k])+fabs(bj[2*k+1]))/rootm[j];
            for(int d=0;d<2;d++)col[d]+=sqrt(C[j])*fabs(bi[2*k+d]);
        }
        row_bound=fmax(row_bound,sqrt(C[i])*row);column_bound=fmax(column_bound,fmax(col[0],col[1])/rootm[i]);
    }
    double lambda=0,previous=0,adjoint=0;
    for(int iter=0;iter<40;iter++){
#pragma omp parallel for schedule(static)
        for(int i=0;i<n;i++){
            double volume_rate=0;
            for(int k=f->face_offset[i];k<f->face_offset[i+1];k++){
                int j=alias[f->neighbor_idx[k]];
                for(int d=0;d<2;d++)volume_rate+=bj[2*k+d]*(z[2*j+d]/rootm[j]-z[2*i+d]/rootm[i]);
            }
            w[i]=C[i]*volume_rate;
        }
#pragma omp parallel for schedule(static)
        for(int i=0;i<n;i++){
            F[2*i]=F[2*i+1]=0;
            for(int k=f->face_offset[i];k<f->face_offset[i+1];k++){
                int j=alias[f->neighbor_idx[k]];
                for(int d=0;d<2;d++)F[2*i+d]+=(w[i]-w[j])*bi[2*k+d]/rootm[i];
            }
        }
        long double dot=0,norm=0,znorm=0,quadratic=0;
        for(int i=0;i<2*n;i++){dot+=(long double)z[i]*F[i];norm+=(long double)F[i]*F[i];znorm+=(long double)z[i]*z[i];}
        for(int i=0;i<n;i++)quadratic+=(long double)w[i]*w[i]/C[i];
        adjoint=fmax(adjoint,(double)(fabsl(dot-quadratic)/fmaxl(1,fabsl(quadratic))));
        previous=lambda;lambda=(double)(dot/znorm);if(!(norm>0)||!isfinite(lambda))return 1;
        double inv=1/sqrt((double)norm);for(int i=0;i<2*n;i++)z[i]=F[i]*inv;
    }
    fprintf(stderr,"[UNIFIED_ACOUSTIC_AUDIT] N=%d min_pair_over_sqrtV=%.17g native_CFL_dt=%.17g positive_part_lambda_est=%.17g lambda_upper_bound=%.17g native_dt_omega_est=%.17g iteration_relative_change=%.17g adjoint_residual=%.17g full_geometric_Hessian_included=0 no_state_update=1\n",n,min_pair_over_h,min_dt,lambda,row_bound*column_bound,min_dt*sqrt(lambda),fabs(lambda-previous)/lambda,adjoint);
    FILE *out=fopen("unified_acoustic_direction.csv","wx");if(!out)return 1;
    fprintf(out,"id,zx,zy\n");for(int i=0;i<n;i++)fprintf(out,"%llu,%.17g,%.17g\n",(unsigned long long)PINDX(p+i),z[2*i],z[2*i+1]);if(fclose(out))return 1;
    free(bi);free(bj);free(diag);free(C);free(rootm);free(z);free(w);free(F);return 0;
}
