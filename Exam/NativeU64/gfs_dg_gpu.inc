/* Read-only native Euler projection gate. No timestep, evolved checkpoint,
 * limiter, mass/weight change, or Euler/Pi split is applied. */
__device__ static void gfs_dg_assess_point(const GfsDgConfig *cfg,
    GfsDgProjectionAudit *a,const GfsDgPolynomial *coarse,double x,double y,
    int kind,double weight,double integral[4]){
    double u[4],uc[4];gfs_dg_evaluate(&a->basis,&a->polynomial,cfg->order,x,y,u);
    gfs_dg_evaluate(&a->basis,coarse,cfg->order,x,y,uc);
    double *m=a->metrics;m[DG_MIN_RHO]=fmin(m[DG_MIN_RHO],u[0]);
    if(!(u[0]>0) || !isfinite(u[0])){a->nonpositive_rho++;return;}
    double p=gfs_dg_pressure(u,cfg->ic.gamma),pc=gfs_dg_pressure(uc,cfg->ic.gamma);
    if(!(p>0) || !isfinite(p))a->nonpositive_p++;
    m[DG_MIN_P]=fmin(m[DG_MIN_P],p);m[DG_MAX_P]=fmax(m[DG_MAX_P],p);
    double error=fabs(p/cfg->ic.pressure-1);
    if(error>m[DG_MAX_PERR]){m[DG_MAX_PERR]=error;m[DG_WORST_X]=x;m[DG_WORST_Y]=y;a->worst_kind=kind;}
    if(uc[0]>0 && isfinite(pc))m[DG_MAX_P_QUAD_DELTA]=fmax(m[DG_MAX_P_QUAD_DELTA],fabs(p-pc)/cfg->ic.pressure);
    if(weight>0){
        m[DG_L1_PERR]+=weight*error;m[DG_MEAN_P]+=weight*p;
        for(int c=0;c<4;c++)integral[c]+=weight*u[c];
    }
}

__global__ static void gfs_dg_project_kernel(ParticleSoA parts,const PolygonOut *polygons,
    int np,GfsDgConfig cfg,GfsDgProjectionAudit *audit,PolygonOut *canonical){
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=np)return;
    GfsDgProjectionAudit a={};a.status=-10;
    double x[GPU_MAX_CORNERS],y[GPU_MAX_CORNERS];
    if(!gfs_moment_support_polygon(parts,polygons+i,i,x,y)){audit[i]=a;return;}
    GfsMomentPolygon polygon={x,y,polygons[i].n_corners,0};
    PolygonOut out_poly=polygons[i];
    for(int k=0;k<polygon.n;k++){out_poly.cx[k]=x[k];out_poly.cy[k]=y[k];}
    canonical[i]=out_poly;
    GfsDgPolynomial coarse={};
    if(!gfs_dg_project(polygon,&cfg,parts.x[i],parts.y[i],&a,&coarse)){audit[i]=a;return;}
    a.metrics[DG_MIN_RHO]=a.metrics[DG_MIN_P]=DBL_MAX;a.metrics[DG_MAX_P]=-DBL_MAX;
    double integrated[4]={0,0,0,0};
    int nq=polygon.n*cfg.ng[2]*cfg.ng[2];
    for(int q=0;q<nq;q++){
        double xx,yy,w;gfs_dg_quad(polygon,&a.basis,&cfg,2,q,&xx,&yy,&w);
        if(w==0)continue;
        gfs_dg_assess_point(&cfg,&a,&coarse,xx,yy,0,w,integrated);
    }
    for(int k=0;k<polygon.n;k++){
        int j=(k+1)%polygon.n;
        gfs_dg_assess_point(&cfg,&a,&coarse,x[k],y[k],1,0,integrated);
        for(int q=0;q<cfg.ng[2];q++){
            double t=cfg.node[2][q];
            gfs_dg_assess_point(&cfg,&a,&coarse,(1-t)*x[k]+t*x[j],(1-t)*y[k]+t*y[j],2,0,integrated);
        }
    }
    a.metrics[DG_L1_PERR]/=a.basis.area;a.metrics[DG_MEAN_P]/=a.basis.area;
    double mass=a.integral[0],energy=a.integral[3];
    double scales[4]={mass,sqrt(2*mass*energy),sqrt(2*mass*energy),energy};
    for(int c=0;c<4;c++)a.metrics[DG_CONSERV_RHO+c]=fabs(integrated[c]-a.integral[c])/scales[c];
    a.metrics[DG_QUAD_MASS_REL]=fabs(a.integral[0]-a.coarse_integral[0])/mass;
    a.metrics[DG_QUAD_PX_REL]=fabs(a.integral[1]-a.coarse_integral[1])/scales[1];
    a.metrics[DG_QUAD_PY_REL]=fabs(a.integral[2]-a.coarse_integral[2])/scales[2];
    a.metrics[DG_QUAD_E_OVER_THERMAL]=fabs(a.integral[3]-a.coarse_integral[3])/
        (cfg.ic.pressure*a.basis.area/(cfg.ic.gamma-1));
    a.status=1;audit[i]=a;
}

extern "C" int gpu_dg_project_initial(GPUContext *ctx,int np,int rank,const GfsDgConfig *cfg,
    GfsDgProjectionAudit *audit){
    GfsDgProjectionAudit *da=NULL;PolygonOut *dp=NULL;
    if(np){
        CUDA_CHECK(cudaMalloc(&da,(size_t)np*sizeof(*da)));CUDA_CHECK(cudaMalloc(&dp,(size_t)np*sizeof(*dp)));
        cudaStream_t s=(cudaStream_t)ctx->stream;
        gfs_dg_project_kernel<<<(np+63)/64,64,0,s>>>(ctx->d_parts,ctx->d_polygons,np,*cfg,da,dp);
        CUDA_CHECK(cudaGetLastError());
        CUDA_CHECK(cudaMemcpyAsync(audit,da,(size_t)np*sizeof(*da),cudaMemcpyDeviceToHost,s));
        CUDA_CHECK(cudaStreamSynchronize(s));
    }
    /* Preserve the exact canonical native polygons for independent CPU
     * quadrature. This is a geometry artifact, not an evolved state. */
    char path[160];snprintf(path,sizeof(path),"dg_polygons_r%d.bin",rank);
    FILE *out=fopen(path,"wx");if(!out){perror(path);exit(1);}
    unsigned long long header[4]={0x4746534447504f31ULL,1,(unsigned long long)np,0};
    fwrite(header,sizeof(header),1,out);
    PolygonOut *poly=(PolygonOut*)malloc((size_t)(np?np:1)*sizeof(*poly));if(!poly)exit(1);
    if(np)CUDA_CHECK(cudaMemcpy(poly,dp,(size_t)np*sizeof(*poly),cudaMemcpyDeviceToHost));
    int bad=0;
    for(int i=0;i<np;i++){
        unsigned long long id=ctx->h_parts.indx[i];int count=audit[i].status!=-10?poly[i].n_corners:0;
        int pair[2]={count,audit[i].status};fwrite(&id,sizeof(id),1,out);fwrite(pair,sizeof(pair),1,out);
        for(int k=0;k<count;k++){double xy[2]={poly[i].cx[k],poly[i].cy[k]};fwrite(xy,sizeof(xy),1,out);}
        if(audit[i].status!=1)bad++;
    }
    if(ferror(out) || fclose(out))exit(1);free(poly);
    if(np){CUDA_CHECK(cudaFree(da));CUDA_CHECK(cudaFree(dp));}
    return bad;
}
