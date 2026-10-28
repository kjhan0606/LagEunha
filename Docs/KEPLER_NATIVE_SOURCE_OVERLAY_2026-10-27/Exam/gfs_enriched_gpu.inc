/* Fresh native CUDA enriched projection only. No state evolution. */
__device__ static void gfs_enrich_assess_point(const GfsEnrichConfig *cfg,
    GfsEnrichAudit *a,double x,double y,double w,double integrated[4],double exact[4]){
    double u[4],v[4];
    if(!gfs_enrich_evaluate(&a->basis,cfg,&a->state,x,y,u) ||
       !gfs_enrich_reference(cfg,a->basis.generator[0]+x,a->basis.generator[1]+y,v,0,0)){
        a->status=-20;return;
    }
    double *m=a->metrics;m[ER_MIN_RHO]=fmin(m[ER_MIN_RHO],u[0]);
    double literal[4];gfs_dg_ic_state(&cfg->field,a->basis.generator[0]+x,a->basis.generator[1]+y,literal);
    m[ER_LITERAL_RHO_REL]=fmax(m[ER_LITERAL_RHO_REL],fabs(v[0]/literal[0]-1));
    if(!(u[0]>0) || !isfinite(u[0])){a->nonpositive_rho++;return;}
    double p=gfs_dg_pressure(u,cfg->field.ic.gamma);
    if(!(p>0) || !isfinite(p))a->nonpositive_p++;
    double err=fabs(p/cfg->field.ic.pressure-1);m[ER_MIN_P]=fmin(m[ER_MIN_P],p);
    m[ER_MAX_PERR]=fmax(m[ER_MAX_PERR],err);
    /* Linear conservative boost of the same state, no reset to reference. */
    double boosted[4]={u[0],u[1]+.013*u[0],u[2]-.007*u[0],
        u[3]+.013*u[1]-.007*u[2]+.5*(.013*.013+.007*.007)*u[0]};
    m[ER_BOOST_PERR]=fmax(m[ER_BOOST_PERR],fabs(gfs_dg_pressure(boosted,cfg->field.ic.gamma)/cfg->field.ic.pressure-1));
    for(int c=0;c<4;c++){
        m[ER_FIELD_ERROR]=fmax(m[ER_FIELD_ERROR],fabs(u[c]-v[c])/a->basis.scale[c]);
        if(w>0){integrated[c]+=w*u[c];exact[c]+=w*v[c];}
    }
    if(w>0)m[ER_L1_PERR]+=w*err;
}

__global__ static void gfs_enrich_project_kernel(ParticleSoA parts,const PolygonOut *polygons,
    int np,const GfsEnrichConfig *cfg,GfsEnrichAudit *audit){
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=np)return;
    GfsEnrichAudit a={};a.status=-10;
    double x[GPU_MAX_CORNERS],y[GPU_MAX_CORNERS];
    if(!gfs_moment_support_polygon(parts,polygons+i,i,x,y)){audit[i]=a;return;}
    GfsMomentPolygon p={x,y,polygons[i].n_corners,0};
    if(!gfs_enrich_build(p,cfg,parts.x[i],parts.y[i],&a.basis,&a.state,a.integral)){
        a.status=a.basis.status;audit[i]=a;return;
    }
    a.status=1;a.metrics[ER_MIN_RHO]=a.metrics[ER_MIN_P]=DBL_MAX;
    double integrated[4]={},exact[4]={};GfsEnrichQuadrature rule;
    if(!gfs_enrich_quadrature_init(p,cfg+1,&a.basis,&rule)){a.status=-11;audit[i]=a;return;}
    for(int q=0;q<gfs_enrich_nquad(p,cfg+1,&rule);q++){
        double xx,yy,w;if(!gfs_enrich_quad(p,cfg+1,&rule,q,&xx,&yy,&w)){a.status=-12;break;}
        if(w>0)gfs_enrich_assess_point(cfg+1,&a,xx,yy,w,integrated,exact);
    }
    for(int k=0;k<p.n;k++){
        int j=(k+1)%p.n;
        gfs_enrich_assess_point(cfg+1,&a,x[k],y[k],0,integrated,exact);
        for(int q=0;q<cfg[1].na;q++){
            double t=cfg[1].an[q];
            gfs_enrich_assess_point(cfg+1,&a,(1-t)*x[k]+t*x[j],(1-t)*y[k]+t*y[j],0,integrated,exact);
        }
    }
    a.metrics[ER_L1_PERR]/=a.basis.geometry.area;
    for(int c=0;c<4;c++){
        double scale=a.basis.geometry.area*a.basis.scale[c];
        a.metrics[ER_CONSERV_RHO+c]=fabs(integrated[c]-a.integral[c])/scale;
        a.metrics[ER_QUAD_RHO+c]=fabs(exact[c]-a.integral[c])/scale;
    }
    a.metrics[ER_QUAD_E_OVER_THERMAL]=fabs(exact[3]-a.integral[3])/
        (cfg->field.ic.pressure*a.basis.geometry.area/(cfg->field.ic.gamma-1));
    audit[i]=a;
}

extern "C" int gpu_enriched_project_initial(GPUContext *ctx,int np,
    const GfsEnrichConfig cfg[2],GfsEnrichAudit *audit){
    if(np==0)return 0;
    GfsEnrichAudit *da=NULL;GfsEnrichConfig *dc=NULL;
    CUDA_CHECK(cudaMalloc(&da,(size_t)np*sizeof(*da)));CUDA_CHECK(cudaMalloc(&dc,2*sizeof(*dc)));
    cudaStream_t stream=(cudaStream_t)ctx->stream;
    CUDA_CHECK(cudaMemcpyAsync(dc,cfg,2*sizeof(*dc),cudaMemcpyHostToDevice,stream));
    gfs_enrich_project_kernel<<<(np+31)/32,32,0,stream>>>(ctx->d_parts,ctx->d_polygons,np,dc,da);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpyAsync(audit,da,(size_t)np*sizeof(*da),cudaMemcpyDeviceToHost,stream));
    CUDA_CHECK(cudaStreamSynchronize(stream));
    CUDA_CHECK(cudaFree(da));CUDA_CHECK(cudaFree(dc));
    int bad=0;for(int i=0;i<np;i++)if(audit[i].status!=1)bad++;return bad;
}
