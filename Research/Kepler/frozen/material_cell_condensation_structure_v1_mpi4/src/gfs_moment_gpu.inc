/* Native CUDA initialization/realizability gate. No evolution yet: this
 * exports full-stride restartable moment states and terminates collectively.
 * The later transport must use gfs_moment_face_rhs/source, never legacy
 * pressure-only force/thermal update on the tagged payload.
 */
static __device__ unsigned long long d_gfs_moment_zero_edges;
static __device__ unsigned long long d_gfs_moment_rebuilt_zero_edges;
__device__ static int gfs_moment_support_polygon(ParticleSoA p,
    const PolygonOut *po,int i,double *x,double *y){
    int n=po->n_corners;if(n<3 || n>GPU_MAX_CORNERS)return 0;
    for(int k=0;k<n;k++){
        int prev=(k+n-1)%n,a=po->edge_neigh[prev],b=po->edge_neigh[k];
        if(a<0 || b<0 || a>=p.n_total || b>=p.n_total || p.w2[i]!=0 || p.w2[a]!=0 || p.w2[b]!=0)return 0;
        double ax=p.x[a]-p.x[i],ay=p.y[a]-p.y[i];
        double bx=p.x[b]-p.x[i],by=p.y[b]-p.y[i];
        double ab=.5*fma(ax,ax,ay*ay),bb=.5*fma(bx,bx,by*by);
        if(!gfs_moment_support_intersection(ax,ay,ab,bx,by,bb,x+k,y+k))return 0;
        double scale=fabs(ab)+fabs(bb)+(fabs(ax)+fabs(bx))*fabs(x[k])+
            (fabs(ay)+fabs(by))*fabs(y[k]);
        double err=fmax(fabs(ax*x[k]+ay*y[k]-ab),fabs(bx*x[k]+by*y[k]-bb));
        if(err>1024*DBL_EPSILON*scale)return 0;
    }
    return 1;
}
__global__ static void gfs_moment_initial_kernel(ParticleSoA p,
    const PolygonOut *polygons,const GfsMomentPayload *input,int np,
    GfsMomentIC ic,GfsMomentInitialAudit *audit){
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=np)return;
    GfsMomentInitialAudit a={};a.payload=input[i];
    GfsMomentPolygon poly={polygons[i].cx,polygons[i].cy,polygons[i].n_corners,1};
    for(int k=0;k<poly.n;k++){
        int j=(k+1)%poly.n;
        if(poly.x[k]==poly.x[j] && poly.y[k]==poly.y[j])atomicAdd(&d_gfs_moment_zero_edges,1ULL);
    }
    double vertices_x[GPU_MAX_CORNERS],vertices_y[GPU_MAX_CORNERS];
    a.valid=-7;
    if(!gfs_moment_support_polygon(p,polygons+i,i,vertices_x,vertices_y)){audit[i]=a;return;}
    poly.x=vertices_x;poly.y=vertices_y;
    for(int k=0;k<poly.n;k++){
        int j=(k+1)%poly.n;
        if(poly.x[k]==poly.x[j] && poly.y[k]==poly.y[j])atomicAdd(&d_gfs_moment_rebuilt_zero_edges,1ULL);
    }
    GfsMomentDensity geom={};
    a.valid=-1;if(!gfs_moment_geometry(poly,&geom)){audit[i]=a;return;}
    a.restart=a.payload.tag==GFS_MOMENT_TAG;
    if(!a.restart){
        gfs_moment_integrate_ic(poly,&geom,&ic,p.x[i],p.y[i],&a.payload.s,&a.thermal);
        GfsMomentPolygon coarse=poly;coarse.refine=0;
        gfs_moment_integrate_ic(coarse,&geom,&ic,p.x[i],p.y[i],&a.coarse,&a.thermal_coarse);
        a.payload.ox=p.x[i];a.payload.oy=p.y[i];
        a.payload.m0=a.payload.s.a[GFM_M];a.payload.mu=a.payload.m0/geom.area;
        a.payload.w=0;a.payload.tag=GFS_MOMENT_TAG;
    }
    /* This gate is W0 and uses the initial generator as the fixed origin.
     * It cannot be used to bless a moving/power-cell restart as implemented. */
    a.valid=-2;if(p.x[i]!=a.payload.ox || p.y[i]!=a.payload.oy ||
       a.payload.w!=0 || !(a.payload.m0>0) || !(a.payload.mu>0) ||
       a.payload.s.a[GFM_M]!=a.payload.m0){audit[i]=a;return;}
    GfsMomentPrimitive r={};a.valid=-3;if(!gfs_moment_primitive(&a.payload.s,&r)){audit[i]=a;return;}
    GfsMomentDensity density={};
    if(!gfs_moment_density_fit(poly,&a.payload.s,0,0,&density)){
        a.valid=-40+density.status;a.moment_residual=density.residual;
        a.iterations=density.iterations;audit[i]=a;return;
    }
    for(int k=0;k<5;k++)a.payload.lambda[k]=density.lambda[k];
    a.min_density=DBL_MAX;
    for(int k=0;k<gfs_moment_nquad(poly);k++){
        double x,y,w,rho,vx,vy;gfs_moment_quad(poly,&geom,k,&x,&y,&w);
        if(w==0)continue;
        double fitted=gfs_moment_density_value(&density,a.payload.s.a[GFM_M],x,y);
        a.valid=-5;if(!(fitted>0) || !isfinite(fitted)){audit[i]=a;return;}
        a.min_density=fmin(a.min_density,fitted);
        if(!a.restart){
            gfs_moment_ic_value(&ic,p.x[i]+x,p.y[i]+y,&rho,&vx,&vy);
            double ex=vx-r.v[0]-r.B[0]*(x-r.c[0])-r.B[1]*(y-r.c[1]);
            double ey=vy-r.v[1]-r.B[2]*(x-r.c[0])-r.B[3]*(y-r.c[1]);
            a.nonaffine_kinetic+=.5*w*rho*(ex*ex+ey*ey);
        }
    }
    for(int k=0;k<poly.n;k++){
        double rho=gfs_moment_density_value(&density,a.payload.s.a[GFM_M],poly.x[k],poly.y[k]);
        a.valid=-6;if(!(rho>0) || !isfinite(rho)){audit[i]=a;return;}
        a.min_density=fmin(a.min_density,rho);
    }
    a.energy_identity_error=a.restart?0:r.U-a.thermal-a.nonaffine_kinetic;
    a.moment_residual=density.residual;a.iterations=density.iterations;
    a.valid=1;audit[i]=a;
}

extern "C" int gpu_moment_initialize(GPUContext *ctx,int np,const GfsMomentIC *ic,
    const GfsMomentPayload *input,GfsMomentInitialAudit *audit){
    if(np<=0)return 0;
    GfsMomentPayload *di=NULL;GfsMomentInitialAudit *da=NULL;
    CUDA_CHECK(cudaMalloc(&di,(size_t)np*sizeof(*di)));
    CUDA_CHECK(cudaMalloc(&da,(size_t)np*sizeof(*da)));
    cudaStream_t stream=(cudaStream_t)ctx->stream;
    unsigned long long zero_edges=0;
    CUDA_CHECK(cudaMemcpyToSymbolAsync(d_gfs_moment_zero_edges,&zero_edges,sizeof(zero_edges),0,cudaMemcpyHostToDevice,stream));
    CUDA_CHECK(cudaMemcpyToSymbolAsync(d_gfs_moment_rebuilt_zero_edges,&zero_edges,sizeof(zero_edges),0,cudaMemcpyHostToDevice,stream));
    CUDA_CHECK(cudaMemcpyAsync(di,input,(size_t)np*sizeof(*di),cudaMemcpyHostToDevice,stream));
    gfs_moment_initial_kernel<<<(np+127)/128,128,0,stream>>>(ctx->d_parts,
        ctx->d_polygons,di,np,*ic,da);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cudaMemcpyAsync(audit,da,(size_t)np*sizeof(*da),cudaMemcpyDeviceToHost,stream));
    CUDA_CHECK(cudaMemcpyFromSymbolAsync(&zero_edges,d_gfs_moment_zero_edges,sizeof(zero_edges),0,cudaMemcpyDeviceToHost,stream));
    unsigned long long rebuilt_zero_edges=0;
    CUDA_CHECK(cudaMemcpyFromSymbolAsync(&rebuilt_zero_edges,d_gfs_moment_rebuilt_zero_edges,sizeof(rebuilt_zero_edges),0,cudaMemcpyDeviceToHost,stream));
    CUDA_CHECK(cudaStreamSynchronize(stream));CUDA_CHECK(cudaFree(di));CUDA_CHECK(cudaFree(da));
    int bad=0;for(int k=0;k<np;k++)if(audit[k].valid!=1)bad++;
    fprintf(stderr,"[MOMENT_GEOMETRY] cells=%d raw_zero_edges=%llu rebuilt_zero_edges=%llu invalid=%d\n",np,zero_edges,rebuilt_zero_edges,bad);
    int shown=0;
    for(int k=0;k<np;k++)if(audit[k].valid==-1 && shown++<8){
        PolygonOut po;CUDA_CHECK(cudaMemcpy(&po,ctx->d_polygons+k,sizeof(po),cudaMemcpyDeviceToHost));
        GfsMomentPolygon poly={po.cx,po.cy,po.n_corners,1};GfsMomentDensity geom={};
        int host_valid=gfs_moment_geometry(poly,&geom);
        fprintf(stderr,"[MOMENT_POLYGON_DIAG] id=%lld nc=%d host_valid=%d area=%.17g centroid=%.17g,%.17g\n",
            ctx->h_parts.indx[k],po.n_corners,host_valid,geom.area,geom.gx,geom.gy);
        for(int e=0;e<po.n_corners;e++)
            fprintf(stderr,"[MOMENT_POLYGON_VERTEX] id=%lld k=%d x=%.17g y=%.17g\n",ctx->h_parts.indx[k],e,po.cx[e],po.cy[e]);
    }
    return bad;
}
