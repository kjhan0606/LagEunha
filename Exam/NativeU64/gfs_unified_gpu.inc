/* Native CUDA adapter for the shared W0 operator.
 * First qualification scope: fully owned one-rank cells and static wall
 * faces. Periodic image synchronization and the new KDK driver are separate
 * integration gates. No old solver path is silently reinterpreted.
 */
#include "../gfs_variational_core.h"
#include "../gfs_velocity_fit.h"
/* ID mapping is independent of the permutation introduced by migration. */
extern "C" const int *gpu_unified_owner_aliases(GPUContext *ctx,int n,int nt,
        double lx,double ly)
{
    static std::vector<int> alias,values;
    static std::vector<long long> keys;
    if(n<1||nt<n){fprintf(stderr,"[UNIFIED_ALIAS] invalid counts\n");exit(99);}
    alias.resize(nt);
    for(int i=0;i<n;i++)alias[i]=i;
    if(nt==n)return alias.data();
    if(!(lx>0)||!(ly>0)){fprintf(stderr,"[UNIFIED_ALIAS] periodic lengths required\n");exit(99);}
    size_t cap=1;while(cap<2*(size_t)n)cap*=2;
    keys.assign(cap,-1);values.resize(cap);
    ParticleSoA &p=ctx->h_parts;
    for(int i=0;i<n;i++){
        long long id=p.indx[i];size_t k=(uint64_t(id)*UINT64_C(11400714819323198485))&(cap-1);
        if(id<0){fprintf(stderr,"[UNIFIED_ALIAS] invalid owned ID\n");exit(99);}
        while(keys[k]>=0){if(keys[k]==id){fprintf(stderr,"[UNIFIED_ALIAS] duplicate owned ID\n");exit(99);}k=(k+1)&(cap-1);}
        keys[k]=id;values[k]=i;
    }
    for(int i=n;i<nt;i++){
        long long id=p.indx[i];size_t k=(uint64_t(id)*UINT64_C(11400714819323198485))&(cap-1);
        while(keys[k]>=0&&keys[k]!=id)k=(k+1)&(cap-1);
        if(keys[k]!=id){fprintf(stderr,"[UNIFIED_ALIAS] padding ID has no local owner\n");exit(99);}
        int j=values[k];double dx=p.x[i]-p.x[j],dy=p.y[i]-p.y[j];
        double tol=128*DBL_EPSILON*(lx+ly+fabs(p.x[i])+fabs(p.y[i])+fabs(p.x[j])+fabs(p.y[j]));
        if(fabs(dx-lx*nearbyint(dx/lx))>tol||fabs(dy-ly*nearbyint(dy/ly))>tol){
            fprintf(stderr,"[UNIFIED_ALIAS] padding is not a periodic translation\n");exit(99);
        }
        alias[i]=j;
    }
    return alias.data();
}
static __global__ void gfs_unified_alias_kernel(ParticleSoA p,const int *owner,
        int n,int nt,int gradients)
{
    int i=n+blockIdx.x*blockDim.x+threadIdx.x;if(i>=nt)return;int j=owner[i];
    if(!gradients){
        p.mass[i]=p.mass[j];p.vx[i]=p.vx[j];p.vy[i]=p.vy[j];
        p.weight_U64[i]=p.weight_U64[j];p.weight_volume64[i]=p.weight_volume64[j];
        p.weight_W[i]=p.weight_W[j];
    }else{
        p.gUxx[i]=p.gUxx[j];p.gUxy[i]=p.gUxy[j];p.gUyx[i]=p.gUyx[j];p.gUyy[i]=p.gUyy[j];
    }
}
static __device__ unsigned long long gfs_unified_fit_failures;

__global__ void gfs_unified_fit_kernel(ParticleSoA p,FaceCSR f,int n,int nt)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;
    GfsVelocityFit fit;gfs_velocity_fit_init(&fit,2);int valid=1;
    for(int k=f.face_offset[i];k<f.face_offset[i+1];k++){
        int j=f.neighbor_idx[k];double r[3]={0},dv[3]={0};
        if(j>=0&&j<nt&&!f.is_ghost[k]){
            r[0]=p.x[j]-p.x[i];r[1]=p.y[j]-p.y[i];
            dv[0]=p.vx[j]-p.vx[i];dv[1]=p.vy[j]-p.vy[i];
        }else if(j<=-1000){
            double sx=f.c2y[k]-f.c1y[k],sy=f.c1x[k]-f.c2x[k],area=hypot(sx,sy);
            if(!(area>0)){valid=0;break;}
            double nx=sx/area,ny=sy/area;
            double distance=.5*((f.c1x[k]+f.c2x[k])*nx+(f.c1y[k]+f.c2y[k])*ny);
            double vn=p.vx[i]*nx+p.vy[i]*ny;
            r[0]=2*distance*nx;r[1]=2*distance*ny;
            dv[0]=-2*vn*nx;dv[1]=-2*vn*ny;
        }else{valid=0;break;}
        if(!gfs_velocity_fit_add(&fit,r,dv)){valid=0;break;}
    }
    double G[3][3];
    if(!valid||!gfs_velocity_fit_solve(&fit,G,NULL)){
        atomicAdd(&gfs_unified_fit_failures,1ull);
        p.gUxx[i]=p.gUxy[i]=p.gUyx[i]=p.gUyy[i]=NAN;return;
    }
    p.gUxx[i]=G[0][0];p.gUxy[i]=G[0][1];p.gUyx[i]=G[1][0];p.gUyy[i]=G[1][1];
}

__global__ void gfs_unified_force_kernel(ParticleSoA p,FaceCSR f,int n,
        double gamma,double courant,int implicit,int nt)
{
    int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;
    double mass=p.mass[i],V=p.weight_volume64[i],U=p.weight_U64[i];
    double rho=mass/V,Pi=(gamma-1)*U/V,ci=sqrt(gamma*Pi/rho),h=sqrt(V);
    double vi[3]={p.vx[i],p.vy[i],0},F[3]={0},vdot=0,heat=0,damping=0,sweep=0;
    double dti=courant*h/ci;
    if(!(mass>0)||!(V>0)||!(U>0)||!isfinite(ci)||!isfinite(p.gUxx[i]))goto invalid;
    for(int k=f.face_offset[i];k<f.face_offset[i+1];k++){
        int j=f.neighbor_idx[k];
        double S[3]={f.c2y[k]-f.c1y[k],f.c1x[k]-f.c2x[k],0};
        double area=hypot(S[0],S[1]);
        double centroid[3]={.5*(f.c1x[k]+f.c2x[k]),.5*(f.c1y[k]+f.c2y[k]),0};
        if(!(area>0))goto invalid;
        if(j<=-1000){
            double nx=S[0]/area,ny=S[1]/area;
            double distance=centroid[0]*nx+centroid[1]*ny;
            double vn=vi[0]*nx+vi[1]*ny;
            double normal_derivative=nx*(p.gUxx[i]*nx+p.gUxy[i]*ny)+ny*(p.gUyx[i]*nx+p.gUyy[i]*ny);
            double r[3]={2*distance*nx,2*distance*ny,0};
            double vg[3]={vi[0]-2*vn*nx,vi[1]-2*vn*ny,0};
            double jump[3]={-2*(vn+distance*normal_derivative)*nx,
                            -2*(vn+distance*normal_derivative)*ny,0};
            GfsVariationalFace b;if(!gfs_variational_face(2,r,centroid,area,&b))goto invalid;
            double div=p.gUxx[i]+p.gUyy[i],curl=p.gUyx[i]-p.gUxy[i];
            double sf=gfs_variational_compression_fraction(div,curl,div,curl);
            double q=gfs_variational_pair_pressure(&b,vi,vg,jump,rho,ci,rho,ci,sf,h,h);
            double rate=-(vi[0]*S[0]+vi[1]*S[1]);vdot+=rate;sweep+=fabs(rate);
            if(!implicit){F[0]-=q*S[0];F[1]-=q*S[1];heat+=q*vn*area;}
            if(q>0&&vn>0)damping+=2*q*area/(vn*mass);
            continue;
        }
        if(j<0||j>=nt||f.is_ghost[k])goto invalid;
        double r[3]={p.x[j]-p.x[i],p.y[j]-p.y[i],0};
        double vj[3]={p.vx[j],p.vy[j],0},jump[3];
        jump[0]=vj[0]-vi[0]-.5*((p.gUxx[i]+p.gUxx[j])*r[0]+(p.gUxy[i]+p.gUxy[j])*r[1]);
        jump[1]=vj[1]-vi[1]-.5*((p.gUyx[i]+p.gUyx[j])*r[0]+(p.gUyy[i]+p.gUyy[j])*r[1]);jump[2]=0;
        double Vj=p.weight_volume64[j],rhoj=p.mass[j]/Vj;
        double Pj=(gamma-1)*p.weight_U64[j]/Vj,cj=sqrt(gamma*Pj/rhoj);
        GfsVariationalFace b;if(!gfs_variational_face(2,r,centroid,area,&b))goto invalid;
        double fp[3];gfs_variational_pressure(&b,Pi,Pj,fp);
        double rate=gfs_variational_sweep(&b,vi,vj);vdot+=rate;sweep+=fabs(rate);
        double sf=gfs_variational_compression_fraction(p.gUxx[i]+p.gUyy[i],p.gUyx[i]-p.gUxy[i],
                    p.gUxx[j]+p.gUyy[j],p.gUyx[j]-p.gUxy[j]);
        double q=gfs_variational_pair_pressure(&b,vi,vj,jump,rho,ci,rhoj,cj,sf,h,sqrt(Vj));
        double relative=(vj[0]-vi[0])*b.normal[0]+(vj[1]-vi[1])*b.normal[1];
        for(int a=0;a<2;a++)F[a]+=fp[a]-(implicit?0:q*area*b.normal[a]);
        if(!implicit)heat+=-.5*q*area*relative;
        if(q>0&&relative<0)damping+=2*q*area/((-relative)*mass);
        dti=fmin(dti,courant*h/(ci+fabs(relative)));
        if(relative<0)dti=fmin(dti,.2*b.separation/(-relative));
    }
    if(!implicit&&damping>0)dti=fmin(dti,courant/damping);
    if(sweep>0)dti=fmin(dti,courant*V/sweep);
    {
        double Udot=-Pi*vdot+heat;
        double Edot=vi[0]*F[0]+vi[1]*F[1]+Udot;
        if(!isfinite(Udot)||!isfinite(Edot)||!(dti>0)||!isfinite(dti))goto invalid;
        p.ax_out[i]=F[0]/mass;p.ay_out[i]=F[1]/mass;
        p.weight_Udot64[i]=Udot;p.weight_Edot64[i]=Edot;p.unified_heat64[i]=heat;
        p.die_out[i]=(float)Edot;p.dK_out[i]=0;p.dt_out[i]=(float)dti;
        p.vsig_max_out[i]=ci;p.vsmoothx_out[i]=p.vsmoothy_out[i]=0;
        p.weight_force[i]=p.weight_Hrow[i]=p.weight_Hcol[i]=0;
        return;
    }
invalid:
    p.ax_out[i]=p.ay_out[i]=p.weight_Udot64[i]=p.weight_Edot64[i]=NAN;
    p.unified_heat64[i]=NAN;p.dt_out[i]=NAN;
    atomicAdd(&gfs_unified_fit_failures,1ull);
}

static double gfs_unified_gpu_launch(GPUContext *ctx,int n,int nt,
        const GPUPhysicsParams *params)
{
    if(nt<n||params->weight_active||
       !params->phase1||params->nu_phys!=0||params->alphavis!=0||
       params->betavis!=0||params->cd_amax!=0||params->hyperv_alpha!=0){
        fprintf(stderr,"[UNIFIED_GPU] requires one-rank U64 W0 with explicit local aliases, one numerical viscosity, no extra stress\n");exit(99);
    }
    if(ctx->weight_geometry_count!=nt||memcmp(ctx->weight_geometry_snapshot,
            ctx->h_parts.weight_W,(size_t)nt*sizeof(double))){
        fprintf(stderr,"[UNIFIED_GPU] immutable geometry/weight snapshot mismatch\n");exit(99);
    }
    for(int i=0;i<nt;i++)if(ctx->h_parts.weight_W[i]!=0){
        fprintf(stderr,"[UNIFIED_GPU] nonzero weight is unsupported\n");exit(99);
    }
    cudaStream_t stream=(cudaStream_t)ctx->stream;
    const int *alias=gpu_unified_owner_aliases(ctx,n,nt,params->domain_x,params->domain_y);
    static int *device_alias=NULL;static int alias_capacity=0;
    if(nt>alias_capacity){if(device_alias)CUDA_CHECK(cudaFree(device_alias));CUDA_CHECK(cudaMalloc(&device_alias,(size_t)nt*sizeof(int)));alias_capacity=nt;}
    CUDA_CHECK(cudaMemcpyAsync(device_alias,alias,(size_t)nt*sizeof(int),cudaMemcpyHostToDevice,stream));
    if(nt>n){gfs_unified_alias_kernel<<<(nt-n+127)/128,128,0,stream>>>(ctx->d_parts,device_alias,n,nt,0);CUDA_CHECK(cudaGetLastError());}
    unsigned long long zero=0,failures=0;
    CUDA_CHECK(cudaMemcpyToSymbolAsync(gfs_unified_fit_failures,&zero,sizeof(zero),0,cudaMemcpyHostToDevice,stream));
    gfs_unified_fit_kernel<<<(n+127)/128,128,0,stream>>>(ctx->d_parts,ctx->d_faces,n,nt);
    if(nt>n){gfs_unified_alias_kernel<<<(nt-n+127)/128,128,0,stream>>>(ctx->d_parts,device_alias,n,nt,1);CUDA_CHECK(cudaGetLastError());}
    CUDA_CHECK(cudaGetLastError());
    const char *implicit_env=getenv("GFS_UNIFIED_IMPLICIT");
    int implicit=implicit_env&&!strcmp(implicit_env,"1");
    gfs_unified_force_kernel<<<(n+127)/128,128,0,stream>>>(ctx->d_parts,ctx->d_faces,n,params->Gamma,params->Courant,implicit,nt);
    CUDA_CHECK(cudaGetLastError());
    CUDA_CHECK(cub::DeviceReduce::Min(ctx->d_cub_temp,ctx->cub_temp_bytes,
        ctx->d_parts.dt_out,(float*)ctx->d_cub_min_out,n,stream));
    CUDA_CHECK(cudaMemcpyAsync(ctx->h_cub_min_out,ctx->d_cub_min_out,sizeof(float),cudaMemcpyDeviceToHost,stream));
    CUDA_CHECK(cudaMemcpyFromSymbolAsync(&failures,gfs_unified_fit_failures,sizeof(failures),0,cudaMemcpyDeviceToHost,stream));
    CUDA_CHECK(cudaStreamSynchronize(stream));
    if(failures){fprintf(stderr,"[UNIFIED_GPU] invalid geometry/fit/RHS cells=%llu\n",failures);exit(99);}
    const char *probe=getenv("GFS_WEIGHT_NATIVE_RHS_ONLY");
    if(probe&&!strcmp(probe,"1")){
        const size_t bytes=(size_t)n*sizeof(double);
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.unified_heat64,ctx->d_parts.unified_heat64,bytes,cudaMemcpyDeviceToHost,stream));
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.ax_out,ctx->d_parts.ax_out,bytes,cudaMemcpyDeviceToHost,stream));
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.ay_out,ctx->d_parts.ay_out,bytes,cudaMemcpyDeviceToHost,stream));
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.dt_out,ctx->d_parts.dt_out,(size_t)n*sizeof(float),cudaMemcpyDeviceToHost,stream));
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.gUxx,ctx->d_parts.gUxx,bytes,cudaMemcpyDeviceToHost,stream));
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.gUxy,ctx->d_parts.gUxy,bytes,cudaMemcpyDeviceToHost,stream));
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.gUyx,ctx->d_parts.gUyx,bytes,cudaMemcpyDeviceToHost,stream));
        CUDA_CHECK(cudaMemcpyAsync(ctx->h_parts.gUyy,ctx->d_parts.gUyy,bytes,cudaMemcpyDeviceToHost,stream));
        CUDA_CHECK(cudaStreamSynchronize(stream));
        struct Row {long long id;double v[8];};
        static_assert(sizeof(Row)==72,"Unified RHS auxiliary record");
        FILE *out=fopen("unified_rhs_aux.bin","wbx");
        if(!out){fprintf(stderr,"[UNIFIED_GPU] cannot create unique RHS audit\n");exit(99);}
        for(int i=0;i<n;i++){
            Row row={ctx->h_parts.indx[i],{ctx->h_parts.unified_heat64[i],ctx->h_parts.dt_out[i],
                ctx->h_parts.gUxx[i],ctx->h_parts.gUxy[i],ctx->h_parts.gUyx[i],ctx->h_parts.gUyy[i],
                ctx->h_parts.ax_out[i],ctx->h_parts.ay_out[i]}};
            if(fwrite(&row,sizeof(row),1,out)!=1){fprintf(stderr,"[UNIFIED_GPU] RHS audit write failed\n");exit(99);}
        }
        if(ferror(out)||fclose(out)){fprintf(stderr,"[UNIFIED_GPU] RHS audit close failed\n");exit(99);}
    }
    return (double)(*(float*)ctx->h_cub_min_out);
}
