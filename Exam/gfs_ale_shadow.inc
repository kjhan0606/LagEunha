/* Read-only native CUDA diagnostic. Each global pair is evaluated once, by
 * the rank owning its lower-ID endpoint. The packet carries one outward
 * flux; the two cell contributions are its exact negative and positive.
 * No atomic updates of hydro state, caches or physical ledgers occur.
 */
#include <stdint.h>
#include <vector>
#include <set>
#include <utility>
#include "gfs_ale_flux.h"

struct GfsAlePair { int i,j,f; };
struct GfsAleRecord { double a[29]; int flags[4]; };
struct GfsAleAffine { int active; double a[10]; };

__global__ static void gfs_ale_shadow_kernel(ParticleSoA p,FaceCSR face,
    const GfsAlePair *pairs,int count,double gamma,const double *native_p,
    const double *native_work,GfsAleAffine affine,GfsAleRecord *records){
    int k=blockIdx.x*blockDim.x+threadIdx.x;if(k>=count)return;
    GfsAlePair pair=pairs[k];int i=pair.i,j=pair.j,f=pair.f;
    GfsAleRecord rec={};rec.flags[0]=1;rec.flags[3]=affine.active;
    double x1=face.c1x[f],y1=face.c1y[f],x2=face.c2x[f],y2=face.c2y[f];
    rec.a[0]=x1;rec.a[1]=y1;rec.a[2]=x2;rec.a[3]=y2;
    double sx=y2-y1,sy=-(x2-x1),area=hypot(sx,sy);
    double dx=p.x[j]-p.x[i],dy=p.y[j]-p.y[i],dist=hypot(dx,dy);
    if(!(area>0) || !(dist>0)){rec.flags[0]=0;records[k]=rec;return;}
    double nx=sx/area,ny=sy/area,ex=dx/dist,ey=dy/dist;
    double dvx=p.vx[j]-p.vx[i],dvy=p.vy[j]-p.vy[i];
    for(int g=0;g<2;g++){
        double xi=(g?1.:-1.)*.57735026918962576451;
        double ox=.5*(x1+x2)+.5*xi*(x2-x1),oy=.5*(y1+y2)+.5*xi*(y2-y1);
        double delta=-(dvx*(ox-.5*dx)+dvy*(oy-.5*dy))/dist;
        double wx=.5*(p.vx[i]+p.vx[j])+delta*ex,wy=.5*(p.vy[i]+p.vy[j])+delta*ey;
        double rx=ox-dx,ry=oy-dy;
        GfsAlePrimitive l={p.den[i]+p.dRhodx[i]*ox+p.dRhody[i]*oy,
            p.pressure[i]+p.dPdx[i]*ox+p.dPdy[i]*oy,
            p.vx[i]+p.gUxx[i]*ox+p.gUxy[i]*oy,
            p.vy[i]+p.gUyx[i]*ox+p.gUyy[i]*oy};
        GfsAlePrimitive r={p.den[j]+p.dRhodx[j]*rx+p.dRhody[j]*ry,
            p.pressure[j]+p.dPdx[j]*rx+p.dPdy[j]*ry,
            p.vx[j]+p.gUxx[j]*rx+p.gUxy[j]*ry,
            p.vy[j]+p.gUyx[j]*rx+p.gUyy[j]*ry};
        /* Match the existing face reconstruction clamps, recording every
         * intervention. This does not change the actual native inputs. */
        if(l.rho<1.e-10){l.rho=1.e-10;rec.flags[1]|=1;}
        if(r.rho<1.e-10){r.rho=1.e-10;rec.flags[1]|=2;}
        if(l.p<1.e-10){l.p=1.e-10;rec.flags[1]|=4;}
        if(r.p<1.e-10){r.p=1.e-10;rec.flags[1]|=8;}
        GfsAleFlux z=gfs_ale_flux(l,r,nx,ny,wx,wy,gamma);
        if(!z.valid)rec.flags[0]=0;rec.flags[2]|=1<<z.mode;
        rec.a[4]+=.5*area*(wx*nx+wy*ny);
        for(int c=0;c<4;c++){
            rec.a[5+c]+=.5*area*z.flux[c];rec.a[9+c]+=.5*area*z.pressure[c];
            rec.a[13+c]+=.5*area*z.advection[c];rec.a[17+c]+=.5*area*z.residual[c];
        }
        if(affine.active){
            double xx=p.x[i]+ox-affine.a[2],yy=p.y[i]+oy-affine.a[3];
            GfsAlePrimitive q={affine.a[0],affine.a[1],
                affine.a[4]+affine.a[6]*xx+affine.a[7]*yy,
                affine.a[5]+affine.a[8]*xx+affine.a[9]*yy};
            GfsAleFlux ref=gfs_ale_flux(q,q,nx,ny,wx,wy,gamma);
            if(!ref.valid)rec.flags[0]=0;
            for(int c=0;c<4;c++)rec.a[25+c]+=.5*area*ref.flux[c];
        }
    }
    rec.a[21]=0;rec.a[22]=native_p[f]*sx;rec.a[23]=native_p[f]*sy;
    rec.a[24]=-native_work[f];records[k]=rec;
}

static void gfs_ale_shadow_capture(GPUContext *ctx,int np,int nt,int rank,int call,
    const GPUPhysicsParams *par,const double *native_p,const double *native_work){
    if(par->av_mode!=5 || !par->use_muscl || !par->phase1 || par->alphavis!=0 ||
       par->betavis!=0 || par->nu_phys!=0 || par->hyperv_alpha!=0 || par->cd_amax!=0){
        fprintf(stderr,"[ALE_SHADOW] unsupported force contract\n");exit(1);
    }
    GfsAleAffine affine={};const char *s=getenv("GFS_ALE_SHADOW_AFFINE");
    if(s && s[0]){
        affine.active=1;
        for(int n=0;n<10;n++){
            char *end=NULL;affine.a[n]=strtod(s,&end);
            if(end==s || !isfinite(affine.a[n]) || (n<9?*end!=',':*end!=0)){
                fprintf(stderr,"[ALE_SHADOW] expected rho,P,cx,cy,vx,vy,Bxx,Bxy,Byx,Byy\n");exit(1);
            }
            s=end+(n<9);
        }
        if(!(affine.a[0]>0) || !(affine.a[1]>0))exit(1);
    }
    gpu_download_face_csr(ctx,np,ctx->d_faces.n_faces_total,
        ctx->h_faces.face_offset,ctx->h_faces.c1x,ctx->h_faces.c1y,
        ctx->h_faces.c2x,ctx->h_faces.c2y,ctx->h_faces.neighbor_idx,
        ctx->h_faces.kp_idx,ctx->h_faces.km_idx,ctx->h_faces.boundary_kind);
#define GFS_ALE_COPY(field,type) CUDA_CHECK(cudaMemcpy(ctx->h_parts.field,ctx->d_parts.field,(size_t)nt*sizeof(type),cudaMemcpyDeviceToHost))
    GFS_ALE_COPY(indx,long long);GFS_ALE_COPY(x,double);GFS_ALE_COPY(y,double);
    GFS_ALE_COPY(vx,double);GFS_ALE_COPY(vy,double);GFS_ALE_COPY(mass,double);
    GFS_ALE_COPY(den,float);GFS_ALE_COPY(pressure,float);GFS_ALE_COPY(ie,float);
    GFS_ALE_COPY(volume,float);GFS_ALE_COPY(w2,float);
    GFS_ALE_COPY(ax_out,double);GFS_ALE_COPY(ay_out,double);GFS_ALE_COPY(die_out,float);
#undef GFS_ALE_COPY
    std::vector<GfsAlePair> pairs;std::set<std::pair<long long,long long> > seen;
    long long skipped_walls=0,duplicates=0;double tiny_length=0;
    for(int i=0;i<np;i++){
        if(ctx->h_parts.w2[i]!=0){fprintf(stderr,"[ALE_SHADOW] requires W0\n");exit(1);}
        for(int f=ctx->h_faces.face_offset[i];f<ctx->h_faces.face_offset[i+1];f++){
            int j=ctx->h_faces.neighbor_idx[f];
            int kind=ctx->h_faces.boundary_kind[f];
            if(kind==GFS_FACE_VACUUM){
                fprintf(stderr,"[ALE_SHADOW] vacuum faces are outside this diagnostic contract\n");exit(1);
            }
            if(j<0 || j>=nt || kind!=GFS_FACE_INTERIOR){skipped_walls++;continue;}
            long long a=ctx->h_parts.indx[i],b=ctx->h_parts.indx[j];
            if(a>=b)continue;
            if(ctx->h_parts.w2[j]!=0){fprintf(stderr,"[ALE_SHADOW] nonzero neighbor W\n");exit(1);}
            if(!seen.insert(std::make_pair(a,b)).second){duplicates++;continue;}
            pairs.push_back({i,j,f});
            double length=hypot(ctx->h_faces.c2x[f]-ctx->h_faces.c1x[f],ctx->h_faces.c2y[f]-ctx->h_faces.c1y[f]);
            if(length<1.e-14)tiny_length+=length;
        }
    }
    if(duplicates || skipped_walls){fprintf(stderr,"[ALE_SHADOW] duplicate/wall pair contract failed\n");exit(1);}
    GfsAlePair *dp=NULL;GfsAleRecord *dr=NULL;
    size_t count=pairs.size();std::vector<GfsAleRecord> records(count);
    if(count){
        CUDA_CHECK(cudaMalloc(&dp,count*sizeof(*dp)));CUDA_CHECK(cudaMalloc(&dr,count*sizeof(*dr)));
        CUDA_CHECK(cudaMemcpy(dp,pairs.data(),count*sizeof(*dp),cudaMemcpyHostToDevice));
        gfs_ale_shadow_kernel<<<(count+255)/256,256,0,(cudaStream_t)ctx->stream>>>(ctx->d_parts,ctx->d_faces,dp,(int)count,par->Gamma,native_p,native_work,affine,dr);
        CUDA_CHECK(cudaGetLastError());CUDA_CHECK(cudaStreamSynchronize((cudaStream_t)ctx->stream));
        CUDA_CHECK(cudaMemcpy(records.data(),dr,count*sizeof(*dr),cudaMemcpyDeviceToHost));
        CUDA_CHECK(cudaFree(dp));CUDA_CHECK(cudaFree(dr));
    }
    char path[160];snprintf(path,sizeof(path),"ale_shadow_r%d_call%d.bin",rank,call);
    FILE *out=fopen(path,"wx");if(!out){perror(path);exit(1);}
    uint64_t header[4]={UINT64_C(0x474653414c453031),1,(uint64_t)count,29};
    fwrite(header,sizeof(header),1,out);long invalid=0,clamped=0;
    for(size_t k=0;k<count;k++){
        long long ids[2]={ctx->h_parts.indx[pairs[k].i],ctx->h_parts.indx[pairs[k].j]};
        fwrite(ids,sizeof(ids),1,out);fwrite(records[k].a,sizeof(double),29,out);
        fwrite(records[k].flags,sizeof(int),4,out);
        if(!records[k].flags[0])invalid++;if(records[k].flags[1])clamped++;
    }
    if(ferror(out) || fclose(out)){fprintf(stderr,"[ALE_SHADOW] packet write failed\n");exit(1);}
    snprintf(path,sizeof(path),"ale_cells_r%d_call%d.tsv",rank,call);out=fopen(path,"wx");if(!out)exit(1);
    fprintf(out,"id\tx\ty\tvx\tvy\tmass\tV\trho\tP\tU\tpxdot_native\tpydot_native\tEdot_native\n");
    for(int i=0;i<np;i++){
        const ParticleSoA *p=&ctx->h_parts;
        fprintf(out,"%lld\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\t%.17g\n",
            p->indx[i],p->x[i],p->y[i],p->vx[i],p->vy[i],p->mass[i],(double)p->volume[i],
            (double)p->den[i],(double)p->pressure[i],(double)p->ie[i],
            p->mass[i]*p->ax_out[i],p->mass[i]*p->ay_out[i],(double)p->die_out[i]);
    }
    if(ferror(out) || fclose(out))exit(1);
    snprintf(path,sizeof(path),"ale_meta_r%d_call%d.json",rank,call);out=fopen(path,"wx");if(!out)exit(1);
    fprintf(out,"{\"rank\":%d,\"call\":%d,\"Gamma\":%.17g,\"pairs\":%zu,\"invalid\":%ld,\"clamped_pairs\":%ld,\"tiny_face_length_sum\":%.17g,\"affine_active\":%d,\"record_bytes\":264,\"data_endian\":\"little\",\"affine\":[",rank,call,par->Gamma,count,invalid,clamped,tiny_length,affine.active);
    for(int k=0;k<10;k++)fprintf(out,"%s%.17g",k?",":"",affine.a[k]);
    fprintf(out,"],\"contract\":\"outward flux lower global ID to upper, opposite cell rates from one packet; full=pressure+advection+explicit_Riemann_residual; no state update\"}\n");
    if(ferror(out) || fclose(out))exit(1);
    fprintf(stderr,"[ALE_SHADOW] rank=%d call=%d pairs=%zu invalid=%ld clamped=%ld state_updates=0\n",rank,call,count,invalid,clamped);
    if(invalid)exit(1);
}
