/* Read-only analytic reference integration on actual saved native power
 * polygons. Positive centroid triangles, with independent order 8/16
 * quadrature controls. No fluid evolution or checkpoint modification. */
#define main archived_fv_evolution_main
#include "fv_native_evolution.c"
#undef main
static void gauss(int n,double *x,double *w){
    for(int i=0;i<(n+1)/2;i++){
        double z=cos(M_PI*(i+.75)/(n+.5)),derivative=0;
        for(int k=0;k<30;k++){
            double p=1,p0=0;
            for(int j=1;j<=n;j++){double p1=p0;p0=p;p=((2*j-1)*z*p0-(j-1)*p1)/j;}
            derivative=n*(z*p-p0)/(z*z-1);double dz=p/derivative;z-=dz;
            if(fabs(dz)<4*DBL_EPSILON)break;
        }
        x[i]=.5*(1-z);x[n-1-i]=.5*(1+z);
        w[i]=w[n-1-i]=1/((1-z*z)*derivative*derivative);
    }
}
static int cell_reference(const GfsFvMesh *m,int i,int problem,int order,double out[5]){
    double node[16],weight[16];gauss(order,node,weight);memset(out,0,5*sizeof(double));
    const GfsFvGeometry *g=m->geometry+i;
    for(int k=0;k<g->faces;k++){
        const GfsFvNativeFace *f=g->face+k;double a[2],b[2];
        for(int d=0;d<2;d++){
            double tangent=d==0?-f->normal[1]:f->normal[0];
            a[d]=f->centroid[d]-.5*f->area*tangent-g->centroid[d];
            b[d]=f->centroid[d]+.5*f->area*tangent-g->centroid[d];
        }
        double determinant=fabs(a[0]*b[1]-a[1]*b[0]);
        for(int u=0;u<order;u++)for(int v=0;v<order;v++){
            double x[3]={0};
            for(int d=0;d<2;d++){
                x[d]=m->cell[i].q[d]+g->centroid[d]+node[u]*a[d]+(1-node[u])*node[v]*b[d];
                x[d]-=floor(x[d]);
            }
            GfsFvPrimitive p;analytic(problem,x,&p);
            double dw=determinant*(1-node[u])*weight[u]*weight[v];
            out[0]+=dw;out[1]+=dw*p.rho;out[2]+=dw*p.rho*p.u[0];out[3]+=dw*p.rho*p.u[1];out[4]+=dw*p.p;
        }
    }
    if(fabs(out[0]/g->V-1)>1e-10||!(out[1]>0))return 1;
    out[2]/=out[1];out[3]/=out[1];out[1]/=out[0];out[4]/=out[0];return 0;
}
int main(int argc,char **argv){
    if(argc!=3)return 2;FILE *in=fopen(argv[1],"rb");FvCheckpoint h;
    if(!in||fread(&h,sizeof(h),1,in)!=1||fclose(in)||h.magic!=magic||h.dim!=2||h.count>1000000||h.side<4||(h.problem!=GRESHO&&h.problem!=COLD_SHEAR))return 2;
    GfsFvMesh m={0};m.dim=h.dim;m.n=h.count;m.gamma=h.gamma;memcpy(m.bc,h.bc,sizeof(m.bc));memcpy(m.lo,h.lo,sizeof(m.lo));memcpy(m.L,h.L,sizeof(m.L));m.cell=calloc(m.n,sizeof(*m.cell));
    double time,boundary[5],initial[5];uint64_t steps;
    if(!m.cell||load_state(argv[1],&m,h.side,h.problem,&time,&steps,boundary,initial)||gfs_fv_geometry(&m))return 3;
    FILE *out=fopen(argv[2],"wx");if(!out)return 4;
    fprintf(out,"id,V_ref,rho_ref,vx_mean_ref,vy_mean_ref,P_mean_ref,vx_node_ref,vy_node_ref,P_node_ref,quadrature_dv,quadrature_dP\n");
    double max_dv=0,max_dp=0;
    for(int i=0;i<m.n;i++){
        double coarse[5],fine[5],x[3];if(cell_reference(&m,i,h.problem,8,coarse)||cell_reference(&m,i,h.problem,16,fine))return 5;
        for(int d=0;d<3;d++)x[d]=m.cell[i].q[d]-floor(m.cell[i].q[d]);GfsFvPrimitive p;analytic(h.problem,x,&p);
        double dv=hypot(coarse[2]-fine[2],coarse[3]-fine[3]),dp=fabs(coarse[4]-fine[4]);max_dv=fmax(max_dv,dv);max_dp=fmax(max_dp,dp);
        fprintf(out,"%llu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",(unsigned long long)m.cell[i].id,fine[0],fine[1],fine[2],fine[3],fine[4],p.u[0],p.u[1],p.p,dv,dp);
    }
    int bad=fclose(out);gfs_fv_geometry_free(&m);free(m.cell);
    printf("{\"time\":%.17g,\"cells\":%llu,\"reference\":\"analytic actual native polygon mean and node\",\"quadrature_orders\":[8,16],\"maximum_mean_velocity_difference\":%.17g,\"maximum_mean_pressure_difference\":%.17g}\n",time,(unsigned long long)h.count,max_dv,max_dp);return bad?1:0;
}
