#ifndef FV_ORIGINAL_CASES_H
#define FV_ORIGINAL_CASES_H
/* Original Hydro1DExam problem definitions, evaluated with the SAME native
 * dimensional FV operator. Open problems use an enlarged moving material
 * patch; assessment is restricted to the original physical interval.
 * Padding convergence is required because the weight solve is screened,
 * not strictly causally local. No generator steering or cell merging. */
typedef struct {
    double lo,hi,cut[2],gamma,end,rho[3],u[3],p[3];
    int regions,lower,upper,sinusoidal;
} FvOriginalCase;
static int original_case(int problem,FvOriginalCase *s){
    *s=(FvOriginalCase){0,1,{.5,1},1.4,0,{1,1,1},{0,0,0},{1,1,1},2,GFS_FV_MOVING_TRANSMIT,GFS_FV_MOVING_TRANSMIT,0};
    switch(problem){
    case BLAST:s->regions=3;s->cut[0]=.1;s->cut[1]=.9;s->p[0]=1000;s->p[1]=.01;s->p[2]=100;s->lower=s->upper=GFS_FV_REFLECT;s->end=.038;break;
    case COLLISION:s->cut[0]=.4;s->rho[0]=5.99924;s->rho[1]=5.99242;s->u[0]=19.5975;s->u[1]=-6.19633;s->p[0]=460.894;s->p[1]=46.095;s->end=.035;break;
    case SHU_OSHER:s->lo=-5;s->hi=5;s->cut[0]=-4;s->rho[0]=3.857143;s->u[0]=2.629369;s->p[0]=10.33333;s->end=1.8;s->sinusoidal=1;break;
    case NOH:s->gamma=5./3;s->regions=1;s->u[0]=-1;s->p[0]=1e-6;s->lower=GFS_FV_REFLECT;s->end=.6;break;
    case LAX:s->rho[0]=.445;s->rho[1]=.5;s->u[0]=.698;s->p[0]=3.528;s->p[1]=.571;s->end=.13;break;
    case DOUBLE_FAN:s->u[0]=-2;s->u[1]=2;s->p[0]=s->p[1]=.4;s->end=.15;break;
    case CONTACT_OPEN:s->rho[1]=.125;s->end=.25;break;
    case MOVING_CONTACT_OPEN:s->rho[1]=.125;s->u[0]=s->u[1]=1;s->end=.25;break;
    default:return 0;
    }return 1;
}
static void original_point(const FvOriginalCase *s,double x,GfsFvPrimitive *p){
    int region=0;while(region<s->regions-1&&x>=s->cut[region])region++;
    *p=(GfsFvPrimitive){s->rho[region],s->p[region],{s->u[region],0,0}};
    if(s->sinusoidal&&region==1)p->rho=1+.2*sin(5*M_PI*x);
}
static int original_mesh(GfsFvMesh *m,int side,const FvOriginalCase *s,int multiplier,int padding[2]){
    if(m->dim!=1||multiplier<1||multiplier>2)return 1;
    double h=(s->hi-s->lo)/side,signal=0;
    for(int k=0;k<s->regions;k++){
        double minrho=s->sinusoidal&&k==1?.8:s->rho[k];
        signal=fmax(signal,fabs(s->u[k])+sqrt(s->gamma*s->p[k]/minrho));
    }
    int guard=multiplier*((int)ceil(signal*s->end/h)+4);
    padding[0]=s->lower==GFS_FV_MOVING_TRANSMIT?guard:0;
    padding[1]=s->upper==GFS_FV_MOVING_TRANSMIT?guard:0;
    m->n=side+padding[0]+padding[1];m->gamma=s->gamma;
    m->lo[0]=s->lo-padding[0]*h;m->L[0]=(m->n)*h;
    m->bc[0]=gfs_fv_boundary_pair(s->lower,s->upper);return 0;
}
static int initialize_original(GfsFvMesh *m,int side,const FvOriginalCase *s,const int padding[2],double *projection){
    double h=(s->hi-s->lo)/side;*projection=0;
    for(int i=0;i<m->n;i++){m->cell[i].q[0]=s->lo+(i-padding[0]+.5)*h;m->cell[i].id=i;}
    if(gfs_fv_geometry(m))return 1;
    for(int i=0;i<m->n;i++){
        GfsFvCell *c=m->cell+i;GfsFvGeometry *g=m->geometry+i;
        double left=c->q[0]+g->centroid[0]-.5*g->V,right=left+g->V,thermal=0,kinetic=0;
        for(int k=0;k<s->regions;k++){
            double a=fmax(left,k?s->cut[k-1]:left),b=fmin(right,k<s->regions-1?s->cut[k]:right);
            if(!(b>a))continue;
            double mass=s->rho[k]*(b-a);
            if(s->sinusoidal&&k==1)mass=(b-a)+.2*(cos(5*M_PI*a)-cos(5*M_PI*b))/(5*M_PI);
            c->Q[0]+=mass;c->Q[1]+=mass*s->u[k];thermal+=(b-a)*s->p[k]/(s->gamma-1);kinetic+=.5*mass*s->u[k]*s->u[k];
        }
        double mean=.5*c->Q[1]*c->Q[1]/c->Q[0];c->Q[4]=thermal+mean;c->m0=c->Q[0];c->mu=c->m0/(g->V*g->V);*projection+=kinetic-mean;
    }return 0;
}
static int original_receipt(const FvOriginalCase *s,int side,int n,const int padding[2],int multiplier){
    FILE *f=fopen("problem_spec.json","wx");if(!f)return 1;
    fprintf(f,"{\"source\":\"Hydro1DExam/laguerre_sod.c::setup_problem\",\"physical_domain\":[%.17g,%.17g],\"physical_resolution\":%d,\"total_cells\":%d,\"padding_cells\":[%d,%d],\"guard_multiplier\":%d,\"gamma\":%.17g,\"original_endpoint\":%.17g,\"cuts\":[%.17g,%.17g],\"regions\":%d,\"rho\":[%.17g,%.17g,%.17g],\"velocity\":[%.17g,%.17g,%.17g],\"pressure\":[%.17g,%.17g,%.17g],\"lower_bc\":%d,\"upper_bc\":%d,\"right_density_wavenumber\":%.17g,\"assessment_scope\":\"original physical interval only; guard padding is validated separately, not assumed causally sufficient\"}\n",s->lo,s->hi,side,n,padding[0],padding[1],multiplier,s->gamma,s->end,s->cut[0],s->cut[1],s->regions,s->rho[0],s->rho[1],s->rho[2],s->u[0],s->u[1],s->u[2],s->p[0],s->p[1],s->p[2],s->lower,s->upper,s->sinusoidal?5*M_PI:0);
    return fclose(f)!=0;
}
#endif
