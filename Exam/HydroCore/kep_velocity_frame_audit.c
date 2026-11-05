/* Local reconstruction audit on ACTUAL native power meshes. Smooth data
 * are centroid samples on shrinking local patches, not evolved PDE norms. */
#define main archived_native_main
#include "fv_native_evolution.c"
#undef main
#include "../gfs_fv_native.c"
typedef struct {double G[3][3],raw[3][3],value[3],mean[3],error,bounds;int unchanged;} Frame;
static void matrix(int d,int rotated,double R[3][3]){
    memset(R,0,9*sizeof(double));for(int a=0;a<3;a++)R[a][a]=1;
    if(!rotated||d==1)return;
    double c=cos(.37),s=sin(.37);R[0][0]=c;R[0][1]=-s;R[1][0]=s;R[1][1]=c;
    if(d==3){double D[3][3]={{cos(.53),0,sin(.53)},{0,1,0},{-sin(.53),0,cos(.53)}},out[3][3]={{0}};
        for(int a=0;a<3;a++)for(int b=0;b<3;b++)for(int k=0;k<3;k++)out[a][b]+=R[a][k]*D[k][b];memcpy(R,out,sizeof(out));}
}
static void smooth(int d,const double x[3],double u[3]){
    memset(u,0,3*sizeof(double));
    u[0]=.2*x[0]*x[0]+.1*x[0]*x[0]*x[0];
    if(d>=2){u[0]+=.15*x[1]*x[1]+.07*x[0]*x[1];u[1]=x[0]+.3*x[1]+.2*x[0]*x[0];}
    if(d==3){u[0]+=.11*x[2]*x[2];u[1]+=.1*x[2]-.05*x[1]*x[2];u[2]=.2*x[0]-.4*x[1]+.5*x[2]+.12*x[2]*x[2];}
}
static int evaluate(int d,int kind,double h,int rotated,int boosted,Frame *out){
    GfsFvMesh m={0};m.dim=d;m.gamma=1.4;m.n=1;
    for(int a=0;a<3;a++){m.L[a]=1;if(a<d){m.n*=7;m.lo[a]=-6*h;m.L[a]=12*h;m.bc[a]=GFS_FV_OUTFLOW;}}
    m.cell=calloc(m.n,sizeof(*m.cell));double R[3][3];matrix(d,rotated,R);
    int center=0,stride=1;for(int a=0;a<d;a++){center+=3*stride;stride*=7;}
    for(int i=0;i<m.n;i++){
        int index=i;double x[3]={0};for(int a=0;a<d;a++){x[a]=(index%7-3)*h;index/=7;}
        for(int a=0;a<d;a++)for(int b=0;b<d;b++)m.cell[i].q[a]+=R[a][b]*x[b];m.cell[i].id=i;
    }
    if(gfs_fv_geometry(&m))return 2;
    double A[3][3]={{.2,.4,.1},{-.3,-.2,.05},{.2,-.1,0}},boost[3]={.7,-.3,.2},rank[3]={.4,.7,.5};
    for(int i=0;i<m.n;i++){
        GfsFvCell *c=m.cell+i;double base_x[3]={0},u[3]={0},world[3]={0};
        for(int b=0;b<d;b++)for(int a=0;a<d;a++)base_x[b]+=R[a][b]*(c->q[a]+m.geometry[i].centroid[a]);
        if(kind==0)for(int a=0;a<d;a++){u[a]=.1*(a+1);for(int b=0;b<d;b++)u[a]+=A[a][b]*base_x[b];}
        if(kind==1){double value=i==center+1?1:i==center-1?-.1:0;for(int a=0;a<d;a++)u[a]=rank[a]*value;}
        if(kind==2)smooth(d,base_x,u);
        if(boosted)for(int a=0;a<d;a++)u[a]+=boost[a];
        for(int a=0;a<d;a++)for(int b=0;b<d;b++)world[a]+=R[a][b]*u[b];
        double V=m.geometry[i].V;c->Q[0]=V;c->Q[4]=V/(m.gamma-1);c->m0=V;c->mu=V/pow(V,2./d);
        for(int a=0;a<d;a++){c->Q[a+1]=V*world[a];c->Q[4]+=.5*V*world[a]*world[a];}
    }
    GfsFvCell *copy=malloc(m.n*sizeof(*copy));memcpy(copy,m.cell,m.n*sizeof(*copy));
    GfsFvWorkspace *w=gfs_fv_workspace_create();allocate(w,m.n);if(reconstruct(w,&m))return 3;
    memset(out,0,sizeof(*out));memcpy(out->raw,w->raw_velocity_gradient[center],sizeof(out->raw));
    for(int a=0;a<d;a++){out->mean[a]=w->primitive[center].u[a];for(int b=0;b<d;b++)out->G[a][b]=w->gradient[center][a+1][b];}
    double probe[3]={.2*h,.15*h,.1*h},world[3]={0},expected[3]={0};
    for(int a=0;a<d;a++)for(int b=0;b<d;b++)world[a]+=R[a][b]*probe[b];
    int owner,hops;double rel[3];if(gfs_fv_point_owner(&m,center,world,&owner,rel,&hops)||owner!=center)return 3;
    for(int a=0;a<d;a++)out->value[a]=physical_velocity(w,&m,owner,a,rel,world);
    if(kind==2){double u[3];smooth(d,probe,u);for(int a=0;a<d;a++)for(int b=0;b<d;b++)expected[a]+=R[a][b]*(u[b]+(boosted?boost[b]:0));for(int a=0;a<d;a++)out->error=hypot(out->error,out->value[a]-expected[a]);}
    for(int a=0;a<d;a++){
        double lo=0,hi=0;
        for(int k=0;k<m.geometry[center].faces;k++){GfsFvPrimitive p;double r[3];neighbor(&m,w,center,m.geometry[center].face+k,&p,r);double dv=p.u[a]-w->primitive[center].u[a];lo=fmin(lo,dv);hi=fmax(hi,dv);}
        for(int k=0;k<m.geometry[center].faces;k++){GfsFvPrimitive p;double r[3];neighbor(&m,w,center,m.geometry[center].face+k,&p,r);double prediction=.5*dot3(w->gradient[center][a+1],r,d);out->bounds=fmax(out->bounds,fmax(0,fmax(prediction-hi,lo-prediction))/(1+fabs(lo)+fabs(hi)));}
    }
    out->unchanged=!memcmp(copy,m.cell,m.n*sizeof(*copy));
    free(copy);gfs_fv_workspace_free(w);gfs_fv_geometry_free(&m);free(m.cell);return 0;
}
static double covariance(int d,const Frame *a,const Frame *b,int raw){
    double R[3][3],maximum=0;matrix(d,1,R);
    for(int i=0;i<d;i++)for(int j=0;j<d;j++){
        double expected=0;for(int p=0;p<d;p++)for(int q=0;q<d;q++)expected+=R[i][p]*(raw?a->raw[p][q]:a->G[p][q])*R[j][q];
        maximum=fmax(maximum,fabs(expected-(raw?b->raw[i][j]:b->G[i][j])));
    }return maximum;
}
static double value_covariance(int d,const Frame *a,const Frame *b){
    double R[3][3],maximum=0;matrix(d,1,R);
    for(int i=0;i<d;i++){double expected=0;for(int j=0;j<d;j++)expected+=R[i][j]*a->value[j];maximum=fmax(maximum,fabs(expected-b->value[i]));}return maximum;
}
static int contact(int d,int moving){
    GfsFvMesh m={0};m.dim=d;m.gamma=1.4;m.n=1;for(int a=0;a<3;a++){m.L[a]=1;if(a<d)m.n*=8;}
    m.cell=calloc(m.n,sizeof(*m.cell));double projection;initialize(&m,8,moving?MOVING_CONTACT:CONTACT,&projection);if(gfs_fv_geometry(&m))return 2;
    GfsFvWorkspace *w=gfs_fv_workspace_create();GfsFvRate *r=calloc(m.n,sizeof(*r));GfsFvRhsStats s;
    int status=gfs_fv_native_rhs(w,&m,r,&s);double equilibrium=0;
    if(!status)for(int i=0;i<m.n;i++){
        for(int c=0;c<5;c++)equilibrium=fmax(equilibrium,fabs(r[i].dQ[c]));
        for(int a=0;a<d;a++)equilibrium=fmax(equilibrium,fabs(r[i].qdot[a]-(moving&&a==0?1:0)));
    }
    int pass=!status&&equilibrium<1e-11&&fabs(s.entropy_balance_defect)<1e-10;
    printf("{\"case\":\"contact_RHS\",\"dim\":%d,\"moving\":%s,\"maximum_equilibrium_defect\":%.17g,\"rhs_status\":%d,\"pass\":%s}\n",d,moving?"true":"false",equilibrium,status,pass?"true":"false");
    free(r);gfs_fv_workspace_free(w);gfs_fv_geometry_free(&m);free(m.cell);return pass?0:1;
}
int main(void){
    int failed=0;
    for(int d=1;d<=3;d++){
        for(int kind=0;kind<2;kind++){
            Frame a,b,c;int status=evaluate(d,kind,.2,0,0,&a)||evaluate(d,kind,.2,1,0,&b)||evaluate(d,kind,.2,0,1,&c);if(status)return 2;
            double raw=covariance(d,&a,&b,1),limited=covariance(d,&a,&b,0),value=value_covariance(d,&a,&b),boost_error=0;
            const double boost[3]={.7,-.3,.2};for(int i=0;i<d;i++){
                boost_error=fmax(boost_error,fabs(c.value[i]-a.value[i]-boost[i]));
                for(int j=0;j<d;j++)boost_error=fmax(boost_error,fabs(c.G[i][j]-a.G[i][j]));
            }
            int pass=raw<1e-12&&limited<1e-12&&value<1e-12&&boost_error<1e-12&&a.bounds<1e-12&&b.bounds<1e-12&&c.bounds<1e-12&&a.unchanged&&b.unchanged&&c.unchanged;
            printf("{\"case\":\"%s\",\"dim\":%d,\"raw_rotation_error\":%.17g,\"limited_rotation_error\":%.17g,\"query_rotation_error\":%.17g,\"boost_error\":%.17g,\"scalar_midpoint_bounds_error\":%.17g,\"pass\":%s}\n",kind?"rank_one_alias":"affine",d,raw,limited,value,boost_error,fmax(a.bounds,fmax(b.bounds,c.bounds)),pass?"true":"false");failed|=!pass;
        }
        Frame a[2],b[2];double rotation[2];
        for(int level=0;level<2;level++){
            double h=level?.1:.2;if(evaluate(d,2,h,0,0,a+level)||evaluate(d,2,h,1,0,b+level))return 2;
            rotation[level]=value_covariance(d,a+level,b+level);
            if(covariance(d,a+level,b+level,1)>1e-12||a[level].bounds>1e-12||b[level].bounds>1e-12)failed=1;
        }
        int pass=a[1].error<a[0].error&&b[1].error<b[0].error&&(rotation[0]<1e-12?rotation[1]<1e-12:rotation[1]<rotation[0]);
        printf("{\"case\":\"smooth_local_patch_refinement\",\"dim\":%d,\"centroid_samples_not_PDE\":true,\"h\":[0.2,0.1],\"unrotated_error\":[%.17g,%.17g],\"rotated_error\":[%.17g,%.17g],\"nonlinear_rotation_difference\":[%.17g,%.17g],\"pass\":%s}\n",d,a[0].error,a[1].error,b[0].error,b[1].error,rotation[0],rotation[1],pass?"true":"false");failed|=!pass;
        for(int moving=0;moving<2;moving++)failed|=contact(d,moving);
    }
    return failed;
}
