/* Frozen Gresho diagnostics only: actual geometry volume work, affine
 * transport moments, and a linear decomposition of pressure-force error.
 * Exact pressure is sampled at generators, not presented as a cell integral.
 * No state, force, timestep, or constitutive coefficient is changed. */
static void gfs_gresho_reference(double x,double y,double *P,double v[3],double G[2][2])
{
    double dx=x-.5,dy=y-.5,r=hypot(dx,dy),a=0,ap=0;
    if(r<=.2){*P=5+12.5*r*r;a=5;}
    else if(r<=.4){*P=9+12.5*r*r-20*r+4*log(5*r);a=2/r-5;ap=-2/(r*r);}
    else *P=3+4*log(2.);
    v[0]=-a*dy;v[1]=a*dx;v[2]=0;
    G[0][0]=r>0?-ap*dx*dy/r:0;G[0][1]=-a-(r>0?ap*dy*dy/r:0);
    G[1][0]=a+(r>0?ap*dx*dx/r:0);G[1][1]=-G[0][0];
}
static int gfs_unified_transport_audit(SimParameters *simpar,GPUContext *ctx,
        const int *alias,double gamma,const double *G)
{
    const char *flag=getenv("GFS_UNIFIED_TRANSPORT_AUDIT");
    if(!flag||strcmp(flag,"1"))return 0;
    const char *ic=getenv("EUNHA_IC");
    if(!gfs_unified_periodic_mode()||!ic||strcmp(ic,"gresho"))return 1;
    int n=VORO_NP(simpar),nt=ctx->h_parts.n_total;
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(simpar);
    FaceCSR *f=&ctx->h_faces;ParticleSoA *a=&ctx->h_parts;
    FILE *out=fopen("unified_transport.csv","wx");if(!out)return 1;
    fprintf(out,"id,x,y,m,V,P,P_exact_node,div_geometric,div_fit,div_geometric_affine_fit,div_geometric_exact_velocity,div_fit_exact_velocity,div_geometric_exact_gradient,Kxx,Kxy,Kyx,Kyy,pressure_moment_error,pressure_gradient_error_x,pressure_gradient_error_y,pressure_noise_force_x,pressure_noise_force_y,cpu_gpu_force_error,thermal_rhs_error\n");
    for(int i=0;i<n;i++){
        double V=p[i].weight_thermal.V,Pi=(gamma-1)*p[i].weight_thermal.U/V;
        double P0,v0[3],G0[2][2],vi[3]={p[i].vx,p[i].vy,0};
        gfs_gresho_reference(p[i].x,p[i].y,&P0,v0,G0);
        double K[2][2]={{0}},M[2][2]={{0}},F[3]={0},projected[3]={0},rate=0,exact_rate=0;
        GfsVelocityFit fit;gfs_velocity_fit_init(&fit,2);
        for(int k=f->face_offset[i];k<f->face_offset[i+1];k++){
            int raw=f->neighbor_idx[k];if(raw<0||raw>=nt||f->is_ghost[k])return 1;
            int j=alias[raw];if(j<0||j>=n)return 1;
            double r[3]={a->x[raw]-p[i].x,a->y[raw]-p[i].y,0};
            double c[3]={.5*(f->c1x[k]+f->c2x[k]),.5*(f->c1y[k]+f->c2y[k]),0};
            double area=hypot(f->c2y[k]-f->c1y[k],f->c1x[k]-f->c2x[k]);GfsVariationalFace b;
            if(!gfs_variational_face(2,r,c,area,&b))return 1;
            double Pj=(gamma-1)*p[j].weight_thermal.U/p[j].weight_thermal.V;
            double vj[3]={p[j].vx,p[j].vy,0},Pj0,vj0[3],Gj0[2][2],dv[3];
            gfs_gresho_reference(p[j].x,p[j].y,&Pj0,vj0,Gj0);
            rate+=gfs_variational_sweep(&b,vi,vj);exact_rate+=gfs_variational_sweep(&b,v0,vj0);
            for(int d=0;d<3;d++)dv[d]=vj0[d]-v0[d];
            if(!gfs_velocity_fit_add_face(&fit,r,dv,area))return 1;
            for(int d=0;d<2;d++){
                F[d]+=(Pi-Pj)*b.bi[d];projected[d]+=(P0-Pj0)*b.bi[d];
                for(int e=0;e<2;e++){K[d][e]+=b.bj[d]*r[e]/V;M[d][e]+=b.bi[d]*r[e]/V;}
            }
        }
        double fit0[3][3];if(!gfs_velocity_fit_solve(&fit,fit0,NULL))return 1;
        double moment_error=0,affine_fit=0,affine_exact=0;
        for(int d=0;d<2;d++)for(int e=0;e<2;e++){
            moment_error=fmax(moment_error,fabs(M[d][e]-(d==e)));
            affine_fit+=K[d][e]*G[4*i+2*d+e];affine_exact+=K[d][e]*G0[d][e];
        }
        double dx=p[i].x-.5,dy=p[i].y-.5,r2=dx*dx+dy*dy,speed2=v0[0]*v0[0]+v0[1]*v0[1];
        double grad[2]={r2>0?speed2*dx/r2:0,r2>0?speed2*dy/r2:0};
        double gpu_error=hypot(F[0]-p[i].mass*p[i].ax,F[1]-p[i].mass*p[i].ay);
        double rhs_error=p[i].weight_thermal.Udot+Pi*rate;
        fprintf(out,"%llu,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
            (unsigned long long)PINDX(p+i),p[i].x,p[i].y,p[i].mass,V,Pi,P0,rate/V,G[4*i]+G[4*i+3],affine_fit,
            exact_rate/V,fit0[0][0]+fit0[1][1],affine_exact,K[0][0],K[0][1],K[1][0],K[1][1],moment_error,
            projected[0]/V+grad[0],projected[1]/V+grad[1],(F[0]-projected[0])/V,(F[1]-projected[1])/V,gpu_error,rhs_error);
    }
    if(ferror(out)||fclose(out))return 1;
    fprintf(stderr,"[UNIFIED_TRANSPORT_AUDIT] native_geometry=1 exact_pressure_is_nodal_sample=1 no_state_update=1\n");return 0;
}
