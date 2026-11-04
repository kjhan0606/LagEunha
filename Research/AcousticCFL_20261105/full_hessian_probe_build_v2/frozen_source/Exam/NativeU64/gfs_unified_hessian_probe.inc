/* Frozen directional derivative of the complete isentropic pressure force.
 * Uses the original native geometry/force wrappers at both displaced states.
 * This diagnostic exits without accepting an evolution step. It does not
 * supply a timestep bound or modify the production integrator. */
typedef struct {
    unsigned long long id;
    treevorostressrk4particletype particle;
    double z[2];
} GfsUnifiedHessianReference;
static int gfs_unified_hessian_compare(const void *a,const void *b){
    unsigned long long x=((const GfsUnifiedHessianReference*)a)->id;
    unsigned long long y=((const GfsUnifiedHessianReference*)b)->id;
    return x<y?-1:x>y;
}
static int gfs_unified_hessian_probe(SimParameters *simpar,double gamma,
        double xmin,double ymin,double xmax,double ymax,double OoA,double courant,
        void (*padding)(SimParameters *,postype),
        Voro2D_point *(*neighbors)(SimParameters *,int,int,int *),
        treevorork4particletype *(*cells)(SimParameters *,int,int,int *),
        void (*linked)(SimParameters *,postype,postype,postype,postype,postype,
                       void (*)(SimParameters *,postype)),
        void (*post)(SimParameters *))
{
    if(!gfs_unified_periodic_mode()||NID(simpar)!=1)return 1;
    int n=VORO_NP(simpar);
    treevorostressrk4particletype *p=(void*)VORORK4_TBP(simpar);
    GfsUnifiedHessianReference *ref=calloc(n,sizeof(*ref));
    double *side[2]={calloc(2*n,sizeof(double)),calloc(2*n,sizeof(double))};
    double *previous=calloc(2*n,sizeof(double));
    if(!ref||!side[0]||!side[1]||!previous)return 1;
    FILE *in=fopen("unified_acoustic_direction.csv","r");char header[128];
    if(!in||!fgets(header,sizeof(header),in))return 1;
    long double norm=0;double max_dx=0,hmin=HUGE_VAL,dt=HUGE_VAL;
    for(int i=0;i<n;i++){
        unsigned long long id;double zx,zy;
        if(fscanf(in,"%llu,%lf,%lf",&id,&zx,&zy)!=3||id!=(unsigned long long)PINDX(p+i))return 1;
        ref[i].id=id;ref[i].particle=p[i];ref[i].z[0]=zx;ref[i].z[1]=zy;
        norm+=(long double)zx*zx+(long double)zy*zy;
        max_dx=fmax(max_dx,hypot(zx,zy)/sqrt(p[i].mass));
        hmin=fmin(hmin,sqrt(p[i].weight_thermal.V));dt=fmin(dt,p[i].dt);
    }
    if(fscanf(in," %127s",header)!=EOF||fclose(in)||!(norm>0)||!(max_dx>0))return 1;
    GPUContext *ctx=gpu_get_context();FaceCSR *face=&ctx->h_faces;
    ParticleSoA *a=&ctx->h_parts;
    const int *alias=gpu_unified_owner_aliases(ctx,n,a->n_total,
        SIMBOX(simpar).x.max-SIMBOX(simpar).x.min,
        SIMBOX(simpar).y.max-SIMBOX(simpar).y.min);
    long double positive=0;
    for(int i=0;i<n;i++){
        double dV=0;
        for(int k=face->face_offset[i];k<face->face_offset[i+1];k++){
            int raw=face->neighbor_idx[k];if(raw<0||face->is_ghost[k])return 1;
            int j=alias[raw];GfsVariationalFace b;
            double r[3]={a->x[raw]-p[i].x,a->y[raw]-p[i].y,0};
            double c[3]={.5*(face->c1x[k]+face->c2x[k]),.5*(face->c1y[k]+face->c2y[k]),0};
            double area=hypot(face->c2y[k]-face->c1y[k],face->c1x[k]-face->c2x[k]);
            if(!gfs_variational_face(2,r,c,area,&b))return 1;
            for(int d=0;d<2;d++)dV+=b.bj[d]*(ref[j].z[d]/sqrt(p[j].mass)-ref[i].z[d]/sqrt(p[i].mass));
        }
        double V=p[i].weight_thermal.V,P=(gamma-1)*p[i].weight_thermal.U/V;
        positive+=(long double)gamma*P/V*dV*dV;
    }
    positive/=norm;
    qsort(ref,n,sizeof(*ref),gfs_unified_hessian_compare);
    for(int i=1;i<n;i++)if(ref[i-1].id==ref[i].id)return 1;
    const double fractions[]={0,.001,.0003,.0001,.00003,.00001,.000003};
    long double E0=0;
    for(int test=0;test<(int)(sizeof(fractions)/sizeof(*fractions));test++){
        double epsilon=fractions[test]*hmin/max_dx;long double energy[2]={0,0};
        for(int sign=0;sign<(test?2:1);sign++){
            p=(void*)VORORK4_TBP(simpar);
            for(int i=0;i<n;i++){
                GfsUnifiedHessianReference key={.id=(unsigned long long)PINDX(p+i)};
                const GfsUnifiedHessianReference *base=bsearch(&key,ref,n,sizeof(*ref),gfs_unified_hessian_compare);
                if(!base)return 1;p[i]=base->particle;
                double step=(sign?-epsilon:epsilon)/sqrt(p[i].mass);
                p[i].x+=step*base->z[0];p[i].y+=step*base->z[1];
            }
            post(simpar);
            updateDenW2Pressure2DBlend(simpar,xmin,ymin,xmax,ymax,gamma,
                padding,neighbors,cells,linked,1);
            if(VORO_NP(simpar)!=n)return 1;p=(void*)VORORK4_TBP(simpar);
            for(int i=0;i<n;i++){
                GfsUnifiedHessianReference key={.id=(unsigned long long)PINDX(p+i)};
                const GfsUnifiedHessianReference *base=bsearch(&key,ref,n,sizeof(*ref),gfs_unified_hessian_compare);
                double U;
                if(!base||!gfs_variational_adiabatic(base->particle.weight_thermal.U,
                        base->particle.weight_thermal.V,p[i].weight_thermal.V,gamma,&U))return 1;
                p[i].weight_thermal.U=U;p[i].ie=(float)U;
                energy[sign]+=(long double)base->particle.weight_thermal.U*
                    powl((long double)base->particle.weight_thermal.V/p[i].weight_thermal.V,gamma-1);
            }
            (void)getAccVoro2DBlend_GPU(simpar,xmin,ymin,xmax,ymax,OoA,courant,gamma,
                padding,neighbors,cells,linked);
            /* Geometry packing precedes the entropy-preserving U update.
             * Refresh that changed state, as the production final kick does. */
            (void)gpu_unified_pressure_refresh(simpar);
            p=(void*)VORORK4_TBP(simpar);
            for(int i=0;i<n;i++){
                GfsUnifiedHessianReference key={.id=(unsigned long long)PINDX(p+i)};
                const GfsUnifiedHessianReference *base=bsearch(&key,ref,n,sizeof(*ref),gfs_unified_hessian_compare);
                if(!base)return 1;int index=(int)(base-ref);
                side[sign][2*index]=p[i].ax;side[sign][2*index+1]=p[i].ay;
            }
        }
        if(!test){E0=energy[0];continue;}
        long double rayleigh=0,change=0,magnitude=0;
        for(int i=0;i<n;i++)for(int d=0;d<2;d++){
            int k=2*i+d;double value=-sqrt(ref[i].particle.mass)*(side[0][k]-side[1][k])/(2*epsilon);
            rayleigh+=(long double)ref[i].z[d]*value;
            change+=(long double)(value-previous[k])*(value-previous[k]);
            magnitude+=(long double)value*value;previous[k]=value;
        }
        rayleigh/=norm;
        long double energy_curvature=(energy[0]+energy[1]-2*E0)/(epsilon*epsilon*norm);
        fprintf(stderr,"[UNIFIED_FULL_HESSIAN] max_displacement_over_hmin=%.9g epsilon=%.17g full_directional_lambda=%.17Lg positive_directional_lambda=%.17Lg geometric_directional_lambda=%.17Lg dt_omega_direction=%.17g vector_relative_change=%.9g energy_curvature=%.17Lg no_accepted_steps=1\n",
            fractions[test],epsilon,rayleigh,positive,rayleigh-positive,
            rayleigh>0?dt*sqrt((double)rayleigh):NAN,
            test>1?sqrt((double)(change/magnitude)):NAN,energy_curvature);
    }
    free(ref);free(side[0]);free(side[1]);free(previous);return 0;
}
