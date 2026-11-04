/* Same local geometric identities in all dimensions. The second layout is a
 * refinement by exact periodic repetition of the first relative site pattern.
 * Affine velocity differences are local patches, not a periodic global shear.
 * This exposes a consistency defect; a nonzero shear result is not hidden by
 * declaring the numerical bookkeeping checks to be a fluid-accuracy pass. */
#define Voro3D_FindVC Voro3D_FindVCStrict
#define NC 256
#define main archived_geometry_audit_main
#include "geometry_audit.c"
#undef main
#include "../gfs_velocity_fit.h"
int main(void){
    int fail=0;
    for(int dim=1;dim<=3;dim++)for(int refinement=1;refinement<=2;refinement++){
        int primitive=dim==1?16:3,side=primitive*refinement,n=1;
        for(int a=0;a<dim;a++)n*=side;if(n>NC)return 2;
        for(int i=0;i<n;i++){
            int q=i,tile=0,stride=1,cell[3]={0};
            for(int a=0;a<dim;a++){cell[a]=q%side;q/=side;tile+=(cell[a]%primitive)*stride;stride*=primitive;}
            for(int a=0;a<3;a++){pos[i][a]=a<dim?(cell[a]+.5+.07*sin(1.3*(tile+1)*(a+2)))/side:0;vel[i][a]=0;}
        }
        double V[NC];tess(dim,n,0,V,1);
        double K[NC][3][3]={{{0}}},M[NC][3][3]={{{0}}};
        for(int k=0;k<nf;k++){
            Face *f=faces+k;GfsVariationalFace b;if(!gfs_variational_face(dim,f->r,f->c,f->area,&b))return 2;
            for(int a=0;a<dim;a++)for(int c=0;c<dim;c++){
                K[f->i][a][c]+=b.bj[a]*f->r[c]/V[f->i];M[f->i][a][c]+=b.bi[a]*f->r[c]/V[f->i];
            }
        }
        double fitted[NC][3][3]={{{0}}};
        for(int i=0;i<n;i++){
            GfsVelocityFit fit;gfs_velocity_fit_init(&fit,dim);
            for(int k=0;k<nf;k++)if(faces[k].i==i){
                double dv[3]={dim>1?faces[k].r[1]:0,0,0};
                if(!gfs_velocity_fit_add_face(&fit,faces[k].r,dv,faces[k].area))return 2;
            }
            if(!gfs_velocity_fit_solve(&fit,fitted[i],NULL))return 2;
        }
        double max_q=0,heat=0,max_fit_div=0;
        for(int k=0;k<nf;k++){
            Face *f=faces+k;GfsVariationalFace b;if(!gfs_variational_face(dim,f->r,f->c,f->area,&b))return 2;
            double vi[3]={0},vj[3]={dim>1?f->r[1]:0,0,0},jump[3];
            for(int a=0;a<3;a++){
                jump[a]=vj[a];
                for(int c=0;c<dim;c++)jump[a]-=.5*(fitted[f->i][a][c]+fitted[f->j][a][c])*f->r[c];
            }
            double di=0,dj=0;for(int a=0;a<dim;a++){di+=fitted[f->i][a][a];dj+=fitted[f->j][a][a];}
            max_fit_div=fmax(max_fit_div,fabs(di));
            double q=gfs_variational_pair_pressure(&b,vi,vj,jump,1,sqrt(1.4),1,sqrt(1.4),
                gfs_variational_compression_fraction(di,dim>1?-1:0,dj,dim>1?-1:0),pow(V[f->i],1./dim),pow(V[f->j],1./dim));
            max_q=fmax(max_q,q);double approach=-dot(vj,b.normal,dim);
            if(approach>0)heat+=.5*q*f->area*approach;
        }
        double pressure_error=0,trace_error=0,rotation_error=0,shear_max=0,shear_L1=0,signed_shear=0;
        for(int i=0;i<n;i++){
            double trace=0;
            for(int a=0;a<dim;a++)for(int c=0;c<dim;c++){
                pressure_error=fmax(pressure_error,fabs(M[i][a][c]-(a==c)));
                rotation_error=fmax(rotation_error,fabs(K[i][a][c]-K[i][c][a]));
                if(a==c)trace+=K[i][a][c];
            }
            trace_error=fmax(trace_error,fabs(trace-dim));
            double shear=dim>1?K[i][0][1]:K[i][0][0]-1;
            shear_max=fmax(shear_max,fabs(shear));shear_L1+=V[i]*fabs(shear);signed_shear+=V[i]*shear;
        }
        int bookkeeping=pressure_error<1e-9&&trace_error<1e-9&&rotation_error<1e-9;
        fail|=!bookkeeping;
        printf("{\"dimension\":%d,\"cells\":%d,\"pattern_refinement\":%d,\"max_linear_pressure_patch_error\":%.17g,\"max_dilation_error\":%.17g,\"max_rotation_error\":%.17g,\"incompressible_shear_geometric_divergence_L1\":%.17g,\"incompressible_shear_geometric_divergence_max\":%.17g,\"global_signed_shear_volume_rate\":%.17g,\"max_fitted_physical_divergence\":%.17g,\"max_direct_viscous_pressure\":%.17g,\"direct_affine_viscous_heat_rate\":%.17g,\"geometry_identity_checks_pass\":%s,\"pointwise_affine_material_transport\":%s}\n",
            dim,n,refinement,pressure_error,trace_error,rotation_error,shear_L1,shear_max,
            signed_shear,max_fit_div,max_q,heat,bookkeeping?"true":"false",shear_max<1e-9?"true":"false");
    }
    return fail;
}
