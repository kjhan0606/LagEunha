/* Actual project geometry in 1D/2D/3D, dense acoustic matrices for independent
 * eigenspectrum checks, and rotation invariance of the block-norm bound. */
#define main archived_geometry_audit_main
#include "geometry_audit.c"
#undef main
#include "../gfs_acoustic_cfl.h"
static void rotate(double x[3],int dim){
    if(dim<2)return;
    double u=cos(.63)*x[0]-sin(.63)*x[1],v=sin(.63)*x[0]+cos(.63)*x[1];x[0]=u;x[1]=v;
    if(dim==3){u=cos(.41)*x[1]-sin(.41)*x[2];v=sin(.41)*x[1]+cos(.41)*x[2];x[1]=u;x[2]=v;}
}
int main(void){
    for(int dim=1;dim<=3;dim++){
        int side=dim==1?16:3,n=1;for(int a=0;a<dim;a++)n*=side;
        double V[NC],m[NC],U[NC],gamma=1.4,rootD[NC],B[NC][NC*3]={{0}};
        for(int i=0;i<n;i++){
            int q=i;
            for(int a=0;a<3;a++){
                int cell=q%side;q/=side;pos[i][a]=a<dim?(cell+.5+.07*sin(1.3*(i+1)*(a+2)))/side:0;vel[i][a]=0;
            }
        }
        tess(dim,n,0,V,1);
        for(int i=0;i<n;i++){
            m[i]=V[i]*(1+.2*sin(.71*i));U[i]=V[i]*(1+.3*cos(.9*i))/(gamma-1);
            rootD[i]=sqrt(gamma*(gamma-1)*U[i])/V[i];
        }
        double bound[2];
        for(int rotated=0;rotated<2;rotated++){
            GfsAcousticBound state[NC];
            for(int i=0;i<n;i++)if(!gfs_acoustic_bound_init(state+i,dim,m[i],V[i],U[i],gamma))return 2;
            for(int k=0;k<nf;k++){
                Face *f=faces+k;GfsVariationalFace b;double r[3],c[3];memcpy(r,f->r,sizeof(r));memcpy(c,f->c,sizeof(c));
                if(rotated){rotate(r,dim);rotate(c,dim);}
                if(!gfs_variational_face(dim,r,c,f->area,&b)||
                   !gfs_acoustic_bound_face(state+f->i,&b,m[f->j],V[f->j],U[f->j],gamma))return 2;
                if(!rotated)for(int a=0;a<dim;a++){
                    B[f->i][dim*f->j+a]+=b.bj[a];B[f->i][dim*f->i+a]-=b.bj[a];
                }
            }
            double row=0,column=0;
            for(int i=0;i<n;i++){gfs_acoustic_bound_finish(state+i);row=fmax(row,state[i].row);column=fmax(column,state[i].column);}
            bound[rotated]=row*column;
        }
        double rotation_error=fabs(bound[1]/bound[0]-1);
        if(rotation_error>1e-12)return 1;
        printf("{\"dimension\":%d,\"cells\":%d,\"positive_lambda_bound\":%.17g,\"rotation_relative_error\":%.17g,\"acoustic_factor\":[",dim,n,bound[0],rotation_error);
        for(int i=0;i<n;i++){
            if(i)printf(",");printf("[");
            for(int k=0;k<dim*n;k++)printf("%s%.17g",k?",":"",rootD[i]*B[i][k]/sqrt(m[k/dim]));
            printf("]");
        }
        printf("]}\n");
    }
    return 0;
}
