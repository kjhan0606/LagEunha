/* Independent interval proof of nodal versus cell-average linear pressure.
 * Uses the existing native variational face/pressure helper. No PDE steps,
 * remeshing, force change, limiter, weight solve or new pressure definition. */
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "../gfs_variational_core.h"

int main(void)
{
    const int primitive=16;
    double mean_errors[2];
    printf("{\"scope\":\"independent actual-interval and native-helper pressure-location proof; zero PDE steps\",\"cases\":[");
    for(int level=0;level<2;level++){
        const int repeats=level?6:3,n=primitive*repeats;
        double x[96],V[96],center[96],offset[96];
        for(int i=0;i<n;i++)x[i]=(i+.5+.07*sin(2.6*(i%primitive+1)))/n;
        for(int i=0;i<n;i++){
            double left=x[(i+n-1)%n]-(i==0),right=x[(i+1)%n]+(i==n-1);
            V[i]=.5*(right-left);
            center[i]=.25*(left+2*x[i]+right);
            offset[i]=center[i]-x[i];
        }
        double defect=0,nodal_max=0,helper_nodal_max=0,helper_average_max=0;
        long double weighted=0,volume=0,analytic_weighted=0;
        int first=(repeats/2)*primitive;
        for(int i=first;i<first+primitive;i++){
            int left=i-1,right=i+1;
            double F_nodal=.5*((2+x[left])-(2+x[right]));
            double F_average=.5*((2+center[left])-(2+center[right]));
            double error=(F_average+V[i])/V[i];
            double analytic=(offset[left]-offset[right])/(2*V[i]);
            defect=fmax(defect,fabs(error-analytic));
            nodal_max=fmax(nodal_max,fabs((F_nodal+V[i])/V[i]));
            double Fhelper[2]={0,0};
            for(int sign=-1;sign<=1;sign+=2){
                int j=i+sign;
                double r[3]={x[j]-x[i],0,0},c[3]={.5*(x[i]+x[j])-x[i],0,0};
                GfsVariationalFace face;
                if(!gfs_variational_face(1,r,c,1,&face))return 2;
                double force[3];
                gfs_variational_pressure(&face,2+x[i],2+x[j],force);Fhelper[0]+=force[0];
                gfs_variational_pressure(&face,2+center[i],2+center[j],force);Fhelper[1]+=force[0];
            }
            helper_nodal_max=fmax(helper_nodal_max,fabs((Fhelper[0]-F_nodal)/V[i]));
            helper_average_max=fmax(helper_average_max,fabs((Fhelper[1]-F_average)/V[i]));
            weighted+=(long double)V[i]*fabs(error);
            analytic_weighted+=(long double)V[i]*fabs(analytic);
            volume+=V[i];
        }
        mean_errors[level]=(double)(weighted/volume);
        if(fmax(fmax(defect,nodal_max),fmax(helper_nodal_max,helper_average_max))>1e-11)return 2;
        printf("%s{\"cells\":%d,\"samples\":%d,\"cell_average_force_density_error_L1\":%.17g,\"closed_form_error_L1\":%.17g,\"closed_form_max_defect\":%.17g,\"nodal_force_density_max_error\":%.17g,\"native_helper_nodal_max_defect\":%.17g,\"native_helper_average_max_defect\":%.17g}",level?",":"",n,primitive,mean_errors[level],(double)(analytic_weighted/volume),defect,nodal_max,helper_nodal_max,helper_average_max);
    }
    double ratio=mean_errors[1]/mean_errors[0];
    if(fabs(ratio-1)>1e-10)return 2;
    printf("],\"refinement_error_ratio\":%.17g,\"accepted_steps\":0,\"interpretation\":\"A cell-average linear pressure differs from its generator sample. The existing B force is exact for the nodal representation; its cell-average force-density defect retains the same amplitude on this repeated irregular pattern. This is not a theorem excluding weak convergence or all Euler discretizations.\"}\n",ratio);
    return 0;
}
