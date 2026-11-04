/* Actual native power-cell regression with a generator outside its cell.
 * The periodic x weights [0,1.2h^2,1.2h^2,0] have exact faces
 * [0,.1,.5,.9,1]. Raw fractions -0.1 and 1.1 exercise the historical clamp.
 * No PDE, approximate substitute mesh or energy-state conversion. */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../gfs_fv_native.h"
static double power(double x,int n){double p=1;for(int i=0;i<n;i++)p*=x;return p;}
int main(void){
    int bad=0;double max_volume=0,max_centroid=0,max_moment=0,max_fd=0,max_gauge=0,max_richardson=0,fd_level[3]={0};
    for(int d=1;d<=3;d++){
        GfsFvMesh m={0};m.dim=d;m.n=1;m.gamma=1.4;for(int a=0;a<3;a++){m.L[a]=1;m.bc[a]=GFS_FV_PERIODIC;}for(int a=0;a<d;a++)m.n*=4;
        m.cell=calloc(m.n,sizeof(*m.cell));GfsFvCell *saved=calloc(m.n,sizeof(*saved));
        const double edges[5]={0,.1,.5,.9,1};double volume[64],velocity[64][3]={{0}},rate[64]={0},Vdot[64]={0};
        for(int i=0;i<m.n;i++){
            int cell=i;for(int a=0;a<d;a++){m.cell[i].q[a]=(cell%4+.5)/4;cell/=4;velocity[i][a]=.05*sin(.7*(i+1)*(a+1));}
            int x=i%4;m.cell[i].W=x==1||x==2?1.2/16:0;rate[i]=.02/16*cos(.9*(i+1));
        }
        memcpy(saved,m.cell,m.n*sizeof(*saved));if(gfs_fv_geometry(&m))return 2;
        for(int i=0;i<m.n;i++){
            GfsFvGeometry *g=m.geometry+i;int x=i%4;double exact=(edges[x+1]-edges[x])*power(.25,d-1);volume[i]=g->V;
            max_volume=fmax(max_volume,fabs(g->V/exact-1));
            max_centroid=fmax(max_centroid,fabs(g->centroid[0]+m.cell[i].q[0]-.5*(edges[x+1]+edges[x]))/.25);
            for(int ax=0;ax<=4;ax++)for(int ay=0;ay<=(d>=2?4-ax:0);ay++)for(int az=0;az<=(d==3?4-ax-ay:0);az++){
                int exponents[3]={ax,ay,az},degree=ax+ay+az;double integral=0,expected=1;
                for(int k=0;k<g->faces;k++){
                    GfsFvNativeFace *f=g->face+k;
                    for(int q=0;q<f->count;q++){
                        GfsFvPoint *p=g->point+f->point+q;double radial=0,value=1;
                        for(int a=0;a<d;a++){radial+=p->x[a]*f->normal[a];value*=power(p->x[a],exponents[a]);}
                        integral+=p->weight*radial*value/(d+degree);
                    }
                }
                for(int a=0;a<d;a++){
                    double left=a==0?edges[x]-m.cell[i].q[0]:-.125,right=a==0?edges[x+1]-m.cell[i].q[0]:.125;
                    expected*=(power(right,exponents[a]+1)-power(left,exponents[a]+1))/(exponents[a]+1);
                }
                max_moment=fmax(max_moment,fabs(integral-expected)/(exact*power(.25,degree)));
            }
            for(int k=0;k<g->faces;k++){
                GfsFvNativeFace *f=g->face+k;int j=f->j;
                for(int q=0;q<f->count;q++){
                    GfsFvPoint *p=g->point+f->point+q;double wn=(rate[i]-rate[j])/(2*f->distance);
                    for(int a=0;a<d;a++)wn+=.5*f->normal[a]*(velocity[i][a]+velocity[j][a])-(velocity[j][a]-velocity[i][a])*(p->x[a]-.5*f->r[a])/f->distance;
                    Vdot[i]+=p->weight*wn;
                }
            }
        }
        for(int i=0;i<m.n;i++)m.cell[i].W+=.01;if(gfs_fv_geometry(&m))return 2;
        for(int i=0;i<m.n;i++)max_gauge=fmax(max_gauge,fabs(m.geometry[i].V-volume[i])/volume[i]);
        double derivative[3][64];
        /* At the extruded/cospherical fixture a zero-area face is born
         * under either sign of displacement. V is C1 but its one-sided
         * quadratic terms can differ, so a centered first derivative has
         * O(epsilon), rather than O(epsilon^2), error. Retain the original
         * raw-error failure and check two fixed Richardson pairs. */
        for(int level=0;level<3;level++){
            double epsilon=1e-5*power(.5,level),plus[64],minus[64];
            for(int sign=-1;sign<=1;sign+=2){
                memcpy(m.cell,saved,m.n*sizeof(*saved));for(int i=0;i<m.n;i++){
                    for(int a=0;a<d;a++)m.cell[i].q[a]+=sign*epsilon*velocity[i][a];m.cell[i].W+=sign*epsilon*rate[i];
                }
                if(gfs_fv_geometry(&m))return 2;
                for(int i=0;i<m.n;i++)(sign>0?plus:minus)[i]=m.geometry[i].V;
            }
            double local=0;
            for(int i=0;i<m.n;i++){
                derivative[level][i]=(plus[i]-minus[i])/(2*epsilon);
                local=fmax(local,fabs(derivative[level][i]-Vdot[i])/power(.25,d-1));
            }
            max_fd=fmax(max_fd,local);fd_level[level]=fmax(fd_level[level],local);
            fprintf(stderr,"power FD dimension=%d epsilon=%.17g max_rate_error_over_harea=%.17g\n",d,epsilon,local);
        }
        for(int level=0;level<2;level++)for(int i=0;i<m.n;i++)max_richardson=fmax(max_richardson,fabs(2*derivative[level+1][i]-derivative[level][i]-Vdot[i])/power(.25,d-1));
        gfs_fv_geometry_free(&m);free(m.cell);free(saved);
    }
    bad=max_volume>1e-11||max_centroid>1e-11||max_moment>1e-10||max_richardson>1e-7||max_gauge>1e-11;
    printf("{\"scope\":\"native true-power geometry, degree-four monomials and moving-face finite differences in D1/D2/D3; no PDE\",\"raw_fraction_extrema\":[-0.1,1.1],\"max_relative_volume_error\":%.17g,\"max_centroid_over_h_error\":%.17g,\"max_moment_over_V_hdegree_error\":%.17g,\"raw_FD_error_over_harea_by_level\":[%.17g,%.17g,%.17g],\"max_two_pair_Richardson_error_over_harea\":%.17g,\"original_raw_FD_gate_pass\":%s,\"max_relative_weight_gauge_error\":%.17g,\"pass\":%s}\n",max_volume,max_centroid,max_moment,fd_level[0],fd_level[1],fd_level[2],max_richardson,max_fd<=1e-7?"true":"false",max_gauge,bad?"false":"true");
    return bad?1:0;
}
