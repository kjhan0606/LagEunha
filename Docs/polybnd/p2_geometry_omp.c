#include <stdint.h>
#include <omp.h>
static int used_threads=0;
int p2_threads(void) { return used_threads; }
void p2_geometry(int64_t n, const int64_t *ids, const double *x,
                 double *area, double *da, double *poly) {
#pragma omp parallel
  {
#pragma omp single
    used_threads=omp_get_num_threads();
#pragma omp for schedule(static)
    for(int64_t i=0;i<n;i++) {
      double a[2],b[2],c[2];
      for(int k=0;k<2;k++) {
        a[k]=x[2*ids[3*i]+k]; b[k]=x[2*ids[3*i+1]+k]; c[k]=x[2*ids[3*i+2]+k];
        poly[6*i+k]=a[k];poly[6*i+2+k]=b[k];poly[6*i+4+k]=c[k];
      }
      area[i]=0.5*((b[0]-a[0])*(c[1]-a[1])-(c[0]-a[0])*(b[1]-a[1]));
      da[6*i]=-0.5*(c[1]-b[1]);da[6*i+1]=0.5*(c[0]-b[0]);
      da[6*i+2]=-0.5*(a[1]-c[1]);da[6*i+3]=0.5*(a[0]-c[0]);
      da[6*i+4]=-0.5*(b[1]-a[1]);da[6*i+5]=0.5*(b[0]-a[0]);
    }
  }
}
