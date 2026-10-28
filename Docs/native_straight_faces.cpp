// Same quadratic roots and ordered positive Gauss panels as the Python oracle.
// Long-double coordinates/speeds retain the audited tiny-face representation.
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>

extern "C" int lag_line_cuts(
    int64_t count, const long double *A, const long double *delta,
    const long double *samples, const double *base_dot, const double *dg,
    const double *position, const double *distance, const double *dz,
    const int64_t *prior_offsets, const double *prior, int64_t capacity,
    double *cuts, int64_t *sizes, int64_t *status) {
  if (capacity < 22 || count < 0) return 1;
  #pragma omp parallel for schedule(static)
  for (int64_t f=0; f<count; ++f) {
    long double relative[5][8];
    const double times[5]={0.,.5,1.,.25,.75};
    for (int t=0;t<5;++t) {
      const long double px=A[2*f]+static_cast<long double>(times[t])*delta[2*f];
      const long double py=A[2*f+1]+static_cast<long double>(times[t])*delta[2*f+1];
      const long double term0=dg[2*f]*(px-position[2*f]);
      const long double term1=dg[2*f+1]*(py-position[2*f+1]);
      const long double dot=(0.L+term0)+term1;
      const long double wn=(base_dot[f]-dot)/distance[f]+dz[f]/(2*distance[f]);
      for (int k=0;k<8;++k) relative[t][k]=samples[40*f+8*t+k]-wn;
    }
    double *out=cuts+capacity*f;int used=2;out[0]=0.;out[1]=1.;status[f]=0;
    for (int k=0;k<8;++k) {
      const long double left=relative[0][k],middle=relative[1][k],right=relative[2][k];
      const long double c2=2*(right+left-2*middle),c1=right-left-c2;
      for (int t=3;t<5;++t) {
        const double x=times[t];
        const long double predicted=(left+x*c1)+(x*x)*c2;
        const long double bound=2e-12+2e-12*std::abs(predicted);
        if (!(std::abs(relative[t][k]-predicted)<=bound)) status[f]=2;
      }
      long double a=c2,b=c1,c=left;
      const long double scale=std::max(std::abs(a),std::max(std::abs(b),std::abs(c)));
      if (scale==0) continue;
      a/=scale;b/=scale;c/=scale;
      long double roots[2];int nr=0;
      if (std::abs(a)<=32*std::numeric_limits<double>::epsilon()) {
        if (b!=0) roots[nr++]=-c/b;
      } else {
        const long double disc=b*b-(4*a)*c;
        if (disc<0) continue;
        const long double q=-.5L*(b+std::copysign(std::sqrt(disc),b));
        if (q==0) roots[nr++]=-b/(2*a);
        else {roots[nr++]=q/a;roots[nr++]=c/q;}
      }
      for (int r=0;r<nr;++r) if (roots[r]>0 && roots[r]<1) {
        if (used>=capacity) {status[f]=3;break;}
        out[used++]=static_cast<double>(roots[r]);
      }
    }
    for (int64_t p=prior_offsets[f];p<prior_offsets[f+1];++p) {
      if (used>=capacity) {status[f]=3;break;}
      out[used++]=prior[p];
    }
    std::sort(out,out+used);sizes[f]=std::unique(out,out+used)-out;
    if (out[0]!=0. || out[sizes[f]-1]!=1.) status[f]=4;
  }
  return 0;
}

extern "C" int lag_line_fill(
    int64_t count,int64_t order,int64_t capacity,const long double *A,
    const long double *delta,const long double *normal,const long double *length,
    const double *cuts,const int64_t *sizes,const int64_t *offsets,
    const double *g,const double *w,double *points,double *normals,double *measure) {
  if (count<0 || order<1) return 1;
  #pragma omp parallel for schedule(static)
  for (int64_t f=0;f<count;++f) {
    int64_t out=offsets[f];
    for (int64_t panel=0;panel<sizes[f]-1;++panel) {
      const double lo=cuts[capacity*f+panel],hi=cuts[capacity*f+panel+1];
      for (int64_t k=0;k<order;++k,++out) {
        const double t=lo+(hi-lo)*g[k],weight=(hi-lo)*w[k];
        points[2*out]=A[2*f]+static_cast<long double>(t)*delta[2*f];
        points[2*out+1]=A[2*f+1]+static_cast<long double>(t)*delta[2*f+1];
        normals[2*out]=normal[2*f];normals[2*out+1]=normal[2*f+1];
        measure[out]=length[f]*weight;
      }
    }
  }
  return 0;
}
