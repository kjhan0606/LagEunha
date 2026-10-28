#include "full_moment_coupled_state.h"
#include <math.h>
#include <float.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

int full_moment_centered_decode(const double a[20],double mass,FullMomentCenteredState *out){
 if(!a||!out||!isfinite(mass)||mass<=0.)return 1;
 for(int k=0;k<20;k++)if(!isfinite(a[k]))return 2;
 if(a[6]!=a[7]||a[16]!=a[17])return 3;
 FullMomentCenteredState s={0};long double m=mass;
 long double q[2]={a[3],a[4]},p[2]={a[9],a[10]},Pi[4];
 for(int i=0;i<2;i++){s.center[i]=q[i]/m;s.velocity[i]=p[i]/m;}
 for(int i=0;i<2;i++)for(int j=0;j<2;j++){
  int k=2*i+j;
  s.S[k]=(long double)a[5+k]-(q[i]*q[j])/m;
  s.C[k]=(long double)a[11+k]-(q[i]*p[j])/m;
  Pi[k]=(long double)a[15+k]-(p[i]*p[j])/m;
 }
 if(!(s.S[0]>0.))return 4;
 long double l00=sqrtl(s.S[0]),l10=s.S[2]/l00,rest=s.S[3]-l10*l10;
 if(!(rest>0.))return 4;
 long double l11=sqrtl(rest),y[4];
 for(int j=0;j<2;j++){y[j]=s.C[j]/l00;y[2+j]=(s.C[2+j]-l10*y[j])/l11;}
 for(int i=0;i<2;i++)for(int j=0;j<2;j++)
  s.R[2*i+j]=Pi[2*i+j]-(y[i]*y[j]+y[2+i]*y[2+j]);
 if(!(s.R[0]>0.))return 5;
 long double r00=sqrtl(s.R[0]),r10=s.R[2]/r00;
 if(!(s.R[3]-r10*r10>0.))return 5;
 s.thermal=(long double)a[19]-.5L*((long double)a[15]+a[18]);
 if(!(s.thermal>0.))return 6;
 for(int k=0;k<4;k++)if(!isfinite(s.S[k])||!isfinite(s.C[k])||!isfinite(s.R[k]))return 2;
 if(!isfinite(s.thermal))return 2;
 *out=s;return 0;
}

int full_moment_coupled_midpoint(int n,const double *base,const double *trial,
 const double *mass,const long double *base_Psi,const long double *trial_Psi,
 double *midpoint,long double *midpoint_Psi,double *midpoint_L,int *failed_cell){
 if(failed_cell)*failed_cell=-1;
 if(n<1||!base||!trial||!mass||!base_Psi||!trial_Psi||!midpoint||!midpoint_Psi||!midpoint_L)return 1;
 if((size_t)n>(SIZE_MAX/sizeof(double)-4)/20)return 1;
 size_t count=(size_t)n*20+4;
 double *state=malloc(count*sizeof(double)),*L=malloc((size_t)n*4*sizeof(double));
 long double *Psi=malloc((size_t)n*sizeof(long double));
 if(!state||!L||!Psi){free(state);free(L);free(Psi);return 7;}
 int status=0;
 for(size_t k=0;k<count;k++){
  state[k]=.5*(base[k]+trial[k]);
  if(!isfinite(base[k])||!isfinite(trial[k])||!isfinite(state[k]))status=2;
 }
 for(int i=0;i<n&&!status;i++){
  FullMomentCenteredState s;
  status=full_moment_centered_decode(state+20*i,mass[i],&s);
  if(!status){
   long double a=sqrtl(s.R[0]),b=s.R[2]/a,c=sqrtl(s.R[3]-b*b);
   L[4*i]=(double)a;L[4*i+1]=0.;L[4*i+2]=(double)b;L[4*i+3]=(double)c;
   Psi[i]=.5L*(base_Psi[i]+trial_Psi[i]);
   if(!isfinite(base_Psi[i])||!isfinite(trial_Psi[i])||!isfinite(Psi[i]))status=2;
   if(!isfinite(L[4*i])||!isfinite(L[4*i+2])||!isfinite(L[4*i+3])||L[4*i]<=0.||L[4*i+3]<=0.)status=5;
  }
  if(status&&failed_cell)*failed_cell=i;
 }
 if(!status){memcpy(midpoint,state,count*sizeof(double));memcpy(midpoint_L,L,(size_t)n*4*sizeof(double));memcpy(midpoint_Psi,Psi,(size_t)n*sizeof(long double));}
 free(state);free(L);free(Psi);return status;
}

static long double norm_difference(const long double *a,const long double *b,int n){
 long double sum=0.;for(int k=0;k<n;k++){long double d=a[k]-b[k];sum+=d*d;}return sqrtl(sum);
}
static void potential_features(const double *a,long double m,const FullMomentPotentialNorm *p,long double *v){
 long double ox=p->origin[0],oy=p->origin[1],s=p->scale,s2=s*s;
 long double qx=a[3],qy=a[4];
 v[0]=m;v[1]=(qx-m*ox)/s;v[2]=(qy-m*oy)/s;
 v[3]=((long double)a[5]-2*ox*qx+m*ox*ox)/s2;
 v[4]=((long double)a[6]-ox*qy-oy*qx+m*ox*oy)/s2;
 v[5]=((long double)a[8]-2*oy*qy+m*oy*oy)/s2;
}
static long double dependent_ratio(const double *a,long double m,long double psi,const FullMomentPotentialNorm *p){
 long double v[6],predicted=0.,absolute=fabsl(psi);potential_features(a,m,p,v);
 for(int k=0;k<6;k++){predicted+=v[k]*p->reference[k];absolute+=fabsl(v[k])*fabsl(p->reference[k]);}
 long double error=fabsl(psi-predicted),bound=64*LDBL_EPSILON*absolute;
 return bound>0.?error/bound:(error==0.?0.:INFINITY);
}
int full_moment_coupled_residual(int n,double gamma,const double *base,
 const double *mass,const double *base_area,const double *trial,const double *mapped,
 const long double *trial_Psi,const long double *mapped_Psi,
 const FullMomentPotentialNorm *potential,FullMomentCoupledResidual *out,int *failed_cell){
 if(failed_cell)*failed_cell=-1;
 if(n<1||!isfinite(gamma)||gamma<=1.||!base||!mass||!base_area||!trial||!mapped||!trial_Psi||!mapped_Psi||!potential||!out)return 1;
 FullMomentCoupledResidual result={0};for(int k=0;k<FM_RES_COUNT;k++)result.cell[k]=-1;
 long double weight_gauge=0.;
 for(int i=0;i<n;i++)weight_gauge+=(long double)mapped[20*i+2]-trial[20*i+2];
 weight_gauge/=n;
 for(int i=0;i<n;i++){
  const double *b=base+20*i,*a=trial+20*i,*c=mapped+20*i;
  FullMomentCenteredState sb,sa,sc;int status;
  if((status=full_moment_centered_decode(b,mass[i],&sb))||
     (status=full_moment_centered_decode(a,mass[i],&sa))||
     (status=full_moment_centered_decode(c,mass[i],&sc))){if(failed_cell)*failed_cell=i;return status;}
  if(!isfinite(base_area[i])||base_area[i]<=0.){if(failed_cell)*failed_cell=i;return 8;}
  long double m=mass[i],area=base_area[i],length=sqrtl(area);
  long double v=hypotl(sb.velocity[0],sb.velocity[1]);
  long double sound=sqrtl((long double)gamma*((long double)gamma-1)*sb.thermal/m);
  if(sound>v)v=sound;
  if(!(v>0.)||!isfinite(v)){if(failed_cell)*failed_cell=i;return 8;}
  long double dx[2]={(long double)a[0],a[1]},dy[2]={(long double)c[0],c[1]};
  long double values[FM_RES_COUNT]={0};
  values[FM_RES_X]=norm_difference(dx,dy,2)/length;
  values[FM_RES_W]=fabsl(((long double)c[2]-a[2])-weight_gauge)/area;
  values[FM_RES_CENTER]=norm_difference(sa.center,sc.center,2)/length;
  values[FM_RES_VELOCITY]=norm_difference(sa.velocity,sc.velocity,2)/v;
  values[FM_RES_S]=norm_difference(sa.S,sc.S,4)/(m*area);
  values[FM_RES_C]=norm_difference(sa.C,sc.C,4)/(m*length*v);
  values[FM_RES_R]=norm_difference(sa.R,sc.R,4)/(m*v*v);
  values[FM_RES_THERMAL]=fabsl(sa.thermal-sc.thermal)/sb.thermal;
  const FullMomentPotentialNorm *p=potential+i;
  if(p->dependent!=0&&p->dependent!=1){if(failed_cell)*failed_cell=i;return 8;}
  if(!isfinite(p->scale)||p->scale<=0.||!isfinite(trial_Psi[i])||!isfinite(mapped_Psi[i])){if(failed_cell)*failed_cell=i;return 8;}
  for(int k=0;k<2;k++)if(!isfinite(p->origin[k])){if(failed_cell)*failed_cell=i;return 8;}
  for(int k=0;k<6;k++)if(!isfinite(p->reference[k])||!isfinite(p->projection[k])){if(failed_cell)*failed_cell=i;return 8;}
  if(p->dependent){
   long double first=dependent_ratio(a,m,trial_Psi[i],p),second=dependent_ratio(c,m,mapped_Psi[i],p);
   values[FM_RES_DEPENDENT_PSI]=first>second?first:second;result.dependent_cells++;
  }else{
   if(!isfinite(p->sigma)||!isfinite(p->rank_error)||p->rank_error<0.||!(p->sigma>p->rank_error)){if(failed_cell)*failed_cell=i;return 8;}
   // Subtract feature DIFFERENCES directly, avoiding a second cancellation
   // of two large reconstructed potential integrals.
   long double qx=(long double)c[3]-a[3],qy=(long double)c[4]-a[4],s=p->scale,s2=s*s;
   long double delta[6]={0.,qx/s,qy/s,
    (((long double)c[5]-a[5])-2*p->origin[0]*qx)/s2,
    (((long double)c[6]-a[6])-p->origin[0]*qy-p->origin[1]*qx)/s2,
    (((long double)c[8]-a[8])-2*p->origin[1]*qy)/s2};
   long double correction=0.;for(int k=1;k<6;k++)correction+=p->projection[k]*delta[k];
   long double remainder=fabsl((mapped_Psi[i]-trial_Psi[i])-correction);
   if(remainder>result.psi_remainder_numerator_max)result.psi_remainder_numerator_max=remainder;
   values[FM_RES_PSI_REMAINDER]=remainder/(m*p->sigma);
  }
  for(int k=0;k<FM_RES_COUNT;k++){
   if(isnan(values[k])){if(failed_cell)*failed_cell=i;return 8;}
   if(values[k]>result.maximum[k]){result.maximum[k]=values[k];result.cell[k]=i;}
  }
 }
 *out=result;return 0;
}
int full_moment_coupled_residual_vector(int n,double gamma,const double *base,
 const double *mass,const double *base_area,const double *trial,const double *mapped,
 const long double *trial_Psi,const long double *mapped_Psi,
 const FullMomentPotentialNorm *potential,long double *out,size_t out_count,
 int *failed_cell){
 if(failed_cell)*failed_cell=-1;
 if(n<1||!isfinite(gamma)||gamma<=1.||!base||!mass||!base_area||!trial||!mapped||
    !trial_Psi||!mapped_Psi||!potential||!out||(size_t)n>SIZE_MAX/21||out_count!=(size_t)n*21)return 1;
 long double weight_gauge=0.;
 for(int i=0;i<n;i++)weight_gauge+=(long double)mapped[20*i+2]-trial[20*i+2];
 weight_gauge/=n;
 size_t cursor=0;
 for(int i=0;i<n;i++){
  const double *b=base+20*i,*a=trial+20*i,*c=mapped+20*i;
  FullMomentCenteredState sb,sa,sc;int status;
  if((status=full_moment_centered_decode(b,mass[i],&sb))||
     (status=full_moment_centered_decode(a,mass[i],&sa))||
     (status=full_moment_centered_decode(c,mass[i],&sc))){if(failed_cell)*failed_cell=i;return status;}
  long double area=base_area[i],m=mass[i];
  if(!isfinite(area)||area<=0.){if(failed_cell)*failed_cell=i;return 8;}
  long double length=sqrtl(area);
  long double v=hypotl(sb.velocity[0],sb.velocity[1]);
  long double sound=sqrtl((long double)gamma*((long double)gamma-1.)*sb.thermal/m);
  if(sound>v)v=sound;
  if(!(v>0.)||!isfinite(v)){if(failed_cell)*failed_cell=i;return 8;}
  out[cursor++]=((long double)c[0]-a[0])/length;
  out[cursor++]=((long double)c[1]-a[1])/length;
  out[cursor++]=(((long double)c[2]-a[2])-weight_gauge)/area;
  out[cursor++]=(sc.center[0]-sa.center[0])/length;
  out[cursor++]=(sc.center[1]-sa.center[1])/length;
  out[cursor++]=(sc.velocity[0]-sa.velocity[0])/v;
  out[cursor++]=(sc.velocity[1]-sa.velocity[1])/v;
  for(int k=0;k<4;k++)out[cursor++]=(sc.S[k]-sa.S[k])/(m*area);
  for(int k=0;k<4;k++)out[cursor++]=(sc.C[k]-sa.C[k])/(m*length*v);
  for(int k=0;k<4;k++)out[cursor++]=(sc.R[k]-sa.R[k])/(m*v*v);
  out[cursor++]=(sc.thermal-sa.thermal)/sb.thermal;
  const FullMomentPotentialNorm *p=potential+i;
  if(!isfinite(trial_Psi[i])||!isfinite(mapped_Psi[i])||p->dependent<0||p->dependent>1){if(failed_cell)*failed_cell=i;return 8;}
  if(p->dependent){
   out[cursor++]=0.;
  }else{
   if(!isfinite(p->scale)||p->scale<=0.||!isfinite(p->sigma)||!isfinite(p->rank_error)||
      p->rank_error<0.||!(p->sigma>p->rank_error)||!isfinite(p->origin[0])||!isfinite(p->origin[1])){
    if(failed_cell)*failed_cell=i;return 8;
   }
   for(int k=0;k<6;k++)if(!isfinite(p->projection[k])){if(failed_cell)*failed_cell=i;return 8;}
   long double qx=(long double)c[3]-a[3],qy=(long double)c[4]-a[4];
   long double s=p->scale,s2=s*s;
   long double delta[6]={0.,qx/s,qy/s,
    (((long double)c[5]-a[5])-2*p->origin[0]*qx)/s2,
    (((long double)c[6]-a[6])-p->origin[0]*qy-p->origin[1]*qx)/s2,
    (((long double)c[8]-a[8])-2*p->origin[1]*qy)/s2};
   long double correction=0.;for(int k=1;k<6;k++)correction+=p->projection[k]*delta[k];
   out[cursor++]=((mapped_Psi[i]-trial_Psi[i])-correction)/(m*p->sigma);
  }
  for(size_t k=cursor-21;k<cursor;k++)if(!isfinite(out[k])){if(failed_cell)*failed_cell=i;return 8;}
 }
 return cursor==out_count?0:9;
}
int full_moment_coupled_converged(const FullMomentCoupledResidual *r){
 if(!r)return 0;
 for(int k=0;k<=FM_RES_R;k++)if(!isfinite(r->maximum[k])||r->maximum[k]>1e-10L)return 0;
 for(int k=FM_RES_THERMAL;k<=FM_RES_PSI_REMAINDER;k++)if(!isfinite(r->maximum[k])||r->maximum[k]>1e-8L)return 0;
 return isfinite(r->maximum[FM_RES_DEPENDENT_PSI])&&r->maximum[FM_RES_DEPENDENT_PSI]<=1.;
}
