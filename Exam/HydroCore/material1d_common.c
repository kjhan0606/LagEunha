/* Fixed-mass 1D counterpart of the current 2D reconstruction candidate.
 * Validation and limitations: Docs/MATERIAL_RT_PORT_AUDIT_2026-11-06.md.
 * This driver does not certify 2D/3D mesh evolution or replace orbital code. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "gfs_exact_riemann_ale.h"
#include "gfs_material_pressure_offset.h"
#include "gfs_material_gradient_project.h"
#include "gfs_material_ordered_jump.h"
#define CAP 16384
#define GFS_FV_MOVING_TRANSMIT 0
#define GFS_FV_REFLECT 1
enum {CONTACT,MOVING_CONTACT,SOD,GRESHO,COLD_SHEAR,ALTERNATING,BLAST,COLLISION,SHU_OSHER,NOH,LAX,DOUBLE_FAN,CONTACT_OPEN,MOVING_CONTACT_OPEN};
#include "material1d_cases.h"
typedef struct{double x[CAP+1],p[CAP],E[CAP],energy_boundary,momentum_boundary;} State;
static double mass[CAP],gamma_=1.4;static int N,wall_left,wall_right,order=2;static GfsFvPrimitive q[CAP];
static int primitive(const State*s){for(int i=0;i<N;i++){double V=s->x[i+1]-s->x[i],Q[5]={mass[i],s->p[i],0,0,s->E[i]};if(!gfs_fv_primitive(1,Q,V,gamma_,q+i))return 0;}return 1;}
static double minmod(double a,double b){return a*b<=0?0:copysign(fmin(fabs(a),fabs(b)),a);}
static GfsFvPrimitive slope[CAP];
static int reconstruct(const State*s){
 memset(slope,0,sizeof(slope));if(order==1)return 1;
 for(int i=1;i<N-1;i++){
  double dl=.5*(s->x[i+1]-s->x[i-1]),dr=.5*(s->x[i+2]-s->x[i]);
  slope[i].rho=minmod((q[i].rho-q[i-1].rho)/dl,(q[i+1].rho-q[i].rho)/dr);
  double lo=fmin(q[i].p,fmin(q[i-1].p,q[i+1].p)),hi=fmax(q[i].p,fmax(q[i-1].p,q[i+1].p));
  double center=.5*(s->x[i]+s->x[i+1]),M=0,B=0,A[6],bounds[6],g0[1],g[1],radius[2];int nc=0,it=0;
  for(int side=0;side<2;side++){
   int j=side?i+1:i-1;double delta[1]={side?dr:-dl},ri[1]={(side?s->x[i+1]:s->x[i])-center},effective[1];radius[side]=ri[0];
   if(!gfs_material_pressure_offset(1,q[i].rho,q[j].rho,ri,delta,effective)||effective[0]==0)return 0;
   double distance=fabs(delta[0]),weight=1/(distance*distance*distance),e=effective[0];M+=weight*e*e;B+=weight*e*(q[j].p-q[i].p);
   A[nc]=-copysign(1,ri[0]);bounds[nc++]=q[i].p/fabs(ri[0]);
   A[nc]=copysign(1,e);bounds[nc++]=(hi-q[i].p)/fabs(e);
   A[nc]=-copysign(1,e);bounds[nc++]=(q[i].p-lo)/fabs(e);
  }
  if(!(M>0))return 0;g0[0]=B/M;if(!mg_project(1,nc,A,bounds,g0,g,&it))return 0;
  double theta=1;
  for(int side=0;side<2;side++){double delta=g[0]*radius[side],margin=4096*DBL_EPSILON*(q[i].p+fabs(delta));if(q[i].p+delta<=margin&&delta<0)theta=fmin(theta,nextafter((q[i].p-margin)/(-delta),0.));}
  slope[i].p=theta*g[0];
  slope[i].u[0]=((q[i].u[0]-q[i-1].u[0])/(dl*dl)+(q[i+1].u[0]-q[i].u[0])/(dr*dr))/(1/dl+1/dr);
 }
 return 1;
}
static GfsFvPrimitive trace(const State*s,int i,int face){
 GfsFvPrimitive out=q[i];double dx=s->x[face]-.5*(s->x[i]+s->x[i+1]);
 out.rho+=dx*slope[i].rho;out.p+=dx*slope[i].p;out.u[0]+=dx*slope[i].u[0];return out;
}
static int rhs(const State*s,State*r,double*stable){
 if(!primitive(s))return 0;
 if(!reconstruct(s))return 0;
 double pf[CAP+1]={0};*stable=1e30;
 for(int f=0;f<=N;f++){
  GfsFvPrimitive l=trace(s,f?f-1:0,f),rr=trace(s,f<N?f:N-1,f);
  if(f>0&&f<N){
   double left=q[f-1].u[0],right=q[f].u[0],il=l.u[0]-left,ir=rr.u[0]-right;
   double jump=right-left,change=ir-il,alpha=1;
   if(!gfs_material_ordered_jump_alpha(jump,change,&alpha))return 0;
   l.u[0]=left+alpha*il;rr.u[0]=right+alpha*ir;
  }
  if(!(l.rho>0&&rr.rho>0&&l.p>0&&rr.p>0))return 0;
  if(wall_left&&f==0){l.u[0]=-rr.u[0];}
  if(wall_right&&f==N){rr.u[0]=-l.u[0];}
  GfsExactStar star;if(!gfs_exact_star(l,rr,l.u[0],rr.u[0],gamma_,&star)||star.vacuum)return 0;
  pf[f]=star.p;r->x[f]=((wall_left&&f==0)||(wall_right&&f==N))?0:star.u;
 }
 r->energy_boundary=pf[0]*r->x[0]-pf[N]*r->x[N];
 r->momentum_boundary=pf[0]-pf[N];
 for(int i=0;i<N;i++){
  r->p[i]=pf[i]-pf[i+1];r->E[i]=pf[i]*r->x[i]-pf[i+1]*r->x[i+1];
  double V=s->x[i+1]-s->x[i],c=sqrt(gamma_*q[i].p/q[i].rho);
  double speed=c+fmax(fabs(r->x[i]-q[i].u[0]),fabs(r->x[i+1]-q[i].u[0]));
  *stable=fmin(*stable,.3*V/speed);
  if(r->x[i]>r->x[i+1])*stable=fmin(*stable,.3*V/(r->x[i]-r->x[i+1]));
 }
 return 1;
}
static void advance(State*out,const State*in,const State*r,double dt){for(int i=0;i<=N;i++)out->x[i]=in->x[i]+dt*r->x[i];for(int i=0;i<N;i++){out->p[i]=in->p[i]+dt*r->p[i];out->E[i]=in->E[i]+dt*r->E[i];}}
static void dump(const char*name,const State*s){FILE*f=fopen(name,"wx");if(!f)exit(8);if(!primitive(s))exit(9);fprintf(f,"# xleft xright mass rho velocity pressure energy\n");for(int i=0;i<N;i++)fprintf(f,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",s->x[i],s->x[i+1],mass[i],q[i].rho,q[i].u[0],q[i].p,s->E[i]);fclose(f);}
int main(int argc,char**argv){
 if(argc!=3)return 2;
 int side=atoi(argv[2]),problem=-1;if(side<4||side>8192)return 2;
 const char*name=argv[1];const char*names[]={"sod","blast","noh","shu","lax","doublefan","collision","contact_open","moving_contact_open"};
 int ids[]={SOD,BLAST,NOH,SHU_OSHER,LAX,DOUBLE_FAN,COLLISION,CONTACT_OPEN,MOVING_CONTACT_OPEN};
 for(int k=0;k<9;k++){if(!strcmp(name,names[k]))problem=ids[k];}if(problem<0)return 2;
 FvOriginalCase sp;
 if(problem==SOD){sp=(FvOriginalCase){0,1,{.5,1},1.4,.2,{1,.125,1},{0,0,0},{1,.1,1},2,GFS_FV_REFLECT,GFS_FV_REFLECT,0};}
 else if(!original_case(problem,&sp))return 2;
 gamma_=sp.gamma;wall_left=sp.lower==GFS_FV_REFLECT;wall_right=sp.upper==GFS_FV_REFLECT;
 int contact=problem==CONTACT_OPEN||problem==MOVING_CONTACT_OPEN;
 double h=(sp.hi-sp.lo)/side,signal=0;for(int k=0;k<sp.regions;k++)signal=fmax(signal,fabs(sp.u[k])+sqrt(sp.gamma*sp.p[k]/(sp.sinusoidal&&k==1?.8:sp.rho[k])));
 int guard=(int)ceil(signal*sp.end/h)+4,padl=wall_left?0:guard,padr=wall_right?0:guard;
 N=side+padl+padr;if(N>CAP)return 2;
 double stop=sp.end;static State s,a,b,k1,k2;
 for(int i=0;i<=N;i++)s.x[i]=sp.lo+(i-padl)*h;
 double E0=0,P0=0,M0=0;
 for(int i=0;i<N;i++){
  double thermal=0;
  for(int k=0;k<sp.regions;k++){
   double lo=fmax(s.x[i],k?sp.cut[k-1]:s.x[i]),hi=fmin(s.x[i+1],k<sp.regions-1?sp.cut[k]:s.x[i+1]);if(!(hi>lo))continue;
   double m=sp.rho[k]*(hi-lo);if(sp.sinusoidal&&k==1)m=(hi-lo)+.2*(cos(5*acos(-1.)*lo)-cos(5*acos(-1.)*hi))/(5*acos(-1.));
   mass[i]+=m;s.p[i]+=m*sp.u[k];thermal+=(hi-lo)*sp.p[k]/(gamma_-1);
  }
  s.E[i]=thermal+.5*s.p[i]*s.p[i]/mass[i];E0+=s.E[i];P0+=s.p[i];M0+=mass[i];
 }
 FILE*spec=fopen("problem_spec.json","wx");if(!spec)return 8;
 fprintf(spec,"{\"case\":\"%s\",\"physical_domain\":[%.17g,%.17g],\"physical_resolution\":%d,\"total_cells\":%d,\"padding_cells\":[%d,%d],\"gamma\":%.17g,\"target\":%.17g,\"lower_reflecting\":%d,\"upper_reflecting\":%d,\"order\":2}\n",name,sp.lo,sp.hi,side,N,padl,padr,gamma_,stop,wall_left,wall_right);fclose(spec);
 dump("initial.dat",&s);double t=0,eb=0,pb=0;int steps=0,retries=0,failed=0;
 while(t<stop&&steps<1000000){double dt,lim;if(!rhs(&s,&k1,&dt)){failed=1;break;}dt=fmin(dt,stop-t);int ok=0;
  for(int retry=0;retry<40;retry++){
   advance(&a,&s,&k1,dt);
   if(rhs(&a,&k2,&lim)&&dt<=lim*(1+1e-12)){
    advance(&b,&a,&k2,dt);for(int i=0;i<=N;i++)b.x[i]=.5*(s.x[i]+b.x[i]);for(int i=0;i<N;i++){b.p[i]=.5*(s.p[i]+b.p[i]);b.E[i]=.5*(s.E[i]+b.E[i]);}
    if(primitive(&b)){ok=1;break;}
   }
   dt*=.5;retries++;
  }
  if(!ok||!(t+dt>t)){failed=1;break;}
  eb+=.5*dt*(k1.energy_boundary+k2.energy_boundary);pb+=.5*dt*(k1.momentum_boundary+k2.momentum_boundary);
  s=b;t+=dt;steps++;
 }
 if(t<stop){failed=1;}
 dump("final.dat",&s);double E=0,P=0,M=0,minP=1e30,minV=1e30,maxContact=0;
 for(int i=0;i<N;i++){E+=s.E[i];P+=s.p[i];M+=mass[i];minP=fmin(minP,q[i].p);minV=fmin(minV,s.x[i+1]-s.x[i]);if(contact)maxContact=fmax(maxContact,fabs(q[i].p-1));}
 printf("{\"case\":\"%s\",\"N\":%d,\"time\":%.17g,\"target\":%.17g,\"failed\":%d,\"steps\":%d,\"retries\":%d,\"mass_change\":%.17g,\"energy_balance_error\":%.17g,\"momentum_balance_error\":%.17g,\"minimum_pressure\":%.17g,\"minimum_volume\":%.17g,\"contact_pressure_error\":%.17g}\n",name,N,t,stop,failed,steps,retries,M-M0,(E-E0-eb)/fmax(1,fabs(E0)),(P-P0-pb)/fmax(1,fabs(P0)),minP,minV,maxContact);
 return failed;
}
