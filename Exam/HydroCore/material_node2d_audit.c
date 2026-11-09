#include <stdio.h>
#include <math.h>
#include "../gfs_material_node2d.h"
static int pair(double rl,double pl,double ul,double rr,double pr,double ur){
 GfsFvPrimitive l={rl,pl,{ul,0,0}},r={rr,pr,{ur,0,0}};GfsExactStar exact;
 if(!gfs_exact_star(l,r,ul,ur,1.4,&exact)||exact.vacuum)return 1;
 GfsMaterialHalfEdge e[4]={{{1,0},.5,rl,pl,{ul,0}},{{-1,0},.5,rr,pr,{ur,0}},{{0,1},.5,1,1,{0,0}},{{0,-1},.5,1,1,{0,0}}};
 double v[2]={.5*(ul+ur),0},res;int it;
 if(!gfs_material_node_solve(4,e,1.4,1,1,v,&it,&res))return 2;
 double p,d;if(!gfs_exact_piston(rl,pl,1.4,ul-v[0],&p,&d))return 3;
 double err=fmax(fabs(v[0]-exact.u)/fmax(1.,fabs(exact.u)),fabs(p-exact.p)/fmax(1.,fabs(exact.p)));
 printf("1D reduction: velocity %.17g pressure %.17g relative difference %.3g iterations %d\n",v[0],p,err,it);
 return err>1e-9;
}
int main(void){
 if(pair(1,1,0,.125,.1,0)||pair(1,1000,0,1,.01,0)||pair(1,.4,-2,1,.4,2)||pair(5.99924,460.894,19.5975,5.99242,46.095,-6.19633))return 1;
 /* Four incident quadrants; two outward half edges per cell. */
 GfsMaterialHalfEdge e[8];
 for(int c=0;c<4;c++)for(int a=0;a<2;a++){
  int k=2*c+a;e[k].n[0]=a==0?((c&1)?-1:1):0;e[k].n[1]=a==1?((c&2)?-1:1):0;
  e[k].length=.5;e[k].rho=1+.2*c;e[k].p=1+.7*c;e[k].u[0]=.1*c;e[k].u[1]=-.2*c;
 }
 double v[2]={0,0},res;int it;if(!gfs_material_node_solve(8,e,1.4,1,1,v,&it,&res))return 2;
 double sum[2]={0,0},power=0,heat=0;
 for(int k=0;k<8;k++){double f[2],de,dh;if(!gfs_material_halfedge_rate(e+k,1.4,v,f,&de,&dh)||dh < -1e-14)return 3;sum[0]+=f[0];sum[1]+=f[1];power+=de;heat+=dh;}
 printf("2D pressure-driven node: velocity=(%.17g,%.17g) force residual %.3g power residual %.3g heat %.17g\n",v[0],v[1],hypot(sum[0],sum[1]),fabs(power),heat);
 if(hypot(sum[0],sum[1])>1e-10||fabs(power)>1e-10)return 4;
 double boost[2]={5,-3},vb[2]={5,-3};for(int k=0;k<8;k++)for(int a=0;a<2;a++)e[k].u[a]+=boost[a];
 if(!gfs_material_node_solve(8,e,1.4,1,1,vb,&it,&res))return 5;
 if(fabs(vb[0]-boost[0]-v[0])>1e-12||fabs(vb[1]-boost[1]-v[1])>1e-12)return 6;
 puts("PASS exact 1D reduction, 2D node momentum/work balance, nonnegative dissipative heat, Galilean shift");return 0;
}
