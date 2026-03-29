/* icc -o mkRtidal mkRtidal.c -lm */
#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>


float getgc(float c){
	float res = 1./(log(1+c) - c/(1+c));
	return res;
}

float getmM(float c,float dist_over_Rv, float rs_over_dist){

	float rs_over_Rv = rs_over_dist*dist_over_Rv;
	float gc = getgc(c);

	float s1 = dist_over_Rv - rs_over_Rv;
	float s2 = dist_over_Rv;
	float mnratio = gc*(rs_over_Rv)*(rs_over_Rv)*(-c/s1/(1.+c*s1) + log(1.+c*s1)/s1/s1);
	float mnratio2 = gc*(rs_over_Rv)*(rs_over_Rv)*(-c/s2/(1.+c*s2) + log(1.+c*s2)/s2/s2);

	mnratio -= mnratio2;

	return mnratio;
}

#define Nrs 128
#define Nm Nrs
#define Nc  32
#define Nd  128


float tidal_den[Nc*Nd*Nm];

float interpole(float *mM, float *rsdist, int np, float massratio){
	int i,j,k;

	if(massratio < mM[0]) return rsdist[0];
	else if(massratio >= mM[np-1]) return rsdist[np-1];
	else {
		for(i=0;i<np-1;i++){
			if(massratio >= mM[i] && massratio < mM[i+1]){
				float res = (rsdist[i+1] - rsdist[i])/(mM[i+1]-mM[i])*(massratio-mM[i]) + rsdist[i];
				return res;
			}
		}
	}
	return 0;
}
/*
float cinit = 1;
float cfinal = 32.;
float dRvinit = -3.;
float dRvfinal = 4;
float rsinit = 0;
float rsfinal = 5;
*/

float nfw_rtidal(float mM, float dRv, float c){
	int i =  (log10(mM) + 4.)/4.*Nm;
	int j =  (log10(dRv) + 3.)/4.*Nd;
	int k =  (c - 1.)/32.*Nc;

	if(i < 0) i = 0;
	else if(i >= Nm) i = Nm-1;

	if(j < 0) j = 0;
	else if(j >= Nd) j = Nd-1;

	if(k < 0) k = 0;
	else if(k >= Nc) k = Nc-1;

	float res = tidal_den[i+Nm*(j+Nd*k)];
	return res;
}

void mkRtidal(void){
	int nrs, nc, nd,nm;

	int i,j,k;
	float c;

	nm = Nm; nrs = Nrs; nc = Nc, nd = Nd;



	for(k=0;k<Nc;k++)
	{
		c = 1. + (32.*k)/Nc;
		float gc = getgc(c);
		for(j=0;j<Nd;j++)
		{
			float dist_over_Rv;
			dist_over_Rv = -3. + 4*j/(float)Nd;
			dist_over_Rv = pow(10., dist_over_Rv);

			float mM[Nrs], rsdist[Nrs];
			int np;

			for(i=0;i<Nrs;i++)
			{
				float rs_over_dist = (float)(i+1) /(float)Nrs * 5.;
				float rs_over_Rv = rs_over_dist * dist_over_Rv;
				float mratio = getmM(c, dist_over_Rv, rs_over_dist);
				if(mratio > 1) break;
				mM[i] = log10(mratio);
				rsdist[i] = log10(rs_over_dist);

//				printf("c= %g d/Rv= %g m/M = %g ::: rs/d= %g \n", c, dist_over_Rv,  mratio,rs_over_dist);
			}
			np = i;

			for(i=0;i<Nm;i++){
				float mr = -4 + 4*(float)i/(float)Nm;
				tidal_den[i + Nm*(j+Nd*k)] = pow(10, interpole(mM, rsdist, np, mr));
			}
		}
	}
	/*
	FILE *wp = fopen("Modelled_Tidal_Radius.dat","w");
	fwrite(&nm, sizeof(int), 1, wp);
	fwrite(&nd, sizeof(int), 1, wp);
	fwrite(&nc, sizeof(int), 1, wp);
	fwrite(tidal_den, sizeof(float), nm*nd*nc,wp);
	fclose(wp);
	*/
}
