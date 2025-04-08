#define icopysign(a) (signbit(a)? -1:1)

#define FortranArray(mx,ny,arr,i,j,k) arr[i+mx*(j+ny*k)]

#define Find2LPT(simpar,type,fx,fy,fz,vr) do{\
	PosType xstart,ystart,zstart;\
	size_t i,mx = 2*(NX(simpar)/2+1);\
	size_t nx = NX(simpar);\
	size_t ny = NY(simpar);\
	size_t nz = NZ(simpar);\
	xstart = SIM_LXMIN(simpar,dm);\
	ystart = SIM_LYMIN(simpar,dm);\
	zstart = SIM_LZMIN(simpar,dm);\
	PosType rnx = nx;\
	PosType rny = ny;\
	PosType rnz = nz;\
	for(i=0;i<type##_NP(simpar);i++){\
		PosType xp,yp,zp,xmin,ymin,zmin;\
		int nearx,neary,nearz,isign,jsign,ksign;\
		int i1,j1,k1,i2,j2,k2,i3,j3,k3;\
		PosType xd1,yd1,zd1;\
		PosType wx1,wy1,wz1,wx2,wy2,wz2,wx3,wy3,wz3;\
		PosType wt1,wt2,wt3,wt4,wt5,wt6,wt7,wt8,wt9,wt10;\
		PosType wt11,wt12,wt13,wt14,wt15,wt16,wt17,wt18,wt19,wt20;\
		PosType wt21,wt22,wt23,wt24,wt25,wt26,wt27;\
		xp = XofP(simpar,type##_BP(simpar)+i)-xstart;\
		yp = YofP(simpar,type##_BP(simpar)+i)-ystart;\
		zp = ZofP(simpar,type##_BP(simpar)+i)-zstart;\
		nearx = rintf(xp);\
		neary = rintf(yp);\
		nearz = rintf(zp);\
		xmin = xp-nearx;\
		ymin = yp-neary;\
		zmin = zp-nearz;\
		isign = icopysign(xmin);\
		jsign = icopysign(ymin);\
		ksign = icopysign(zmin);\
		i1 = nearx;\
		j1 = neary;\
		k1 = nearz;\
		i2 = nearx + isign;\
		j2 = neary + jsign;\
		k2 = nearz + ksign;\
		i3 = nearx-isign;\
		j3 = neary-jsign;\
		k3 = nearz-ksign;\
		xd1 = fabs(xmin);\
		yd1 = fabs(ymin);\
		zd1 = fabs(zmin);\
		wx1 = 0.75-xd1*xd1;\
		wy1 = 0.75-yd1*yd1;\
		wz1 = 0.75-zd1*zd1;\
		wx2 = 0.5*(0.5+xd1)*(0.5+xd1);\
		wy2 = 0.5*(0.5+yd1)*(0.5+yd1);\
		wz2 = 0.5*(0.5+zd1)*(0.5+zd1);\
		wx3 = 0.5*(0.5-xd1)*(0.5-xd1);\
		wy3 = 0.5*(0.5-yd1)*(0.5-yd1);\
		wz3 = 0.5*(0.5-zd1)*(0.5-zd1);\
		wt1 = wx1*wy1*wz1;\
		wt2 = wx2*wy1*wz1;\
		wt3 = wx3*wy1*wz1;\
		wt4 = wx1*wy2*wz1;\
		wt5 = wx2*wy2*wz1;\
		wt6 = wx3*wy2*wz1;\
		wt7 = wx1*wy3*wz1;\
		wt8 = wx2*wy3*wz1;\
		wt9 = wx3*wy3*wz1;\
		wt10 = wx1*wy1*wz2;\
		wt11 = wx2*wy1*wz2;\
		wt12 = wx3*wy1*wz2;\
		wt13 = wx1*wy2*wz2;\
		wt14 = wx2*wy2*wz2;\
		wt15 = wx3*wy2*wz2;\
		wt16 = wx1*wy3*wz2;\
		wt17 = wx2*wy3*wz2;\
		wt18 = wx3*wy3*wz2;\
		wt19 = wx1*wy1*wz3;\
		wt20 = wx2*wy1*wz3;\
		wt21 = wx3*wy1*wz3;\
		wt22 = wx1*wy2*wz3;\
		wt23 = wx2*wy2*wz3;\
		wt24 = wx3*wy2*wz3;\
		wt25 = wx1*wy3*wz3;\
		wt26 = wx2*wy3*wz3;\
		wt27 = wx3*wy3*wz3;\
		type##_BP(simpar)[i].vr##x += (wt1 *FortranArray(mx,ny,fx,i1,j1,k1)+wt2 *FortranArray(mx,ny,fx,i2,j1,k1)+wt3 *FortranArray(mx,ny,fx,i3,j1,k1)+\
									wt4 *FortranArray(mx,ny,fx,i1,j2,k1)+wt5 *FortranArray(mx,ny,fx,i2,j2,k1)+wt6 *FortranArray(mx,ny,fx,i3,j2,k1)+\
									wt7 *FortranArray(mx,ny,fx,i1,j3,k1)+wt8 *FortranArray(mx,ny,fx,i2,j3,k1)+wt9 *FortranArray(mx,ny,fx,i3,j3,k1)+\
									wt10*FortranArray(mx,ny,fx,i1,j1,k2)+wt11*FortranArray(mx,ny,fx,i2,j1,k2)+wt12*FortranArray(mx,ny,fx,i3,j1,k2)+\
									wt13*FortranArray(mx,ny,fx,i1,j2,k2)+wt14*FortranArray(mx,ny,fx,i2,j2,k2)+wt15*FortranArray(mx,ny,fx,i3,j2,k2)+\
									wt16*FortranArray(mx,ny,fx,i1,j3,k2)+wt17*FortranArray(mx,ny,fx,i2,j3,k2)+wt18*FortranArray(mx,ny,fx,i3,j3,k2)+\
									wt19*FortranArray(mx,ny,fx,i1,j1,k3)+wt20*FortranArray(mx,ny,fx,i2,j1,k3)+wt21*FortranArray(mx,ny,fx,i3,j1,k3)+\
									wt22*FortranArray(mx,ny,fx,i1,j2,k3)+wt23*FortranArray(mx,ny,fx,i2,j2,k3)+wt24*FortranArray(mx,ny,fx,i3,j2,k3)+\
									wt25*FortranArray(mx,ny,fx,i1,j3,k3)+wt26*FortranArray(mx,ny,fx,i2,j3,k3)+wt27*FortranArray(mx,ny,fx,i3,j3,k3));\
		type##_BP(simpar)[i].vr##y += (wt1 *FortranArray(mx,ny,fy,i1,j1,k1)+wt2 *FortranArray(mx,ny,fy,i2,j1,k1)+wt3 *FortranArray(mx,ny,fy,i3,j1,k1)+\
									wt4 *FortranArray(mx,ny,fy,i1,j2,k1)+wt5 *FortranArray(mx,ny,fy,i2,j2,k1)+wt6 *FortranArray(mx,ny,fy,i3,j2,k1)+\
									wt7 *FortranArray(mx,ny,fy,i1,j3,k1)+wt8 *FortranArray(mx,ny,fy,i2,j3,k1)+wt9 *FortranArray(mx,ny,fy,i3,j3,k1)+\
									wt10*FortranArray(mx,ny,fy,i1,j1,k2)+wt11*FortranArray(mx,ny,fy,i2,j1,k2)+wt12*FortranArray(mx,ny,fy,i3,j1,k2)+\
									wt13*FortranArray(mx,ny,fy,i1,j2,k2)+wt14*FortranArray(mx,ny,fy,i2,j2,k2)+wt15*FortranArray(mx,ny,fy,i3,j2,k2)+\
									wt16*FortranArray(mx,ny,fy,i1,j3,k2)+wt17*FortranArray(mx,ny,fy,i2,j3,k2)+wt18*FortranArray(mx,ny,fy,i3,j3,k2)+\
									wt19*FortranArray(mx,ny,fy,i1,j1,k3)+wt20*FortranArray(mx,ny,fy,i2,j1,k3)+wt21*FortranArray(mx,ny,fy,i3,j1,k3)+\
									wt22*FortranArray(mx,ny,fy,i1,j2,k3)+wt23*FortranArray(mx,ny,fy,i2,j2,k3)+wt24*FortranArray(mx,ny,fy,i3,j2,k3)+\
									wt25*FortranArray(mx,ny,fy,i1,j3,k3)+wt26*FortranArray(mx,ny,fy,i2,j3,k3)+wt27*FortranArray(mx,ny,fy,i3,j3,k3));\
		type##_BP(simpar)[i].vr##z += (wt1 *FortranArray(mx,ny,fz,i1,j1,k1)+wt2 *FortranArray(mx,ny,fz,i2,j1,k1)+wt3 *FortranArray(mx,ny,fz,i3,j1,k1)+\
									wt4 *FortranArray(mx,ny,fz,i1,j2,k1)+wt5 *FortranArray(mx,ny,fz,i2,j2,k1)+wt6 *FortranArray(mx,ny,fz,i3,j2,k1)+\
									wt7 *FortranArray(mx,ny,fz,i1,j3,k1)+wt8 *FortranArray(mx,ny,fz,i2,j3,k1)+wt9 *FortranArray(mx,ny,fz,i3,j3,k1)+\
									wt10*FortranArray(mx,ny,fz,i1,j1,k2)+wt11*FortranArray(mx,ny,fz,i2,j1,k2)+wt12*FortranArray(mx,ny,fz,i3,j1,k2)+\
									wt13*FortranArray(mx,ny,fz,i1,j2,k2)+wt14*FortranArray(mx,ny,fz,i2,j2,k2)+wt15*FortranArray(mx,ny,fz,i3,j2,k2)+\
									wt16*FortranArray(mx,ny,fz,i1,j3,k2)+wt17*FortranArray(mx,ny,fz,i2,j3,k2)+wt18*FortranArray(mx,ny,fz,i3,j3,k2)+\
									wt19*FortranArray(mx,ny,fz,i1,j1,k3)+wt20*FortranArray(mx,ny,fz,i2,j1,k3)+wt21*FortranArray(mx,ny,fz,i3,j1,k3)+\
									wt22*FortranArray(mx,ny,fz,i1,j2,k3)+wt23*FortranArray(mx,ny,fz,i2,j2,k3)+wt24*FortranArray(mx,ny,fz,i3,j2,k3)+\
									wt25*FortranArray(mx,ny,fz,i1,j3,k3)+wt26*FortranArray(mx,ny,fz,i2,j3,k3)+wt27*FortranArray(mx,ny,fz,i3,j3,k3));\
	}\
} while(0)

#define PeriodicPosition(simpar,type) do {\
	size_t i;\
	PosType rnx,rny,rnz; rnx = NX(simpar); rny = NY(simpar); rnz=NZ(simpar);\
	for(i=0;i<type##_NP(simpar);i++){\
		type##_BP(simpar)[i].x = fmod( type##_BP(simpar)[i].x+rnx,rnx);\
		type##_BP(simpar)[i].y = fmod( type##_BP(simpar)[i].y+rny,rny);\
		type##_BP(simpar)[i].z = fmod( type##_BP(simpar)[i].z+rnz,rnz);\
	}\
	printf("P%d has %s %g %g %g %g %g %g\n",MYID(simpar),#type,type##_BP(simpar)[0].x,type##_BP(simpar)[0].y,type##_BP(simpar)[0].z,\
			type##_BP(simpar)[0].vx,type##_BP(simpar)[0].vy,type##_BP(simpar)[0].vz);\
} while(0)

#define Initialize_Vel(simpar,type) do{\
	size_t ii;\
	for(ii=0;ii<type##_NP(simpar);ii++){\
		type##_BP(simpar)[ii].vx = type##_BP(simpar)[ii].vy = type##_BP(simpar)[ii].vz = 0;\
	}\
} while(0)
