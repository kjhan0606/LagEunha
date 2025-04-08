void TwoLPTmain(SimParameters *);
void Zeldovichmain(SimParameters *);
void Set_Particle_Zeld(SimParameters *,  GridInfo *, GridInfo *, GridInfo *, enum mtype );
void Set_Particle_2LPT_Vel(SimParameters *, GridInfo *, GridInfo *, GridInfo *, enum mtype);
void Set_Particle_2LPT_Pos(SimParameters *, GridInfo *, GridInfo *, GridInfo *, enum mtype);
void PMSeedForce(SimParameters *, DenType *, DenType *, DenType *, DenType *, DenType *, DenType *, enum mtype);
void expectDMandSphNumParticle(SimParameters *);

#define FortranArray(mx,ny,arr,i,j,k) arr[i+mx*(j+ny*k)]

#define GlacialFind2LPT(simpar,type,fx,fy,fz,vr) do{\
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
	DEBUGPRINT("P%d has %s %g %g %g %g %g %g\n",MYID(simpar),#type,type##_BP(simpar)[0].x,type##_BP(simpar)[0].y,type##_BP(simpar)[0].z,\
			type##_BP(simpar)[0].vx,type##_BP(simpar)[0].vy,type##_BP(simpar)[0].vz);\
} while(0)

#define Initialize_Vel(simpar,type) do{\
	size_t ii;\
	for(ii=0;ii<type##_NP(simpar);ii++){\
		type##_BP(simpar)[ii].vx = type##_BP(simpar)[ii].vy = type##_BP(simpar)[ii].vz = 0;\
	}\
} while(0)

#define TwoLPT_GetVR(simpar,ffr,den1,den2,amp1,amp2,ixyz,type) do{\
	ffr = (GridInfo*)malloc(sizeof(GridInfo));\
	GridInfo *fr = (GridInfo*)malloc(sizeof(DenType)*densize + sizeof(GridInfo));\
	DenType *dfr = (DenType*)(fr+1);\
	ptrdiff_t mx = 2*(NX(simpar)/2+1);\
	size_t ii,ij,ik;\
	for(ik=0;ik<LOCAL_NZ(simpar);ik++) for(ij=0;ij<NY(simpar);ij++) for(ii=0;ii<mx;ii++){\
		dfr[ii+mx*(ij+NY(simpar)*ik)] = amp1*den1[ii+mx*(ij+NY(simpar)*ik)]+amp2*den2[ii+mx*(ij+NY(simpar)*ik)];\
	}\
	if(1){\
		double maxval = 0.;\
		for(ik=0;ik<local_nz;ik++) for(ij=0;ij<ny;ij++) for(ii=0;ii<nx;ii++){\
			maxval = MAX(maxval, dfr[ii+mx*(ij+ny*ik)]);\
		}\
		DEBUGPRINT("-P%d has fftw maxval of composite = %g ::: %g %g\n",MYID(simpar), maxval, amp1, amp2);\
	}\
	if(0){\
		void dumpphi(SimParameters *, float*);\
		dumpphi(simpar,dfr);\
	}\
	int nbuff = 0;\
	void selfdiffphi(SimParameters *, DenType *, enum dimension);\
	selfdiffphi(simpar, dfr,ixyz);\
	if(0){\
		void dumpden(SimParameters *, float*, ptrdiff_t);\
		dumpden(simpar,dfr,NX(simpar));\
	}\
	if(1){\
		double maxval = 0.;\
		for(ik=0;ik<local_nz;ik++) for(ij=0;ij<ny;ij++) for(ii=0;ii<nx;ii++){\
			maxval = MAX(maxval, dfr[ii+nx*(ij+ny*ik)]);\
		}\
		DEBUGPRINT("+P%d has fftw maxval of composite = %g ::: %g %g ::: %g -- %g\n",\
				MYID(simpar), maxval, amp1, amp2,SIM_LXMIN(simpar,type), SIM_LXMAX(simpar,type));\
	}\
	fr->ix=fr->iy = 0;\
	fr->jx=NX(simpar) - 1; fr->jy = NY(simpar) - 1;\
	fr->iz = LOCAL_Z_START(simpar); fr->jz = fr->iz + LOCAL_NZ(simpar) - 1;\
	ffr->nx = fr->nx = NX(simpar); ffr->ny = fr->ny = NY(simpar); ffr->nz = fr->nz = NZ(simpar);\
	ffr->ix = rint(SIM_LXMIN(simpar,type)) - nbuff;\
	ffr->iy = rint(SIM_LYMIN(simpar,type)) - nbuff;\
	ffr->iz = rint(SIM_LZMIN(simpar,type)) - nbuff;\
	ffr->jx = rint(SIM_LXMAX(simpar,type)) - 1 + nbuff;\
	ffr->jy = rint(SIM_LYMAX(simpar,type)) - 1 + nbuff;\
	ffr->jz = rint(SIM_LZMAX(simpar,type)) - 1 + nbuff;\
	void getingridnpix(GridInfo *);\
	getingridnpix(fr);\
	void *getoutgridnpix(GridInfo *, DenType );\
	ffr = getoutgridnpix(ffr,0.);\
	gmigrate(fr, ffr, DM_DDINFO(simpar), NDDINFO(simpar));\
	/*\
	maxval = 0.;\
	dfr = (DenType*)(ffr+1);\
	for(ik=0;ik<ffr->jz-ffr->iz+1;ik++) for(ij=0;ij<ffr->jy-ffr->iy+1;ij++) for(ii=0;ii<ffr->jx-ffr->ix+1;ii++){\
		maxval = MAX(maxval, dfr[ii+(ffr->jx-ffr->ix+1)*(ij+(ffr->jy-ffr->iy+1)*ik)]);\
	}\
	DEBUGPRINT("P%d has rms maxval of composite = %g ::: %g %g\n",MYID(simpar), maxval, amp1, amp2);\
	*/\
} while(0)


#define Zeld_GetVR(simpar,ffr,den,amp,type) do{\
	ffr = (GridInfo*)malloc(sizeof(GridInfo));\
	GridInfo *fr = (GridInfo*)malloc(sizeof(DenType)*densize + sizeof(GridInfo));\
	DenType *dfr = (DenType*)(fr + 1);\
	ptrdiff_t mx = 2*(NX(simpar)/2+1);\
	size_t ii,ij,ik;\
	for(ik=0;ik<LOCAL_NZ(simpar);ik++) for(ij=0;ij<NY(simpar);ij++) for(ii=0;ii<NX(simpar);ii++){\
		dfr[ii+NX(simpar)*(ij+NY(simpar)*ik)] = amp*den[ii+mx*(ij+NY(simpar)*ik)];\
	}\
	int nbuff = 4;\
	fr->ix=fr->iy = 0;\
	fr->jx=NX(simpar) - 1; fr->jy = NY(simpar) - 1;\
	fr->iz = LOCAL_Z_START(simpar); fr->jz = fr->iz + LOCAL_NZ(simpar) - 1;\
	ffr->nx = fr->nx = NX(simpar); ffr->ny = fr->ny = NY(simpar); ffr->nz = fr->nz = NZ(simpar);\
	ffr->ix = rint(SIM_LXMIN(simpar,type)) - nbuff;\
	ffr->iy = rint(SIM_LYMIN(simpar,type)) - nbuff;\
	ffr->iz = rint(SIM_LZMIN(simpar,type)) - nbuff;\
	ffr->jx = rint(SIM_LXMAX(simpar,type)) - 1 + nbuff;\
	ffr->jy = rint(SIM_LYMAX(simpar,type)) - 1 + nbuff;\
	ffr->jz = rint(SIM_LZMAX(simpar,type)) - 1 + nbuff;\
	void getingridnpix(GridInfo *);\
	getingridnpix(fr);\
	void *getoutgridnpix(GridInfo *, DenType );\
	ffr = getoutgridnpix(ffr,0.);\
	gmigrate(fr, ffr, DM_DDINFO(simpar), NDDINFO(simpar));\
} while(0)


typedef struct VelOnly{
	float vx,vy,vz;
} VelOnly;


#ifdef XYZDBL
#define TwoLPTPosAssign(simpar, TYPE,type) do{\
	size_t np = TYPE##_NP(simpar);\
	TYPE##_BP(simpar) = (type##particletype*)realloc(TYPE##_BP(simpar), sizeof(type##particletype)*np);\
	VelOnly *vel = (VelOnly*) TYPE##_BP(simpar);\
	type##particletype tmp,*ptl = TYPE##_BP(simpar);\
	for(i=np-1;i>=0;i--){\
		tmp.vx = vel[i].vx;\
		tmp.vy = vel[i].vy;\
		tmp.vz = vel[i].vz;\
		ptl[i] = tmp;\
	}\
	size_t mp = 0;\
	for(k=iz;k<=jz;k++){\
		size_t zi = k-iz;\
		for(j=iy;j<=jy;j++){\
			size_t yi = j - iy;\
			for(i=ix;i<=jx;i++){\
				size_t xi = i - ix;\
				size_t ipos = xi +mx*(yi+my*zi);\
				ptl[mp].x = ffx[ipos];\
				ptl[mp].y = ffy[ipos];\
				ptl[mp].z = ffz[ipos];\
				ptl[mp].u4if.indx = i + NX(simpar)*(j+NY(simpar)*k);\
				mp++;\
			}\
		}\
	}\
	PosType maxdist,tmaxdist,dist;\
	maxdist = 0;\
	for(i=0;i<mp;i++){\
		dist = sqrt(ptl[i].x*ptl[i].x + ptl[i].y*ptl[i].y + ptl[i].z*ptl[i].z);\
		maxdist = MAX(dist, maxdist);\
	}\
	MPI_Reduce( &maxdist, &tmaxdist, 1, MPI_POSTYPE, MPI_MAX, 0, MPI_COMM(simpar));\
	if(MYID(simpar)==0) printf("####################################\n\nThe maximum displace with 2LPT is %g\n\n####################################\n", tmaxdist);\
} while(0) 
#define ZeldPosAssign(type,simpar,vel) do{\
		vel = (VelOnly*)realloc(vel, sizeof(type)*np);\
		type tmp,*ptl = (type*)vel;\
		for(i=np-1;i>=0;i--){\
			tmp.vx = vel[i].vx;\
			tmp.vy = vel[i].vy;\
			tmp.vz = vel[i].vz;\
			ptl[i] = tmp;\
		}\
		size_t mp = 0;\
		for(k=iz;k<jz;k++){\
			for(j=iy;j<jy;j++){\
				for(i=ix;i<jx;i++){\
					ptl[mp].x = pfact*ptl[mp].vx;\
					ptl[mp].y = pfact*ptl[mp].vy;\
					ptl[mp].z = pfact*ptl[mp].vz;\
					ptl[mp].u4if.indx = i + NX(simpar)*(j+NY(simpar)*k);\
					mp++;\
				}\
			}\
		}\
		PosType maxdist,tmaxdist,dist;\
		maxdist = 0;\
		for(i=0;i<mp;i++){\
			dist = sqrt(ptl[i].x*ptl[i].x + ptl[i].y*ptl[i].y + ptl[i].z*ptl[i].z);\
			maxdist = MAX(dist, maxdist);\
		}\
		MPI_Reduce( &maxdist, &tmaxdist, 1, MPI_POSTYPE, MPI_MAX, 0, MPI_COMM(simpar));\
		if(MYID(simpar)==0) printf("####################################\n\nThe maximum displace with Zeld is %g\n\n####################################\n", tmaxdist);\
	} while(0) 
#else
#define TwoLPTPosAssign(simpar, TYPE,type) do{\
	PosType xboxsize = NX(simpar);\
	PosType yboxsize = NY(simpar);\
	PosType zboxsize = NZ(simpar);\
	size_t np = TYPE##_NP(simpar);\
	TYPE##_BP(simpar) = (type##particletype*)realloc(TYPE##_BP(simpar), sizeof(type##particletype)*np);\
	VelOnly *vel = (VelOnly*) TYPE##_BP(simpar);\
	type##particletype tmp,*ptl = TYPE##_BP(simpar);\
	for(i=np-1;i>=0;i--){\
		tmp.vx = vel[i].vx;\
		tmp.vy = vel[i].vy;\
		tmp.vz = vel[i].vz;\
		ptl[i] = tmp;\
	}\
	PosType maxdist,tmaxdist,dist; maxdist = 0;\
	size_t mp = 0;\
	for(k=iz;k<=jz;k++){\
		size_t zi = k-iz;\
		for(j=iy;j<=jy;j++){\
			size_t yi = j - iy;\
			for(i=ix;i<=jx;i++){\
				size_t xi = i - ix;\
				size_t ipos = xi +mx*(yi+my*zi);\
				ptl[mp].x = fmod(i+ffx[ipos] + xboxsize, zboxsize);\
				ptl[mp].y = fmod(j+ffy[ipos] + yboxsize, yboxsize);\
				ptl[mp].z = fmod(k+ffz[ipos] + zboxsize, zboxsize);\
				ptl[mp].u4if.indx = i + NX(simpar)*(j+NY(simpar)*k);\
				dist = sqrt(ffx[ipos]*ffx[ipos] + ffy[ipos]*ffy[ipos] + ffz[ipos]*ffz[ipos]);\
				maxdist = MAX(dist, maxdist);\
				mp++;\
			}\
		}\
	}\
	MPI_Reduce( &maxdist, &tmaxdist, 1, MPI_POSTYPE, MPI_MAX, 0, MPI_COMM(simpar));\
	if(MYID(simpar)==0) printf("####################################\n\nThe maximum displace with 2LPT is %g\n\n####################################\n", tmaxdist);\
} while(0) 
#define ZeldPosAssign(type,simpar,vel) do{\
		PosType xboxsize = NX(simpar);\
		PosType yboxsize = NY(simpar);\
		PosType zboxsize = NZ(simpar);\
		vel = (VelOnly*)realloc(vel, sizeof(type)*np);\
		type tmp,*ptl = (type*)vel;\
		for(i=np-1;i>=0;i--){\
			tmp.vx = vel[i].vx;\
			tmp.vy = vel[i].vy;\
			tmp.vz = vel[i].vz;\
			ptl[i] = tmp;\
		}\
		PosType maxdist,tmaxdist,dist;maxdist = 0;\
		size_t mp = 0;\
		for(k=iz;k<=jz;k++){\
			for(j=iy;j<=jy;j++){\
				for(i=ix;i<=jx;i++){\
					ptl[mp].x = fmod(i+pfact*ptl[mp].vx + xboxsize, xboxsize);\
					ptl[mp].y = fmod(j+pfact*ptl[mp].vy + yboxsize, yboxsize);\
					ptl[mp].z = fmod(k+pfact*ptl[mp].vz + zboxsize, zboxsize);\
					ptl[mp].u4if.indx = i + NX(simpar)*(j+NY(simpar)*k);\
					dist = pfact*sqrt(ptl[mp].vx*ptl[mp].vx + ptl[mp].vy*ptl[mp].vy + ptl[mp].vz*ptl[mp].vz);\
					maxdist = MAX(dist, maxdist);\
					mp++;\
				}\
			}\
		}\
		MPI_Reduce( &maxdist, &tmaxdist, 1, MPI_POSTYPE, MPI_MAX, 0, MPI_COMM(simpar));\
		if(MYID(simpar)==0) printf("####################################\n\nThe maximum displace with Zeld is %g\n\n####################################\n", tmaxdist);\
	} while(0) 

#endif
