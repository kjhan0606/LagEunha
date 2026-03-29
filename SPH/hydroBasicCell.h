
/*
#define CellWidth 4.L
*/


#define LinkParticles(simpar, TYPE, type,BasicCell,mx,my,xmin,ymin,zmin, CellWidth) do{\
	tree##type##particletype *bp = TYPE##_TBP(simpar);\
	for(i=0;i<TYPE##_NP(simpar);i++){\
		ptrdiff_t ix,iy,iz;\
		ix = (XofP(simpar,bp+i)-xmin)/CellWidth;\
		iy = (YofP(simpar,bp+i)-ymin)/CellWidth;\
		iz = (ZofP(simpar,bp+i)-zmin)/CellWidth;\
		ptrdiff_t ipos = ix+mx*(iy+my*iz);\
		if(GetSphTsubPower(bp+i) > BasicCell[ipos].calflag)\
			BasicCell[ipos].calflag = GetSphTsubPower(bp+i);\
		if(GetTsubPower(bp+i) > BasicCell[ipos].calflag)\
			BasicCell[ipos].calflag = GetTsubPower(bp+i);\
		linkedlisttype *tmp = BasicCell[ipos].link;\
		BasicCell[ipos].link = (linkedlisttype*)(bp+i);\
		BasicCell[ipos].nmem ++;\
		(bp+i)->next = tmp;\
	}\
	bp = TYPE##_TBPP(simpar);\
	for(i=0;i<TYPE##_NPP(simpar);i++){\
		ptrdiff_t ix,iy,iz;\
		ix = (XofP(simpar,bp+i)-xmin)/CellWidth;\
		iy = (YofP(simpar,bp+i)-ymin)/CellWidth;\
		iz = (ZofP(simpar,bp+i)-zmin)/CellWidth;\
		ptrdiff_t ipos = ix+mx*(iy+my*iz);\
		linkedlisttype *tmp = BasicCell[ipos].link;\
		BasicCell[ipos].link = (linkedlisttype*)(bp+i);\
		BasicCell[ipos].nmem ++;\
		(bp+i)->next = tmp;\
	}\
}while(0)


#define PaddingHydroParticles(simpar) do{\
	SimBoxRange box;\
	box.x.min= box.y.min= box.z.min= 0;\
	box.x.max= NX(simpar); box.y.max= NY(simpar); box.z.max= NZ(simpar);\
	float width = SPH_INTERACTIONSPHERE(simpar);\
	ppadding(SPH_TBP(simpar), SPH_NP(simpar), (void**)(&SPH_TBPP(simpar)), &SPH_NP(simpar), SPH_DDINFO(simpar), NDDINFO(simpar), \
			box, width, &GRIDINFO(simpar), 3);\
	ppadding(STAR_TBP(simpar), STAR_NP(simpar), (void**)(&STAR_TBPP(simpar)), &STAR_NP(simpar), STAR_DDINFO(simpar), NDDINFO(simpar), \
			box, width, &GRIDINFO(simpar), 3);\
	ppadding(AGN_TBP(simpar), AGN_NP(simpar), (void**)(&AGN_TBPP(simpar)), &AGN_NP(simpar), AGN_DDINFO(simpar), NDDINFO(simpar), \
			box, width, &GRIDINFO(simpar), 3);\
}while(0)

#define ComPaddingHydroParticles(simpar) do{\
	SimBoxRange box;\
	box.x.min= box.y.min= box.z.min= 0;\
	box.x.max= NX(simpar); box.y.max= NY(simpar); box.z.max= NZ(simpar);\
	float width = SPH_INTERACTIONSPHERE(simpar);\
	size_t i;\
	for(i=0;i<SPH_NP(simpar);i++) SPH_TBP(simpar)[i].bp = (SPH_TBP(simpar)+i);\
	ppadding(SPH_TBP(simpar), SPH_NP(simpar), (void**)(&SPH_TBPP(simpar)), &SPH_NP(simpar), SPH_DDINFO(simpar), NDDINFO(simpar), \
			box, width, &GRIDINFO(simpar), 3);\
	ppadding(STAR_TBP(simpar), STAR_NP(simpar), (void**)(&STAR_TBPP(simpar)), &STAR_NP(simpar), STAR_DDINFO(simpar), NDDINFO(simpar), \
			box, width, &GRIDINFO(simpar), 3);\
	ppadding(AGN_TBP(simpar), AGN_NP(simpar), (void**)(&AGN_TBPP(simpar)), &AGN_NP(simpar), AGN_DDINFO(simpar), NDDINFO(simpar), \
			box, width, &GRIDINFO(simpar), 3);\
}while(0)



void HydroBasicLinkedList(SimParameters *);
size_t SPHBuildLinkedList(SimParameters *simpar);
