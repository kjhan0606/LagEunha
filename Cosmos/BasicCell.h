#define CellWidth 4.L

void BuildLinkedList(SimParameters *,ptrdiff_t, ptrdiff_t, ptrdiff_t, PosType, PosType, PosType, PosType);

#define LinkParticles(simpar, TYPE, type,BasicCell,mx,my,xmin,ymin,zmin) do{\
	tree##type##particletype *bp = TYPE##_TBP(simpar);\
	ptrdiff_t _i;\
	for(_i=0;_i<TYPE##_NP(simpar);_i++){\
		ptrdiff_t ix,iy,iz;\
		ix = (XofP(simpar,bp+_i)-xmin)/CellWidth;\
		iy = (YofP(simpar,bp+_i)-ymin)/CellWidth;\
		iz = (ZofP(simpar,bp+_i)-zmin)/CellWidth;\
		ptrdiff_t ipos = ix+mx*(iy+my*iz);\
		linkedlisttype *tmp = BasicCell[ipos].link;\
		BasicCell[ipos].link = (linkedlisttype*)(bp+_i);\
		BasicCell[ipos].nmem ++;\
		(bp+_i)->next = tmp;\
	}\
	bp = TYPE##_TBPP(simpar);\
	for(_i=0;_i<TYPE##_NPP(simpar);_i++){\
		ptrdiff_t ix,iy,iz;\
		ix = (bp[_i].x-xmin)/CellWidth;\
		iy = (bp[_i].y-ymin)/CellWidth;\
		iz = (bp[_i].z-zmin)/CellWidth;\
		ptrdiff_t ipos = ix+mx*(iy+my*iz);\
		linkedlisttype *tmp = BasicCell[ipos].link;\
		BasicCell[ipos].link = (linkedlisttype*)(bp+_i);\
		BasicCell[ipos].nmem ++;\
		(bp+_i)->next = tmp;\
	}\
}while(0)

#define OMP_LinkParticles(simpar, TYPE, type,BasicCell,mx,my,xmin,ymin,zmin,izs,izf) do{\
	tree##type##particletype *bp = TYPE##_TBP(simpar);\
	ptrdiff_t _i;\
	for(_i=0;_i<TYPE##_NP(simpar);_i++){\
		ptrdiff_t ix,iy,iz;\
		iz = (ZofP(simpar,bp+_i)-zmin)/CellWidth;\
		if(iz >=izs && iz < izf){\
			ix = (XofP(simpar,bp+_i)-xmin)/CellWidth;\
			iy = (YofP(simpar,bp+_i)-ymin)/CellWidth;\
			ptrdiff_t ipos = ix+mx*(iy+my*iz);\
			linkedlisttype *tmp = BasicCell[ipos].link;\
			BasicCell[ipos].link = (linkedlisttype*)(bp+_i);\
			BasicCell[ipos].nmem ++;\
			(bp+_i)->next = tmp;\
		}\
	}\
	bp = TYPE##_TBPP(simpar);\
	for(_i=0;_i<TYPE##_NPP(simpar);_i++){\
		ptrdiff_t ix,iy,iz;\
		iz = (bp[_i].z-zmin)/CellWidth;\
		if(iz>=izs && iz < izf){\
			ix = (bp[_i].x-xmin)/CellWidth;\
			iy = (bp[_i].y-ymin)/CellWidth;\
			ptrdiff_t ipos = ix+mx*(iy+my*iz);\
			linkedlisttype *tmp = BasicCell[ipos].link;\
			BasicCell[ipos].link = (linkedlisttype*)(bp+_i);\
			BasicCell[ipos].nmem ++;\
			(bp+_i)->next = tmp;\
		}\
	}\
}while(0)

