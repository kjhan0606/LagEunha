/*                          */
#define p8overpi (2.546479089)
#define p16overpi (5.0929581789)
#define p48overpi (15.2788745)
#define m48overpi (-15.2788745)
#define m24overpi (-7.639437268)

#define mumin 0.59212
#define mumax 1.23009

enum time {ONCE=1, NEVER=0};


/*
 * float cputime0[10],cputime1[10];
 * float gettime();
 * */
#define min(a,b) ((a)<(b) ? (a):(b))
#define max(a,b) ((a)>(b) ? (a):(b))

/*
#define CellWidth 4.L 
#define invCellWidth 0.25L 
*/

#define rsphere 4.

/***********************/
#define nstride 1
#define twonstridep1 3L
/***********************/


#define kB 1.380650424E-16L
/*
#define mH 1.66053878283E-24L
*/

#define mH 1.673575E-24L
#define kBovermH 0.83144726174176E8L
#define mHoverkB 1.20272210399148E-8L
#define Mpccgs 3.085677488036136E24L
#define P_G 6.673E-8L
#define not4 3.9715 /* mass ratio of He/H */
#define H0 100.E5L 
#define Myr 3.153600E13L
#define onesolarmass 1.989E33L

#define LinkedTreeSphType(a) ((treesphparticletype*)((a)->bp))
#define LinkedTreeBPSphType(a) ((treebpsphparticletype*)((a)->bp))
#define LinkedTreeStarType(a) ((treestarparticletype*)((a)->bp))
#define LinkedTreeDmType(a) ((treedmparticletype*)((a)->bp))

#define GetTemp(p) (p->Temp)
#define GetCs (p) sqrtf(GetTemp(p)/GetMu(p)*kBovermH);




/*
#define StarCellLinkedList(BasicCell,NBasicCell,ptype,localzstart) do{\
    long i,ix,iy,iz,mpixel;\
    long mxm1,mym1,mzm1;\
	linkedlisttype *p,*tp;\
    for(i=0;i<mx*my*mz;i++){\
        BasicCell[i] = NULL;\
        NBasicCell[i] = 0;\
    }\
    mxm1 = mx - 1;\
    mym1 = my - 1;\
    mzm1 = mz - 1;\
	if(simpar.ptype.np>0){\
		for(i=0;i<simpar.ptype.np;i++){\
			tp = (linkedlisttype*) (simpar.##ptype##.u.tbp+i);\
			ix = (long)(XofP(tp)*invCellWidth);\
			iy = (long)(YofP(tp)*invCellWidth);\
			iz = (long)((ZofP(tp)-localzstart)*invCellWidth);\
			ix = min(ix,mxm1); iy = min(iy,mym1); iz = min(iz,mzm1);\
			mpixel = ix+mx*(iy+my*iz);\
			NBasicCell[mpixel] ++;\
			p = BasicCell[mpixel];\
			BasicCell[mpixel] = tp;\
			tp->next = p;\
		}\
	}\
} while(0)
*/


#define HydroCellLinkedList(BasicCell,NBasicCell,ptype,localzstart) do{\
    long i,ix,iy,iz,mpixel;\
    long mxm1,mym1,mzm1;\
	linkedlisttype *p,*tp;\
    for(i=0;i<mx*my*mz;i++){\
        BasicCell[i].link = NULL;\
        BasicCell[i].nmem = 0;\
        BasicCell[i].sphcalcul = BasicCell[i].starcalcul = BasicCell[i].agncalcul = BasicCell[i].dmcalcul = 0;\
    }\
    mxm1 = mx - 1;\
    mym1 = my - 1;\
    mzm1 = mz - 1;\
	if(simpar.ptype.np>0){\
		for(i=0;i<simpar.ptype.np;i++){\
			tp = (linkedlisttype*) (simpar.ptype.u.tbp+i);\
			ix = (long)(XofP(tp)*invCellWidth);\
			iy = (long)(YofP(tp)*invCellWidth);\
			iz = (long)((ZofP(tp)-localzstart)*invCellWidth);\
			ix = min(ix,mxm1); iy = min(iy,mym1); iz = min(iz,mzm1);\
			mpixel = ix+mx*(iy+my*iz);\
			NBasicCell[mpixel] ++;\
			if(GetSphTsubPower(tp) > calflag[mpixel]) calflag[mpixel]=GetSphTsubPower(tp);\
			if(GetTsubPower(tp) > calflag[mpixel]) calflag[mpixel]=GetTsubPower(tp);\
			p = BasicCell[mpixel];\
			BasicCell[mpixel] = tp;\
			tp->next = p;\
		}\
	}\
	for(i=0;i<mx*my*mz;i++){\
		if(IsNowOneOfTwoSteps(nowTsubdiv,calflag[i],calflag[i],maxTsubpower)) {\
			calflag[i] = 1;\
		}\
		else {\
			calflag[i] = 0;\
		}\
	}\
} while(0)

#define DetermineMxMyMz(nx,ny,zheight){\
	mx = (long)ceil((double)nx/(double)CellWidth);\
	my = (long)ceil((double)ny/(double)CellWidth);\
	mz = (long)ceil((double)zheight/(double)CellWidth);\
}



float getmu(SimParameters *, float,float);
float GetMu(treesphparticletype*);


#define InitAccel_SphTsub(type){\
	long i;\
	for(i=0;i<simpar.type.np;i++){\
		(simpar.type.u.tbp+i)->ax=0;\
		(simpar.type.u.tbp+i)->ay=0;\
		(simpar.type.u.tbp+i)->az=0;\
		SetSphTsubPower(simpar.type.u.tbp+i,0);\
		SetSphNextTsubPower(simpar.type.u.tbp+i,0);\
	}\
}
#define InitAccel(type){\
	long i;\
	for(i=0;i<simpar.type.np;i++){\
		(simpar.type.u.tbp+i)->ax=0;\
		(simpar.type.u.tbp+i)->ay=0;\
		(simpar.type.u.tbp+i)->az=0;\
		SetSphNextTsubPower(simpar.type.u.tbp+i,0);\
	}\
}
#define InitSphTsub(type){\
	long i;\
	for(i=0;i<simpar.type.np;i++){\
		SetSphTsubPower(simpar.type.u.tbp+i,0);\
		SetSphNextTsubPower(simpar.type.u.tbp+i,0);\
	}\
}


#define FreeBoundaryArray(a){\
	long i;\
	for(i=0;i<twonstridep1;i++) Free(a[i]);\
}

#define TimeStep_limiter(imax,bpneighbor,Num_neighbor) do {\
	int inew = (imax-SPH_TIMESTEPLIMITER(simpar)); /* This is the power */\
	for(i=0;i<Num_neighbor;i++){\
		treesphparticletype *pj =  LinkedTreeSphType(bpneighbor[i]);\
		if(inew > GetSphTsubPower(pj) && inew > GetSphNextTsubPower(pj)){\
			SetSphNextTsubPower(pj,inew);\
		}\
	}\
} while(0)

#define Fixed_TimeStep_limiter(simpar, imax,ip,pi) do {\
	int inew = (imax-SPH_TIMESTEPLIMITER(simpar)); /* This is the power */\
	for(j=ip;neighborlist[j].pi == pi && j<Nneighborlist;j++){\
		long long nj = neighborlist[j].njnow;\
		treesphparticletype *pj;\
		if(neighborlist[j].pid == MYID(simpar)) pj = (treesphparticletype*) neighborlist[j].pj;\
		else DEBUGPRINT2("Something wrong here\n");\
		if(inew > GetSphTsubPower(pj) && inew > GetSphNextTsubPower(pj)){\
			SetSphNextTsubPower(pj,inew);\
		}\
	}\
} while(0)






float GetT(SimParameters *, treesphparticletype *);
float getT(SimParameters *, treesphparticletype *);
float GetTempGiveMu(SimParameters *, treesphparticletype *);

float rho_w(float , float );
float drhodr_w(float ,float );
void DestroyHydroLinkedCell(SimParameters *);

void flagboundaryghostparticles(SimParameters *);
void Init_dAs(SimParameters *);
size_t SPHDump2Pos(SimParameters *, size_t ,size_t ,size_t , particle *);
