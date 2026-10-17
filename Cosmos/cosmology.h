#define icopysign(a) (signbit(a)? -1:1)
void SetEnvGridInfo(SimParameters *);
void TreeAllParticleMigrate(SimParameters *);
void AllParticleMigrate(SimParameters *);
void readGrafic(SimParameters *);
void readGraficFile(SimParameters *);



#define DefineProtoTypeFunctions( type) \
    int(type##particletypexcompare)(GridInfo *,const void *, const void *);\
    int(type##particletypeycompare)(GridInfo *,const void *, const void *);\
    int(type##particletypezcompare)(GridInfo *,const void *, const void *);\
    int(type##particletypexpinner)(GridInfo *,const void *, const void *, const void *, PosType, int, Range *);\
    int(type##particletypeypinner)(GridInfo *,const void *, const void *, const void *, PosType, int, Range *);\
    int(type##particletypezpinner)(GridInfo *,const void *, const void *, const void *, PosType, int, Range *);\
    int(type##particletypemuladd)(GridInfo *,const void *, const void *, float, char, int);\
    int(type##particletypeInSideBox)(GridInfo *,const void *, BoxMinMax *, PosType *, SimBoxRange *, int, int,\
			enum PadPtl);\
    int(type##particletypedetermineEdgePtl)(GridInfo *,const void *, SimBoxRange *, PosType *);\
    void(type##particletypeshift2pos)(GridInfo *,const void *);\
    int(type##particletypeDDRule3D)(GridInfo *, const void *, ptrdiff_t, MPI_Comm)




#define CopyDDFuncDDInfoFromDM(simpar,TYPE, type) do{\
	int i,nddinfo = NDDINFO(simpar);\
	TYPE##_DDFUNC(simpar).xcompare = type##particletypexcompare;\
	TYPE##_DDFUNC(simpar).ycompare = type##particletypeycompare;\
	TYPE##_DDFUNC(simpar).zcompare = type##particletypezcompare;\
	TYPE##_DDFUNC(simpar).xpinner = type##particletypexpinner;\
	TYPE##_DDFUNC(simpar).ypinner = type##particletypeypinner;\
	TYPE##_DDFUNC(simpar).zpinner = type##particletypezpinner;\
	TYPE##_DDFUNC(simpar).muladd = type##particletypemuladd;\
	TYPE##_DDFUNC(simpar).insidebox = type##particletypeInSideBox;\
	TYPE##_DDFUNC(simpar).edgeptl = type##particletypedetermineEdgePtl;\
	TYPE##_DDFUNC(simpar).shift2pos = type##particletypeshift2pos;\
	TYPE##_DDFUNC(simpar).divdir = type##particletypeDDRule3D;\
	TYPE##_DDFUNC(simpar).xyzchip = DM_DDFUNC(simpar).xyzchip;\
	for(i=0;i<NDDINFO(simpar);i++){\
		TYPE##_DDINFO(simpar)[i].xcompare = TYPE##_DDFUNC(simpar).xcompare;\
		TYPE##_DDINFO(simpar)[i].ycompare = TYPE##_DDFUNC(simpar).ycompare;\
		TYPE##_DDINFO(simpar)[i].zcompare = TYPE##_DDFUNC(simpar).zcompare;\
		TYPE##_DDINFO(simpar)[i].xpinner = TYPE##_DDFUNC(simpar).xpinner;\
		TYPE##_DDINFO(simpar)[i].ypinner = TYPE##_DDFUNC(simpar).ypinner;\
		TYPE##_DDINFO(simpar)[i].zpinner = TYPE##_DDFUNC(simpar).zpinner;\
		TYPE##_DDINFO(simpar)[i].muladd = TYPE##_DDFUNC(simpar).muladd;\
		TYPE##_DDINFO(simpar)[i].insidebox = TYPE##_DDFUNC(simpar).insidebox;\
		TYPE##_DDINFO(simpar)[i].edgeptl = TYPE##_DDFUNC(simpar).edgeptl;\
		TYPE##_DDINFO(simpar)[i].shift2pos = TYPE##_DDFUNC(simpar).shift2pos;\
		TYPE##_DDINFO(simpar)[i].xyzchip = DM_DDINFO(simpar)[i].xyzchip;\
		TYPE##_DDINFO(simpar)[i].com = DM_DDINFO(simpar)[i].com;\
		TYPE##_DDINFO(simpar)[i].n_size = sizeof(type##particletype);\
		TYPE##_DDINFO(simpar)[i].nid = DM_DDINFO(simpar)[i].nid;\
		TYPE##_DDINFO(simpar)[i].myid = DM_DDINFO(simpar)[i].myid;\
		TYPE##_DDINFO(simpar)[i].nsubgroup = DM_DDINFO(simpar)[i].nsubgroup;\
		TYPE##_DDINFO(simpar)[i].idirection = DM_DDINFO(simpar)[i].idirection;\
		TYPE##_DDINFO(simpar)[i].subgroupsize = DM_DDINFO(simpar)[i].subgroupsize;\
		TYPE##_DDINFO(simpar)[i].subgroupid = DM_DDINFO(simpar)[i].subgroupid;\
		TYPE##_DDINFO(simpar)[i].sublocalgrid = DM_DDINFO(simpar)[i].sublocalgrid;\
		TYPE##_DDINFO(simpar)[i].memoffset.xyz.x = offsetof(type##particletype, x);\
		TYPE##_DDINFO(simpar)[i].memoffset.xyz.y = offsetof(type##particletype, y);\
		TYPE##_DDINFO(simpar)[i].memoffset.xyz.z = offsetof(type##particletype, z);\
		TYPE##_DDINFO(simpar)[i].npivot = DM_DDINFO(simpar)[i].npivot;\
		TYPE##_DDINFO(simpar)[i].lgroup = DM_DDINFO(simpar)[i].lgroup;\
		TYPE##_DDINFO(simpar)[i].pivot = (void*)malloc(TYPE##_DDINFO(simpar)[i].npivot*TYPE##_DDINFO(simpar)[i].n_size);\
		int j;\
		for(j=0;j<TYPE##_DDINFO(simpar)[i].npivot;j++){\
			((type##particletype*)(TYPE##_DDINFO(simpar)[i].pivot)+j)->x = \
				((dmparticletype*)(DM_DDINFO(simpar)[i].pivot)+j)->x;\
			((type##particletype*)(TYPE##_DDINFO(simpar)[i].pivot)+j)->y = \
				((dmparticletype*)(DM_DDINFO(simpar)[i].pivot)+j)->y;\
			((type##particletype*)(TYPE##_DDINFO(simpar)[i].pivot)+j)->z = \
				((dmparticletype*)(DM_DDINFO(simpar)[i].pivot)+j)->z;\
			((type##particletype*)(TYPE##_DDINFO(simpar)[i].pivot)+j)->u4if= \
				((dmparticletype*)(DM_DDINFO(simpar)[i].pivot)+j)->u4if;\
		}\
	}\
} while(0)



#define EunhaParticleFuncs(ptype) \
	GenerateMemberCompare(xcompare, ptype,x,PosType);\
	GenerateMemberCompare(ycompare, ptype,y,PosType);\
	GenerateMemberCompare(zcompare, ptype,z,PosType);\
	GenerateInsideRange(xpinner, ptype,x,PosType);\
	GenerateInsideRange(ypinner, ptype,y,PosType);\
	GenerateInsideRange(zpinner, ptype,z,PosType);\
	GenerateShift2Pos(Shift2Pos, ptype);\



#define InSideOrNot(type, a, lbox, width, simbox, mflag, pflag, def_ptl) int type##InSideBox(\
		GridInfo *gridinfo, const void *a, BoxMinMax *lbox, PosType *width, \
		SimBoxRange *simbox, int mflag , int pflag, enum PadPtl def_ptl){ \
	type *aa = (type*)a;\
	PosType aax,aay,aaz;\
	int i,j,k;\
	if(pflag == 13) {\
		i = j = k = 0;\
	}\
	else {\
		i = pflag%3-1;\
		j = (pflag%9)/3 -1;\
		k = (pflag/9)-1;\
	}\
	if(def_ptl == SimPtl){\
		aax = xPos(gridinfo, aa);\
		aay = yPos(gridinfo, aa);\
		aaz = zPos(gridinfo, aa);\
	}\
	else {\
		aax = XofBp(aa);\
		aay = YofBp(aa);\
		aaz = ZofBp(aa);\
	}\
	aax += i*(simbox->x).max;\
	aay += j*(simbox->y).max;\
	aaz += k*(simbox->z).max;\
	if( (aax >= (lbox->xmin - (*width)) && aax <(lbox->xmax + (*width))) &&\
			(aay >=(lbox->ymin - (*width)) && aay <(lbox->ymax + (*width))) &&\
			(aaz >=(lbox->zmin - (*width)) && aaz <(lbox->zmax + (*width)))) {\
		if(mflag) {\
			aa->x = aax;\
			aa->y = aay;\
			aa->z = aaz;\
			CHANGEINDX(aa,0L);\
		}\
		return 1;\
	}\
	else return 0;\
}


#define determineEdgePtl(type,a, simbox, width) \
	int type##determineEdgePtl(GridInfo *gridinfo, const void *a, SimBoxRange *simbox,\
		PosType *width){\
	type *aa = (type *)a ;\
	PosType aax = xPos(gridinfo, aa);\
	PosType aay = yPos(gridinfo, aa);\
	PosType aaz = zPos(gridinfo, aa);\
	if(aax<=*width) return 1;\
	else if ( (simbox->x).max - aax < *width) return 1;\
	if(aay<=*width) return 1;\
	else if ( (simbox->y).max - aay < *width) return 1;\
	if(aaz<=*width) return 1;\
	else if ( (simbox->z).max - aaz < *width) return 1;\
	return 0;\
}

#define Dmuladd(type, a, b, fact, xyz, ifalg) int type##muladd(GridInfo *gridinfo, const void *a, const void *b, float fact,\
		char xyz, int iflag){\
	PosType aa;\
	float *av;\
	PosType bb;\
	switch (xyz) {\
		case 'x':\
			av = &(((type*)a)->x);\
			if(iflag<0) bb = ( (type*)b)->x;\
			else bb = xPos(gridinfo, (type*)b);\
			break;\
		case 'y': \
			av = &(((type*)a)->y);\
			if(iflag<0) bb = ( (type*)b)->y;\
			else bb = yPos(gridinfo, (type*)b);\
			break;\
		default: \
			av = &(((type*)a)->z);\
			if(iflag<0) bb = ( (type*)b)->z;\
			else bb = zPos(gridinfo, (type*)b);\
	}\
	aa = (PosType)(*av);\
	if(iflag){\
		if(isnan(aa) || isnan(bb)) return 0;\
		else {\
			(*av) += bb*fact;\
			return 1;\
		}\
	}\
	else {\
		Range *bbb = (Range*)b;\
		(*av) += bbb->max*fact;\
	}\
	return 1;\
}





#define DDRule3D(type, a, nmem, com) \
	int type##DDRule3D(GridInfo *gridinfo,const void *a, ptrdiff_t nmem, MPI_Comm com){\
	ptrdiff_t _i;\
	double xmean,xstd;\
	double ymean,ystd;\
	double zmean,zstd;\
	double totalval;\
	ptrdiff_t totalnum;\
	type *aa = (type *)a;\
	xmean = xstd = 0;\
	ymean = ystd = 0;\
	zmean = zstd = 0;\
	for(_i=0;_i<nmem;_i++){\
		type *aai = aa + _i;\
		PosType aax = xPos(gridinfo, aai);\
		PosType aay = yPos(gridinfo, aai);\
		PosType aaz = zPos(gridinfo, aai);\
		xmean += aax;\
		ymean += aay;\
		zmean += aaz;\
		xstd += aax*aax;\
		ystd += aay*aay;\
		zstd += aaz*aaz;\
	}\
	MPI_Reduce(&xmean,&totalval,1, MPI_DOUBLE, MPI_SUM, 0, com);\
	MPI_Bcast(&totalval, 1, MPI_DOUBLE, 0, com);\
	\
	xmean = totalval;\
	MPI_Reduce(&ymean,&totalval,1, MPI_DOUBLE, MPI_SUM, 0, com);\
	MPI_Bcast(&totalval, 1, MPI_DOUBLE, 0, com);\
	ymean = totalval;\
\
	MPI_Reduce(&zmean,&totalval,1, MPI_DOUBLE, MPI_SUM, 0, com);\
	MPI_Bcast(&totalval, 1, MPI_DOUBLE, 0, com);\
	zmean = totalval;\
\
	MPI_Reduce(&xstd,&totalval,1, MPI_DOUBLE, MPI_SUM, 0, com);\
	MPI_Bcast(&totalval, 1, MPI_DOUBLE, 0, com);\
	xstd = totalval;\
\
	MPI_Reduce(&ystd,&totalval,1, MPI_DOUBLE, MPI_SUM, 0, com);\
	MPI_Bcast(&totalval, 1, MPI_DOUBLE, 0, com);\
	ystd = totalval;\
\
	MPI_Reduce(&zstd,&totalval,1, MPI_DOUBLE, MPI_SUM, 0, com);\
	MPI_Bcast(&totalval, 1, MPI_DOUBLE, 0, com);\
	zstd = totalval;\
\
	MPI_Reduce(&nmem,&totalnum,1, MPI_INT64_T, MPI_SUM, 0, com);\
	MPI_Bcast(&totalnum, 1, MPI_INT64_T, 0, com);\
	xmean = xmean/totalnum;\
	ymean = ymean/totalnum;\
	zmean = zmean/totalnum;\
	xstd = xstd/totalnum - xmean*xmean;\
	ystd = ystd/totalnum - ymean*ymean;\
	zstd = zstd/totalnum - zmean*zmean;\
	if(xstd >= ystd && xstd >= zstd) return 0;\
	else if(zstd >= ystd && zstd >= xstd) return 2;\
	else return 1;\
}


int RunCosmos(SimParameters *, int);



#define debugparticles(simpar) do { \
	PosType xmin,ymin,zmin,xmax,ymax,zmax; \
	xmin = ymin = zmin = 1.E20; \
	xmax = ymax = zmax = -1.E20; \
	ptrdiff_t i,j,k; \
	for(i=0;i<DM_NP(simpar);i++){ \
		xmin = MIN(xmin, XofP(simpar, DM_BP(simpar)+i)); \
		xmax = MAX(xmax, XofP(simpar, DM_BP(simpar)+i)); \
		ymin = MIN(ymin, YofP(simpar, DM_BP(simpar)+i)); \
		ymax = MAX(ymax, YofP(simpar, DM_BP(simpar)+i)); \
		zmin = MIN(zmin, ZofP(simpar, DM_BP(simpar)+i)); \
		zmax = MAX(zmax, ZofP(simpar, DM_BP(simpar)+i)); \
	} \
	DEBUGPRINT("P%d has %ld pregion %g %g : %g %g : %g %g\n",MYID(simpar),DM_NP(simpar),xmin,xmax, ymin,ymax,zmin,zmax); \
	DEBUGPRINT("P%d has %p region %g %g : %g %g : %g %g\n",MYID(simpar), DM_BP(simpar),SIM_LXMIN(simpar,dm), \
			SIM_LXMAX(simpar,dm), SIM_LYMIN(simpar,dm), SIM_LYMAX(simpar,dm), \
			SIM_LZMIN(simpar,dm), SIM_LZMAX(simpar,dm));\
	DEBUGPRINT("P%d has  p value %g %g %g \n",MYID(simpar),XofP(simpar,DM_BP(simpar)),\
			YofP(simpar, DM_BP(simpar)), ZofP(simpar, DM_BP(simpar))) ; \
}while(0)
