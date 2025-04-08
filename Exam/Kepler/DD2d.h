#define icopysign(a) (signbit(a)? -1:1)
void SetEnvGridInfo(SimParameters *);
void TreeAllParticleMigrate(SimParameters *);
void AllParticleMigrate(SimParameters *);
void readGrafic(SimParameters *);
void readGraficFile(SimParameters *);



#define DefineProtoTypeFunctions2D(type) \
    int(type##particletypexcompare2D_kp)(GridInfo *,const void *, const void *);\
    int(type##particletypeycompare2D_kp)(GridInfo *,const void *, const void *);\
    int(type##particletypexpinner2D)(GridInfo *,const void *, const void *, const void *, PosType, int, Range *);\
    int(type##particletypeypinner2D)(GridInfo *,const void *, const void *, const void *, PosType, int, Range *);\
    int(type##particletypemuladd2D)(GridInfo *,const void *, const void *, float, char, int);\
    int(type##particletypeInSideBox2D)(GridInfo *,const void *, BoxMinMax *, PosType *, SimBoxRange *, int, int,\
			enum PadPtl);\
    int(type##particletypedetermineEdgePtl2D)(GridInfo *,const void *, SimBoxRange *, PosType *);\
    void(type##particletypeshift2pos2D)(GridInfo *,const void *);\
    int(type##particletypeDDRule2D)(GridInfo *, const void *, ptrdiff_t, MPI_Comm)







#define EunhaParticleFuncs2D(vorotype) \
	GenerateInsideRange(xpinner2D, vorotype,x,PosType);\
	GenerateInsideRange(ypinner2D, vorotype,y,PosType);\
//	GenerateMemberCompare(xcompare, vorotype,x,PosType);\
//	GenerateMemberCompare(ycompare, vorotype,y,PosType);\
//	GenerateShift2Pos(Shift2Pos, vorotype);



#define InSideOrNot2D(type, a, lbox, width, simbox, mflag, pflag, def_ptl) \
	int type##InSideBox2D(\
		GridInfo *gridinfo, const void *a, BoxMinMax *lbox, PosType *width, \
		SimBoxRange *simbox, int mflag , int pflag, enum PadPtl def_ptl){ \
	type *aa = (type*)a;\
	PosType aax,aay,aaz;\
	int i,j,k;\
	if(pflag == 4) {\
		i = j = 0;\
	}\
	else {\
		i = pflag%3-1;\
		j = pflag/3 -1;\
	}\
		aax = XofBp(aa);\
		aay = YofBp(aa);\
	aax += i*(simbox->x).max;\
	aay += j*(simbox->y).max;\
		/*\
	DEBUGPRINT("+inside box values are %g <= %g < %g : %g <= %g < %g :: %d %d\n", \
			(lbox->xmin - (*width)), aax, (lbox->xmax + (*width)),\
			(lbox->ymin - (*width)), aay, (lbox->ymax + (*width)), mflag, pflag);\
			*/\
	if( (aax >= (lbox->xmin - (*width)) && aax <(lbox->xmax + (*width))) &&\
			(aay >=(lbox->ymin - (*width)) && aay <(lbox->ymax + (*width)))){\
		if(mflag) {\
			aa->x = aax;\
			aa->y = aay;\
			/*\
			CHANGEINDX(aa,0L);\
			*/\
		}\
		return 1;\
	}\
	else return 0;\
}


#define determineEdgePtl2D(type,a, simbox, width) \
	int type##determineEdgePtl2D(GridInfo *gridinfo, const void *a, SimBoxRange *simbox,\
		PosType *width){\
	type *aa = (type *)a ;\
	PosType aax = aa->x;\
	PosType aay = aa->y;\
	if(aax<=*width) return 1;\
	else if ( (simbox->x).max - aax < *width) return 1;\
	if(aay<=*width) return 1;\
	else if ( (simbox->y).max - aay < *width) return 1;\
	return 0;\
}

#define Dmuladd2D(type, a, b, fact, xyz, ifalg) int type##muladd2D(\
		GridInfo *gridinfo, const void *a, const void *b, float fact,\
		char xyz, int iflag){\
	float *av;\
	PosType bb;\
	switch (xyz) {\
		case 'x':\
			av = &(((type*)a)->x);\
			bb = ((type*)b)->x;\
			break;\
		default: \
			av = &(((type*)a)->y);\
			bb = ((type*)b)->y;\
	}\
	if(iflag){\
		if(isnan(*av) || isnan(bb)) return 0;\
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

#define DShift2Pos2D(shift2pos,vorotype) \
	void vorotype##shift2pos(GridInfo *gridinfo,const void *a){\
		vorotype *aa = (vorotype*)a;\
		/*\
		aa->u4if.indx = -1;\
		*/\
	}







#define DDRule2D(type, a, nmem, com) \
	int type##DDRule2D(GridInfo *gridinfo,const void *a, ptrdiff_t nmem, MPI_Comm com){\
	ptrdiff_t _i;\
	double xmean,xstd;\
	double ymean,ystd;\
	double totalval;\
	ptrdiff_t totalnum;\
	type *aa = (type *)a;\
	xmean = xstd = 0;\
	ymean = ystd = 0;\
	for(_i=0;_i<nmem;_i++){\
		type *aai = aa + _i;\
		/*\
		PosType aax = xPos(gridinfo, aai);\
		PosType aay = yPos(gridinfo, aai);\
		*/\
		PosType aax = aai->x;\
		PosType aay = aai->y;\
		xmean += aax;\
		ymean += aay;\
		xstd += aax*aax;\
		ystd += aay*aay;\
	}\
	MPI_Reduce(&xmean,&totalval,1, MPI_DOUBLE, MPI_SUM, 0, com);\
	MPI_Bcast(&totalval, 1, MPI_DOUBLE, 0, com);\
	\
	xmean = totalval;\
	MPI_Reduce(&ymean,&totalval,1, MPI_DOUBLE, MPI_SUM, 0, com);\
	MPI_Bcast(&totalval, 1, MPI_DOUBLE, 0, com);\
	ymean = totalval;\
\
	MPI_Reduce(&xstd,&totalval,1, MPI_DOUBLE, MPI_SUM, 0, com);\
	MPI_Bcast(&totalval, 1, MPI_DOUBLE, 0, com);\
	xstd = totalval;\
\
	MPI_Reduce(&ystd,&totalval,1, MPI_DOUBLE, MPI_SUM, 0, com);\
	MPI_Bcast(&totalval, 1, MPI_DOUBLE, 0, com);\
	ystd = totalval;\
\
	MPI_Reduce(&nmem,&totalnum,1, MPI_INT64_T, MPI_SUM, 0, com);\
	MPI_Bcast(&totalnum, 1, MPI_INT64_T, 0, com);\
	xmean = xmean/totalnum;\
	ymean = ymean/totalnum;\
	xstd = xstd/totalnum - xmean*xmean;\
	ystd = ystd/totalnum - ymean*ymean;\
	if(xstd >= ystd) return 0;\
	else return 1;\
}




#define debugparticles(simpar) do { \
	PosType xmin,ymin,zmin,xmax,ymax,zmax; \
	xmin = ymin = zmin = 1.E20; \
	xmax = ymax = zmax = -1.E20; \
	ptrdiff_t i,j,k; \
	for(i=0;i<DM_NP(simpar);i++){ \
		xmin = MIN(xmin, XOFP(simpar, DM_BP(simpar)+i)); \
		xmax = MAX(xmax, XOFP(simpar, DM_BP(simpar)+i)); \
		ymin = MIN(ymin, YOFP(simpar, DM_BP(simpar)+i)); \
		ymax = MAX(ymax, YOFP(simpar, DM_BP(simpar)+i)); \
		zmin = MIN(zmin, ZOFP(simpar, DM_BP(simpar)+i)); \
		zmax = MAX(zmax, ZOFP(simpar, DM_BP(simpar)+i)); \
	} \
	DEBUGPRINT("P%d has %ld pregion %g %g : %g %g : %g %g\n",MYID(simpar),DM_NP(simpar),xmin,xmax, ymin,ymax,zmin,zmax); \
	DEBUGPRINT("P%d has %p region %g %g : %g %g : %g %g\n",MYID(simpar), DM_BP(simpar),SIM_LXMIN(simpar,dm), \
			SIM_LXMAX(simpar,dm), SIM_LYMIN(simpar,dm), SIM_LYMAX(simpar,dm), \
			SIM_LZMIN(simpar,dm), SIM_LZMAX(simpar,dm));\
	DEBUGPRINT("P%d has  p value %g %g %g \n",MYID(simpar),XOFP(simpar,DM_BP(simpar)),\
			YOFP(simpar, DM_BP(simpar)), ZOFP(simpar, DM_BP(simpar))) ; \
}while(0)


#define Dcompare2D(compare,vorotype,direction) \
	int vorotype##direction##compare##2D_kp(GridInfo *gridinfo,const void *a, const void *b){\
		postype _a,_b;\
		_a = (((vorotype*)a)->direction);\
		_b = (((vorotype*)b)->direction);\
		if(_a>_b) return 1;\
		else if(_a<_b) return -1;\
		else return 0;\
	}



#define MakeAllDDFuncDDinfo2D(type)  \
	DDRule2D(type##particletype, a, nmem, com);\
	DDRule2D(tree##type##particletype, a, nmem, com);\
	InSideOrNot2D(type##particletype, a, lbox, width, simbox, mflag, pfflag, padptl);\
	InSideOrNot2D(tree##type##particletype, a, lbox, width, simbox, mflag, pfflag, padptl);\
	Dmuladd2D(type##particletype, a, b, fact, xyz, iflag);\
	Dmuladd2D(tree##type##particletype, a, b, fact, xyz, iflag);\
	DShift2Pos2D(shift2pos2D,type##particletype);\
	DShift2Pos2D(shift2pos2D,tree##type##particletype);\
	Dcompare2D(compare,tree##type##particletype,x);\
	Dcompare2D(compare,tree##type##particletype,y);\
	Dcompare2D(compare,type##particletype,x);\
	Dcompare2D(compare,type##particletype,y);\
	EunhaParticleFuncs2D(type##particletype);\
	EunhaParticleFuncs2D(tree##type##particletype);\
	determineEdgePtl2D(type##particletype,a,simbox,width);\
	determineEdgePtl2D(tree##type##particletype,a,simbox,width);



#define CopyDDFuncDDInfoFromTYPE1(simpar,TYPE1,type1, TYPE, type) do{\
    int i,nddinfo = NDDINFO(simpar);\
    TYPE##_DDFUNC(simpar).xcompare = type##particletypexcompare2D_kp;\
    TYPE##_DDFUNC(simpar).ycompare = type##particletypeycompare2D_kp;\
    TYPE##_DDFUNC(simpar).xpinner = type##particletypexpinner2D;\
    TYPE##_DDFUNC(simpar).ypinner = type##particletypeypinner2D;\
    TYPE##_DDFUNC(simpar).muladd = type##particletypemuladd2D;\
    TYPE##_DDFUNC(simpar).insidebox = type##particletypeInSideBox2D;\
    TYPE##_DDFUNC(simpar).edgeptl = type##particletypedetermineEdgePtl2D;\
    TYPE##_DDFUNC(simpar).shift2pos = type##particletypeshift2pos2D;\
    TYPE##_DDFUNC(simpar).divdir = type##particletypeDDRule2D;\
    TYPE##_DDFUNC(simpar).xyzchip = TYPE1##_DDFUNC(simpar).xyzchip;\
    for(i=0;i<NDDINFO(simpar);i++){\
        TYPE##_DDINFO(simpar)[i].xcompare = TYPE##_DDFUNC(simpar).xcompare;\
        TYPE##_DDINFO(simpar)[i].ycompare = TYPE##_DDFUNC(simpar).ycompare;\
        TYPE##_DDINFO(simpar)[i].xpinner = TYPE##_DDFUNC(simpar).xpinner;\
        TYPE##_DDINFO(simpar)[i].ypinner = TYPE##_DDFUNC(simpar).ypinner;\
        TYPE##_DDINFO(simpar)[i].muladd = TYPE##_DDFUNC(simpar).muladd;\
        TYPE##_DDINFO(simpar)[i].insidebox = TYPE##_DDFUNC(simpar).insidebox;\
        TYPE##_DDINFO(simpar)[i].edgeptl = TYPE##_DDFUNC(simpar).edgeptl;\
        TYPE##_DDINFO(simpar)[i].shift2pos = TYPE##_DDFUNC(simpar).shift2pos;\
        TYPE##_DDINFO(simpar)[i].xyzchip = TYPE1##_DDINFO(simpar)[i].xyzchip;\
        TYPE##_DDINFO(simpar)[i].com = TYPE1##_DDINFO(simpar)[i].com;\
        TYPE##_DDINFO(simpar)[i].n_size = sizeof(type##particletype);\
        TYPE##_DDINFO(simpar)[i].nid = TYPE1##_DDINFO(simpar)[i].nid;\
        TYPE##_DDINFO(simpar)[i].myid = TYPE1##_DDINFO(simpar)[i].myid;\
        TYPE##_DDINFO(simpar)[i].nsubgroup = TYPE1##_DDINFO(simpar)[i].nsubgroup;\
        TYPE##_DDINFO(simpar)[i].idirection = TYPE1##_DDINFO(simpar)[i].idirection;\
        TYPE##_DDINFO(simpar)[i].subgroupsize = TYPE1##_DDINFO(simpar)[i].subgroupsize;\
        TYPE##_DDINFO(simpar)[i].subgroupid = TYPE1##_DDINFO(simpar)[i].subgroupid;\
        TYPE##_DDINFO(simpar)[i].sublocalgrid = TYPE1##_DDINFO(simpar)[i].sublocalgrid;\
        TYPE##_DDINFO(simpar)[i].memoffset.xyz.x = offsetof(type##particletype, x);\
        TYPE##_DDINFO(simpar)[i].memoffset.xyz.y = offsetof(type##particletype, y);\
        TYPE##_DDINFO(simpar)[i].memoffset.xyz.z = offsetof(type##particletype, z);\
        TYPE##_DDINFO(simpar)[i].npivot = TYPE1##_DDINFO(simpar)[i].npivot;\
        TYPE##_DDINFO(simpar)[i].lgroup = TYPE1##_DDINFO(simpar)[i].lgroup;\
        TYPE##_DDINFO(simpar)[i].pivot = (void*)malloc(TYPE##_DDINFO(simpar)[i].npivot*TYPE##_DDINFO(simpar)[i].n_size);\
        int j;\
        for(j=0;j<TYPE##_DDINFO(simpar)[i].npivot;j++){\
            ((type##particletype*)(TYPE##_DDINFO(simpar)[i].pivot)+j)->x = \
                ((type1##particletype*)(TYPE1##_DDINFO(simpar)[i].pivot)+j)->x;\
            ((type##particletype*)(TYPE##_DDINFO(simpar)[i].pivot)+j)->y = \
                ((type1##particletype*)(TYPE1##_DDINFO(simpar)[i].pivot)+j)->y;\
            ((type##particletype*)(TYPE##_DDINFO(simpar)[i].pivot)+j)->z = \
                ((type1##particletype*)(TYPE1##_DDINFO(simpar)[i].pivot)+j)->z;\
            ((type##particletype*)(TYPE##_DDINFO(simpar)[i].pivot)+j)->u4if= \
                ((type1##particletype*)(TYPE1##_DDINFO(simpar)[i].pivot)+j)->u4if;\
        }\
    }\
} while(0)
