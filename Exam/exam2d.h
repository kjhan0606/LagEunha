#define XOFP(simpar,bp) (bp->x)
#define YOFP(simpar,bp) (bp->y)
/* These are only for the header file for proto type declaration of functions */
#define DefineProtoTypeFunctions2D(type) \
    int(type##particletypexcompare2D)(GridInfo *,const void *, const void *);\
    int(type##particletypeycompare2D)(GridInfo *,const void *, const void *);\
    int(type##particletypexpinner2D)(GridInfo *,const void *, const void *, const void *, PosType, int, Range *);\
    int(type##particletypeypinner2D)(GridInfo *,const void *, const void *, const void *, PosType, int, Range *);\
    int(type##particletypemuladd2D)(GridInfo *,const void *, const void *, float, char, int);\
    int(type##particletypeInSideBox2D)(GridInfo *,const void *, BoxMinMax *, PosType *, SimBoxRange *, int, int,\
			enum PadPtl);\
    int(type##particletypedetermineEdgePtl2D)(GridInfo *,const void *, SimBoxRange *, PosType *);\
    void(type##particletypeshift2pos2D)(GridInfo *,const void *);\
    int(type##particletypeDDRule2D)(GridInfo *, const void *, ptrdiff_t, MPI_Comm)

DefineProtoTypeFunctions2D(vorork4);
DefineProtoTypeFunctions2D(treevorork4);

#define makeLink4DDFunc2D(type,ddfunc,xcompare2d,ycompare2d,xpinner2d,ypinner2d,dmuladd2d,\
        divdir2d,InSideBox2d,Edgeptl2d,Dshift2pos2d,ddinfo) do{\
    ddfunc.xcompare = type##xcompare2d;\
    ddfunc.ycompare = type##ycompare2d;\
    ddfunc.xpinner = type##xpinner2d;\
    ddfunc.ypinner = type##ypinner2d;\
    ddfunc.muladd = type##dmuladd2d;\
    ddfunc.divdir = type##divdir2d;\
    ddfunc.insidebox = type##InSideBox2d;\
    ddfunc.edgeptl = type##Edgeptl2d;\
    ddfunc.shift2pos = type##Dshift2pos2d;\
    ddfunc.xyzchip = 'y';\
} while(0)

#define buildDDInfo2D(simpar) do{\
	int nprime;\
    PrimeNumber prime[100];\
    nprime = getprimenumber(NID(simpar), prime);\
    MakeDoDeInfo2D(NID(simpar), prime, vorork4particletype, x,y,\
            VORORK4_DDINFO(simpar), NDDINFO(simpar));\
    MakeDoDeInfo2D(NID(simpar), prime, treevorork4particletype, x,y,\
            TVORORK4_DDINFO(simpar), NDDINFO(simpar));\
    makeLink4DDFunc2D(vorork4particletype, VORORK4_DDFUNC(simpar), \
			xcompare2D,ycompare2D,\
            xpinner2D, ypinner2D, muladd2D, DDRule2D, InSideBox2D,\
            determineEdgePtl2D, shift2pos2D, VORORK4_DDINFO(simpar));\
    makeLink4DDFunc2D(treevorork4particletype, TVORORK4_DDFUNC(simpar), \
			xcompare2D,ycompare2D,\
            xpinner2D, ypinner2D, muladd2D, DDRule2D, InSideBox2D,\
            determineEdgePtl2D, shift2pos2D, TVORORK4_DDINFO(simpar));\
}while(0)

#define startRkSDD2D(simpar, model) do {\
	DEBUGPRINT("P%d: "#model" is now starting buildDDInfo2D(simpar)\n", MYID(simpar));\
	buildDDInfo2D(simpar); \
    NDDINFO(simpar) = 0; \
    SIMBOX(simpar).x.min = model##_XMIN(simpar); \
    SIMBOX(simpar).x.max = model##_XMAX(simpar); \
    SIMBOX(simpar).y.min = model##_YMIN(simpar); \
    SIMBOX(simpar).y.max = model##_YMAX(simpar); \
	DEBUGPRINT("P%d: "#model" has x/y min and max = %g %g %g %g\n", MYID(simpar),model##_XMIN(simpar),\
			model##_YMIN(simpar), model##_XMAX(simpar), model##_YMAX(simpar));\
    SIMBOX(simpar).z.min = 0; \
    SIMBOX(simpar).z.max = 0; \
	buildSimpleRMS(simpar, SIMBOX(simpar), sizeof(vorork4particletype), \
			&VORORK4_DDFUNC(simpar),\
            VORORK4_DDINFO(simpar), MPI_COMM(simpar));\
	extractLocalDomainVolume(simpar);\
    copyDDFuncDDInfoFromTYPE1(simpar, VORORK4, vorork4, TVORORK4, treevorork4);\
} while(0)




#define icopysign(a) (signbit(a)? -1:1)
void SetEnvGridInfo(SimParameters *);
void TreeAllParticleMigrate(SimParameters *);
void AllParticleMigrate(SimParameters *);
void readGrafic(SimParameters *);
void readGraficFile(SimParameters *);
void migrateTreeVoroParticles(SimParameters *);
void migrateVoroParticles(SimParameters *);
void buildSimpleRMS(SimParameters *, SimBoxRange ,  size_t , DoDeFunc *, DoDeInfo *, MPI_Comm );
void extractLocalDomainVolume(SimParameters *);
void paddingVoroParticles(SimParameters *, postype );
void paddingTreeVoroParticles(SimParameters *, postype );
/*
void updateDenW2Pressure2D(SimParameters *, postype , postype , postype , postype , postype ,
		void (*)(SimParameters *, postype));
		*/







#define generateFunc4Shift2Pos(Shift2Pos,ptype) void ptype##shift2pos(GridInfo *gridinfo,const void *a){\
    ptype *aa = (ptype*)a;\
    aa->x = xPos(gridinfo, aa);\
    aa->y = yPos(gridinfo, aa);\
    aa->u4if.indx = 0;\
}
#define generateFunc4MemberCompare(compare,ptype,member,memtype) int ptype##compare(GridInfo *gridinfo,const void *a, const void *b){\
    memtype _a,_b;\
    ptype *aa = (ptype*)a;\
    ptype *bb = (ptype*)b;\
    _a = member##Pos(gridinfo, aa);\
    _b = member##Pos(gridinfo, bb);\
    if(_a>_b) return 1;\
    else if(_a == _b) return 0;\
    else return -1;\
}
#define generateFunc4InsideRange(dist,ptype,member,memtype) int ptype##dist(GridInfo *gridinfo, \
        const void *a, const void *l, const void *u, PosType width, int iperiod, Range *L){\
    memtype _a,_u,_l;\
    ptype *aa = (ptype*)a;\
    _a = member##Pos(gridinfo, aa);\
    ptype *ll = (ptype*)l;\
    ptype *uu = (ptype*)u;\
    _l = member##Pos(gridinfo,ll) - width + iperiod*(L->max);\
    _u = member##Pos(gridinfo,uu) + width + iperiod*(L->max);\
    if(_a >= _l && _a <_u) return 1;\
    else return 0;\
}





#define EunhaParticleFuncs2D(vorotype) \
	generateFunc4InsideRange(xpinner2D, vorotype,x,PosType);\
	generateFunc4InsideRange(ypinner2D, vorotype,y,PosType);\
	generateFunc4MemberCompare(xcompare, vorotype,x,PosType);\
	generateFunc4MemberCompare(ycompare, vorotype,y,PosType);\
	generateFunc4Shift2Pos(Shift2Pos, vorotype);



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




/*
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
*/


#define Dcompare2D(compare,vorotype,direction) \
	int vorotype##direction##compare##2D(GridInfo *gridinfo,const void *a, const void *b){\
		postype _a,_b;\
		_a = (((vorotype*)a)->direction);\
		_b = (((vorotype*)b)->direction);\
		if(_a>_b) return 1;\
		else if(_a<_b) return -1;\
		else return 0;\
	}



#define generateAllDDFunc2D(type)  \
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


#define copyDDFuncDDInfoFromTYPE1(simpar,TYPE1,type1, TYPE, type) do{\
    int i,nddinfo = NDDINFO(simpar);\
    TYPE##_DDFUNC(simpar).xcompare = type##particletypexcompare2D;\
    TYPE##_DDFUNC(simpar).ycompare = type##particletypeycompare2D;\
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
