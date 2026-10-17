#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include<omp.h>
#include "eunha.h"
#include "eunha.sub.h"
//#include "cosmology.h"
#include "mpirks.h"
#include "voro.h"
#include "rt.h"
#include "DD2d.h"
#include "timerutil.h"
#include "mpiaux.h"


float cputime0[1024];
float cputime1[1024];



/* Declare prototypes of functions for RkS */

MakeAllDDFuncDDinfo(voro);
DefineProtoTypeFunctions2D(voro);
DefineProtoTypeFunctions2D(treevoro);


void RT_BuildDDInfo(SimParameters *simpar){
	int nprime;
    PrimeNumber prime[100];
    nprime = getprimenumber(NID(simpar), prime);

    MakeDoDeInfo2D(NID(simpar), prime, voroparticletype, x,y,
			VORO_DDINFO(simpar), NDDINFO(simpar));


    MakeDoDeInfo2D(NID(simpar), prime, treevoroparticletype, x,y,
			TVORO_DDINFO(simpar), NDDINFO(simpar));

    MakeDoDeFunc2D(voroparticletype, VORO_DDFUNC(simpar), xcompare2D_rt,ycompare2D_rt, 
			xpinner, ypinner, muladd2D, DDRule2D, InSideBox2D, 
			determineEdgePtl, shift2pos2D, VORO_DDINFO(simpar));

    MakeDoDeFunc2D(treevoroparticletype, TVORO_DDFUNC(simpar), xcompare2D_rt,ycompare2D_rt, 
			xpinner, ypinner, muladd2D, DDRule2D, InSideBox2D, 
			determineEdgePtl, shift2pos2D, TVORO_DDINFO(simpar));
}

void RT_AllParticleMigrate(SimParameters *simpar){
    MPI_Allreduce(&VORO_NP(simpar), &VORO_TNP(simpar), 1, MPI_PTRDIFF_T, MPI_SUM, MPI_COMM(simpar));
    if(VORO_TNP(simpar)>0){
        pmigrate((void**)(&VORO_BP(simpar)), &VORO_NP(simpar), 
				VORO_DDINFO(simpar), &GRIDINFO(simpar));
    }
    NPSUM(simpar) = VORO_NP(simpar);
}



void RT_TreeAllParticleMigrate(SimParameters *simpar){
    MPI_Allreduce(&VORO_NP(simpar), &VORO_TNP(simpar), 1, MPI_PTRDIFF_T, 
			MPI_SUM, MPI_COMM(simpar));
    DEBUGPRINT("preMigration: P%d has %ld with total np= %ld \n", MYID(simpar), 
			VORO_NP(simpar), VORO_TNP(simpar));
//  HAMB(simpar);
    if(VORO_TNP(simpar)>0){
        pmigrate((void**)(&VORO_TBP(simpar)), &VORO_NP(simpar), 
				TVORO_DDINFO(simpar), &GRIDINFO(simpar));
    }

    NPSUM(simpar) = VORO_TNP(simpar);
    DEBUGPRINT("postMigration: P%d has %ld \n", MYID(simpar), VORO_NP(simpar));
}



void RT_AllParticlePadding(SimParameters *simpar, float width){
	int ndim = 2;
	size_t i;
	VORO_BPP(simpar) = NULL;
	VORO_NPAD(simpar) = 0;
    ppadding(VORO_BP(simpar), VORO_NP(simpar), (void**)(&VORO_BPP(simpar)), 
			&VORO_NPAD(simpar), VORO_DDINFO(simpar), NDDINFO(simpar), 
			COS_SIMBOX(simpar), width, &GRIDINFO(simpar),
			ndim);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<VORO_NPAD(simpar);i++) { 
		CLEAR_FLAG(VORO_BPP(simpar)+i); 
		SET_FLAG(VORO_BPP(simpar)+i, BoundaryGhostflag); 
		SET_FLAG(VORO_BPP(simpar)+i, VOROflag);
	} 
}

void RT_TreeAllParticlePadding(SimParameters *simpar, float width){
	int ndim = 2;
	size_t i;
	VORO_TBPP(simpar) = NULL;
	VORO_NPAD(simpar) = 0;
	DEBUGPRINT("P%d has cos simbox %g %g %g : %g %g %g\n", MYID(simpar), 
			COS_SIMBOX(simpar).x.min,COS_SIMBOX(simpar).y.min,COS_SIMBOX(simpar).z.min,
			COS_SIMBOX(simpar).x.max,COS_SIMBOX(simpar).y.max,COS_SIMBOX(simpar).z.max);

    ppadding(VORO_TBP(simpar), VORO_NP(simpar), (void**)(&VORO_TBPP(simpar)), &VORO_NPAD(simpar),
            TVORO_DDINFO(simpar), NDDINFO(simpar), COS_SIMBOX(simpar), width, &GRIDINFO(simpar),
			ndim);

	DEBUGPRINT("P%d has npad = %ld \n", MYID(simpar), VORO_NPAD(simpar));
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<VORO_NPAD(simpar);i++) { 
		CLEAR_FLAG(VORO_TBPP(simpar)+i); 
		SET_FLAG(VORO_TBPP(simpar)+i, BoundaryGhostflag); 
		SET_FLAG(VORO_TBPP(simpar)+i, VOROflag);
	} 
}



/*
int getdirection2d(SimBoxRange box){
    PosType xwidth = box.x.max-box.x.min;
    PosType ywidth = box.y.max-box.y.min;
    if(ywidth >= xwidth){
        return 1;
    }
	return 0;
}
*/






void RT_ExtractLocalDomainVolume(SimParameters *simpar){
    DoDeInfo *ddinfo = VORO_DDINFO(simpar);
    int i,nddinfo = NDDINFO(simpar);
    SimBoxRange box = SIMBOX(simpar);

    ddinfo->lgroup.xyz.xmin = box.x.min;
    ddinfo->lgroup.xyz.ymin = box.y.min;
    ddinfo->lgroup.xyz.zmin = box.z.min;
    ddinfo->lgroup.xyz.wmin = box.w.min;

    ddinfo->lgroup.xyz.xmax = box.x.max;
    ddinfo->lgroup.xyz.ymax = box.y.max;
    ddinfo->lgroup.xyz.zmax = box.z.max;
    ddinfo->lgroup.xyz.wmax = box.w.max;


	DEBUGPRINT("box range:, %lg %lg %lg : %lg %lg %lg\n", box.x.min, box.y.min, box.z.min,
            box.x.max, box.y.max, box.z.max);
    for(i=0;i<nddinfo;i++){
        /* Just copy the upper level box info and patch the splitting in the target dimension */
        if(i!=0) ddinfo[i].lgroup.r = ddinfo[i-1].lgroup.r;
        if(ddinfo[i].subgroupid==0) {
            voroparticletype  *bp = (voroparticletype*)(ddinfo[i].pivot)+ddinfo[i].subgroupid;
            if(ddinfo[i].xyzchip == 'x') ddinfo[i].lgroup.r.rmax[0] = XOFP(simpar,bp);
            else  ddinfo[i].lgroup.r.rmax[1] = YOFP(simpar,bp);
#ifdef DEBUG
            DEBUGPRINT("-P%d: Now p%d rmin= %g %g : rmax= %g %g: %c\n", ddinfo[i].myid,i,
					ddinfo[i].lgroup.r.rmin[0], ddinfo[i].lgroup.r.rmin[1],
					ddinfo[i].lgroup.r.rmax[0], ddinfo[i].lgroup.r.rmax[1],
                    ddinfo[i].xyzchip);
			DEBUGPRINT("--P%d %g %g %g in %c direction\n", ddinfo[i].myid,bp->x,bp->y, bp->z, ddinfo[i].xyzchip);
#endif
        }
        else if(ddinfo[i].subgroupid != ddinfo[i].nsubgroup-1){
            voroparticletype  *bp1 = (voroparticletype*)(ddinfo[i].pivot)+ddinfo[i].subgroupid;
            voroparticletype  *bp2 = (voroparticletype*)(ddinfo[i].pivot)+ddinfo[i].subgroupid-1;
            if(ddinfo[i].xyzchip == 'x') {
                ddinfo[i].lgroup.r.rmax[0] = XOFP(simpar, bp1);
                ddinfo[i].lgroup.r.rmin[0] = XOFP(simpar, bp2);
            }
            else {
                ddinfo[i].lgroup.r.rmin[1] = YOFP(simpar, bp2);
                ddinfo[i].lgroup.r.rmax[1] = YOFP(simpar, bp1);
            }
#ifdef DEBUG
            DEBUGPRINT("Now p%d rmax %g %g %g :: %g %g %g : %c\n",i,
                    XOFP(simpar,bp2), YOFP(simpar,bp2), ZofP(simpar, bp2),
                    XOFP(simpar,bp1), YOFP(simpar,bp1), ZofP(simpar, bp1), ddinfo[i].xyzchip);
#endif
        }
        else {
            voroparticletype  *bp = (voroparticletype*)(ddinfo[i].pivot)+ddinfo[i].subgroupid-1;
            if(ddinfo[i].xyzchip == 'x') ddinfo[i].lgroup.r.rmin[0] = XOFP(simpar, bp);
            else ddinfo[i].lgroup.r.rmin[1] = YOFP(simpar, bp);
#ifdef DEBUG
			DEBUGPRINT("-Now p%d rmin= %g %g : rmax= %g %g: %c\n",i,
                    ddinfo[i].lgroup.r.rmin[0], ddinfo[i].lgroup.r.rmin[1],
                    ddinfo[i].lgroup.r.rmax[0], ddinfo[i].lgroup.r.rmax[1],
                    ddinfo[i].xyzchip);
#endif
        }
    }
}

void RT_BuildSimpleRMS(SimParameters *simpar, SimBoxRange box,  size_t n_size, DoDeFunc *ddfunc, DoDeInfo *ddinfo, MPI_Comm Comm){
    MPI_Status status;
    int src,dest;
    ptrdiff_t nrecv,nsend;
    ptrdiff_t i,j,k;
    int myid,nid;
    int igroup,subgroupsize,subgroupid;
    ptrdiff_t nowmem=0;


    SimBoxRange sbox = box;

    MPI_Comm_rank(Comm,&myid);
    MPI_Comm_size(Comm,&nid);

    int nsubgroup = getNextPrimeNumber(nid);
    if(nsubgroup == 1) {
        return;
    }
    NDDINFO(simpar) ++;
    int npivot = nsubgroup-1;
    voroparticletype pivot[npivot];

    int idir = getdirection2d(box);
    if(idir==0) ddfunc->xyzchip = 'x';
    else ddfunc->xyzchip = 'y';

#ifdef DEBUG
    DEBUGPRINT("P%d divided along with %c\n",myid, ddfunc->xyzchip);
#endif
    subgroupsize = nid/nsubgroup;
    subgroupid = myid/subgroupsize;
    char xyzchip = ddfunc->xyzchip;


    /*
    void HAMB(SimParameters *);
    HAMB(simpar);
    */

    for(i=0;i<npivot;i++){
        pivot[i].x = pivot[i].y = pivot[i].z = 0;

        if(idir ==0) {
            pivot[i].x = (box.x.max-box.x.min)/(npivot+1.L) *(i+1) + box.x.min;
            DEBUGPRINT("P%d divided along %c with %g\n",myid, ddfunc->xyzchip, pivot[i].x);
        }
        else {
            pivot[i].y = (box.y.max-box.y.min)/(npivot+1.L) *(i+1) + box.y.min;
            DEBUGPRINT("P%d divided along %c with %g\n",myid, ddfunc->xyzchip, pivot[i].y);
        }
#ifdef XYZDBL
        CHANGEINDX( (pivot+i), 0L);
#endif
    }

    if(idir==0) {
        sbox.x.min = (box.x.max-box.x.min)/(npivot+1.L) * subgroupid + box.x.min;
        sbox.x.max = (box.x.max-box.x.min)/(npivot+1.L) * (subgroupid+1) + box.x.min;
    }
    else {
        sbox.y.min = (box.y.max-box.y.min)/(npivot+1.L) * subgroupid + box.y.min;
        sbox.y.max = (box.y.max-box.y.min)/(npivot+1.L) * (subgroupid+1) + box.y.min;
    }
    {
        ddinfo->com = Comm;
        ddinfo->myid = myid;
        ddinfo->nid = nid;
        ddinfo->nsubgroup = nsubgroup;
        ddinfo->npivot = npivot;
        ddinfo->idirection = idir;
        ddinfo->subgroupsize = subgroupsize;
        ddinfo->subgroupid = subgroupid;
        ddinfo->xyzchip = ddfunc->xyzchip;
        ddinfo->xcompare = ddfunc->xcompare;
        ddinfo->ycompare = ddfunc->ycompare;
        ddinfo->zcompare = ddfunc->zcompare;
        ddinfo->wcompare = ddfunc->wcompare;

        ddinfo->xpinner = ddfunc->xpinner;
        ddinfo->ypinner = ddfunc->ypinner;
        ddinfo->zpinner = ddfunc->zpinner;
        ddinfo->muladd = ddfunc->muladd;
        ddinfo->wdist = ddfunc->wdist;

        ddinfo->insidebox = ddfunc->insidebox;
        ddinfo->edgeptl = ddfunc->edgeptl;
        ddinfo->n_size = n_size;
        memmove(ddinfo->pivot, pivot, npivot*n_size);
    }
    int (*compare)(GridInfo *, const void *, const void *);
    if(xyzchip=='x') compare = ddfunc->xcompare;
    else compare = ddfunc->ycompare;
    {
        MPI_Comm newcom;
        int key = myid % subgroupsize;
        MPI_Comm_split(Comm, subgroupid, key, &newcom);
        RT_BuildSimpleRMS(simpar, sbox,  n_size, ddfunc, ddinfo+1,newcom);
        if(subgroupsize==1) MPI_Comm_free(&newcom);
    }
    return;
}



void StartRTRkSDD(SimParameters *simpar){
	int nprime;
	PrimeNumber prime[100];
	nprime = getprimenumber(NID(simpar), prime);

    RT_BuildDDInfo(simpar);
    NDDINFO(simpar) = 0;
	/*
    KH_SIMBOX(simpar).x.min = KH_SIMBOX(simpar).y.min = KH_SIMBOX(simpar).z.min = 0;
    KH_SIMBOX(simpar).x.max = NX(simpar);
    KH_SIMBOX(simpar).y.max = NY(simpar);
    KH_SIMBOX(simpar).z.max = NZ(simpar);
	*/
    SIMBOX(simpar).x.min = RT_XMIN(simpar);
    SIMBOX(simpar).x.max = RT_XMAX(simpar);
    SIMBOX(simpar).y.min = RT_YMIN(simpar);
    SIMBOX(simpar).y.max = RT_YMAX(simpar);
    SIMBOX(simpar).z.min = 0;
    SIMBOX(simpar).z.max = 0;
	DEBUGPRINT("P%d has initial volume %g %g : %g %g\n", MYID(simpar), RT_XMIN(simpar),
			RT_XMAX(simpar), RT_YMIN(simpar), RT_YMAX(simpar));

    RT_BuildSimpleRMS(simpar, RT_SIMBOX(simpar), sizeof(voroparticletype), &VORO_DDFUNC(simpar),
            VORO_DDINFO(simpar), MPI_COMM(simpar));


    RT_ExtractLocalDomainVolume(simpar);
    /*
    ExtractLocalDomainVolume(DM_DDINFO(simpar), NDDINFO(simpar),RT_SIMBOX(simpar));
    */

//    BroadCastDDFuncAndDDInfo(simpar);
	CopyDDFuncDDInfoFromTYPE1(simpar, VORO, voro, TVORO, treevoro);
}
void  rt_DomainDecomp2D(SimParameters *simpar, int iflag){
    ptrdiff_t np, tnp;

    np = VORO_NP(simpar);
    MPI_Reduce(&np, &tnp, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM(simpar));
    MPI_Bcast(&tnp, 1, MPI_INT64_T, 0, MPI_COMM(simpar));
    VORO_TNP(simpar) = tnp;


    float percent = (double)np/(double)tnp*NID(simpar)*100.L;
    float minper, maxper;
    MPI_Reduce(&percent, &minper, 1, MPI_FLOAT, MPI_MIN, 0, MPI_COMM(simpar));
    MPI_Reduce(&percent, &maxper, 1, MPI_FLOAT, MPI_MAX, 0, MPI_COMM(simpar));
    MPI_Bcast(&minper, 1, MPI_FLOAT, 0, MPI_COMM(simpar));
    MPI_Bcast(&maxper, 1, MPI_FLOAT, 0, MPI_COMM(simpar));
    if(MYID(simpar)==0) DEBUGPRINT("The min/max percentage in DD is %g %g\n",minper, maxper);

    if(iflag ==0 &&(minper>=95 && maxper < 105)) return;

    TIMER_START(10);

#ifdef DEBUG
    if(MYID(simpar)==0) DEBUGPRINT("P%d has %ld in DD\n",MYID(simpar), VORO_NP(simpar));
#endif
    mpirks( (void**)(&VORO_BP(simpar)),&VORO_NP(simpar), sizeof(voroparticletype), &VORO_DDFUNC(simpar), VORO_DDINFO(simpar),COM(simpar),
            &GRIDINFO(simpar), oldway, OldCom);

    void SetEntireSimBox(SimParameters *);
    SetEntireSimBox(simpar);

    /*
    void UpdateDdinfoBox(SimParameters *); UpdateDdinfoBox(simpar);
    */

    RT_ExtractLocalDomainVolume(simpar);
    /*
    ExtractLocalDomainVolume(DM_DDINFO(simpar), NDDINFO(simpar),RT_SIMBOX(simpar));
    */

    BroadCastDDFuncAndDDInfo(simpar);


    if(VORO_TNP(simpar)>0){
        pmigrate((void**)(&VORO_BP(simpar)), &VORO_NP(simpar), VORO_DDINFO(simpar), &GRIDINFO(simpar));
        MPI_Allreduce(&VORO_NP(simpar), &VORO_TNP(simpar), 1, MPI_PTRDIFF_T, MPI_SUM, MPI_COMM(simpar));
    }
    TIMER_STOP(10);
    MPI_Barrier(MPI_COMM(simpar));
    if(MYID(simpar) == 0) printf("%d CPU(Domain Decomp) = %lf\n",STEPCOUNT(simpar),ELAPSED_TIME(10));
#ifdef DEBUG
    DEBUGPRINT("P%d has %ld voronoi particles\n",MYID(simpar), VORO_NP(simpar));
#endif
}

