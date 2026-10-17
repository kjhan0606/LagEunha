#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include "eunha.h"
//#include "mpirms.h"
#include "cosmology.h"
#include "cosmosic.h"
#include "eunha.sub.h"
#include "timerutil.h"
#include "domaindecomp.h"

static float cputime0[16];
static float cputime1[16];

void BuildCosmosDDInfo(SimParameters *);
void StartCosmosRkSDD(SimParameters *simpar){ 
	BuildCosmosDDInfo(simpar); 
	NDDINFO(simpar) = 0; 
	COS_SIMBOX(simpar).x.min = COS_SIMBOX(simpar).y.min = COS_SIMBOX(simpar).z.min = 0;
	COS_SIMBOX(simpar).x.max = NX(simpar); 
	COS_SIMBOX(simpar).y.max = NY(simpar); 
	COS_SIMBOX(simpar).z.max = NZ(simpar);

	BuildSimpleRMS(simpar, COS_SIMBOX(simpar), sizeof(dmparticletype), &DM_DDFUNC(simpar), 
			DM_DDINFO(simpar), MPI_COMM(simpar)); 


	DM_ExtractLocalDomainVolume(simpar);
	/*
	ExtractLocalDomainVolume(DM_DDINFO(simpar), NDDINFO(simpar),COS_SIMBOX(simpar));
	*/

	BroadCastDDFuncAndDDInfo(simpar);
}

void  DomainDecomp(SimParameters *simpar, int iflag){
	ptrdiff_t np, tnp;

	np = DM_NP(simpar);
	MPI_Reduce(&np, &tnp, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM(simpar));
	MPI_Bcast(&tnp, 1, MPI_INT64_T, 0, MPI_COMM(simpar));
	DM_TNP(simpar) = tnp;

#ifndef GOTPM
	np = SPH_NP(simpar);
	MPI_Reduce(&np, &tnp, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM(simpar));
	MPI_Bcast(&tnp, 1, MPI_INT64_T, 0, MPI_COMM(simpar));
	SPH_TNP(simpar) = tnp;

	np = STAR_NP(simpar);
	MPI_Reduce(&np, &tnp, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM(simpar));
	MPI_Bcast(&tnp, 1, MPI_INT64_T, 0, MPI_COMM(simpar));
	STAR_TNP(simpar) = tnp;

	np = AGN_NP(simpar);
	MPI_Reduce(&np, &tnp, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM(simpar));
	MPI_Bcast(&tnp, 1, MPI_INT64_T, 0, MPI_COMM(simpar));
	AGN_TNP(simpar) = tnp;
#endif


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
	if(MYID(simpar)==0) DEBUGPRINT("P%d has %ld %ld %ld %ld in DD\n",MYID(simpar), DM_TNP(simpar), SPH_TNP(simpar), STAR_TNP(simpar), AGN_TNP(simpar));
#endif
	mpirks( (void**)(&DM_BP(simpar)),&DM_NP(simpar), sizeof(dmparticletype), &DM_DDFUNC(simpar), DM_DDINFO(simpar),COM(simpar), 
			&GRIDINFO(simpar), oldway, OldCom);

	/*
	COS_SIMBOX(simpar).x.min  = 0;
	COS_SIMBOX(simpar).y.min  = 0;
	COS_SIMBOX(simpar).z.min  = 0;
	COS_SIMBOX(simpar).x.max  = NX(simpar);
	COS_SIMBOX(simpar).y.max  = NY(simpar);
	COS_SIMBOX(simpar).z.max  = NZ(simpar);
	*/

	void SetEntireSimBox(SimParameters *); 
	SetEntireSimBox(simpar);

	/*
	void UpdateDdinfoBox(SimParameters *); UpdateDdinfoBox(simpar);
	*/

	DM_ExtractLocalDomainVolume(simpar);
	/*
	ExtractLocalDomainVolume(DM_DDINFO(simpar), NDDINFO(simpar),COS_SIMBOX(simpar));
	*/

	BroadCastDDFuncAndDDInfo(simpar);


	if(SPH_TNP(simpar)>0){ 
		pmigrate((void**)(&SPH_BP(simpar)), &SPH_NP(simpar), SPH_DDINFO(simpar), &GRIDINFO(simpar)); 
		MPI_Allreduce(&SPH_NP(simpar), &SPH_TNP(simpar), 1, MPI_PTRDIFF_T, MPI_SUM, MPI_COMM(simpar)); 
	} 
	if(VORO_TNP(simpar)>0){ 
		pmigrate((void**)(&VORO_BP(simpar)), &VORO_NP(simpar), VORO_DDINFO(simpar), &GRIDINFO(simpar)); 
		MPI_Allreduce(&VORO_NP(simpar), &VORO_TNP(simpar), 1, MPI_PTRDIFF_T, MPI_SUM, MPI_COMM(simpar)); 
	} 
	if(STAR_TNP(simpar)>0){ 
		pmigrate((void**)(&STAR_BP(simpar)), &STAR_NP(simpar), STAR_DDINFO(simpar), &GRIDINFO(simpar)); 
		MPI_Allreduce(&STAR_NP(simpar), &STAR_TNP(simpar), 1, MPI_PTRDIFF_T, MPI_SUM, MPI_COMM(simpar)); 
	} 
	if(AGN_TNP(simpar)>0){ 
		pmigrate((void**)(&AGN_BP(simpar)), &AGN_NP(simpar), AGN_DDINFO(simpar), &GRIDINFO(simpar)); 
		MPI_Allreduce(&AGN_NP(simpar), &AGN_TNP(simpar), 1, MPI_PTRDIFF_T, MPI_SUM, MPI_COMM(simpar)); 
	}
	TIMER_STOP(10);
	MPI_Barrier(MPI_COMM(simpar));
	if(MYID(simpar) == 0) printf("%d CPU(Domain Decomp) = %lf\n",STEPCOUNT(simpar),ELAPSED_TIME(10));
#ifdef DEBUG
	DEBUGPRINT("P%d has %ld dm particles\n",MYID(simpar), DM_NP(simpar));
#endif
}


void DM_ExtractLocalDomainVolume(SimParameters *simpar){
	DoDeInfo *ddinfo = DM_DDINFO(simpar);
	int i,nddinfo = NDDINFO(simpar);
	SimBoxRange box = COS_SIMBOX(simpar);

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
			dmparticletype  *bp = (dmparticletype*)(ddinfo[i].pivot)+ddinfo[i].subgroupid;
			if(ddinfo[i].xyzchip == 'x') ddinfo[i].lgroup.r.rmax[0] = XofP(simpar,bp);
			else if (ddinfo[i].xyzchip == 'y') ddinfo[i].lgroup.r.rmax[1] = YofP(simpar,bp);
			else if (ddinfo[i].xyzchip == 'z') ddinfo[i].lgroup.r.rmax[2] = ZofP(simpar,bp);
#ifdef DEBUG
			DEBUGPRINT("-Now p%d rmin %g %g %g : rmax %g %g %g %c\n",i, 
					ddinfo[i].lgroup.r.rmin[0], ddinfo[i].lgroup.r.rmin[1],ddinfo[i].lgroup.r.rmin[2],
					ddinfo[i].lgroup.r.rmax[0], ddinfo[i].lgroup.r.rmax[1],ddinfo[i].lgroup.r.rmax[2],
					ddinfo[i].xyzchip);
#endif
		}
		else if(ddinfo[i].subgroupid != ddinfo[i].nsubgroup-1){
			dmparticletype  *bp1 = (dmparticletype*)(ddinfo[i].pivot)+ddinfo[i].subgroupid;
			dmparticletype  *bp2 = (dmparticletype*)(ddinfo[i].pivot)+ddinfo[i].subgroupid-1;
			if(ddinfo[i].xyzchip == 'x') {
				ddinfo[i].lgroup.r.rmax[0] = XofP(simpar, bp1);
				ddinfo[i].lgroup.r.rmin[0] = XofP(simpar, bp2);
			}
			else if (ddinfo[i].xyzchip == 'y') {
				ddinfo[i].lgroup.r.rmin[1] = YofP(simpar, bp2);
				ddinfo[i].lgroup.r.rmax[1] = YofP(simpar, bp1);
			}
			else if  (ddinfo[i].xyzchip == 'z') {
				ddinfo[i].lgroup.r.rmin[2] = ZofP(simpar, bp2);
				ddinfo[i].lgroup.r.rmax[2] = ZofP(simpar, bp1);
			}
#ifdef DEBUG
			DEBUGPRINT("Now p%d rmin %g %g %g :: rmax %g %g %g : %c\n",i,
					ddinfo[i].lgroup.r.rmin[0], ddinfo[i].lgroup.r.rmin[1],ddinfo[i].lgroup.r.rmin[2],
					ddinfo[i].lgroup.r.rmax[0], ddinfo[i].lgroup.r.rmax[1],ddinfo[i].lgroup.r.rmax[2],
					ddinfo[i].xyzchip);
#endif
		}
		else {
			dmparticletype  *bp = (dmparticletype*)(ddinfo[i].pivot)+ddinfo[i].subgroupid-1;
			if(ddinfo[i].xyzchip == 'x') ddinfo[i].lgroup.r.rmin[0] = XofP(simpar, bp);
			else if (ddinfo[i].xyzchip == 'y') ddinfo[i].lgroup.r.rmin[1] = YofP(simpar, bp);
			else if (ddinfo[i].xyzchip == 'z') ddinfo[i].lgroup.r.rmin[2] = ZofP(simpar, bp);
#ifdef DEBUG
			DEBUGPRINT("+Now p%d rmin %g %g %g : rmax= %g %g %g %c\n",i, 
					ddinfo[i].lgroup.r.rmin[0], ddinfo[i].lgroup.r.rmin[1],ddinfo[i].lgroup.r.rmin[2],
					ddinfo[i].lgroup.r.rmax[0], ddinfo[i].lgroup.r.rmax[1],ddinfo[i].lgroup.r.rmax[2],
					ddinfo[i].xyzchip);
#endif
		}
	}
}
