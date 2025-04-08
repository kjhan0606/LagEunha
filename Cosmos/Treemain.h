/*
#define Initialize_Particle_Accel(simpar) do{\
	Init_Accel(simpar, DM, dm);\
	Init_Accel(simpar, SPH, sph);\
	Init_Accel(simpar, STAR, star);\
	Init_Accel(simpar, AGN, agn);\
}while(0)
*/


#define HubbleFlowForNbody(simpar, evolfactor) do{\
	HubbleFlowForNbodyStep(simpar, DM, dm, evolfactor);\
	HubbleFlowForNbodyStep(simpar, SPH, sph, evolfactor);\
	HubbleFlowForNbodyStep(simpar, STAR, star, evolfactor);\
	HubbleFlowForNbodyStep(simpar, AGN, agn, evolfactor);\
}while(0)


void TreeMain(SimParameters *, DeterminedEvolFact *);
