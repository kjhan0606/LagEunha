#define NSave (1<<18)
#define GROUPID(a,b) ((a)/(b))
#define RANKINGROUP(a,b) ( (a)%(b) )

#define RESET_TYPE_FOFFLAG(simpar, TYPE) do{\
	size_t i,np = TYPE##_NP(simpar);\
	for(i=0;i<np;i++){\
		UNSET_P_FLAG(simpar, TYPE, i,FoFflag);\
	}\
}while(0)


#define RESET_FOFFLAG(simpar) do{\
	RESET_TYPE_FOFFLAG(simpar, DM);\
	RESET_TYPE_FOFFLAG(simpar, SPH);\
	RESET_TYPE_FOFFLAG(simpar, STAR);\
	RESET_TYPE_FOFFLAG(simpar, AGN);\
}while(0)



void write_zminmax_header(ptrdiff_t, SimParameters *, FILE *, double);
void PreFoF(SimParameters *);
