/*
 * This is the header file for the Exam source directories.
 * Last update (27/08/2025)
*/

/*
#define getw2forHydroParticle(simpar, bp) \
	GAS_Kappa(simpar)*GAS_Kappa(simpar)*GAS_dMean(simpar)*GAS_dMean(simpar)\
	*pow((bp)->pressure*GAS_invw2Scale(simpar), GAS_w2Power(simpar))
	*/
postype getw2forHydroParticle(SimParameters *, treevorork4particletype *, postype);


void ex2d_findCentroid(TStruct *);
void ex3d_findCentroid(TStruct *);
void ex2d_findCellSize(TStruct *);
void ex3d_findCellSize(TStruct *);
void ex2d_idivision(TStruct *, TStruct *);
void ex3d_idivision(TStruct *, TStruct *);
void det3d_dpqRK4(SimParameters *, void (*)(SimParameters *, postype));
void det2d_dpqRK4(SimParameters *, void (*)(SimParameters *, postype));
void det3d_dpq(SimParameters *, void (*)(SimParameters *, postype));
void det2d_dpq(SimParameters *, void (*)(SimParameters *, postype));
int nearest2dOpen(void *, TStruct *, postype );
int nearest3dOpen(void *, TStruct *, postype );

typedef HydroTreeLinkedCell CellType;
void mkLinkedList2D_oExam(SimParameters *, postype , 
		postype , postype , postype , postype , 
		void (*)(SimParameters *, postype)); 
void mkLinkedList2D_rt(SimParameters *, postype , 
		postype , postype , postype , postype , 
		void (*)(SimParameters *, postype)); 

int periodicity(int , int , int );


treevorork4particletype *findCellRk4BP2D(SimParameters *,int , int , int *);

Voro2D_point *searchCellRk4Neighbors2D(SimParameters *, int , int , int *);

double getAccVoro2D(SimParameters *, postype , postype , 
		postype , postype , postype ,  postype, postype,
		void (*)(SimParameters *, postype),
		Voro2D_point *(*)(SimParameters *, int , int, int *),
		treevorork4particletype *(*)(SimParameters *, int , int , int *), 
		void (*)(SimParameters *, postype , 
			postype , postype , postype , postype, void (*)(SimParameters *, postype))
		);

double exam2d_vph_int_rt(
		SimParameters *, 
		void (*)(SimParameters *, postype),
		double (*)(SimParameters *, postype, postype, postype),
		Voro2D_point *(*)(SimParameters *, int , int, int *),
		treevorork4particletype *(*)(SimParameters *, int , int , int *),
		int (*)(treevorork4particletype*, postype, postype),
		void (*)(SimParameters *, postype , 
			postype , postype , postype , postype, void (*)(SimParameters *, postype))
		);
double exam2d_vph_int(
		SimParameters *, 
		void (*)(SimParameters *, postype),
		double (*)(SimParameters *, postype, postype, postype),
		Voro2D_point *(*)(SimParameters *, int , int, int *),
		treevorork4particletype *(*)(SimParameters *, int , int , int *),
		int (*)(treevorork4particletype*, postype, postype),
		void (*)(SimParameters *, postype , 
			postype , postype , postype , postype, void (*)(SimParameters *, postype))
		);
double exam2d_vph_rk4_int(
		SimParameters *, 
		void (*)(SimParameters *, postype),
		double (*)(SimParameters *, postype, postype, postype),
		Voro2D_point *(*)(SimParameters *, int , int, int *),
		treevorork4particletype *(*)(SimParameters *, int , int , int *),
		int (*)(treevorork4particletype*, postype, postype),
		void (*)(SimParameters *, postype , 
			postype , postype , postype , postype, void (*)(SimParameters *, postype))
		);
double exam2d_vph_rk4_int_rt(
		SimParameters *, 
		void (*)(SimParameters *, postype),
		double (*)(SimParameters *, postype, postype, postype),
		Voro2D_point *(*)(SimParameters *, int , int, int *),
		treevorork4particletype *(*)(SimParameters *, int , int , int *),
		int (*)(treevorork4particletype*, postype, postype),
		void (*)(SimParameters *, postype , 
			postype , postype , postype , postype, void (*)(SimParameters *, postype))
		);

double exam2d_vph_rk4( SimParameters *, void (*)(SimParameters *, postype),
		double (*)(SimParameters *, postype, postype, postype),
		Voro2D_point *(*)(SimParameters *, int , int, int *),
		treevorork4particletype *(*)(SimParameters *, int , int , int *),
		void (*)(SimParameters *, postype , 
			postype , postype , postype , postype, void (*)(SimParameters *, postype))
		);
void exam2dUpdateVol( SimParameters *, void (*)(SimParameters *, postype),
		Voro2D_point *(*)(SimParameters *, int , int, int *),
		treevorork4particletype *(*)(SimParameters *, int , int , int *),
		void (*)(SimParameters *, postype , 
			postype , postype , postype , postype, void (*)(SimParameters *, postype))
		);
void updateDenW2Pressure2D(SimParameters *, postype , postype , postype , postype ,
		 postype , 
		 void (*)(SimParameters *, postype),
		 Voro2D_point *(*)(SimParameters *, int, int, int *), 
		 treevorork4particletype *(*)(SimParameters *, int , int , int *),
		 void (*)(SimParameters *, postype ,
             postype , postype , postype , postype, void (*)(SimParameters *, postype)),
		 postype
		 );

/* NS stress + HLLC + CD10 two-tier blending functions */
void updateDenW2Pressure2DBlend(SimParameters *, postype , postype , postype , postype ,
		postype ,
		void (*)(SimParameters *, postype),
		Voro2D_point *(*)(SimParameters *, int, int, int *),
		treevorork4particletype *(*)(SimParameters *, int , int , int *),
		void (*)(SimParameters *, postype ,
			postype , postype , postype , postype, void (*)(SimParameters *, postype)),
		postype
		);
double getAccVoro2DBlend(SimParameters *, postype , postype ,
		postype , postype , postype , postype, postype,
		void (*)(SimParameters *, postype),
		Voro2D_point *(*)(SimParameters *, int , int, int *),
		treevorork4particletype *(*)(SimParameters *, int , int , int *),
		void (*)(SimParameters *, postype ,
			postype , postype , postype , postype, void (*)(SimParameters *, postype))
		);
double exam2d_vph_rk4_int_blend(
		SimParameters *,
		void (*)(SimParameters *, postype),
		double (*)(SimParameters *, postype, postype, postype),
		Voro2D_point *(*)(SimParameters *, int , int, int *),
		treevorork4particletype *(*)(SimParameters *, int , int , int *),
		int (*)(treevorork4particletype*, postype, postype),
		void (*)(SimParameters *, postype ,
			postype , postype , postype , postype, void (*)(SimParameters *, postype))
		);

void dumpRk4Data2D(SimParameters *, int , postype , postype );
void readRk4Data2D(SimParameters *, postype *, postype *,int );

void exam2d_centroidShift(
		SimParameters *,
		void (*)(SimParameters *, postype),
		Voro2D_point *(*)(SimParameters *, int, int, int *),
		treevorork4particletype *(*)(SimParameters *, int , int , int *),
		void (*)(SimParameters *, postype , 
			postype , postype , postype , postype, void (*)(SimParameters *, postype))
		);

postype ex2d_dist( void *, void *);
postype ex3d_dist( void *, void *);


postype w2in2D(SimParameters *, treevorork4particletype *);


/*
double exam2d_vph_kickV(
		SimParameters *, 
		void (*)(SimParameters *, postype),
		Voro2D_point *(*)(SimParameters *, int, int, int *),
		treevorork4particletype *(*)(SimParameters *, int , int , int *)
		);
void exam2d_vph_dragE(
		SimParameters *,
		void (*)(SimParameters *, postype),
		Voro2D_point *(*)(SimParameters *, int, int, int *),
		treevorork4particletype *(*)(SimParameters *, int , int , int *)
		);
void exam2d_vph_centroidShift(
		SimParameters *,
		void (*)(SimParameters *, postype),
		Voro2D_point *(*)(SimParameters *, int, int, int *),
		treevorork4particletype *(*)(SimParameters *, int , int , int *)
		);
		*/
