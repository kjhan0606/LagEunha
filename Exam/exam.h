void ex2d_findCentroid(TStruct *);
void ex3d_findCentroid(TStruct *);
void ex2d_findCellSize(TStruct *);
void ex3d_findCellSize(TStruct *);
void ex2d_idivision(TStruct *, TStruct *);
void ex3d_idivision(TStruct *, TStruct *);
void det3d_dpq(SimParameters *, void (*)(SimParameters *, postype));
void det2d_dpq(SimParameters *, void (*)(SimParameters *, postype));

typedef HydroTreeLinkedCell CellType;
void mkLinkedList2D(SimParameters *, postype , 
		postype , postype , postype , postype , 
		void (*)(SimParameters *, postype)); 

int periodicity(int , int , int );


treevorork4particletype *findCellBP2D(SimParameters *,int , int , int *);

Voro2D_point *searchCellNeighbors2D(SimParameters *, int , int , int *);

double getAccVoro2D(SimParameters *, postype , postype , 
		postype , postype , postype ,  postype, postype,
		void (*)(SimParameters *, postype),
		Voro2D_point *(*)(SimParameters *, int , int, int *),
		treevorork4particletype *(*)(SimParameters *, int , int , int *)
		);

double exam2d_vph_rk4( SimParameters *, void (*)(SimParameters *, postype),
		postype (*)(SimParameters *, postype, postype, postype),
		Voro2D_point *(*)(SimParameters *, int , int, int *),
		treevorork4particletype *(*)(SimParameters *, int , int , int *)
		);
void exam2d_findVol( SimParameters *, void (*)(SimParameters *, postype),
		Voro2D_point *(*)(SimParameters *, int , int, int *),
		treevorork4particletype *(*)(SimParameters *, int , int , int *)
		);
void updateDenW2Pressure2D(SimParameters *, postype , postype , postype , postype ,
		 postype , 
		 void (*)(SimParameters *, postype),
		 Voro2D_point *(*)(SimParameters *, int, int, int *), 
		 treevorork4particletype *(*)(SimParameters *, int , int , int *)
		 );
