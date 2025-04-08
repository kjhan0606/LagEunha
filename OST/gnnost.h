/*      */
/*
#define GENERAL_TYPE unsigned int type: 1
#define GENERALTPtlPOINTER GENERAL_TYPE;void *sibling
*/
/*      */
typedef struct nearestneighbor{
	PosType dist2;
	void *bp;
} Neighbor;

typedef struct gpu_nearestneighbor{
	PosType dist2;
	int ibp;
} GPU_Neighbor;


void Make_GNN_Tree( TStruct *, size_t , 
		void (*)(TStruct *, TStruct *), 
		void (*)(TStruct *), void (*)(TStruct *), int 
		);
PosType find_GNearest( void *, TStruct *, 
		int (*)(void *, TStruct *, PosType ),
		PosType (*)(void *, void *)
		);
