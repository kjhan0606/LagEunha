typedef struct Voro2dGasParticle{
	indxflag u4if;
	int iregion;
	postype x,y,den,mass,pressure,vx,vy,ax,ay,energy,te,ie,ke,die,dke,dte,csound,dt;
	postype volume,minsize2;
	struct Voro2dGasParticle *next;
} Voro2dGasParticle;

typedef struct Voro3D_GasParticle{
	indxflag u4if;
	int iregion;
	postype x,y,z,den,mass,pressure,vx,vy,vz,ax,ay,az,energy,te,ie,ke,die,dke,dte,csound,dt;
	postype volume,minsize2;
	struct Voro3D_GasParticle *next;
} Voro3D_GasParticle;


postype findSDen3D(Voro3D_Vertex *, int , Voro3D_point *, int *, int , Voro3D_GasParticle *);
/* This returns the average pressure on the edge according to cbp (cf. Lilek & Peric (1995) )*/
postype getAvgPressure2D(int , Voro2D_Corner *, Voro2dGasParticle *, Voro2D_point *);
void get2dAreaAvgNeighorPressure(
        treevorork4particletype *,
        Voro2D_Corner *,
        Voro2D_point *,
        treevorork4particletype *);

Voro2D_point get2dUpqradRk4(treevorork4particletype *, treevorork4particletype *, postype );
Voro3D_point get3dUpqradRk4(treevorork4particletype *, treevorork4particletype *, postype );
Voro2D_point get2dUpqrad(treevoroparticletype *, treevoroparticletype *, postype );
Voro3D_point get3dUpqrad(treevoroparticletype *, treevoroparticletype *, postype );

#define Voro_Pressure2D(pi,pj,tmp,neighwork,simpar,MODEL) \
	( (0.5-MODEL##_OA(simpar)/3.)*(pj + pi) + MODEL##_OA(simpar)/3* \
	  (((treevoroparticletype*) (neighwork[(tmp->upperlink)->upperrelated].bp))->pressure +\
	   ((treevoroparticletype*) (neighwork[(tmp           )->lowerrelated].bp))->pressure))

#define getPressure2D(pi,pj,tmp,neighwork,OoA) (\
		(0.5-OoA/3.)*(pj + pi) + OoA/3* \
		( ((treevoroparticletype*) (neighwork[(tmp->upperlink)->upperrelated].bp))->pressure +\
		 ((treevoroparticletype*) (neighwork[(tmp           )->lowerrelated].bp))->pressure ) \
	   )


#define Voro_Velocity2D(ai,aj,tmp,neighwork,simpar,MODEL) \
	({Voro2D_point res; \
	 res.x=(0.5-MODEL##_OA(simpar)/3.)*(ai->vx + aj->vx) + MODEL##_OA(simpar)/3* \
	  (((treevorork4particletype*) (neighwork[(tmp->upperlink)->upperrelated].vx))+\
	   ((treevorork4particletype*) (neighwork[(tmp           )->lowerrelated].vx))) \
	 res.y=(0.5-MODEL##_OA(simpar)/3.)*(ai->vy + aj->vy) + MODEL##_OA(simpar)/3* \
	  (((treevorork4particletype*) (neighwork[(tmp->upperlink)->upperrelated].vy))+\
	   ((treevorork4particletype*) (neighwork[(tmp           )->lowerrelated].vy)));\
	   res;}\
	 )

#define VoroRK4_Pressure2D(pi,pj,tmp,neighwork,simpar,MODEL) \
	( (0.5-MODEL##_OA(simpar)/3.)*(pj + pi) + MODEL##_OA(simpar)/3* \
	  (((treevorork4particletype*) (neighwork[(tmp->upperlink)->upperrelated].bp))->pressure +\
	   ((treevorork4particletype*) (neighwork[(tmp           )->lowerrelated].bp))->pressure))

#define getPressure2DRK4(pi,pj,tmp,neighwork,OoA) (\
		(0.5-OoA/3.)*(pj + pi) + OoA/3* \
		( ((treevorork4particletype*) (neighwork[(tmp->upperlink)->upperrelated].bp))->pressure +\
		 ((treevorork4particletype*) (neighwork[(tmp           )->lowerrelated].bp))->pressure ) \
	   )

#define VoroRK4_Velocity2D(ai,aj,tmp,neighwork,simpar,MODEL) \
	({Voro2D_point res; \
	 res.x=(0.5-MODEL##_OA(simpar)/3.)*(ai->vx + aj->vx) + MODEL##_OA(simpar)/3* \
	  (((treevorork4particletype*) (neighwork[(tmp->upperlink)->upperrelated].vx))+\
	   ((treevorork4particletype*) (neighwork[(tmp           )->lowerrelated].vx))) \
	 res.y=(0.5-MODEL##_OA(simpar)/3.)*(ai->vy + aj->vy) + MODEL##_OA(simpar)/3* \
	  (((treevorork4particletype*) (neighwork[(tmp->upperlink)->upperrelated].vy))+\
	   ((treevorork4particletype*) (neighwork[(tmp           )->lowerrelated].vy)));\
	   res;}\
	 )

/* M(n,m) macros for treevorostressrk4particletype */
#define VoroStressRK4_Pressure2D(pi,pj,tmp,neighwork,simpar,MODEL) \
	( (0.5-MODEL##_OA(simpar)/3.)*(pj + pi) + MODEL##_OA(simpar)/3* \
	  (((treevorostressrk4particletype*) (neighwork[(tmp->upperlink)->upperrelated].bp))->pressure +\
	   ((treevorostressrk4particletype*) (neighwork[(tmp           )->lowerrelated].bp))->pressure))

/* M(n,m) interpolation of a scalar field for stress particles:
   val_i, val_j = scalar at particles i and j
   val_k+, val_k- = scalar at the two stencil neighbors */
#define VoroStressRK4_Scalar2D(vi,vj,tmp,neighwork,simpar,MODEL,field) \
	( (0.5-MODEL##_OA(simpar)/3.)*(vj + vi) + MODEL##_OA(simpar)/3* \
	  (((treevorostressrk4particletype*) (neighwork[(tmp->upperlink)->upperrelated].bp))->field +\
	   ((treevorostressrk4particletype*) (neighwork[(tmp           )->lowerrelated].bp))->field))
