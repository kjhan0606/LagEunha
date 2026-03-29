#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#ifdef VORO_MAIN
#include "../indxflag.h"
#endif
//#include "vorowrapper.h"



#define MIN(a,b) ( (a)<(b)?(a):(b))
#define MAX(a,b) ( (a)>(b)?(a):(b))


/*
#define EPS (1.e-11L)
#define EPS2_InOut (1.e-6L)
*/

#define EPS (1.e-9)
#define EPS2_InOut (1.e-6)

#ifdef VORO_SNGL
typedef float postype;
#define MPI_POSTYPE MPI_FLOAT
#else
typedef double postype;
#define MPI_POSTYPE MPI_DOUBLE
#endif

#define FMOD(a,b) fmod(a,b)

#define Half (0.5L)
#define None -1

#define Outside 1
#define Onborder -1
#define Inside -1

#define HalfLife 0
#define Active 1
#define Inactive -1

#define Yes 1
#define No -1


typedef struct Voro3D_tensor{
	postype xx,xy,xz,yx,yy,yz,zx,zy,zz;
} Voro3D_tensor;


typedef struct Voro3D_point{
	postype x,y,z,dist2,csound,pressure;
	postype w2,drad2,w2ceil; // Laguerre square radius 
	int indx,type;
	struct Voro3D_point *sibling;
	indxflag u4if;
	void *bp;
} Voro3D_point;


typedef struct Voro2D_point{
	postype x,y,dist2,csound,pressure;
	postype w2,drad2,w2ceil; // Laguerre square radius 
	int indx,type;
	struct Voro2D_point *sibling;
	indxflag u4if;
	void *bp;
} Voro2D_point;

typedef struct Voro3D_Vertex{
	int indx, status,related[3];
	struct Voro3D_Vertex *link[3];
	int considered[3];
	int ilink, facecount[3];
	postype x,y,z;
}Voro3D_Vertex;

typedef struct Voro3D_face{
	int indx;
	Voro3D_Vertex LeftVertex,RightVertex;
	Voro3D_point neighbor;
}Voro3D_face;


typedef struct Voro3D_Vertex2{
	Voro3D_Vertex *vertex,*link;
	int ilink;
} Voro3D_Vertex2;


typedef struct Voro2D_Corner{
	/* lowerrelated: the neighbor id for the lower link
	 * upperrelated: the neighbor id for the upper link */
	int indx, status,lowerrelated,upperrelated;
	postype x,y;
	/* upperlink: a pointer to the upper voronoi corner */
	/* lowerlink: a pointer to the lower voronoi corner */
	struct Voro2D_Corner *lowerlink,*upperlink;
	Voro2D_point *neighlowerlink, *neighupperlink;
}Voro2D_Corner;



typedef struct Voro3D_HessianParticle{
    int indx;
    postype x,y,z,mass,den;
	Voro3D_point denGrad;
	Voro3D_tensor denHessian;
    postype volume;
    struct Voro3D_HessianParticle *next;
} Voro3D_HessianParticle;


int Voro3D_CreateOneVertex(Voro3D_Vertex *, 
		Voro3D_Vertex *,
		Voro3D_Vertex *, 
		int , 
		Voro3D_point *, 
		int , 
		postype );
int sort2Ddist(const void *, const void *);
int sort2Ddrad2(const void *, const void *);
int sort3Ddrad2(const void *, const void *);
postype get3Dw2pCeil(Voro3D_Vertex *, Voro3D_point *, int );
postype get2Dw2pCeil(Voro2D_Corner *, Voro2D_point *);
Voro2D_Corner Voro2D_FindIntersection(Voro2D_Corner *, Voro2D_Corner *, Voro2D_point *, postype , postype *);
Voro3D_Vertex Voro3D_FindIntersection(
		Voro3D_Vertex *, 
		Voro3D_Vertex *, 
		Voro3D_point *,
		postype );
int Voro2D_FindVC(Voro2D_point *, Voro2D_point *, Voro2D_point *, int , Voro2D_Corner *, int , postype);
int Voro2D_InitializeCorner(Voro2D_Corner *,int , postype );
int Voro2D_CutOrStay(Voro2D_point *, Voro2D_Corner *, postype);
Voro2D_Corner *FindingLowerLimit(Voro2D_point *, Voro2D_Corner *, postype);
Voro2D_Corner *FindingUpperLimit(Voro2D_point *, Voro2D_Corner *, postype);
void InactivateLowerAndUpper(Voro2D_Corner *, Voro2D_Corner *);
int CreateNewCorner(Voro2D_Corner *, Voro2D_Corner *, Voro2D_Corner *, int , Voro2D_point *, int, postype );
int Voro2D_Trim_Corner(Voro2D_Corner *, int , postype *);
postype Area2DPolygon(Voro2D_Corner *, int );

int sort3Ddist(const void *a, const void *b);
int Voro3D_CutOrStay(Voro3D_point *, Voro3D_Vertex *, postype, postype *);

int Voro3D_FindVC(Voro3D_point *, Voro3D_point *, Voro3D_point *, int , Voro3D_Vertex *, int , postype, int , int *);
int Voro3D_InitializeVertex(Voro3D_Vertex *,int , postype );
Voro3D_Vertex2 FindNearestBoundaryVertex(Voro3D_point *, Voro3D_Vertex *, postype);
int InactivateOutsideVertex(Voro3D_point *, Voro3D_Vertex *, int , postype);
int CreateNewVertices(Voro3D_point *, int , Voro3D_Vertex2 *, Voro3D_Vertex *, int , postype);
int Voro3D_Trim_Vertex(Voro3D_Vertex *, int, postype * , int *);
int Voro3D_FaceExtract(Voro3D_Vertex *, int ,int , Voro3D_point *);
void voro3D_EulerRot(Voro3D_point *, int, postype , postype , postype );
Voro3D_point voro3D_norm(Voro3D_point *);
Voro2D_point voro2D_norm(Voro2D_point *);

Voro2D_Corner Centroid2DPolygon(Voro2D_Corner *);
int Voro2D_LineExtract(Voro2D_Corner *, int ,Voro2D_point *);
postype getWeight3D(postype , postype , postype , postype );

postype Voro3D_Volume_Polyhedron(Voro3D_Vertex *,int );
postype Voro3D_pyramid_volume(Voro3D_Vertex *, Voro3D_Vertex *);
Voro3D_point Voro3D_norm_polygon(Voro3D_Vertex *, Voro3D_Vertex *);

Voro3D_point Voro3D_denGrad_Polyhedron(Voro3D_Vertex *,int , int , Voro3D_HessianParticle *, Voro3D_point *);
Voro3D_tensor Voro3D_denHessian_Polyhedron(Voro3D_Vertex *,int , int , Voro3D_HessianParticle *, Voro3D_point *);

Voro3D_point findPolyhedronCentroid(Voro3D_Vertex *, int );
void findLine2D(Voro2D_point *, Voro2D_point *, postype , Voro2D_point *, postype *);
Voro2D_point calculateIntersection2D(postype , postype , postype , postype );
void findNewCorner(Voro2D_Corner *, Voro2D_point *, Voro2D_point *);

Voro3D_point calculateIntersection3D( Voro3D_point *, postype , Voro3D_point *, postype , Voro3D_point *, postype );
void findNewVertex(Voro3D_Vertex *, Voro3D_point *, Voro3D_point *);

void findPlane3D(Voro3D_point *, Voro3D_point *, postype , Voro3D_point *, postype *);


postype getAvgPressureOnSurface3D(postype , postype , Voro3D_point *, Voro3D_Vertex *, Voro3D_Vertex *);




#define gIDMakingFace3D(neighbor, vertex,idir) (neighbor[vertex->related[(idir+2)%3]].indx)


#define Vec3DInnP(a,b) dotProduct3D(a,b)


#define dotProduct2D(a,b) ((a)->x*(b)->x + (a)->y*(b)->y)

#define Vec2DInnP(a,b) dotProduct2D(a,b)

#define Vec2DDotP(a,b) dotProduct2D(a,b)

#define EunhaVec2DSub(a,b) ({Voro2D_point _xx; _xx.x = (a)->x-(b)->x; \
		_xx.y = (a)->y-(b)->y; _xx.indx = (a)->u4if.indx; _xx;})

#define Vec2DSub(a,b) ({Voro2D_point _xx; _xx.x = (a)->x-(b)->x; \
		_xx.y = (a)->y-(b)->y; _xx.indx = (a)->indx; _xx;})

#define Vec2DSubCorner(a,b) ({Voro2D_Corner _xx; _xx.x = (a)->x-(b)->x; \
		_xx.y = (a)->y-(b)->y; _xx.indx = (a)->indx; _xx;})

#define Vec2DMulAdd(t,a,b) ({Voro2D_Corner _xx; _xx.x = t*(a)->x+(b)->x;\
		_xx.y = t*(a)->y+(b)->y; _xx;})

#define Vec2DMean(a,b) ({Voro2D_point _xx; _xx.x = 0.5*((a)->x+(b)->x);\
        _xx.y = 0.5*((a)->y+(b)->y); _xx;})


#define Vec3DDotP(a,b) dotProduct3D(a,b)

#define dotProduct3D(a,b) ((a)->x*(b)->x + (a)->y*(b)->y + (a)->z*(b)->z)

#define Vec3DNorm(a) ({postype amp;{ amp=sqrt(Vec3DDotP(a,a)); (a)->x = (a)->x/amp;\
		(a)->y /= amp; (a)->z /= amp;})

#define Vec3DAmp(a) ({postype amp;{ amp=sqrt(Vec3DDotP(a,a)); amp;})

#define Vec3DSub(a,b) ({Voro3D_point _xx; _xx.x = (a)->x-(b)->x;\
		_xx.y = (a)->y-(b)->y; _xx.z=(a)->z-(b)->z; _xx.indx = (a)->indx; _xx;})

#define Vec3DMulAdd(t,a,b) ({Voro3D_Vertex _xx; _xx.x = t*(a)->x+(b)->x;\
		_xx.y = t*(a)->y+(b)->y; _xx.z=t*(a)->z+(b)->z; _xx;})

#define AddPoint3D(a,b) ({Voro3D_point _xx; _xx.x = (a)->x+(b)->x;\
		_xx.y = (a)->y+(b)->y; _xx.z=(a)->z+(b)->z; _xx;})

#define Vec3DCross(a,b) ({Voro3D_point _xx; _xx.x = (a)->y*(b)->z - (a)->z*(b)->y;\
		_xx.y = (a)->z*(b)->x - (a)->x*(b)->z; _xx.z=(a)->x*(b)->y - (a)->y*(b)->x; _xx;})

#define Vec3DMean(a,b) ({Voro3D_point _xx; _xx.x = 0.5*((a)->x+(b)->x);\
        _xx.y = 0.5*((a)->y+(b)->y); _xx.z=0.5*((a)->z+(b)->z); _xx;})

#define Vec2DCross(a,b) ({Voro2D_point _xx; _xx.x = 0; _xx.y=0; _xx.z=(a)->x*(b)->y - (a)->y*(b)->x; _xx;})

#define DumpTo3DPoint(a,b) do{\
		(a)->x = (b)->x; (a)->y = (b)->y; (a)->z = (b)->z;\
	} while(0)

#define Vec2DLength(a,b) sqrt(((a)->x-(b)->x)*((a)->x-(b)->x) + ((a)->y-(b)->y)*((a)->y-(b)->y))

#define Vec2DLineCenter(a,b) ( { Voro2D_point _xx; _xx.x = 0.5*((a)->x + (b)->x);\
		_xx.y = 0.5*((a)->y + (b)->y);_xx; })

#define getNextPolygonVertex(now,next) ( { int _i=0; while(next->link[_i]!=now) { _i = (_i+1)%3; } _i = (_i-1+3)%3; Voro3D_Vertex *nnext = next->link[_i]; nnext;})


#define getwfrac(sqw1,sqw2,sqd12) (0.5 +0.5*((sqw1-sqw2)/sqd12))

#define initStreeUpq(ai) do{\
	Stress *stress=ai->stress;\
	stress->gUxx = stress->gUxy = stress->gUxz = 0;\
	stress->gUyx = stress->gUyy = stress->gUyz = 0;\
	stress->gUzx = stress->gUzy = stress->gUzz = 0;\
	stress->gU33 = 0;\
} while(0)

#define getStressUpq(ai,aj,dS) do {\
	Stress *stress = ai->stress;\
	stress->gUxx += Half*(aj->vx-ai->vx)*dS.x;\
	stress->gUxy += Half*(aj->vx-ai->vx)*dS.y;\
	stress->gUxz += Half*(aj->vx-ai->vx)*dS.z;\
	stress->gUxy += Half*(aj->vy-ai->vy)*dS.x;\
	stress->gUyy += Half*(aj->vy-ai->vy)*dS.y;\
	stress->gUyz += Half*(aj->vy-ai->vy)*dS.z;\
	stress->gUxz += Half*(aj->vz-ai->vz)*dS.x;\
	stress->gUyz += Half*(aj->vz-ai->vz)*dS.y;\
	stress->gUzz += Half*(aj->vz-ai->vz)*dS.z;\
	stress->gU33 += Half*(aj->vx-ai->vx)*dS.x;\
	stress->gU33 += Half*(aj->vy-ai->vy)*dS.y;\
	stress->gU33 += Half*(aj->vz-ai->vz)*dS.z;\
} while(0)

#define getNeighborPressure(jbp,ibp) ( IS_FLAG_ON(jbp,Wallflag) ?  (ibp->pressure): (jbp->pressure))
