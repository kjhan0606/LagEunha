/**
 * @file laguerre.h
 * @brief Header for Laguerre Diagram (Power Diagram / Weighted Voronoi) Construction
 * 
 * This library implements Laguerre tessellation in 2D and 3D.
 * A Laguerre diagram partitions space based on "power distance":
 *   power_distance(point, site) = ||point - site||^2 - weight^2
 * 
 * The radical hyperplane between two weighted sites divides space where
 * power distances to both sites are equal.
 */

#ifndef LAGUERRE_H
#define LAGUERRE_H

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>

/* Stub definition for indxflag if not available */
#ifdef LAGUERRE_MAIN
  #ifdef HAS_INDXFLAG
    #include "../indxflag.h"
  #else
    typedef struct { int indx; int flag; } indxflag;
  #endif
#else
  typedef struct { int indx; int flag; } indxflag;
#endif

/*============================================================================
 * NUMERIC CONSTANTS AND MACROS
 *============================================================================*/

#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

#define EPSILON             (1.0e-9)   /**< Numerical tolerance */
#define BOUNDARY_TOLERANCE  (1.0e-6)   /**< Inside/outside classification tolerance */
#define FMOD(a, b)          fmod(a, b)
#define HALF                (0.5L)

/*============================================================================
 * PRECISION CONFIGURATION
 *============================================================================*/

#ifdef LAGUERRE_SINGLE_PRECISION
    typedef float  real_t;
    #define MPI_REAL_T MPI_FLOAT
#else
    typedef double real_t;
    #define MPI_REAL_T MPI_DOUBLE
#endif

/*============================================================================
 * STATUS AND CLASSIFICATION CONSTANTS
 *============================================================================*/

#define RELATION_NONE     -1   /**< No relationship */
#define POINT_INSIDE      -1   /**< Point inside cell */
#define POINT_BOUNDARY    -1   /**< Point on boundary */
#define POINT_OUTSIDE      1   /**< Point outside cell */
#define STATUS_INACTIVE   -1   /**< Vertex/corner inactive */
#define STATUS_HALFLIFE    0   /**< Transition state */
#define STATUS_ACTIVE      1   /**< Vertex/corner active */
#define RESULT_NO         -1
#define RESULT_YES         1

/*============================================================================
 * 3D DATA STRUCTURES
 *============================================================================*/

/** @brief 3x3 tensor for Hessian matrix */
typedef struct Tensor3D {
    real_t xx, xy, xz;
    real_t yx, yy, yz;
    real_t zx, zy, zz;
} Tensor3D;

/**
 * @brief 3D point/site for Laguerre diagram
 */
typedef struct Point3D {
    real_t x, y, z;           /**< Coordinates */
    real_t squaredDist;       /**< Squared distance from reference */
    real_t soundSpeed;        /**< Speed of sound (CFD) */
    real_t pressure;          /**< Pressure value */
    real_t weightSquared;     /**< w^2: Laguerre weight squared */
    real_t deltaRadiusSq;     /**< For priority sorting */
    real_t weightCeiling;     /**< Upper bound for neighbor search */
    int    index;             /**< Unique identifier */
    int    type;              /**< Classification */
    struct Point3D *sibling;  /**< Linked list pointer */
    indxflag indexFlags;
    void   *backPointer;      /**< User data */
} Point3D;

/**
 * @brief 3D vertex of Laguerre cell (polyhedron corner)
 */
typedef struct Vertex3D {
    int index;
    int status;
    int neighborSites[3];         /**< 3 defining neighbor indices (-n = boundary) */
    struct Vertex3D *adjacent[3]; /**< 3 connected vertices */
    int visited[3];
    int currentLink;
    int faceVertexCount[3];
    real_t x, y, z;
} Vertex3D;

/** @brief 3D face structure */
typedef struct Face3D {
    int index;
    Vertex3D leftVertex, rightVertex;
    Point3D neighbor;
} Face3D;

/** @brief Vertex traversal reference */
typedef struct VertexRef3D {
    Vertex3D *vertex;
    Vertex3D *linkedVertex;
    int linkIndex;
} VertexRef3D;

/** @brief Particle for Hessian computation */
typedef struct HessianParticle3D {
    int index;
    real_t x, y, z, mass, density;
    Point3D densityGradient;
    Tensor3D densityHessian;
    real_t volume;
    struct HessianParticle3D *next;
} HessianParticle3D;

/*============================================================================
 * 2D DATA STRUCTURES
 *============================================================================*/

/** @brief 2D point/site */
typedef struct Point2D {
    real_t x, y;
    real_t squaredDist;
    real_t soundSpeed;
    real_t pressure;
    real_t weightSquared;
    real_t deltaRadiusSq;
    real_t weightCeiling;
    int index, type;
    struct Point2D *sibling;
    indxflag indexFlags;
    void *backPointer;
} Point2D;

/** @brief 2D corner (polygon vertex) */
typedef struct Corner2D {
    int index;
    int status;
    int lowerNeighborIdx;
    int upperNeighborIdx;
    real_t x, y;
    struct Corner2D *lowerLink;
    struct Corner2D *upperLink;
    Point2D *lowerNeighborPtr;
    Point2D *upperNeighborPtr;
} Corner2D;

/*============================================================================
 * VECTOR OPERATION MACROS
 *============================================================================*/

#define dot3D(a, b) ((a)->x*(b)->x + (a)->y*(b)->y + (a)->z*(b)->z)
#define dot2D(a, b) ((a)->x*(b)->x + (a)->y*(b)->y)

#define subtract3D(a, b) ({ \
    Point3D _r = {0}; _r.x=(a)->x-(b)->x; _r.y=(a)->y-(b)->y; _r.z=(a)->z-(b)->z; \
    _r.index=(a)->index; _r; })

#define subtract2D(a, b) ({ \
    Point2D _r = {0}; _r.x=(a)->x-(b)->x; _r.y=(a)->y-(b)->y; _r.index=(a)->index; _r; })

#define subtractCorner2D(a, b) ({ \
    Corner2D _r = {0}; _r.x=(a)->x-(b)->x; _r.y=(a)->y-(b)->y; _r.index=(a)->index; _r; })

#define multiplyAdd3D(t, a, b) ({ \
    Vertex3D _r = {0}; _r.x=(t)*(a)->x+(b)->x; _r.y=(t)*(a)->y+(b)->y; _r.z=(t)*(a)->z+(b)->z; _r; })

#define multiplyAdd2D(t, a, b) ({ \
    Corner2D _r = {0}; _r.x=(t)*(a)->x+(b)->x; _r.y=(t)*(a)->y+(b)->y; _r; })

#define cross3D(a, b) ({ \
    Point3D _r = {0}; _r.x=(a)->y*(b)->z-(a)->z*(b)->y; _r.y=(a)->z*(b)->x-(a)->x*(b)->z; \
    _r.z=(a)->x*(b)->y-(a)->y*(b)->x; _r; })

#define add3D(a, b) ({ \
    Point3D _r = {0}; _r.x=(a)->x+(b)->x; _r.y=(a)->y+(b)->y; _r.z=(a)->z+(b)->z; _r; })

#define midpoint3D(a, b) ({ \
    Point3D _r = {0}; _r.x=0.5*((a)->x+(b)->x); _r.y=0.5*((a)->y+(b)->y); _r.z=0.5*((a)->z+(b)->z); _r; })

#define midpoint2D(a, b) ({ \
    Point2D _r = {0}; _r.x=0.5*((a)->x+(b)->x); _r.y=0.5*((a)->y+(b)->y); _r; })

#define copyPoint3D(dest, src) do { (dest)->x=(src)->x; (dest)->y=(src)->y; (dest)->z=(src)->z; } while(0)

#define distance2D(a, b) sqrt(((a)->x-(b)->x)*((a)->x-(b)->x)+((a)->y-(b)->y)*((a)->y-(b)->y))

/** @brief Weight fraction: t = 0.5 + 0.5*(w1^2-w2^2)/d^2 */
#define computeWeightFraction(w1Sq, w2Sq, distSq) (0.5 + 0.5*((w1Sq)-(w2Sq))/(distSq))

#define getNextFaceVertex(current, next) ({ \
    int _i=0; while((next)->adjacent[_i]!=(current)){_i=(_i+1)%3;} _i=(_i-1+3)%3; \
    Vertex3D *_n=(next)->adjacent[_i]; _n; })

#define getFaceId(neighbors, vertex, dir) ((neighbors)[(vertex)->neighborSites[((dir)+2)%3]].index)

/*============================================================================
 * FUNCTION DECLARATIONS
 *============================================================================*/

/* Transformation */
void applyEulerRotation3D(Point3D *points, int count, real_t theta, real_t phi, real_t psi);
Point3D computeSurfaceNormal3D(Point3D *vertices);
Point2D computeNormal2D(Point2D *direction);

/* Initialization */
int initBoundingBox2D(Corner2D *corners, int maxCorners, real_t boxSize);
int initBoundingBox3D(Vertex3D *vertices, int maxVertices, real_t boxSize);

/* Classification */
int classifyCornerPosition2D(Point2D *neighbor, Corner2D *corner, real_t weightFrac);
int classifyVertexPosition3D(Point3D *neighbor, Vertex3D *vertex, real_t weightFrac, real_t *offset);

/* Intersection */
Corner2D findEdgePlaneIntersection2D(Corner2D *a, Corner2D *b, Point2D *neighbor, real_t wf, real_t *t);
Vertex3D findEdgePlaneIntersection3D(Vertex3D *a, Vertex3D *b, Point3D *neighbor, real_t wf);
Point2D solveLinearSystem2D(real_t s1, real_t i1, real_t s2, real_t i2);
Point3D solveLinearSystem3D(Point3D *n1, real_t d1, Point3D *n2, real_t d2, Point3D *n3, real_t d3);

/* Cell Construction */
int constructLaguerreCell2D(Point2D *ctr, Point2D *in, Point2D *nb, int np, Corner2D *c, int mp, real_t box);
int constructLaguerreCell3D(Point3D *ctr, Point3D *in, Point3D *nb, int np, Vertex3D *v, int mp, real_t box, int sh, int *w);

/* Compaction */
int compactActiveCorners2D(Corner2D *corners, int count, real_t *maxDistSq);
int compactActiveVertices3D(Vertex3D *vertices, int count, real_t *maxDistSq, int *map);
void deactivateCornersBetween(Corner2D *lower, Corner2D *upper);

/* Geometry */
real_t computePolygonArea2D(Corner2D *start, int numCorners);
Corner2D computePolygonCentroid2D(Corner2D *start);
real_t computePolyhedronVolume3D(Vertex3D *vertices, int numVertices);
real_t computePyramidVolume3D(Vertex3D *start, Vertex3D *next);
Point3D computePolyhedronCentroid3D(Vertex3D *vertices, int numVertices);

/* Plane/Line Definition */
void computeRadicalPlane3D(Point3D *s1, Point3D *s2, real_t t, Point3D *n, real_t *d);
void computeRadicalLine2D(Point2D *s1, Point2D *s2, real_t t, Point2D *n, real_t *d);
void updateVertexFromDelaunay(Vertex3D *v, Point3D *neighbors, Point3D *center);

/* Weight Ceiling */
real_t computeMaxWeightCeiling2D(Corner2D *start, Point2D *neighbors);
real_t computeMaxWeightCeiling3D(Vertex3D *v, Point3D *nb, int nv);

/* Sorting */
int compareByDistanceSquared2D(const void *a, const void *b);
int compareByDistanceSquared3D(const void *a, const void *b);
int compareByDeltaRadius2D(const void *a, const void *b);
int compareByDeltaRadius3D(const void *a, const void *b);

/* Utility */
real_t computeLaguerreWeight3D(real_t x, real_t y, real_t z, real_t h);
int extractFaceVertices3D(Vertex3D *v, int nv, int fid, Point3D *out);
int extractPolygonVertices2D(Corner2D *c, int nc, Point2D *out);

/*============================================================================
 * BACKWARD COMPATIBILITY ALIASES
 *============================================================================*/

typedef real_t postype;
typedef Point3D Voro3D_point;
typedef Point2D Voro2D_point;
typedef Vertex3D Voro3D_Vertex;
typedef Corner2D Voro2D_Corner;
typedef Tensor3D Voro3D_tensor;
typedef VertexRef3D Voro3D_Vertex2;
typedef Face3D Voro3D_face;
typedef HessianParticle3D Voro3D_HessianParticle;

#define EPS EPSILON
#define EPS2_InOut BOUNDARY_TOLERANCE
#define None RELATION_NONE
#define Outside POINT_OUTSIDE
#define Inside POINT_INSIDE
#define Onborder POINT_BOUNDARY
#define Active STATUS_ACTIVE
#define Inactive STATUS_INACTIVE
#define HalfLife STATUS_HALFLIFE
#define Yes RESULT_YES
#define No RESULT_NO

#define Vec3DInnP(a,b) dot3D(a,b)
#define Vec2DInnP(a,b) dot2D(a,b)
#define dotProduct3D(a,b) dot3D(a,b)
#define dotProduct2D(a,b) dot2D(a,b)
#define Vec3DSub(a,b) subtract3D(a,b)
#define Vec2DSub(a,b) subtract2D(a,b)
#define Vec3DMulAdd(t,a,b) multiplyAdd3D(t,a,b)
#define Vec2DMulAdd(t,a,b) multiplyAdd2D(t,a,b)
#define Vec3DCross(a,b) cross3D(a,b)
#define AddPoint3D(a,b) add3D(a,b)
#define Vec3DMean(a,b) midpoint3D(a,b)
#define Vec2DMean(a,b) midpoint2D(a,b)
#define DumpTo3DPoint(a,b) copyPoint3D(a,b)
#define Vec2DLength(a,b) distance2D(a,b)
#define getwfrac(w1,w2,d) computeWeightFraction(w1,w2,d)
#define getNextPolygonVertex(c,n) getNextFaceVertex(c,n)
#define gIDMakingFace3D(nb,v,d) getFaceId(nb,v,d)

#define voro3D_EulerRot applyEulerRotation3D
#define voro3D_norm computeSurfaceNormal3D
#define voro2D_norm computeNormal2D
#define Voro2D_InitializeCorner initBoundingBox2D
#define Voro3D_InitializeVertex initBoundingBox3D
#define Voro2D_CutOrStay classifyCornerPosition2D
#define Voro3D_CutOrStay classifyVertexPosition3D
#define Voro2D_FindIntersection findEdgePlaneIntersection2D
#define Voro3D_FindIntersection findEdgePlaneIntersection3D
#define calculateIntersection2D solveLinearSystem2D
#define calculateIntersection3D solveLinearSystem3D
#define Voro2D_FindVC constructLaguerreCell2D
#define Voro3D_FindVC constructLaguerreCell3D
#define Voro2D_Trim_Corner compactActiveCorners2D
#define Voro3D_Trim_Vertex compactActiveVertices3D
#define InactivateLowerAndUpper deactivateCornersBetween
#define Area2DPolygon computePolygonArea2D
#define Centroid2DPolygon computePolygonCentroid2D
#define Voro3D_Volume_Polyhedron computePolyhedronVolume3D
#define Voro3D_pyramid_volume computePyramidVolume3D
#define findPolyhedronCentroid computePolyhedronCentroid3D
#define findPlane3D computeRadicalPlane3D
#define findLine2D computeRadicalLine2D
#define findNewVertex updateVertexFromDelaunay
#define get2Dw2pCeil computeMaxWeightCeiling2D
#define get3Dw2pCeil computeMaxWeightCeiling3D
#define sort2Ddist compareByDistanceSquared2D
#define sort3Ddist compareByDistanceSquared3D
#define sort2Ddrad2 compareByDeltaRadius2D
#define sort3Ddrad2 compareByDeltaRadius3D
#define getWeight3D computeLaguerreWeight3D
#define Voro3D_FaceExtract extractFaceVertices3D
#define Voro2D_LineExtract extractPolygonVertices2D

#endif /* LAGUERRE_H */
