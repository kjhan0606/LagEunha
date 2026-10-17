typedef struct Ray3D{
	Voro3D_point r,k;
} Ray3D;

typedef struct Ray2D{
	Voro2D_point r,k;
} Ray2D;

typedef struct BundleRays2D{
	int nrays;
	Ray2D *ray;
} BundleRays2D;


// prototypes of raytracing routines
postype getS2D(Voro2D_point *, Voro2D_point *, Voro2D_point *, Voro2D_point *);
postype getS3D(Voro3D_point *, Voro3D_point *, Voro3D_point *, Voro3D_point *);
