#include "eunha.h"
#include "voro.h"
#include "raytrace.h"

postype getS3D(Voro3D_point *p, Voro3D_point *q, Voro3D_point *r, Voro3D_point *k){
	Voro3D_point pp = Vec3DMean(q,p);
	Voro3D_point n = Vec3DSub(q,p);
	Voro3D_point ppmr = Vec3DSub(&pp,r);
	postype den = Vec3DDotP(&n,&ppmr);
	postype nom = Vec3DDotP(&n,k);
	postype res;
	if(nom !=0) res = den/nom;
	else {
		res = INFINITY;
	}
	return res;
}
postype getS2D(Voro2D_point *p, Voro2D_point *q, Voro2D_point *r, Voro2D_point *k){
	Voro2D_point pp = Vec2DMean(q,p);
	Voro2D_point n = Vec2DSub(q,p);
	Voro2D_point ppmr = Vec2DSub(&pp,r);
	postype den = Vec2DDotP(&n,&ppmr);
	postype nom = Vec2DDotP(&n,k);
	postype res;
	if(nom !=0) res = den/nom;
	else {
		res = INFINITY;
	}
	return res;
}
