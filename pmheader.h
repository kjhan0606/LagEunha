

#ifdef XYZDBL
typedef double PosType;
#else
typedef float PosType;
#endif




#ifdef XYZDBL
#	ifdef INDEX
#		define XofP(p) fmod((double)((p)->x)+(double)(((p)->indx%simpar.mx))+simpar.nx,simpar.nx)
#		define YofP(p) fmod((double)((p)->y)+(double)((((p)->indx%simpar.mxmy)/simpar.mx))+simpar.ny,simpar.ny)
#		define ZofP(p) fmod((double)((p)->z)+(double)(((p)->indx/simpar.mxmy))+simpar.nz,simpar.nz)
#		define Xoftype(type,p) fmod((double)(((type*)p)->x)+(double)((((type*)p)->indx%simpar.mx))+simpar.nx,simpar.nx)
#		define Yoftype(type,p) fmod((double)(((type*)p)->y)+(double)(((((type*)p)->indx%simpar.mxmy)/simpar.mx))+simpar.ny,simpar.ny)
#		define Zoftype(type,p) fmod((double)(((type*)p)->z)+(double)((((type*)p)->indx/simpar.mxmy))+simpar.nz,simpar.nz)
#	else
#		error the pair definition of XYZDBL and !INDEX  is not allowed in the current version.
#	endif
#else
#	define XofP(p) ((p)->x)
#	define YofP(p) ((p)->y)
#	define ZofP(p) ((p)->z)
#	define Xoftype(type,p) (((type*)p)->x)
#	define Yoftype(type,p) (((type*)p)->y)
#	define Zoftype(type,p) (((type*)p)->z)
#endif
#define MofP(p,lvsimpar) (lvsimpar.mass)


