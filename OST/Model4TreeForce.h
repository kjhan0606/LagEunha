#define NSPLINE (1<<12)

typedef double FORCEARRAYTYPE;

typedef struct DIFFFORCE{
	FORCEARRAYTYPE diff[3];
	FORCEARRAYTYPE slope[3];
} DIFFFORCE;




#define RMAX 6.L
#define RMAXSQ 36.L
#define LOG10RMAX (0.7781512)
#define LOG10RMIN (-3.0)
#define LOG10RMAXmLOG10RMIN (3.7781512)
#define invLOG10RMAXmLOG10RMIN (0.264680)


/*

#ifdef Model4TreeForce
#define EXTERN
#else
#define EXTERN extern
#endif
EXTERN FORCEARRAYTYPE diff[NSPLINE][3], slope[NSPLINE][3];
EXTERN float ran2nran, invran2nran;
#undef EXTERN
*/




#define forcecorrectdiff(a,b) diff[a][b]
#define forcecorrectslope(a,b) slope[a][b]

void i_force_spline(SimParameters *);
void i_potent_spline(SimParameters *);
void GetCellGrav(float , float *, float *, float *);

float GetPointGrav(float);
