#define G 6.672E-8L
#define motherrank 0
#define pc 3.08567802E18L
#define pi 3.1415926535L
#define pow2(A) ((A)*(A))
#define pow3(A) ((A)*(A)*(A))
#define MAX(A,B) ((A)>(B) ? (A) : (B))
#define MIN(A,B) ((A)<(B) ? (A) : (B))
#define MAX_PEAK_NUM 100000
#define MIN_HALO_SEED_NUM 10
#define MAX_SURROUNDING_NUM 100000
#define NSTEP 5
#define BOUND_ITER 4
/*                                        */
/*
#define SNGL_CPU_PTL_NUM_MAX 150000
*/
#define SNGL_CPU_PTL_NUM_MAX 100000000
#define SERIAL_MAXIMUM 10000

#ifdef OLD_DEPTH
#define NDEPTH 20000000
#endif

#define REPEAT_POTENT_NUM 4
#define NODE_SIZE 19
/* #define MAX_BOX_MEMORY 91125000 *//* 450^3 grid points->318 */
#define MAX_BOX_MEMORY 421875000  /* 750^3 grid points->318 */
/*#define MAX_BOX_MEMORY 343000000 */ /* 750^3 grid points->318 */

/*
#define TIDAL_STRIDE_NUM 2000
*/
#define TIDAL_STRIDE_NUM 1000
#define NO_EXIT -999999
#define FOF_LINK 0.05
/* ?? ??Àº ?????? ?Á·Î¼????? ???? serial ?Ï°? tidal radius ??
 * ?????Ï±? À§??  mpeak ?? Á¶?? ???Ì´?.  */
#define SCATTER_TIDAL_CAL 10000000
/*                                     */
/*  tag numbers and status number for communications      */
#define READY 1
#define ParallelPotential 88
#define INIT_TOTAL_PARTICLE_BCAST 1234
#define ParallelTidal 444
#define PARALLEL_FOF_HFIND 7654
#define FoF_Report 7655
#define Tidal_Report 204
#define FREE_IN_TIDAL 777
#define FREE_INITIAL_TIDAL 888
#define TOTAL_PARTICLE_BCAST 10000
#define WRITING 3
#define NP_TAG 17
#define R_TAG 32
#define VR_TAG 83
#define INDX_TAG 24
#define PT2H_TAG 55
#define MPEAK_TAG 96
#define DEN_TAG 100
#define DEN_TAG1 101
#define DEN_TAG2 102
#define DEN_TAG3 103
#define DEN_TAG4 104
#define DEN_TAG5 105
#define DEN_TAG6 106
#define DEN_TAG7 107
#define PEAK_DEN_TAG 108
/*                                     */
#define INIT_NP 1000000

#define OUTPUTFILE "HFIND.DAT"
//int NODE_NUM_SIZE,NODE_SIZE_STRIDE;
static int NODE_BUFFER=1;
typedef int lint;
/*
typedef struct _basicparticletype{
	float x,y,z;
	float vx,vy,vz;
	indxtype indx;
}basicparticletype;
*/

struct PtlPos{
	double x,y,z;
};
struct HaloPos{
	double x,y,z;
};
struct PtlVel{
	float vx;float vy;float vz;
};
struct HaloVel{
	double vx,vy,vz;
};
void final_out_file();
void initial_out_file();
#define NUM_MASS 10000
