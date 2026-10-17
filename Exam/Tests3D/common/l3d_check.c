/* l3d_check: read a LAG3DV1 file with the C reader and print totals.
 *   gcc -O2 -Wall -o l3d_check l3d_check.c -lm
 *   ./l3d_check ic.l3d [ic2.l3d ...]
 *   ./l3d_check -w copy.l3d ic.l3d   (also rewrites ic.l3d with the C writer)
 * Used by selftest.sh to confirm the Python writer and the C reader agree,
 * and by the run scripts as a sanity check on the IC before a job starts. */
#include <math.h>
#include "lag3d_ic.h"

int main(int argc, char **argv){
	int bad = 0, a0 = 1;
	const char *wpath = NULL;
	if(argc > 2 && strcmp(argv[1], "-w") == 0){ wpath = argv[2]; a0 = 3; }
	for(int a=a0;a<argc;a++){
		lag3d_data d;
		int rc = lag3d_read(argv[a], &d);
		if(rc){ fprintf(stderr, "%s: read error %d\n", argv[a], rc); bad = 1; continue; }
		double M=0, px=0, py=0, pz=0, ek=0, ei=0, umin=1e300, mmin=1e300;
		long nout = 0, nnan = 0;
		for(int64_t i=0;i<d.h.np;i++){
			double m = d.mass[i];
			M += m; px += m*d.vx[i]; py += m*d.vy[i]; pz += m*d.vz[i];
			ek += 0.5*m*(d.vx[i]*d.vx[i]+d.vy[i]*d.vy[i]+d.vz[i]*d.vz[i]);
			ei += m*d.u[i];
			if(d.u[i] < umin) umin = d.u[i];
			if(m < mmin) mmin = m;
			if(d.x[i]<d.h.box[0]||d.x[i]>=d.h.box[1]||d.y[i]<d.h.box[2]||d.y[i]>=d.h.box[3]
			   ||d.z[i]<d.h.box[4]||d.z[i]>=d.h.box[5]) nout++;
			if(!isfinite(d.x[i]+d.y[i]+d.z[i]+d.vx[i]+d.vy[i]+d.vz[i]+m+d.u[i])) nnan++;
		}
		printf("%s: np=%lld t=%g gamma=%.6g G=%g eps=%g box=[%g,%g]x[%g,%g]x[%g,%g] per=%d%d%d\n",
			argv[a], (long long)d.h.np, d.h.time, d.h.gamma, d.h.G, d.h.softening,
			d.h.box[0], d.h.box[1], d.h.box[2], d.h.box[3], d.h.box[4], d.h.box[5],
			d.h.periodic[0], d.h.periodic[1], d.h.periodic[2]);
		printf("  M=%.12g P=(%.3e,%.3e,%.3e) Ekin=%.12g Eint=%.12g umin=%.3e mmin=%.3e out_of_box=%ld nonfinite=%ld\n",
			M, px, py, pz, ek, ei, umin, mmin, nout, nnan);
		if(nout || nnan || umin < 0 || mmin <= 0) bad = 1;
		if(wpath && a == a0 && lag3d_write(wpath, &d)){ fprintf(stderr, "write %s failed\n", wpath); bad = 1; }
		lag3d_free(&d);
	}
	return bad;
}
