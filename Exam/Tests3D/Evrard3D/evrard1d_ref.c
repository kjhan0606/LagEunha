/* evrard1d_ref.c -- independent 1D spherical reference for the Evrard (1988)
 * collapse: Lagrangian staggered-mesh hydro (von Neumann-Richtmyer
 * artificial viscosity, KDK leapfrog, compatible energy update solved by
 * fixed-point iteration) plus exact spherical self-gravity g = G M(<r)/r^2.
 *
 * Setup (Evrard 1988; Springel 2010; Hopkins 2015):
 *   M = R = G = 1, rho = M/(2 pi R^2 r) for r < R, u = 0.05 G M/R, v = 0,
 *   gamma = 5/3, vacuum outside. E_pot(0) = -2/3, E_th(0) = 0.05.
 *
 * Output:
 *   <prefix>_energy.txt   t Ekin Eth Epot Etot
 *   <prefix>_prof_t<T>.txt  r rho v P  at the requested times (default 0.8)
 * This is NOT published reference data. It is a 1D spherically symmetric
 * solution computed here, and selftest.sh checks it against the published
 * t=0.8 profile of HydroCode1D (Vandenbroucke; the SWIFT reference file
 * evrardCollapse3D_exact.txt, fetched by getReference.sh).
 *
 *   gcc -O2 -o evrard1d_ref evrard1d_ref.c -lm
 *   ./evrard1d_ref [N=2000] [tend=3] [prefix=evrard1d] [t_prof ...]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define G 1.0
#define GAMMA (5.0/3.0)
#define CQ 2.0     /* quadratic viscosity coefficient */
#define CL 0.25    /* linear viscosity coefficient */
#define CFL 0.25

static int N;
static double *r, *u, *mn, *dm, *e, *rho, *p, *q, *Mi, *vol;

static void zone_state(void){
	for(int j=1;j<=N;j++){
		vol[j] = 4.0*M_PI/3.0*(r[j]*r[j]*r[j]-r[j-1]*r[j-1]*r[j-1]);
		rho[j] = dm[j]/vol[j];
		p[j] = (GAMMA-1.0)*rho[j]*e[j];
	}
}
static void viscosity(const double *uu){
	for(int j=1;j<=N;j++){
		double du = uu[j]-uu[j-1];
		if(du < 0){
			double c = sqrt(GAMMA*p[j]/rho[j]);
			q[j] = rho[j]*(CQ*du*du + CL*c*(-du));
		} else q[j] = 0;
	}
}
static void accel(double *a){
	for(int i=1;i<=N;i++){
		double pin = p[i]+q[i];
		double pout = (i<N) ? p[i+1]+q[i+1] : 0.0;
		a[i] = 4.0*M_PI*r[i]*r[i]*(pin-pout)/mn[i] - G*Mi[i]/(r[i]*r[i]);
	}
	a[0] = 0;
}
static void energies(double *ek, double *eth, double *ep){
	double k=0, th=0, w=0;
	for(int i=1;i<=N;i++) k += 0.5*mn[i]*u[i]*u[i];
	for(int j=1;j<=N;j++){
		th += dm[j]*e[j];
		/* uniform-density shell: -G int M(r) dm/r, M(r)=M_{j-1}+rho*4pi/3 (r^3-ra^3) */
		double ra=r[j-1], rb=r[j], M0=Mi[j-1], rh=rho[j];
		/* integrand dm = 4 pi rh r^2 dr; exact integral of (M0 + a (r^3-ra^3)) 4pi rh r dr */
		double a = 4.0*M_PI/3.0*rh;
		double I = 4.0*M_PI*rh*( (M0 - a*ra*ra*ra)*0.5*(rb*rb-ra*ra) + a*0.2*(pow(rb,5)-pow(ra,5)) );
		w -= G*I;
	}
	*ek=k; *eth=th; *ep=w;
}
static void dump_profile(const char *prefix, double t){
	char fn[256];
	snprintf(fn, sizeof fn, "%s_prof_t%.3f.txt", prefix, t);
	FILE *fp = fopen(fn, "w");
	fprintf(fp, "# evrard1d_ref N=%d t=%.6f\n# r_zone_centre rho v P\n", N, t);
	for(int j=1;j<=N;j++){
		double rc = cbrt(0.5*(r[j]*r[j]*r[j]+r[j-1]*r[j-1]*r[j-1]));
		fprintf(fp, "%.8e %.8e %.8e %.8e\n", rc, rho[j], 0.5*(u[j]+u[j-1]), p[j]);
	}
	fclose(fp);
}

int main(int argc, char **argv){
	N = (argc>1) ? atoi(argv[1]) : 2000;
	double tend = (argc>2) ? atof(argv[2]) : 3.0;
	const char *prefix = (argc>3) ? argv[3] : "evrard1d";
	int nprof = (argc>4) ? argc-4 : 1;
	double tprof[64]; int done[64];
	if(argc>4) for(int k=0;k<nprof && k<64;k++) tprof[k]=atof(argv[4+k]); else tprof[0]=0.8;
	memset(done, 0, sizeof done);

	r = calloc(N+2, sizeof(double)); u = calloc(N+2, sizeof(double));
	mn = calloc(N+2, sizeof(double)); dm = calloc(N+2, sizeof(double));
	e = calloc(N+2, sizeof(double)); rho = calloc(N+2, sizeof(double));
	p = calloc(N+2, sizeof(double)); q = calloc(N+2, sizeof(double));
	Mi = calloc(N+2, sizeof(double)); vol = calloc(N+2, sizeof(double));
	double *a = calloc(N+2, sizeof(double)), *uh = calloc(N+2, sizeof(double));
	double *e0 = calloc(N+2, sizeof(double)), *pold = calloc(N+2, sizeof(double));
	double *rold = calloc(N+2, sizeof(double)), *unew = calloc(N+2, sizeof(double));

	/* uniform radial spacing, M(<r) = r^2 */
	for(int i=0;i<=N;i++) r[i] = (double)i/N;
	for(int j=1;j<=N;j++){ dm[j] = r[j]*r[j]-r[j-1]*r[j-1]; e[j] = 0.05; }
	Mi[0] = 0;
	for(int i=1;i<=N;i++) Mi[i] = Mi[i-1]+dm[i];
	for(int i=1;i<=N;i++) mn[i] = 0.5*(dm[i] + (i<N ? dm[i+1] : 0.0));
	zone_state();
	for(int j=1;j<=N;j++) q[j]=0;
	accel(a);

	char fn[256];
	snprintf(fn, sizeof fn, "%s_energy.txt", prefix);
	FILE *fe = fopen(fn, "w");
	fprintf(fe, "# evrard1d_ref N=%d  (1D spherical Lagrangian VNR + exact self-gravity)\n# t Ekin Eth Epot Etot\n", N);
	double t = 0, ek, eth, ep, E0;
	energies(&ek, &eth, &ep); E0 = ek+eth+ep;
	fprintf(fe, "%.6e %.10e %.10e %.10e %.10e\n", t, ek, eth, ep, E0);
	double nextout = 0.01;
	long step = 0;
	while(t < tend){
		/* time step */
		double dt = 1e30;
		for(int j=1;j<=N;j++){
			double dr = r[j]-r[j-1];
			double c = sqrt(GAMMA*p[j]/rho[j] + 2.0*q[j]/rho[j]);
			double du = fabs(u[j]-u[j-1]);
			double d1 = dr/(c + du + 1e-300);
			if(d1 < dt) dt = d1;
		}
		dt *= CFL;
		for(int i=1;i<=N;i++){
			double g = G*Mi[i]/(r[i]*r[i]);
			double dg = 0.1*sqrt(r[i]/g);
			if(dg < dt) dt = dg;
		}
		double tnext = tend;
		for(int k=0;k<nprof;k++) if(!done[k] && tprof[k] > t && tprof[k] < tnext) tnext = tprof[k];
		if(nextout > t && nextout < tnext) tnext = nextout;
		if(t+dt > tnext) dt = tnext - t;

		/* Compatible KDK (Caramana et al. 1998 style): the zone energy change
		 * equals minus the work its pressure does on its two nodes with the
		 * same forces and the step-averaged node velocity, so the hydro part
		 * conserves Ekin+Eth to the fixed-point tolerance. Gravity is plain
		 * KDK. */
		for(int i=0;i<=N;i++){ uh[i] = u[i] + 0.5*dt*a[i]; rold[i] = r[i]; }
		uh[0] = 0;
		for(int j=1;j<=N;j++){ e0[j] = e[j]; pold[j] = p[j] + q[j]; }
		for(int i=1;i<=N;i++) r[i] += dt*uh[i];
		for(int i=0;i<=N;i++) unew[i] = uh[i] + 0.5*dt*a[i];
		zone_state();
		viscosity(uh);
		for(int it=0; it<8; it++){
			double maxrel = 0;
			for(int j=1;j<=N;j++){
				double ub_o = 0.5*(u[j]+unew[j]), ub_i = 0.5*(u[j-1]+unew[j-1]);
				double Ao0 = 4*M_PI*rold[j]*rold[j], Ai0 = 4*M_PI*rold[j-1]*rold[j-1];
				double Ao1 = 4*M_PI*r[j]*r[j],       Ai1 = 4*M_PI*r[j-1]*r[j-1];
				double P1 = p[j] + q[j];
				double w = pold[j]*(Ao0*ub_o - Ai0*ub_i) + P1*(Ao1*ub_o - Ai1*ub_i);
				double enew = e0[j] - 0.5*dt*w/dm[j];
				if(enew < 1e-12) enew = 1e-12;
				double rel = fabs(enew-e[j])/(fabs(enew)+1e-300);
				if(rel > maxrel) maxrel = rel;
				e[j] = enew;
			}
			zone_state();
			accel(a);
			for(int i=1;i<=N;i++) unew[i] = uh[i] + 0.5*dt*a[i];
			if(maxrel < 1e-13) break;
		}
		/* forces at the start of the next step use this p and q, and the
		 * energy update above used exactly A(r_old)*(p_old+q_old) for the
		 * first half kick, so the ledger closes. */
		for(int i=1;i<=N;i++) u[i] = unew[i];
		t += dt; step++;

		for(int k=0;k<nprof;k++) if(!done[k] && fabs(t-tprof[k]) < 1e-12){ dump_profile(prefix, t); done[k]=1; }
		if(fabs(t-nextout) < 1e-12 || t >= tend){
			energies(&ek, &eth, &ep);
			fprintf(fe, "%.6e %.10e %.10e %.10e %.10e\n", t, ek, eth, ep, ek+eth+ep);
			nextout += 0.01;
		}
	}
	fclose(fe);
	energies(&ek, &eth, &ep);
	printf("evrard1d_ref N=%d t=%.3f steps=%ld E0=%.8f E=%.8f dE/|E0|=%.3e\n",
	       N, t, step, E0, ek+eth+ep, (ek+eth+ep-E0)/fabs(E0));
	return 0;
}
