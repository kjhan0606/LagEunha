/* lag3d_ic.h -- reader/writer for the LAG3DV1 IC/snapshot format used by
 * Exam/Tests3D (layout documented in common/lag3d_io.py and README.md).
 *
 * Header-only, plain C99, no dependency on eunha.h, so a future 3D driver
 * (and the box-side checker common/l3d_check.c) can include it unchanged.
 * u[] is SPECIFIC internal energy; a driver that follows the 2D convention
 * must set ie = mass*u when it loads the file.
 */
#ifndef LAG3D_IC_H
#define LAG3D_IC_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#define LAG3D_MAGIC "LAG3DV1"
#define LAG3D_HEADER_SIZE 128
#define LAG3D_FLAG_POT 1
#define LAG3D_FLAG_VOL 2

typedef struct {
	int32_t version, flags;
	int64_t np;
	double time, gamma;
	double box[6];          /* xmin xmax ymin ymax zmin zmax */
	int32_t periodic[3];
	double G, softening;
} lag3d_header;

typedef struct {
	lag3d_header h;
	double *x, *y, *z, *vx, *vy, *vz, *mass, *u, *rho, *vol, *pot;
	int64_t *id;
} lag3d_data;

static inline int lag3d_rd(void *p, size_t sz, size_t n, FILE *fp){
	return fread(p, sz, n, fp) == n ? 0 : -1;
}

/* Returns 0 on success. Caller frees with lag3d_free(). */
static inline int lag3d_read(const char *path, lag3d_data *d){
	memset(d, 0, sizeof(*d));
	FILE *fp = fopen(path, "rb");
	if(!fp) return -1;
	char magic[8];
	int32_t reserved;
	double dres;
	int err = 0;
	err |= lag3d_rd(magic, 1, 8, fp);
	if(err || memcmp(magic, LAG3D_MAGIC, 8) != 0){ fclose(fp); return -2; }
	err |= lag3d_rd(&d->h.version, 4, 1, fp);
	err |= lag3d_rd(&d->h.flags, 4, 1, fp);
	err |= lag3d_rd(&d->h.np, 8, 1, fp);
	err |= lag3d_rd(&d->h.time, 8, 1, fp);
	err |= lag3d_rd(&d->h.gamma, 8, 1, fp);
	err |= lag3d_rd(d->h.box, 8, 6, fp);
	err |= lag3d_rd(d->h.periodic, 4, 3, fp);
	err |= lag3d_rd(&reserved, 4, 1, fp);
	err |= lag3d_rd(&d->h.G, 8, 1, fp);
	err |= lag3d_rd(&d->h.softening, 8, 1, fp);
	err |= lag3d_rd(&dres, 8, 1, fp);
	if(err || d->h.np < 0){ fclose(fp); return -3; }
	size_t n = (size_t)d->h.np;
	double **f[10] = {&d->x,&d->y,&d->z,&d->vx,&d->vy,&d->vz,&d->mass,&d->u,&d->rho,&d->vol};
	for(int k=0;k<10 && !err;k++){
		*f[k] = (double*)malloc(sizeof(double)*(n?n:1));
		err |= (*f[k]==NULL) ? -1 : lag3d_rd(*f[k], 8, n, fp);
	}
	if(!err && (d->h.flags & LAG3D_FLAG_POT)){
		d->pot = (double*)malloc(sizeof(double)*(n?n:1));
		err |= lag3d_rd(d->pot, 8, n, fp);
	}
	if(!err){
		d->id = (int64_t*)malloc(sizeof(int64_t)*(n?n:1));
		err |= lag3d_rd(d->id, 8, n, fp);
	}
	fclose(fp);
	return err ? -4 : 0;
}

static inline void lag3d_free(lag3d_data *d){
	free(d->x); free(d->y); free(d->z); free(d->vx); free(d->vy); free(d->vz);
	free(d->mass); free(d->u); free(d->rho); free(d->vol); free(d->pot); free(d->id);
	memset(d, 0, sizeof(*d));
}

/* Writer, for a 3D driver's snapshots. Same layout as the reader. */
static inline int lag3d_write(const char *path, const lag3d_data *d){
	FILE *fp = fopen(path, "wb");
	if(!fp) return -1;
	int32_t reserved = 0;
	double dres = 0;
	size_t n = (size_t)d->h.np;
	fwrite(LAG3D_MAGIC "", 1, 8, fp);
	fwrite(&d->h.version, 4, 1, fp);
	fwrite(&d->h.flags, 4, 1, fp);
	fwrite(&d->h.np, 8, 1, fp);
	fwrite(&d->h.time, 8, 1, fp);
	fwrite(&d->h.gamma, 8, 1, fp);
	fwrite(d->h.box, 8, 6, fp);
	fwrite(d->h.periodic, 4, 3, fp);
	fwrite(&reserved, 4, 1, fp);
	fwrite(&d->h.G, 8, 1, fp);
	fwrite(&d->h.softening, 8, 1, fp);
	fwrite(&dres, 8, 1, fp);
	const double *f[10] = {d->x,d->y,d->z,d->vx,d->vy,d->vz,d->mass,d->u,d->rho,d->vol};
	for(int k=0;k<10;k++) fwrite(f[k], 8, n, fp);
	if(d->h.flags & LAG3D_FLAG_POT) fwrite(d->pot, 8, n, fp);
	fwrite(d->id, 8, n, fp);
	int bad = ferror(fp);
	fclose(fp);
	return bad ? -2 : 0;
}
#endif
