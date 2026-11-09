#ifndef GFS_MATERIAL_MPI_H
#define GFS_MATERIAL_MPI_H
#include <stddef.h>
#ifdef __cplusplus
extern "C" {
#endif
/* Active only inside the distributed material operator. Owned rows come
 * from the existing k-section DD, never from ID modulo rank. */
int gfs_material_mpi_rank(void);
void *gfs_material_mpi_view_begin(void *simulation,void **ticket);
void gfs_material_mpi_view_end(void *ticket,int allow_boundary_events);
int gfs_material_mpi_rows(int n,const long long *ids,unsigned char *owned);
void gfs_material_mpi_bcast(void *data,size_t bytes);
void gfs_material_mpi_merge(void *rows,size_t stride,int n,const unsigned char *owned);
void gfs_material_mpi_maxima(double *values,int n);
void gfs_material_mpi_counters(unsigned long long *values,int n);
#ifdef __cplusplus
}
#endif
#endif
