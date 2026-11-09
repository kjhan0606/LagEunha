#ifndef GFS_NATIVE_CHECKPOINT_IO_H
#define GFS_NATIVE_CHECKPOINT_IO_H
/* Include after eunha.h. Shared native snapshot implementation lives in kh.c;
 * KH wrappers retain their historical names and format. */
void gfs_native_checkpoint_outdata(SimParameters *,int,postype,postype,const char *);
void gfs_native_checkpoint_readdata(SimParameters *,postype *,postype *,int,const char *);
#endif
